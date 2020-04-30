using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Jobs;
using Unity.Collections;
using System;

public class WorldSimComponent : MonoBehaviour
{
	private const int _simStateCount = 3;
	private const int _tempStateCount = 2;

	[Header("Init")]
	public int Seed;
	public int Subdivisions = 5;
	public TextAsset WorldGenAsset;
	public TextAsset WorldDataAsset;
	public float TicksPerRealSecond = 1;

	[Header("Simulation Features")]
	public SimSettings SimSettings = new SimSettings()
	{
		MakeAirIncompressible = true,
		MakeWaterIncompressible = true,
		WaterSurfaceFlowEnabled = true,

		Condensation = true,
		Evaporation = true,
		Flora = true,
		Freezing = true,
		IceMelting = true,
		Plankton = true,
		Precipitation = true,
		SoilRespiration = true,
		GroundWater = true,

		ConductionAirFlora = true,
		ConductionAirIce = true,
		ConductionAirTerrain = true,
		ConductionAirWater = true,
		ConductionFloraTerrain = true,
		ConductionIceFlora = true,
		ConductionIceTerrain = true,
		ConductionIceWater = true,
		ConductionWaterTerrain = true,

		IncompressibilityIterations = 20,

		CheckForDegeneracy = false,
		CollectGlobals = false,
		CollectOverlay = false,
		LogState = false,
		LogStateIndex = 0,
		
	};


	public WorldData WorldData;
	[HideInInspector] public StaticState StaticState;

	[HideInInspector] public Icosphere Icosphere;
	[HideInInspector] public bool CollectOverlay = false;
	[HideInInspector] public float TimeScale = 0;

	public int CellCount { get; private set; }
	public ref SimState LastSimState { get { return ref _simStates[_lastSimState]; } }
	public ref TempState LastTempState { get { return ref _tempStates[_lastTempState]; } }
	public float TimeTillTick { get; private set; }
	public float InverseCellCount { get; private set; }

	[SerializeField]
	private WorldSim _worldSim;
	private WorldGenData _worldGenData = new WorldGenData();
	private SimState[] _simStates;
	private TempState[] _tempStates;
	private int _lastSimState;
	private int _lastTempState;
	private int _activeSimState;
	private int _activeTempState;
	private JobHandle _tickJobHandle;

	public delegate void TickEventHandler();
	public event TickEventHandler OnTick;

	public void Awake()
    {
		TimeScale = 0;
		WorldData = JsonUtility.FromJson<WorldData>(WorldDataAsset.text);
		WorldData.Init();

		Icosphere = new Icosphere(Subdivisions);
		CellCount = Icosphere.Vertices.Length;
		InverseCellCount = 1.0f / CellCount;

		_worldSim = new WorldSim(CellCount, ref WorldData);

		TimeTillTick = 0.00001f;
		_activeSimState = 0;
		_lastSimState = 0;
		_activeTempState = 0;
		_lastTempState = 0;
		_simStates = new SimState[_simStateCount];
		for (int i = 0; i < _simStateCount; i++)
		{
			_simStates[i] = new SimState();
			_simStates[i].Init(CellCount, ref WorldData);
		}
		_tempStates = new TempState[_tempStateCount];
		for (int i = 0; i < _tempStateCount; i++)
		{
			_tempStates[i] = new TempState();
			_tempStates[i].Init(CellCount, ref WorldData);
		}

		_worldGenData = JsonUtility.FromJson<WorldGenData>(WorldGenAsset.text);
		WorldGen.Generate(Seed, _worldGenData, Icosphere, ref WorldData, ref StaticState, ref _simStates[0], ref _tempStates[0]);

		OnTick?.Invoke();
	}

	public void OnDestroy()
	{
		_tickJobHandle.Complete();
		_worldSim.Dispose(ref WorldData);
		StaticState.Dispose();
		Icosphere.Dispose();
		for (int i = 0; i < _simStateCount; i++)
		{
			_simStates[i].Dispose();
		}
		for (int i = 0; i < _tempStateCount; i++)
		{
			_tempStates[i].Dispose(ref WorldData);
		}
	}

	private int _ticksToAdvance = 0;
	public void Update()
	{
		if (TimeTillTick > -1)
		{
			TimeTillTick -= Time.deltaTime * TimeScale;
		}
		if (TimeTillTick <= 0)
		{
			while (TimeTillTick <= 0)
			{
				TimeTillTick += TicksPerRealSecond;
				_ticksToAdvance++;
			}

			_tickJobHandle.Complete();

			if (SimSettings.LogState)
			{
				CellInfo.PrintState("State", SimSettings.LogStateIndex, ref StaticState, ref _simStates[_activeSimState], ref WorldData, new List<string>());
				CellInfo.PrintDependentState("Dependent Vars", SimSettings.LogStateIndex, ref _tempStates[_activeTempState], ref WorldData);
			}


			bool degenerate = false;
			if (SimSettings.CheckForDegeneracy)
			{
				SortedSet<int> degenIndices = new SortedSet<int>();
				List<string> degenVarNames = new List<string>();
				if (SimTest.CheckForDegeneracy(
					CellCount, 
					ref _simStates[_activeSimState], 
					ref _simStates[_lastTempState], 
					ref StaticState, 
					ref _tempStates[_activeTempState], 
					ref WorldData,
					ref degenIndices,
					ref degenVarNames))
				{
					foreach (var i in degenIndices)
					{
						CellInfo.PrintState("Degenerate", i, ref StaticState, ref _simStates[_activeSimState], ref WorldData, degenVarNames);
						CellInfo.PrintDependentState("Dependent Vars", i, ref _tempStates[_activeTempState], ref WorldData);
						CellInfo.PrintState("Last State", i, ref StaticState, ref _simStates[_lastTempState], ref WorldData, new List<string>());
					}
					TimeScale = 0;
					Debug.Break();
				}
			}

			_lastSimState = _activeSimState;
			_lastTempState = _activeTempState;
			OnTick?.Invoke();


			//if (degenerate)
			//{
			//	TimeScale = 0;
			//}

		}
	}

	public void LateUpdate()
	{
		if (_ticksToAdvance > 0)
		{
			
			_tickJobHandle = _worldSim.Tick(
				_ticksToAdvance,
				_simStates,
				_tempStates,
				ref StaticState,
				ref WorldData,
				ref SimSettings,
				ref _activeSimState,
				ref _activeTempState);
			JobHandle.ScheduleBatchedJobs();
			_ticksToAdvance = 0;
		}
	}


	public delegate void EditSimStateFunc(ref SimState lastState, ref SimState nextState);
	public void Edit(EditSimStateFunc editSimState)
	{
		_tickJobHandle.Complete();
		_lastSimState = _activeSimState;
		_lastTempState = _activeTempState;

		ref var lastState = ref _simStates[_activeSimState];
		_activeSimState = (_activeSimState + 1) % _simStateCount;
		ref var nextState = ref _simStates[_activeSimState];

		editSimState(ref lastState, ref nextState);

		TempState.Update(_worldSim.SimJob, ref nextState, ref _tempStates[_lastTempState], ref WorldData, ref StaticState, default(JobHandle)).Complete();
		OnTick?.Invoke();
	}

	public void StepTime()
	{
		TimeScale = 0;
		TimeTillTick = 0;
	}


}
