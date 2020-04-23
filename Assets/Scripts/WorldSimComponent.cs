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
		
	};


	public WorldData WorldData;
	[HideInInspector] public StaticState StaticState;

	[HideInInspector] public Icosphere Icosphere;
	[HideInInspector] public bool CollectOverlay = false;
	[HideInInspector] public float TimeScale = 0;

	public int CellCount { get; private set; }
	public ref SimState ActiveSimState { get { return ref _simStates[_activeSimState]; } }
	public ref TempState TempState { get { return ref _tempState; } }
	public float TimeTillTick { get; private set; }
	public float InverseCellCount { get; private set; }

	[SerializeField]
	private WorldSim _worldSim;
	private WorldGenData _worldGenData = new WorldGenData();
	private SimState[] _simStates;
	private TempState _tempState;
	private int _activeSimState;
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
		_simStates = new SimState[_simStateCount];
		for (int i = 0; i < _simStateCount; i++)
		{
			_simStates[i] = new SimState();
			_simStates[i].Init(CellCount, ref WorldData);
		}

		_worldGenData = JsonUtility.FromJson<WorldGenData>(WorldGenAsset.text);
		WorldGen.Generate(Seed, _worldGenData, Icosphere, ref WorldData, ref StaticState, ref _simStates[0], ref _tempState);

		OnTick?.Invoke();
	}

	public void OnDestroy()
	{
		_tickJobHandle.Complete();
		_worldSim.Dispose(ref WorldData);
		StaticState.Dispose();
		Icosphere.Dispose();
		_tempState.Dispose(ref WorldData);
		for (int i = 0; i < _simStateCount; i++)
		{
			_simStates[i].Dispose();
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
			// TODO: make there be 3 sim states and don't overwrite the one that we are using for display
			_tickJobHandle = _worldSim.Tick(
				_simStates,
				_simStateCount,
				_ticksToAdvance,
				ref _tempState,
				ref StaticState,
				ref WorldData,
				ref SimSettings,
				ref _activeSimState);
			JobHandle.ScheduleBatchedJobs();
			_ticksToAdvance = 0;
		}
	}


	public delegate void EditSimStateFunc(ref SimState lastState, ref SimState nextState);
	public void Edit(EditSimStateFunc editSimState)
	{
		ref var lastState = ref _simStates[_activeSimState];
		_activeSimState = (_activeSimState + 1) % _simStateCount;
		ref var nextState = ref _simStates[_activeSimState];

		var tempStateJobHandle = _tempState.Clear(StaticState.Count, ref WorldData, default(JobHandle));
		tempStateJobHandle.Complete();

		editSimState(ref lastState, ref nextState);

		List<NativeArray<float>> tempArrays = new List<NativeArray<float>>();
		TempState.Update(_worldSim.SimJob, ref nextState, ref TempState, ref WorldData, ref StaticState, default(JobHandle)).Complete();
		foreach (var i in tempArrays)
		{
			i.Dispose();
		}

		OnTick?.Invoke();
	}

	public void StepTime()
	{
		TimeScale = 0;
		TimeTillTick = 0;
	}


}
