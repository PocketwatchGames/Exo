using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Jobs;
using Unity.Collections;
using System;

public class WorldSimComponent : MonoBehaviour
{
	private const int _simStateCount = 2;

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
	public ref DisplayState DisplayState { get { return ref _displayState; } }
	public float TimeTillTick { get; private set; }
	public float InverseCellCount { get; private set; }

	[SerializeField]
	private WorldSim _worldSim;
	private WorldGenData _worldGenData = new WorldGenData();
	private SimState[] _simStates;
	private TempState _tempState;
	private DisplayState _displayState;
	private int _activeSimState;
	private Action _prepNextFrame;

	public delegate void TickEventHandler(JobHandle tickJobHandle);
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

		UpdateDisplayIncomplete(CellCount, ref DisplayState, ref ActiveSimState, ref TempState, ref StaticState, ref WorldData, ref SimSettings);
	}

	public void OnDestroy()
	{
		_worldSim.Dispose(ref WorldData);
		StaticState.Dispose();
		Icosphere.Dispose();
		_tempState.Dispose(ref WorldData, _prepNextFrame);
		_displayState.Dispose();
		for (int i = 0; i < _simStateCount; i++)
		{
			_simStates[i].Dispose();
		}
	}

	public void Update()
	{
		if (TimeTillTick > -1)
		{
			TimeTillTick -= Time.deltaTime * TimeScale;
		}
		if (TimeTillTick <= 0)
		{
			int iterations = 0;
			while (TimeTillTick <= 0)
			{
				TimeTillTick += TicksPerRealSecond;
				iterations++;
			}
			Tick(ref _simStates[_activeSimState], iterations);

		}
	}

	private void Tick(ref SimState state, int ticksToAdvance)
	{
		bool degen;
		JobHandle tickJobHandle;
		_worldSim.Tick(
			out tickJobHandle,
			out degen,
			_simStates,
			_simStateCount,
			ticksToAdvance, 
			ref _tempState, 
			ref _displayState, 
			ref StaticState, 
			ref WorldData, 
			ref SimSettings,
			ref _activeSimState,
			ref _prepNextFrame);

		OnTick?.Invoke(tickJobHandle);

		if (degen)
		{
			TimeScale = 0;
		}
	}

	public delegate void EditSimStateFunc(ref SimState lastState, ref SimState nextState);
	public void Edit(EditSimStateFunc editSimState)
	{
		ref var lastState = ref _simStates[_activeSimState];
		_activeSimState = (_activeSimState + 1) % _simStateCount;
		ref var nextState = ref _simStates[_activeSimState];

		_prepNextFrame?.Invoke();

		editSimState(ref lastState, ref nextState);

		List<NativeArray<float>> tempArrays = new List<NativeArray<float>>();
		TempState.Update(_worldSim.SimJob, ref nextState, ref TempState, ref WorldData, default(JobHandle), tempArrays).Complete();
		foreach (var i in tempArrays)
		{
			i.Dispose();
		}

		UpdateDisplayIncomplete(CellCount, ref DisplayState, ref nextState, ref TempState, ref StaticState, ref WorldData, ref SimSettings);

		OnTick?.Invoke(default(JobHandle));
	}

	public void StepTime()
	{
		TimeScale = 0;
		TimeTillTick = 0;
	}


	static public void UpdateDisplayIncomplete(int cellCount, ref DisplayState displayState, ref SimState activeSimState, ref TempState dependentState, ref StaticState staticState, ref WorldData worldData, ref SimSettings settings)
	{
		displayState.Dispose();

		displayState = new DisplayState();
		displayState.Init(cellCount, ref worldData);

		DisplayState.Update(ref displayState, ref displayState, ref worldData, ref dependentState, ref activeSimState, ref staticState, ref settings).Complete();
	}

}
