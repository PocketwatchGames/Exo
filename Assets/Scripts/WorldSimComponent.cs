using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Jobs;
using Unity.Collections;

public class WorldSimComponent : MonoBehaviour
{

	public int Seed;
	public int Subdivisions = 5;
	public bool CheckForDegeneracy;
	public bool CollectGlobals = false;
	public bool LogState;
	public int LogStateIndex;
	public TextAsset WorldGenAsset;
	public TextAsset WorldDataAsset;
	public WorldData WorldData;
	public StaticState StaticState;
	public Icosphere Icosphere;
	public float TimeScale;

	public int CellCount { get; private set; }
	public ref SimState ActiveSimState { get { return ref _simStates[_activeSimState]; } }
	public ref DependentState DependentState { get { return ref _dependentState; } }
	public ref DisplayState DisplayState { get { return ref _displayState; } }
	public float TimeTillTick { get; private set; }

	public float InverseCellCount { get; private set; }

	private WorldGenData _worldGenData = new WorldGenData();

	[SerializeField]
	private WorldSim _worldSim;
	private const int _simStateCount = 2;
	private SimState[] _simStates;
	private DependentState _dependentState;
	private DisplayState _displayState;
	private int _activeSimState;
	private float _ticksPerSecond = 1;

	public delegate void TickEventHandler();
	public event TickEventHandler OnTick;

	public void Awake()
    {
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
			_simStates[i].Init(CellCount, WorldData.AirLayers, WorldData.WaterLayers);
		}

		_worldGenData = JsonUtility.FromJson<WorldGenData>(WorldGenAsset.text);
		WorldGen.Generate(Seed, _worldGenData, Icosphere, ref WorldData, ref StaticState, ref _simStates[0], ref _dependentState);

		UpdateDisplayIncomplete(CellCount, ref DisplayState, ref ActiveSimState, ref DependentState, ref WorldData);
	}

	public void OnDestroy()
	{
		_worldSim.Dispose();
		StaticState.Dispose();
		Icosphere.Dispose();
		_dependentState.Dispose();
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
				TimeTillTick += _ticksPerSecond;
				iterations++;
			}
			Tick(ref _simStates[_activeSimState], iterations);

		}
	}

	private void Tick(ref SimState state, int ticksToAdvance)
	{
		bool degen = _worldSim.Tick(
			_simStates,
			_simStateCount,
			ticksToAdvance, 
			ref _dependentState, 
			ref _displayState, 
			ref StaticState, 
			ref WorldData, 
			ref _activeSimState,
			CheckForDegeneracy,
			LogState, 
			LogStateIndex,
			CollectGlobals, 
			CollectGlobals);

		OnTick?.Invoke();

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

		editSimState(ref lastState, ref nextState);

		UpdateDisplayIncomplete(CellCount, ref DisplayState, ref nextState, ref DependentState, ref WorldData);

		OnTick?.Invoke();
	}

	public void StepTime()
	{
		TimeScale = 0;
		TimeTillTick = 0;
	}


	static public void UpdateDisplayIncomplete(int cellCount, ref DisplayState displayState, ref SimState activeSimState, ref DependentState dependentState, ref WorldData worldData)
	{
		displayState.Dispose();

		displayState = new DisplayState();
		displayState.Init(cellCount, worldData.AirLayers, worldData.WaterLayers, worldData.LayerCount);

		var tempAdvectionDestination = new NativeArray<BarycentricValueVertical>(cellCount, Allocator.TempJob);
		var tempVelocity = new NativeArray<float3>(cellCount, Allocator.TempJob);
		var tempCloudMass = new NativeArray<float>(cellCount, Allocator.TempJob);
		JobHandle initDisplayHandle = default(JobHandle);
		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			var initDisplayJob = new InitDisplayAirLayerJob()
			{
				DisplayPressure = displayState.Pressure[i],
				DisplayPressureGradientForce = displayState.PressureGradientForce[i],
				DisplayCondensationCloud = displayState.CondensationCloud,
				DisplayCondensationGround = displayState.CondensationGround,
				Enthalpy = displayState.EnthalpyAir[i],
				WindVertical = displayState.WindVertical[i],
				DustCoverage = displayState.DustMass,
				CarbonDioxidePercent = displayState.CarbonDioxidePercent[i],

				CarbonDioxide = activeSimState.AirCarbon[i],
				Gravity = activeSimState.PlanetState.Gravity,
				AirTemperaturePotential = activeSimState.AirTemperaturePotential[i],
				AirPressure = dependentState.AirPressure[i],
				LayerMiddle = dependentState.LayerMiddle[i],
				PressureGradientForce = tempVelocity,
				CondensationCloud = tempCloudMass,
				CondensationGround = tempCloudMass,
				AirMass = dependentState.AirMass[i],
				VaporMass = activeSimState.AirVapor[i],
				DustMass = activeSimState.Dust[i],
				AdvectionDestination = tempAdvectionDestination,
			};
			initDisplayHandle = JobHandle.CombineDependencies(initDisplayHandle, initDisplayJob.Schedule(cellCount, 1, initDisplayHandle));
		}
		initDisplayHandle.Complete();
		tempAdvectionDestination.Dispose();
		tempVelocity.Dispose();
		tempCloudMass.Dispose();
	}

}
