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

		_worldSim = new WorldSim(CellCount, WorldData.AirLayers, WorldData.WaterLayers);

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

		_displayState = new DisplayState();
		// TODO: get rid of this constant
		const int totalLayersHack = 15;
		_displayState.Init(CellCount, WorldData.AirLayers, WorldData.WaterLayers, totalLayersHack);

		var tempAdvectionDestination = new NativeArray<BarycentricValueVertical>(CellCount, Allocator.TempJob);
		var tempVelocity = new NativeArray<float3>(CellCount, Allocator.TempJob);
		var tempCloudMass = new NativeArray<float>(CellCount, Allocator.TempJob);
		JobHandle initDisplayHandle = default(JobHandle);
		for (int i = 1; i < WorldData.AirLayers - 1; i++)
		{
			var initDisplayJob = new InitDisplayAirLayerJob()
			{
				DisplayPressure = _displayState.Pressure[i],
				DisplayPressureGradientForce = _displayState.PressureGradientForce[i],
				DisplayCondensationCloud = _displayState.CondensationCloud,
				DisplayCondensationGround = _displayState.CondensationGround,
				Enthalpy = _displayState.EnthalpyAir[i],
				WindVertical = _displayState.WindVertical[i],
				DustCoverage = _displayState.DustMass,
				CarbonDioxidePercent = _displayState.CarbonDioxidePercent[i],

				CarbonDioxide = ActiveSimState.AirCarbonDioxide[i],
				Gravity = ActiveSimState.PlanetState.Gravity,
				AirTemperaturePotential = ActiveSimState.AirTemperaturePotential[i],
				AirPressure = DependentState.AirPressure[i],
				LayerMiddle = DependentState.LayerMiddle[i],
				PressureGradientForce = tempVelocity,
				CondensationCloud = tempCloudMass, 
				CondensationGround = tempCloudMass, 
				AirMass = _dependentState.AirMass[i],
				VaporMass = ActiveSimState.AirVapor[i],
				DustMass = ActiveSimState.Dust[i],
				AdvectionDestination = tempAdvectionDestination,
			};
			initDisplayHandle = JobHandle.CombineDependencies(initDisplayHandle, initDisplayJob.Schedule(CellCount, 1, initDisplayHandle));
		}
		initDisplayHandle.Complete();
		tempAdvectionDestination.Dispose();
		tempVelocity.Dispose();
		tempCloudMass.Dispose();

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
		bool degen = _worldSim.Tick(_simStates, _simStateCount, ticksToAdvance, ref _dependentState, ref _displayState, ref StaticState, ref WorldData, ref _activeSimState, CheckForDegeneracy, LogState, LogStateIndex, CollectGlobals);

		OnTick?.Invoke();

		if (degen)
		{
			TimeScale = 0;
		}
	}


	public void StepTime()
	{
		TimeScale = 0;
		TimeTillTick = 0;
	}


}
