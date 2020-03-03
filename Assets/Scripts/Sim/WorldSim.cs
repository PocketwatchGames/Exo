//#define EnergyTerrainJobDebug
//#define ConductionWaterBottomJobDebug
//#define ConductionWaterTerrainJobDebug
//#define SolarRadiationAbsorbedTerrainJobDebug
//#define EnergyAirJobDebug
//#define EnergyWaterSurfaceJobDebug
//#define DiffusionAirJobDebug
#define AdvectionAirJobDebug
//#define PressureGradientForceAirJobDebug
#define AdvectionCloudJobDebug
//#define WaterFrictionJobDebug
//#define WaterDensityGradientForceJobDebug
//#define StateChangeJobDebug
//#define EnergyIceJobDebug
//#define StateChangeAirLayerJobDebug
//#define FluxCloudJobDebug
//#define GetVectorDestCoordsJobDebug
//#define DiffusionCloudJobDebug

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using Unity.Burst;
using UnityEngine;

[Serializable]
public class WorldSim {

	public JobHelper SolarRadiationInJob;
	public JobHelper UpdateTerrainJob;
	public JobHelper EmissivityAirJob;
	public JobHelper EmissivityWaterJob;
	public JobHelper EmissivityTerrainJob;
	public JobHelper SolarRadiationAbsorbedAirJob;
	public JobHelper SolarRadiationAbsorbedIceJob;
	public JobHelper SolarRadiationAbsorbedWaterJob;
	public JobHelper SolarRadiationAbsorbedTerrainJob;
	public JobHelper ThermalOutCloudJob;
	public JobHelper ThermalOutAirJob;
	public JobHelper ThermalOutIceJob;
	public JobHelper ThermalOutWaterJob;
	public JobHelper ThermalOutTerrainJob;
	public JobHelper ThermalInUpAirJob;
	public JobHelper ThermalInUpIceJob;
	public JobHelper ThermalInUpWaterJob;
	public JobHelper ThermalInDownAirJob;
	public JobHelper ThermalInDownIceJob;
	public JobHelper ThermalInDownWaterJob;
	public JobHelper ThermalInDownTerrainJob;
	public JobHelper PressureGradientForceAirJob;
	public JobHelper PressureGradientForceCloudJob;
	public JobHelper WaterDensityGradientForceJob;
	public JobHelper WindFrictionJob;
	public JobHelper WaterFrictionJob;
	public JobHelper ConductionAirIceJob;
	public JobHelper ConductionAirWaterJob;
	public JobHelper ConductionAirTerrainJob;
	public JobHelper ConductionIceWaterJob;
	public JobHelper ConductionIceTerrainJob;
	public JobHelper ConductionWaterTerrainJob;
	public JobHelper EnergyTerrainJob;
	public JobHelper FluxCloudJob;
	public JobHelper FluxAirJob;
	public JobHelper FluxWaterJob;
	public JobHelper FluxIceJob;
	public JobHelper EnergyIceJob;
	public JobHelper EnergyWaterJob;
	public JobHelper EnergyAirJob;
	public JobHelper StateChangeJob;
	public JobHelper StateChangeAirLayerJob;
	public JobHelper DiffusionAirJob;
	public JobHelper DiffusionWaterJob;
	public JobHelper DiffusionCloudJob;
	public JobHelper AdvectionAirJob;
	public JobHelper AdvectionWaterJob;
	public JobHelper AdvectionCloudJob;
	public JobHelper ApplyAdvectionAirJob;
	public JobHelper ApplyAdvectionWaterJob;
	public JobHelper ApplyAdvectionCloudJob;
	public JobHelper UpdateDependentWaterLayerJob;
	public JobHelper UpdateDependentAirLayerJob;
	public JobHelper UpdateDependentJob;
	public JobHelper GetVectorDestCoordsJob;

	private int _cellCount;
	private int _batchCount = 100;
	private const int _terrainLayers = 2;

	private int _airLayers;
	private int _waterLayers;
	private int _layerCount;
	private int _terrainLayer;
	private int _waterLayer0;
	private int _iceLayer;
	private int _airLayer0;
	private int _surfaceWaterLayer;

	private NativeArray<float> solarRadiation;
	private NativeArray<float> waterSlopeAlbedo;
	private NativeArray<float>[] emissivity;
	private NativeArray<float>[] solarRadiationIn;
	private NativeArray<DiffusionAir>[] diffusionAir;
	private NativeArray<DiffusionWater>[] diffusionWater;
	private NativeArray<DiffusionAir>[] advectionAir;
	private NativeArray<DiffusionWater>[] advectionWater;
	private NativeArray<DiffusionCloud> diffusionCloud;
	private NativeArray<DiffusionCloud> advectionCloud;
	private NativeArray<BarycentricValue> destinationCloud;
	private NativeArray<BarycentricValue>[] destinationAir;
	private NativeArray<BarycentricValue>[] destinationWater;
	private NativeArray<float>[] latentHeat;
	private NativeArray<float3>[] pressureGradientForce;
	private NativeArray<float3> pressureGradientForceCloud;
	private NativeArray<float> windFriction;
	private NativeArray<float3>[] waterForce;
	private NativeArray<float> conductionSurfaceAirIce;
	private NativeArray<float> conductionSurfaceAirWater;
	private NativeArray<float> conductionSurfaceAirTerrain;
	private NativeArray<float> conductionUpperAirIce;
	private NativeArray<float> conductionUpperAirWater;
	private NativeArray<float> conductionUpperAirTerrain;
	private NativeArray<float> conductionIceWater;
	private NativeArray<float> conductionIceTerrain;
	private NativeArray<float>[] conductionWaterTerrain;
	private NativeArray<float>[] fluxEnergyWater;
	private NativeArray<float>[] fluxEnergyAir;
	private NativeArray<float>[] evaporationMass;
	private NativeArray<float>[] frozenMass;
	private NativeArray<float>[] frozenTemperature;
	private NativeArray<float> fluxEnergyIce;

	private NativeArray<float> displaySolarRadiation;

	public WorldSim(int cellCount, int airLayers, int waterLayers)
	{
		_cellCount = cellCount;
		_airLayers = airLayers;
		_waterLayers = waterLayers;
		_layerCount = _airLayers + _waterLayers + _terrainLayers;
		_terrainLayer = 0;
		_waterLayer0 = _terrainLayer + 1;
		_iceLayer = _waterLayer0 + _waterLayers;
		_airLayer0 = _iceLayer + 1;
		_surfaceWaterLayer = _waterLayers - 2;

		#region Job Initialization
		DiffusionAirJob = new JobHelper(_cellCount);
		DiffusionWaterJob = new JobHelper(_cellCount);
		DiffusionCloudJob = new JobHelper(_cellCount);
		SolarRadiationInJob = new JobHelper(_cellCount);
		UpdateTerrainJob = new JobHelper(_cellCount);
		EmissivityAirJob = new JobHelper(_cellCount);
		EmissivityWaterJob = new JobHelper(_cellCount);
		EmissivityTerrainJob = new JobHelper(_cellCount);
		SolarRadiationAbsorbedAirJob = new JobHelper(_cellCount);
		SolarRadiationAbsorbedIceJob = new JobHelper(_cellCount);
		SolarRadiationAbsorbedWaterJob = new JobHelper(_cellCount);
		SolarRadiationAbsorbedTerrainJob = new JobHelper(_cellCount);
		ThermalOutCloudJob = new JobHelper(_cellCount);
		ThermalOutAirJob = new JobHelper(_cellCount);
		ThermalOutIceJob = new JobHelper(_cellCount);
		ThermalOutWaterJob = new JobHelper(_cellCount);
		ThermalOutTerrainJob = new JobHelper(_cellCount);
		ThermalInUpAirJob = new JobHelper(_cellCount);
		ThermalInUpIceJob = new JobHelper(_cellCount);
		ThermalInUpWaterJob = new JobHelper(_cellCount);
		ThermalInDownAirJob = new JobHelper(_cellCount);
		ThermalInDownIceJob = new JobHelper(_cellCount);
		ThermalInDownWaterJob = new JobHelper(_cellCount);
		ThermalInDownTerrainJob = new JobHelper(_cellCount);
		PressureGradientForceAirJob = new JobHelper(_cellCount);
		PressureGradientForceCloudJob = new JobHelper(_cellCount);
		WaterDensityGradientForceJob = new JobHelper(_cellCount);
		WindFrictionJob = new JobHelper(_cellCount);
		WaterFrictionJob = new JobHelper(_cellCount);
		ConductionAirIceJob = new JobHelper(_cellCount);
		ConductionAirWaterJob = new JobHelper(_cellCount);
		ConductionAirTerrainJob = new JobHelper(_cellCount);
		ConductionIceWaterJob = new JobHelper(_cellCount);
		ConductionIceTerrainJob = new JobHelper(_cellCount);
		ConductionWaterTerrainJob = new JobHelper(_cellCount);
		EnergyTerrainJob = new JobHelper(_cellCount);
		FluxCloudJob = new JobHelper(_cellCount);
		FluxWaterJob = new JobHelper(_cellCount);
		FluxAirJob = new JobHelper(_cellCount);
		FluxIceJob = new JobHelper(_cellCount);
		EnergyIceJob = new JobHelper(_cellCount);
		EnergyWaterJob = new JobHelper(_cellCount);
		EnergyAirJob = new JobHelper(_cellCount);
		StateChangeJob = new JobHelper(_cellCount);
		StateChangeAirLayerJob = new JobHelper(_cellCount);
		AdvectionAirJob = new JobHelper(_cellCount);
		AdvectionWaterJob = new JobHelper(_cellCount);
		AdvectionCloudJob = new JobHelper(_cellCount);
		ApplyAdvectionAirJob = new JobHelper(_cellCount);
		ApplyAdvectionWaterJob = new JobHelper(_cellCount);
		ApplyAdvectionCloudJob = new JobHelper(_cellCount);
		UpdateDependentWaterLayerJob = new JobHelper(_cellCount);
		UpdateDependentAirLayerJob = new JobHelper(_cellCount);
		UpdateDependentJob = new JobHelper(_cellCount);
		GetVectorDestCoordsJob = new JobHelper(_cellCount);

#if SolarRadiationInJobDebug
		SolarRadiationInJob.Async = false;
#endif

#if UpdateTerrainJobDebug
		UpdateTerrainJob.Async = false;
#endif

#if EmissivityAirJobDebug
		EmissivityAirJob.Async = false;
#endif

#if EmissivityWaterJobDebug
		EmissivityWaterJob.Async = false;
#endif

#if EmissivityTerrainJobDebug
		EmissivityTerrainJob.Async = false;
#endif

#if SolarRadiationAbsorbedAirJobDebug
		SolarRadiationAbsorbedAirJob.Async = false;
#endif

#if SolarRadiationAbsorbedIceJobDebug
		SolarRadiationAbsorbedIceJob.Async = false;
#endif

#if SolarRadiationAbsorbedWaterJobDebug
		SolarRadiationAbsorbedWaterJob.Async = false;
#endif

#if SolarRadiationAbsorbedTerrainJobDebug
		SolarRadiationAbsorbedTerrainJob.Async = false;
#endif

#if ThermalOutCloudJobDebug
		ThermalOutCloudJob.Async = false;
#endif

#if ThermalOutAirJobDebug
		ThermalOutAirJob.Async = false;
#endif

#if ThermalOutIceJobDebug
		ThermalOutIceJob.Async = false;
#endif

#if ThermalOutWaterJobDebug
		ThermalOutWaterJob.Async = false;
#endif

#if ThermalOutTerrainJobDebug
		ThermalOutTerrainJob.Async = false;
#endif

#if ThermalInUpAirJobDebug
		ThermalInUpAirJob.Async = false;
#endif

#if ThermalInUpIceJobDebug
		ThermalInUpIceJob.Async = false;
#endif

#if ThermalInUpWaterJobDebug
		ThermalInUpWaterJob.Async = false;
#endif

#if ThermalInDownAirJobDebug
		ThermalInDownAirJob.Async = false;
#endif

#if ThermalInDownIceJobDebug
		ThermalInDownIceJob.Async = false;
#endif

#if ThermalInDownWaterJobDebug
		ThermalInDownWaterJob.Async = false;
#endif

#if ThermalInDownTerrainJobDebug
		ThermalInDownTerrainJob.Async = false;
#endif

#if PressureGradientForceAirJobDebug
		PressureGradientForceAirJob.Async = false;
#endif

#if PressureGradientForceCloudJobDebug
		PressureGradientForceCloudJob.Async = false;
#endif

#if WindFrictionJobDebug
		WindFrictionJob.Async = false;
#endif

#if WaterFrictionJobDebug
		WaterFrictionJob.Async = false;
#endif

#if WaterDensityGradientForceJobDebug
		WaterDensityGradientForceJob.Async = false;
#endif

#if DiffusionAirJobDebug
		DiffusionAirJob.Async = false;
#endif

#if DiffusionWaterJobDebug
		DiffusionWaterJob.Async = false;
#endif

#if DiffusionCloudJobDebug
		DiffusionCloudJob.Async = false;
#endif

#if AdvectionAirJobDebug
		AdvectionAirJob.Async = false;
#endif

#if AdvectionWaterJobDebug
		AdvectionWaterJob.Async = false;
#endif

#if AdvectionCloudJobDebug
		AdvectionCloudJob.Async = false;
#endif

#if GetVectorDestCoordsJobDebug
		GetVectorDestCoordsJob.Async = false;
#endif


#if ConductionAirIceJobDebug
		ConductionAirIceJob.Async = false;
#endif

#if ConductionAirWaterJobDebug
		ConductionAirWaterJob.Async = false;
#endif

#if ConductionAirTerrainJobDebug
		ConductionAirTerrainJob.Async = false;
#endif

#if ConductionIceWaterJobDebug
		ConductionIceWaterJob.Async = false;
#endif

#if ConductionIceTerrainJobDebug
		ConductionIceTerrainJob.Async = false;
#endif

#if ConductionWaterTerrainJobDebug
		ConductionWaterTerrainJob.Async = false;
#endif

#if EnergyTerrainJobDebug
		EnergyTerrainJob.Async = false;
#endif

#if EnergyCloudJobDebug
		EnergyCloudJob.Async = false;
#endif

#if EnergyIceJobDebug
		EnergyIceJob.Async = false;
#endif

#if EnergyWaterJobDebug
		EnergyWaterJob.Async = false;
#endif

#if EnergyWaterSurfaceJobDebug
		EnergyWaterSurfaceJob.Async = false;
#endif

#if EnergyAirJobDebug
		EnergyAirJob.Async = false;
#endif

#if StateChangeJobDebug
		StateChangeJob.Async = false;
#endif

#if StateChangeAirLayerJobDebug
		StateChangeAirLayerJob.Async = false;
#endif

#if UpdateDependentWaterLayerJobDebug
		UpdateDependentWaterLayerJob.Async = false;
#endif

#if UpdateDependentAirLayerJobDebug
		UpdateDependentAirLayerJob.Async = false;
#endif

#if UpdateDependentJobDebug
		UpdateDependentJob.Async = false;
#endif
#if FluxAirJobDebug
		FluxAirJob.Async = false;
#endif
#if FluxCloudJobDebug
		FluxCloudJob.Async = false;
#endif
#if FluxIceJobDebug
		FluxIceJob.Async = false;
#endif
#if FluxWaterJobDebug
		FluxWaterJob.Async = false;
#endif



		#endregion

		solarRadiation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		waterSlopeAlbedo = new NativeArray<float>(_cellCount, Allocator.Persistent);
		emissivity = new NativeArray<float>[_layerCount];
		solarRadiationIn = new NativeArray<float>[_layerCount];
		latentHeat = new NativeArray<float>[_layerCount];
		for (int i=0;i<_layerCount;i++)
		{
			latentHeat[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			emissivity[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			solarRadiationIn[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
		}
		diffusionAir = new NativeArray<DiffusionAir>[_airLayers];
		advectionAir = new NativeArray<DiffusionAir>[_airLayers];
		destinationAir = new NativeArray<BarycentricValue>[_airLayers];
		pressureGradientForce = new NativeArray<float3>[_airLayers];
		fluxEnergyAir = new NativeArray<float>[_airLayers];
		for (int i = 0; i < _airLayers; i++)
		{
			fluxEnergyAir[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			diffusionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
			advectionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
			destinationAir[i] = new NativeArray<BarycentricValue>(_cellCount, Allocator.Persistent);
			pressureGradientForce[i] = new NativeArray<float3>(_cellCount, Allocator.Persistent);
		}
		diffusionWater = new NativeArray<DiffusionWater>[_waterLayers];
		advectionWater = new NativeArray<DiffusionWater>[_waterLayers];
		destinationWater = new NativeArray<BarycentricValue>[_waterLayers];
		conductionWaterTerrain = new NativeArray<float>[_waterLayers];
		fluxEnergyWater = new NativeArray<float>[_waterLayers];
		evaporationMass = new NativeArray<float>[_waterLayers];
		frozenMass = new NativeArray<float>[_waterLayers];
		frozenTemperature = new NativeArray<float>[_waterLayers];
		waterForce = new NativeArray<float3>[_waterLayers];
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
			advectionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
			destinationWater[i] = new NativeArray<BarycentricValue>(_cellCount, Allocator.Persistent);
			conductionWaterTerrain[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			fluxEnergyWater[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			evaporationMass[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			frozenMass[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			frozenTemperature[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			waterForce[i] = new NativeArray<float3>(_cellCount, Allocator.Persistent);
		}
		fluxEnergyIce = new NativeArray<float>(_cellCount, Allocator.Persistent);
		pressureGradientForceCloud = new NativeArray<float3>(_cellCount, Allocator.Persistent);
		windFriction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		diffusionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		advectionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		destinationCloud = new NativeArray<BarycentricValue>(_cellCount, Allocator.Persistent);
		conductionSurfaceAirIce = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionSurfaceAirWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionSurfaceAirTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionUpperAirIce = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionUpperAirWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionUpperAirTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionIceWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionIceTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		displaySolarRadiation = new NativeArray<float>(_cellCount, Allocator.Persistent);

	}

	public void Dispose()
	{
		solarRadiation.Dispose();
		waterSlopeAlbedo.Dispose();
		for (int i=0;i<_layerCount;i++)
		{
			latentHeat[i].Dispose();
			emissivity[i].Dispose();
			solarRadiationIn[i].Dispose();
		}
		for (int i = 0; i < _airLayers; i++)
		{
			fluxEnergyAir[i].Dispose();
			diffusionAir[i].Dispose();
			advectionAir[i].Dispose();
			destinationAir[i].Dispose();
			pressureGradientForce[i].Dispose();
		}
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i].Dispose();
			advectionWater[i].Dispose();
			destinationWater[i].Dispose();
			conductionWaterTerrain[i].Dispose();
			fluxEnergyWater[i].Dispose();
			evaporationMass[i].Dispose();
			frozenMass[i].Dispose();
			frozenTemperature[i].Dispose();
			waterForce[i].Dispose();
		}
		pressureGradientForceCloud.Dispose();
		fluxEnergyIce.Dispose();
		windFriction.Dispose();
		diffusionCloud.Dispose();
		advectionCloud.Dispose();
		destinationCloud.Dispose();
		conductionSurfaceAirIce.Dispose();
		conductionSurfaceAirWater.Dispose();
		conductionSurfaceAirTerrain.Dispose();
		conductionUpperAirIce.Dispose();
		conductionUpperAirWater.Dispose();
		conductionUpperAirTerrain.Dispose();
		conductionIceWater.Dispose();
		conductionIceTerrain.Dispose();


		displaySolarRadiation.Dispose();

	}

	public void Tick(SimState[] states, int stateCount, int ticksToAdvance, ref DependentState dependent, ref DisplayState display, ref StaticState staticState, ref WorldData worldData, ref int curStateIndex, bool checkForDegeneracy, bool logState, int logStateIndex)
	{
		JobHandle lastJobHandle = default(JobHandle);
		for (int tick = 0; tick < ticksToAdvance; tick++)
		{
			#region Init Time step

			bool updateDisplay = tick == ticksToAdvance - 1;

			ref var lastState = ref states[curStateIndex];
			curStateIndex = (curStateIndex + 1) % stateCount;
			ref var nextState = ref states[curStateIndex];

			var jobHandleDependencies = new List<NativeList<JobHandle>>();

			var thermalRadiationDeltaIceBottom = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var thermalRadiationDeltaIceTop = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var thermalRadiationDeltaSurfaceWater = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var thermalRadiationDelta = new NativeArray<float>[_layerCount];
			var thermalRadiationTransmittedUp = new NativeArray<float>[_layerCount];
			var thermalRadiationTransmittedDown = new NativeArray<float>[_layerCount];
			var windowRadiationTransmittedUp = new NativeArray<float>[_layerCount];
			var windowRadiationTransmittedDown = new NativeArray<float>[_layerCount];
			var condensationGroundMass = new NativeArray<float>[_layerCount];
			var condensationCloudMass = new NativeArray<float>[_layerCount];
			var solarReflected = new NativeArray<float>[_layerCount];
			var conductionWaterTerrainTotal = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var cloudEvaporationMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var dropletDelta = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var precipitationMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var precipitationTemperature = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var iceMeltedMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var airMassTotal = new NativeArray<float>(staticState.StratosphereMass, Allocator.TempJob);
			var waterDepthTotal = new NativeArray<float>(_cellCount, Allocator.TempJob);
			for (int i = 0; i < _layerCount; i++)
			{
				solarReflected[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				thermalRadiationDelta[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				thermalRadiationTransmittedUp[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				thermalRadiationTransmittedDown[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				windowRadiationTransmittedUp[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				windowRadiationTransmittedDown[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				condensationGroundMass[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				condensationCloudMass[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
			}
			#endregion

			#region Update Planetary Globals

			nextState.PlanetState = lastState.PlanetState;
			nextState.PlanetState.Ticks++;

			// TODO: update
			float distanceToSun = math.length(lastState.PlanetState.Position);
			float angleToSun = math.atan2(lastState.PlanetState.Position.z, lastState.PlanetState.Position.x);
			angleToSun += lastState.PlanetState.OrbitSpeed;
			nextState.PlanetState.Position = new float3(math.cos(angleToSun), 0, math.sin(angleToSun)) * distanceToSun;
			nextState.PlanetState.Rotation = new float3(lastState.PlanetState.Rotation.x, Mathf.Repeat(lastState.PlanetState.Rotation.y + lastState.PlanetState.SpinSpeed * worldData.SecondsPerTick, math.PI * 2), 0);

			float coriolisTerm = 2 * lastState.PlanetState.SpinSpeed;
			#endregion

			#region Init Solar Radiation Per Cell
			var solarInJobHandle = SolarRadiationInJob.Run(new SolarRadiationJob()
			{
				SolarRadiation = solarRadiation,
				DisplaySolarRadiation = displaySolarRadiation,
				WaterSlopeAlbedo = waterSlopeAlbedo,

				SphericalPosition = staticState.SphericalPosition,
				IncomingSolarRadiation = lastState.PlanetState.SolarRadiation * worldData.SecondsPerTick,
				PlanetRotation = quaternion.Euler(lastState.PlanetState.Rotation),
				SunToPlanetDir = math.normalize(lastState.PlanetState.Position),
			}, lastJobHandle);

			var updateTerrainJobHandle = UpdateTerrainJob.Run(new UpdateTerrainJob()
			{
				Terrain = nextState.Terrain,

				LastTerrain = lastState.Terrain,
			}, lastJobHandle);
			#endregion

			// Calculate emissivity per cell
			// TODO: combine cloud emissivity with cloud conduction
			// TODO: we can probably combine this step with the thermal radiation step
			#region Emissivity Per Cell

			JobHandle[] emissivityJobHandles = new JobHandle[_layerCount];
			for (int j = 1; j < _airLayers-1; j++)
			{
				int layerIndex = _airLayer0 + j;
				emissivityJobHandles[layerIndex] = EmissivityAirJob.Run(new EmissivityAirJob()
				{
					Emissivity = emissivity[layerIndex],
					AirMass = dependent.AirMass[j],
					VaporMass = lastState.AirVapor[j],
				});
			}

			for (int j = 1; j < _waterLayers-1; j++)
			{
				int layerIndex = _waterLayer0 + j;
				emissivityJobHandles[j + _waterLayer0] = EmissivityWaterJob.Run(new EmissivityWaterJob()
				{
					Emissivity = emissivity[layerIndex],
					WaterMass = lastState.WaterMass[j],
					SaltMass = lastState.SaltMass[j]
				});
			}
			emissivityJobHandles[_terrainLayer] = EmissivityTerrainJob.Run(new EmissivityTerrainJob()
			{
				Emissivity = emissivity[_terrainLayer],
				Terrain = lastState.Terrain,
				VegetationCoverage = dependent.VegetationCoverage
			});
			#endregion

			// Follow the solar radiation down from the top of the atmosphere to ther terrain, and absorb some as it passes through each layer
			#region Solar Radiation Absorbed
			// process each vertical layer in order

			// atmosphere
			JobHandle[] solarInJobHandles = new JobHandle[_layerCount];
			for (int j = _airLayer0 + _airLayers - 2; j > _airLayer0; j--)
			{
				int airLayerIndex = j - _airLayer0;
				solarInJobHandles[j] = solarInJobHandle = SolarRadiationAbsorbedAirJob.Run(new SolarRadiationAbsorbedAirJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[j],
					SolarRadiationIncoming = solarRadiation,
					SolarRadiationReflected = solarReflected[j],
					AirMass = dependent.AirMass[airLayerIndex],
					VaporMass = lastState.AirVapor[airLayerIndex],
					CloudMass = lastState.CloudMass,
					CloudDropletMass = lastState.CloudDropletMass,
					CloudElevation = dependent.CloudElevation,
					CloudCoverage = dependent.CloudCoverage,
					DewPoint = dependent.DewPoint,
					WaterSlopeAlbedo = waterSlopeAlbedo,
					SolarReflectivityAir = worldData.SolarReflectivityAir,
					SolarAbsorptivityAir = worldData.SolarAbsorptivityAir,
					SolarAbsorptivityWaterVapor = worldData.SolarAbsorptivityWaterVapor,
					SolarAbsorptivityCloud = worldData.SolarAbsorptivityCloud,
					LayerElevation = dependent.LayerElevation[airLayerIndex],
					LayerHeight = dependent.LayerHeight[airLayerIndex],
					LayerIndex = airLayerIndex,
					worldData = worldData
				}, solarInJobHandle);
			}

			
			// ice
			solarInJobHandles[_iceLayer] = solarInJobHandle = SolarRadiationAbsorbedIceJob.Run(new SolarRadiationAbsorbedIceJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_iceLayer],
				SolarRadiationIncoming = solarRadiation,
				SolarRadiationReflected = solarReflected[_iceLayer],
				AlbedoIce = WorldData.AlbedoIce,
				IceCoverage = dependent.IceCoverage
			}, solarInJobHandle);
			
			for (int j = _surfaceWaterLayer; j >= 1; j--)
			{
				int layerIndex = _waterLayer0 + j;
				solarInJobHandles[layerIndex] = solarInJobHandle = SolarRadiationAbsorbedWaterJob.Run(new SolarRadiationAbsorbedWaterJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[layerIndex],
					SolarRadiationIncoming = solarRadiation,
					SolarRadiationReflected = solarReflected[layerIndex],
					WaterCoverage = dependent.WaterCoverage[layerIndex - _waterLayer0],
					WaterSlopeAlbedo = waterSlopeAlbedo,
				}, solarInJobHandle);
			}

			solarInJobHandles[_terrainLayer] = solarInJobHandle = SolarRadiationAbsorbedTerrainJob.Run(new SolarRadiationAbsorbedTerrainJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_terrainLayer],
				SolarRadiationIncoming = solarRadiation,
				SolarRadiationReflected = solarReflected[_terrainLayer],
				VegetationCoverage = dependent.VegetationCoverage,
				worldData = worldData,
				LastTerrain = lastState.Terrain,
			}, solarInJobHandle);
			#endregion

			// Calculate how much thermal radition is being emitted out of each layer
			#region Thermal Radiation
			JobHandle[] thermalOutJobHandles = new JobHandle[_layerCount];

			// ICE
			thermalOutJobHandles[_iceLayer] = ThermalOutIceJob.Run(new ThermalEnergyRadiatedConstantEmissivityJob()
			{
				ThermalRadiationDelta = thermalRadiationDelta[_iceLayer],
				ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[_iceLayer],
				ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[_iceLayer],
				WindowRadiationTransmittedUp = windowRadiationTransmittedUp[_iceLayer],
				WindowRadiationTransmittedDown = windowRadiationTransmittedDown[_iceLayer],

				PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
				Emissivity = WorldData.EmissivityIce,
				Energy = dependent.IceEnergy,
				Temperature = lastState.IceTemperature,
				SurfaceArea = dependent.IceCoverage,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);



			// TERRAIN
			thermalOutJobHandles[_terrainLayer] = ThermalOutTerrainJob.Run(new ThermalEnergyRadiatedTerrainJob()
			{
				ThermalRadiationDelta = thermalRadiationDelta[_terrainLayer],
				ThermalRadiationTransmitted = thermalRadiationTransmittedUp[_terrainLayer],
				WindowRadiationTransmitted = windowRadiationTransmittedUp[_terrainLayer],

				PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
				Emissivity = emissivity[_terrainLayer],
				Temperature = lastState.TerrainTemperature,
				SecondsPerTick = worldData.SecondsPerTick
			}, emissivityJobHandles[_terrainLayer]);


			// ATMOSPHERE
			for (int j = 1; j < _airLayers-1; j++)
			{
				int layer = _airLayer0 + j;
				thermalOutJobHandles[layer] = ThermalOutAirJob.Run(new ThermalEnergyRadiatedJob()
				{
					ThermalRadiationDelta = thermalRadiationDelta[layer],
					ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[layer],
					ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[layer],
					WindowRadiationTransmittedUp = windowRadiationTransmittedUp[layer],
					WindowRadiationTransmittedDown = windowRadiationTransmittedDown[layer],

					PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
					Energy = dependent.AirPotentialEnergy[j],
					Emissivity = emissivity[layer],
					Temperature = lastState.AirTemperature[j],
					SecondsPerTick = worldData.SecondsPerTick
				}, emissivityJobHandles[layer]);
			}

			// WATER
			for (int j = 1; j < _waterLayers - 1; j++)
			{
				int layer = _waterLayer0 + j;
				thermalOutJobHandles[layer] = ThermalOutWaterJob.Run(new ThermalEnergyRadiatedJob()
				{
					ThermalRadiationDelta = thermalRadiationDelta[layer],
					ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[layer],
					ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[layer],
					WindowRadiationTransmittedUp = windowRadiationTransmittedUp[layer],
					WindowRadiationTransmittedDown = windowRadiationTransmittedDown[layer],

					PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
					Emissivity = emissivity[layer],
					Energy = dependent.WaterPotentialEnergy[j],
					Temperature = lastState.WaterTemperature[j],
					SecondsPerTick = worldData.SecondsPerTick
				}, emissivityJobHandles[layer]);
			}
			#endregion


			// Thermal radiation travels upwards, partially reflecting downwards (clouds), partially absorbed, and partially lost to space
			#region Thermal Radiation Absorbed Up
			// Start at bottom water layer and go up, then go back down
			JobHandle[] thermalInUpJobHandles = new JobHandle[_layerCount];
			JobHandle[] thermalInDownJobHandles = new JobHandle[_layerCount];
			NativeArray<float> atmosphericWindowUp = new NativeArray<float>(_cellCount, Allocator.TempJob);
			NativeArray<float> atmosphericWindowDown = new NativeArray<float>(_cellCount, Allocator.TempJob);

			// transmit up from land
			for (int j = 1; j < _layerCount; j++)
			{

				if (j > _airLayer0 && j < _airLayer0 + _airLayers - 1)
				{
					int airLayer = j - _airLayer0;
					int downIndex = airLayer == 1 ? _iceLayer : (j - 1);
					thermalInUpJobHandles[j] = ThermalInUpAirJob.Run(new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[downIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[downIndex],
						CloudSurfaceArea = dependent.CloudCoverage,
						CloudElevation = dependent.CloudElevation,
						CloudMass = lastState.CloudMass,
						LayerElevation = dependent.LayerElevation[airLayer],
						LayerHeight = dependent.LayerHeight[airLayer],
						AirMass = dependent.AirMass[airLayer],
						VaporMass = lastState.AirVapor[airLayer],
						AirAbsorptivity = worldData.ThermalAbsorptivityAir,
						VaporAbsorptivity = worldData.ThermalAbsorptivityWaterVapor,
						WaterAbsorptivity = worldData.ThermalAbsorptivityWaterLiquid,
						CarbonAbsorptivity = worldData.ThermalAbsorptivityCarbonDioxide,
						LayerIndex = airLayer,
						FromTop = false,
					}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex]));
				}
				else if (j == _iceLayer)
				{
					int downIndex = _waterLayer0 + _surfaceWaterLayer;
					thermalInUpJobHandles[j] = ThermalInUpIceJob.Run(new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDeltaIceBottom,
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[downIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[downIndex],
						Coverage = dependent.IceCoverage,

					}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex]));
				}
				else if (j > _waterLayer0 && j < _waterLayer0 + _waterLayers - 1)
				{
					int waterLayerIndex = j - _waterLayer0;
					int downIndex = (waterLayerIndex == 1) ? _terrainLayer : (j - 1);
					thermalInUpJobHandles[j] = ThermalInUpWaterJob.Run(new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDeltaSurfaceWater,
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[downIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[downIndex],
						Coverage = dependent.WaterCoverage[waterLayerIndex],
					}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex], thermalOutJobHandles[downIndex]));
				}
			}

			var thermalInUpJobHandlesCombined = default(JobHandle);
			for (int j=0;j<_layerCount;j++)
			{
				thermalInUpJobHandlesCombined = JobHandle.CombineDependencies(thermalInUpJobHandlesCombined, thermalInUpJobHandles[j]);
			}
			#endregion

			// Thermal radiation is absorbed travelling downwards, collecting and then eventually hitting the earth (back radiation)
			// TODO: we need to include the top layer of atmosphere here, since we calculate cloud absorption as part of the air layer step
			#region Thermal Radiation Absorbed Down

			// transmit down from top of atmosphere			
			for (int j = _layerCount - 2; j >= 0; j--)
			{

				if (j == _terrainLayer)
				{
					// TERRAIN
					int upIndex = _waterLayer0 + 1;
					var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
					thermalInDownJobHandles[j] = ThermalInDownTerrainJob.Run(new ThermalEnergyAbsorbedTerrainJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationDelta[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[upIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[upIndex],
					}, thermalInDependenciesHandle);
				}
				else if (j == _iceLayer)
				{
					// ICE
					int upIndex = _airLayer0 + 1;
					var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
					thermalInDownJobHandles[j] = ThermalInDownIceJob.Run(new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDeltaIceTop,
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[upIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[upIndex],
						Coverage = dependent.IceCoverage,
					}, thermalInDependenciesHandle);
				}
				else if (j > _waterLayer0 && j <= _surfaceWaterLayer)
				{
					// WATER
					int waterLayerIndex = j - _waterLayer0;
					int upIndex = (waterLayerIndex == _surfaceWaterLayer) ? _iceLayer : (j + 1);
					var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
					thermalInDownJobHandles[j] = ThermalInDownWaterJob.Run(new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[upIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[upIndex],
						Coverage = dependent.WaterCoverage[waterLayerIndex],
					}, thermalInDependenciesHandle);
				}
				else if (j > _airLayer0 && j < _airLayer0 + _airLayers - 1)
				{
					int airLayer = j - _airLayer0;
					int upIndex = j + 1;
					var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
					thermalInDownJobHandles[j] = ThermalInDownAirJob.Run(new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[upIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[upIndex],
						CloudSurfaceArea = dependent.CloudCoverage,
						CloudElevation = dependent.CloudElevation,
						LayerElevation = dependent.LayerElevation[airLayer],
						LayerHeight = dependent.LayerHeight[airLayer],
						AirMass = dependent.AirMass[airLayer],
						CloudMass = lastState.CloudMass,
						VaporMass = lastState.AirVapor[airLayer],
						AirAbsorptivity = worldData.ThermalAbsorptivityAir,
						VaporAbsorptivity = worldData.ThermalAbsorptivityWaterVapor,
						WaterAbsorptivity = worldData.ThermalAbsorptivityWaterLiquid,
						CarbonAbsorptivity = worldData.ThermalAbsorptivityCarbonDioxide,
						LayerIndex = airLayer,
						FromTop = true,
					}, thermalInDependenciesHandle);
				}
			}
			#endregion

			// Buoyancy, Updrafts, and mixing occur across air layers and water layers
			// TODO: add an empty air layer on top and bottom so we can calculate up/down diffusion in a single step 
			// Temperature and trace elements diffuse into neighboring horizontal cells based on a diffusion constant
			// Air, Water, Cloud
			#region Pressure Gradient Force

			JobHandle[] pgfJobHandles = new JobHandle[_airLayers];
			for (int j = 1; j < _airLayers - 1; j++)
			{
				pgfJobHandles[j] = PressureGradientForceAirJob.Run(new PressureGradientForceAirJob()
				{
					Delta = pressureGradientForce[j],

					Pressure = dependent.AirPressure[j],
					AirMass = dependent.AirMass[j],
					Temperature = lastState.AirTemperature[j],
					VaporMass = lastState.AirVapor[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					Neighbors = staticState.Neighbors,
					Positions = staticState.SphericalPosition,
					PlanetRadius = staticState.PlanetRadius,
					Gravity = lastState.PlanetState.Gravity,
					UpTemperature = lastState.AirTemperature[j + 1],
					UpHumidity = lastState.AirVapor[j + 1],
					UpAirMass = dependent.AirMass[j + 1],
					UpLayerElevation = dependent.LayerElevation[j + 1],
					UpLayerHeight = dependent.LayerHeight[j + 1],
					DownTemperature = lastState.AirTemperature[j - 1],
					DownHumidity = lastState.AirVapor[j - 1],
					DownAirMass = dependent.AirMass[j - 1],
					DownLayerElevation = dependent.LayerElevation[j - 1],
					DownLayerHeight = dependent.LayerHeight[j - 1],
					IsTop = j == _airLayers - 2,
					IsBottom = j == 1

				}, lastJobHandle);
			}
			var windFrictionJobHandle = WindFrictionJob.Run(new WindFrictionJob()
			{
				Friction = windFriction,
				IceCoverage = dependent.IceCoverage,
				WaterCoverage = dependent.WaterCoverage[_surfaceWaterLayer],
				VegetationCoverage = dependent.VegetationCoverage,
				Terrain = lastState.Terrain,
				IceFriction = worldData.WindIceFriction,
				TerrainFrictionMin = worldData.WindTerrainFrictionMin,
				TerrainFrictionMax = worldData.WindTerrainFrictionMax,
				VegetationFriction = worldData.WindVegetationFriction,
				WaterFriction = worldData.WindWaterFriction,
				MaxTerrainRoughness = worldData.MaxTerrainRoughnessForWindFriction
			}, lastJobHandle);

			JobHandle pgfCloudJobHandle = default(JobHandle);
			for (int j = 1; j < _airLayers - 1; j++)
			{

				pgfCloudJobHandle = JobHandle.CombineDependencies(pgfCloudJobHandle, PressureGradientForceCloudJob.Run(new PressureGradientForceCloudJob()
				{
					Force = pressureGradientForceCloud,

					CloudElevation = dependent.CloudElevation,
					LayerForce = pressureGradientForce[j],
					LayerFriction = windFriction,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					DownForce = pressureGradientForce[j-1],
					DownLayerElevation = dependent.LayerElevation[j-1],
					UpForce = pressureGradientForce[j + 1],
					UpLayerElevation = dependent.LayerElevation[j + 1],
					UpLayerHeight = dependent.LayerHeight[j + 1],
					Velocity = lastState.CloudVelocity,
					LayerFrictionMultiplier = (j == 1) ? 1 : 0,
					UpFrictionMultiplier = ((j + 1) == 1) ? 1 : 0,
					DownFrictionMultiplier = ((j - 1) == 1) ? 1 : 0,
					IsTop = j == _airLayers - 2,
					IsBottom = j == 1,
					SecondsPerTick = worldData.SecondsPerTick
				}, JobHandle.CombineDependencies(pgfCloudJobHandle, windFrictionJobHandle, JobHandle.CombineDependencies(pgfJobHandles[j], pgfJobHandles[j - 1], pgfJobHandles[j + 1]))));
			}

			JobHandle[] waterDensityForceJobHandles = new JobHandle[_airLayers];
			for (int j = 1; j < _waterLayers - 1; j++)
			{
				waterDensityForceJobHandles[j] = WaterDensityGradientForceJob.Run(new WaterDensityGradientForceJob()
				{
					Force = waterForce[j],

					Positions = staticState.SphericalPosition,
					Neighbors = staticState.Neighbors,
					WaterDensity = dependent.WaterDensity[j],
					WaterPressure = dependent.WaterPressure[j],
					LayerDepth = dependent.WaterLayerDepth[j],
					LayerHeight = dependent.WaterLayerHeight[j],
					UpWaterDensity = dependent.WaterDensity[j+1],
					UpWaterPressure = dependent.WaterPressure[j+1],
					UpLayerDepth = dependent.WaterLayerDepth[j + 1],
					UpLayerHeight = dependent.WaterLayerHeight[j + 1],
					DownWaterDensity = dependent.WaterDensity[j - 1],
					DownWaterPressure = dependent.WaterPressure[j - 1],
					DownLayerDepth = dependent.WaterLayerDepth[j - 1],
					DownLayerHeight = dependent.WaterLayerHeight[j - 1],
					Gravity = lastState.PlanetState.Gravity,
					InverseCellDiameter = staticState.InverseCellDiameter,

				}, lastJobHandle);

			}

			var waterFrictionJobHandle = WaterFrictionJob.Run(new WaterFrictionForceJob()
			{
				Force = waterForce[_surfaceWaterLayer],

				Position = staticState.SphericalPosition,
				Current = lastState.WaterVelocity[_surfaceWaterLayer],
				WindUp = lastState.Wind[1],
				WindDown = lastState.WaterVelocity[_surfaceWaterLayer - 1],
				LayerHeight = dependent.WaterLayerHeight[_surfaceWaterLayer],
				CoriolisMultiplier = staticState.CoriolisMultiplier,
				FrictionCoefficientUp = worldData.WindToWaterCurrentFrictionCoefficient,
				FrictionCoefficientDown = 0, // TODO: do we want to add a frictional force between layers of water?
				CoriolisTerm = coriolisTerm,
				WaterSurfaceFrictionDepth = worldData.WaterSurfaceFrictionDepth,
			}, waterDensityForceJobHandles[_surfaceWaterLayer]);


			#endregion


			// Conduction is calculated for each Surface that might touch another surface
			// Air to Cloud, Air to Ice, Air to Water, Air to Terrain, Ice to Water, Ice to Terrain, Water to Terrain
			#region Conduction
			// air to ice
			var conductionAirIceJobHandle = ConductionAirIceJob.Run(new ConductionAirIceJob()
			{
				EnergyDelta = conductionSurfaceAirIce,
				TemperatureA = dependent.SurfaceAirTemperature,
				TemperatureB = lastState.IceTemperature,
				EnergyB = dependent.IceEnergy,
				ConductionCoefficient = WorldData.ConductivityAirIce,
				Coverage = dependent.IceCoverage,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// air to water
			var conductionAirWaterJobHandle = ConductionAirWaterJob.Run(new ConductionAirWaterJob()
			{
				EnergyDelta = conductionSurfaceAirWater,
				TemperatureA = dependent.SurfaceAirTemperature,
				TemperatureB = lastState.WaterTemperature[_surfaceWaterLayer],
				EnergyA = dependent.IceEnergy,
				EnergyB = dependent.WaterPotentialEnergy[_surfaceWaterLayer],
				ConductionCoefficient = WorldData.ConductivityAirWater,
				CoverageIce = dependent.IceCoverage,
				CoverageWater = dependent.WaterCoverage[_surfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// air to terrain
			var conductionAirTerrainJobHandle = ConductionAirTerrainJob.Run(new ConductionAirTerrainJob()
			{
				EnergyDelta = conductionSurfaceAirTerrain,
				TemperatureA = dependent.SurfaceAirTemperature,
				TemperatureB = lastState.TerrainTemperature,
				ConductionCoefficient = WorldData.ConductivityAirTerrain,
				CoverageIce = dependent.IceCoverage,
				CoverageWater = dependent.WaterCoverage[_surfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// ice to water
			var conductionIceWaterJobHandle = ConductionIceWaterJob.Run(new ConductionIceWaterJob()
			{
				EnergyDelta = conductionIceWater,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.WaterTemperature[_surfaceWaterLayer],
				EnergyA = dependent.IceEnergy,
				EnergyB = dependent.WaterPotentialEnergy[_surfaceWaterLayer],
				ConductionCoefficient = WorldData.ConductivityIceWater,
				CoverageA = dependent.IceCoverage,
				CoverageB = dependent.WaterCoverage[_surfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// ice to terrain
			var conductionIceTerrainJobHandle = ConductionIceTerrainJob.Run(new ConductionIceTerrainJob()
			{
				EnergyDelta = conductionIceTerrain,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.TerrainTemperature,
				EnergyA = dependent.IceEnergy,
				ConductionCoefficient = WorldData.ConductivityIceTerrain,
				CoverageIce = dependent.IceCoverage,
				CoverageWater = dependent.WaterCoverage[_surfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// water to terrain
			JobHandle conductionWaterTerrainJobHandle = default(JobHandle);
			{
				conductionWaterTerrainJobHandle = JobHandle.CombineDependencies(conductionWaterTerrainJobHandle, ConductionWaterTerrainJob.Run(new ConductionWaterBottomJob()
				{
					EnergyDelta = conductionWaterTerrain[1],
					EnergyDeltaWaterTotal = conductionWaterTerrainTotal,
					TemperatureA = lastState.WaterTemperature[1],
					TemperatureB = lastState.TerrainTemperature,
					EnergyA = dependent.WaterPotentialEnergy[1],
					ConductionCoefficient = WorldData.ConductivityWaterTerrain,
					Coverage = dependent.WaterCoverage[1],
					SecondsPerTick = worldData.SecondsPerTick
				}, lastJobHandle));
			}
			for (int i = 2; i < _waterLayers-1; i++) {
				conductionWaterTerrainJobHandle = JobHandle.CombineDependencies(conductionWaterTerrainJobHandle, ConductionWaterTerrainJob.Run(new ConductionWaterTerrainJob()
				{
					EnergyDelta = conductionWaterTerrain[i],
					EnergyDeltaWaterTotal = conductionWaterTerrainTotal,
					TemperatureA = lastState.WaterTemperature[i],
					TemperatureB = lastState.TerrainTemperature,
					EnergyA = dependent.WaterPotentialEnergy[i],
					ConductionCoefficient = WorldData.ConductivityWaterTerrain,
					Coverage = dependent.WaterCoverage[i],
					CoverageBelow = dependent.WaterCoverage[i - 1],
					SecondsPerTick = worldData.SecondsPerTick
				}, conductionWaterTerrainJobHandle));
			}

			#endregion

			#region State Change (Flux)

			var energyJobHandles = new NativeList<JobHandle>(Allocator.TempJob);

			// surface water
			JobHandle[] fluxJobs = new JobHandle[_layerCount];
			int surfaceWaterLayerIndex = _waterLayer0 + _surfaceWaterLayer;

			var fluxWaterDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[surfaceWaterLayerIndex],
				thermalOutJobHandles[surfaceWaterLayerIndex],
				thermalInDownJobHandles[surfaceWaterLayerIndex],
				thermalInUpJobHandles[surfaceWaterLayerIndex],
				conductionAirWaterJobHandle,
				conductionIceWaterJobHandle,
				conductionWaterTerrainJobHandle,
				waterDensityForceJobHandles[_surfaceWaterLayer],
				waterFrictionJobHandle,
			};
			fluxJobs[_waterLayer0+_surfaceWaterLayer] = FluxWaterJob.Run(new FluxWaterJob()
			{
				EvaporatedWaterMass = evaporationMass[_surfaceWaterLayer],
				FrozenMass = frozenMass[_surfaceWaterLayer],
				FrozenTemperature = frozenTemperature[_surfaceWaterLayer],
				Energy = fluxEnergyWater[_surfaceWaterLayer],

				LastMass = lastState.WaterMass[_surfaceWaterLayer],
				LastSaltMass = lastState.SaltMass[_surfaceWaterLayer],
				LastTemperature = lastState.WaterTemperature[_surfaceWaterLayer],
				SolarRadiationIn = solarRadiationIn[surfaceWaterLayerIndex],
				ThermalRadiationDeltaTop = thermalRadiationDelta[surfaceWaterLayerIndex],
				ThermalRadiationDeltaBottom = thermalRadiationDeltaSurfaceWater,
				ConductionEnergyAir = conductionSurfaceAirWater,
				ConductionEnergyIce = conductionIceWater,
				ConductionEnergyTerrain = conductionWaterTerrain[_surfaceWaterLayer],
				RelativeHumidity = dependent.AirHumidityRelative[1],
				IceCoverage = dependent.IceCoverage,
				WaterCoverage = dependent.WaterCoverage[_surfaceWaterLayer],
				Salinity = dependent.WaterSalinity[_surfaceWaterLayer],
				EvaporationRate = worldData.EvaporationRate,
				EvapTemperatureMax = worldData.EvapMaxTemperature,
				EvapTemperatureMin = worldData.EvapMinTemperature,
				WaterHeatingDepth = worldData.WaterHeatingDepth,
				FreezePointReductionPerSalinity = worldData.FreezePointReductionPerSalinity,
			}, JobHandle.CombineDependencies(fluxWaterDependencies));


			// CLOUD
			var cloudEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				pgfCloudJobHandle,
				windFrictionJobHandle,
			};
			jobHandleDependencies.Add(cloudEnergyJobHandleDependencies);
			var fluxCloudJobHandle = FluxCloudJob.Run(new FluxCloudJob()
			{
				CloudMass = nextState.CloudMass,
				DropletMass = nextState.CloudDropletMass,
				Velocity = nextState.CloudVelocity,
				EvaporationMass = cloudEvaporationMass,
				PrecipitationMass = precipitationMass,
				PrecipitationTemperature = precipitationTemperature,

				SurfaceAirTemperature = lastState.AirTemperature[1],
				SurfaceSaltMass = lastState.SaltMass[_surfaceWaterLayer],
				LastCloudMass = lastState.CloudMass,
				LastVelocity = lastState.CloudVelocity,
				LastDropletMass = lastState.CloudDropletMass,
				CloudElevation = dependent.CloudElevation,
				DewPoint = dependent.DewPoint,
				AirMassCloud = dependent.AirMassCloud,
				WaterVaporCloud = dependent.AirVaporCloud,
				AirPressureCloud = dependent.AirPressureCloud,
				RelativeHumidityCloud = dependent.AirHumidityRelativeCloud,
				Position = staticState.SphericalPosition,
				Gravity = lastState.PlanetState.Gravity,
				RainDropDragCoefficient = worldData.rainDropDragCoefficient,
				RainDropMaxSize = worldData.rainDropMaxSize,
				RainDropMinSize = worldData.rainDropMinSize,
				SecondsPerTick = worldData.SecondsPerTick,
				CloudDissapationRateDryAir = worldData.CloudDissapationRateDryAir,
				CloudDissapationRateWind = worldData.CloudDissapationRateWind,
				InverseCellDiameter = staticState.InverseCellDiameter,
				SurfaceElevation = dependent.SurfaceElevation,
				TicksPerSecond = worldData.TicksPerSecond,
				WindFriction = windFriction,
				WindFrictionMultiplier = 0,
				CoriolisMultiplier = staticState.CoriolisMultiplier,
				CoriolisTerm = coriolisTerm,
				PressureGradientForce = pressureGradientForceCloud,
				EvaporationRate = worldData.EvaporationRate,
				EvapTemperatureMax = worldData.EvapMaxTemperature,
				EvapTemperatureMin = worldData.EvapMinTemperature,
			}, JobHandle.CombineDependencies(cloudEnergyJobHandleDependencies));
			energyJobHandles.Add(fluxCloudJobHandle);

			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layerIndex = _airLayer0 + j;
				var airDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
				};
				if (j == 1)
				{
					airDependencies.Add(conductionAirWaterJobHandle);
					airDependencies.Add(conductionAirIceJobHandle);
					airDependencies.Add(conductionAirTerrainJobHandle);
				}
				jobHandleDependencies.Add(airDependencies);
				fluxJobs[j + _airLayer0] = FluxAirJob.Run(new FluxAirJob()
				{
					FluxEnergy = fluxEnergyAir[j],
					CondensationCloudMass = condensationCloudMass[j],
					CondensationGroundMass = condensationGroundMass[j],
					LastTemperature = lastState.AirTemperature[j],
					LastVapor = lastState.AirVapor[j],
					AirMass = dependent.AirMass[j],
					ConductionEnergyWater = j == 1 ? conductionSurfaceAirWater : conductionUpperAirWater,
					ConductionEnergyIce = j == 1 ? conductionSurfaceAirIce : conductionUpperAirIce,
					ConductionEnergyTerrain = j == 1 ? conductionSurfaceAirTerrain : conductionUpperAirTerrain,
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationDelta = thermalRadiationDelta[layerIndex],
					CloudElevation = dependent.CloudElevation,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					DewPointZero = worldData.DewPointZero,
					InverseDewPointTemperatureRange = worldData.inverseDewPointTemperatureRange,
					WaterVaporMassToAirMassAtDewPoint = worldData.WaterVaporMassToAirMassAtDewPoint,
				}, JobHandle.CombineDependencies(airDependencies));
			}

			var fluxIceJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_iceLayer],
				thermalOutJobHandles[_iceLayer],
				thermalInDownJobHandles[_iceLayer],
				thermalInUpJobHandles[_iceLayer],
				conductionAirIceJobHandle,
				conductionIceWaterJobHandle,
				conductionIceTerrainJobHandle,
				fluxCloudJobHandle,
			};
			jobHandleDependencies.Add(fluxIceJobHandleDependencies);
			fluxJobs[_iceLayer] = FluxIceJob.Run(new FluxIceJob()
			{
				Energy = fluxEnergyIce,
				MeltedMass = iceMeltedMass,
				LastMass = lastState.IceMass,
				SolarRadiationIn = solarRadiationIn[_iceLayer],
				ThermalRadiationDeltaBottom = thermalRadiationDeltaIceBottom,
				ThermalRadiationDeltaTop = thermalRadiationDeltaIceTop,
				ConductionEnergyAir = conductionSurfaceAirIce,
				ConductionEnergyTerrain = conductionIceTerrain,
				ConductionEnergyWater = conductionIceWater,
				LastTemperature = lastState.IceTemperature,
				IceHeatingDepth = worldData.IceHeatingDepth,
			}, JobHandle.CombineDependencies(fluxIceJobHandleDependencies));

			#endregion

			#region State Changes - Evaporation, Condensation, Melting, Rainfall

			JobHandle stateChangeJobHandle = default(JobHandle);
			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layerIndex = _airLayer0 + j;

				var fluxAirJobsHandle = new NativeList<JobHandle>(Allocator.TempJob)
				{
					fluxCloudJobHandle,
					fluxJobs[_waterLayer0 + _surfaceWaterLayer],
					fluxJobs[layerIndex],
					stateChangeJobHandle
				};
				jobHandleDependencies.Add(fluxAirJobsHandle);

				stateChangeJobHandle = StateChangeAirLayerJob.Run(new StateChangeAirLayerJob()
				{
					AirTemperature = nextState.AirTemperature[j],
					VaporMass = nextState.AirVapor[j],
					CloudDropletMass = nextState.CloudDropletMass,
					CloudEvaporationMass = cloudEvaporationMass,
					CloudMass = nextState.CloudMass,
					SurfaceWaterMass = nextState.WaterMass[_surfaceWaterLayer],
					SurfaceWaterTemperature = nextState.WaterTemperature[_surfaceWaterLayer],

					SurfaceSaltMass = nextState.SaltMass[_surfaceWaterLayer],
					CloudCondensationMass = condensationCloudMass[j],
					GroundCondensationMass = condensationGroundMass[j],
					AirMass = dependent.AirMass[j],
					LastAirTemperature = lastState.AirTemperature[j],
					LastVaporMass = lastState.AirVapor[j],
					LayerIndex = layerIndex,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],

				}, JobHandle.CombineDependencies(fluxAirJobsHandle));
			}

			var fluxJobsHandle = new NativeList<JobHandle>(Allocator.TempJob)
			{
				fluxCloudJobHandle,
				fluxJobs[_iceLayer],
				fluxJobs[_waterLayer0 + _surfaceWaterLayer],
				fluxJobs[_airLayer0 + 1],
				stateChangeJobHandle
			};
			jobHandleDependencies.Add(fluxJobsHandle);
			stateChangeJobHandle = StateChangeJob.Run(new StateChangeJob()
			{
				IceTemperature = nextState.IceTemperature,
				IceMass = nextState.IceMass,
				SurfaceWaterTemperature = nextState.WaterTemperature[_surfaceWaterLayer],
				SurfaceWaterMass = nextState.WaterMass[_surfaceWaterLayer],
				SurfaceAirTemperature = nextState.AirTemperature[1],
				SurfaceAirVapor = nextState.AirVapor[1],

				LastIceMass = lastState.IceMass,
				LastIceTemperature = lastState.IceTemperature,
				SurfaceAirMass = dependent.AirMass[1],
				SurfaceSaltMass = lastState.SaltMass[_surfaceWaterLayer],
				WaterEvaporatedMass = evaporationMass[_surfaceWaterLayer],
				WaterFrozenMass = frozenMass[_surfaceWaterLayer],
				WaterFrozenTemperature = frozenTemperature[_surfaceWaterLayer],
				IceMeltedMass = iceMeltedMass,
				PrecipitationTemperature = precipitationTemperature,
				PrecipitationMass = precipitationMass,
			}, JobHandle.CombineDependencies(fluxJobsHandle));


			#endregion


			#region COMBINE ADVECTION, DIFFUSION, SOLAR, THERMAL DELTA, and State changes within each layer

			var terrainEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_terrainLayer],
				thermalOutJobHandles[_terrainLayer],
				thermalInDownJobHandles[_terrainLayer],
				thermalInUpJobHandles[_terrainLayer],
				conductionAirTerrainJobHandle,
				conductionIceTerrainJobHandle,
				conductionWaterTerrainJobHandle,
			};
			jobHandleDependencies.Add(terrainEnergyJobHandleDependencies);
			energyJobHandles.Add(EnergyTerrainJob.Run(new EnergyTerrainJob()
			{
				Temperature = nextState.TerrainTemperature,
				LastTemperature = lastState.TerrainTemperature,
				Terrain = lastState.Terrain,
				SolarRadiationIn = solarRadiationIn[_terrainLayer],
				ThermalRadiationDelta = thermalRadiationDelta[_terrainLayer],
				ConductionEnergyAir = conductionSurfaceAirTerrain,
				ConductionEnergyIce = conductionIceTerrain,
				ConductionEnergyWater = conductionWaterTerrainTotal,
				VegetationCoverage = dependent.VegetationCoverage,
				GeothermalEnergy = nextState.PlanetState.GeothermalHeat * worldData.SecondsPerTick,
				HeatingDepth = worldData.SoilHeatDepth
			}, JobHandle.CombineDependencies(terrainEnergyJobHandleDependencies)));

			var energyIceHandle = EnergyIceJob.Run(new EnergyIceJob()
			{
				Temperature = nextState.IceTemperature,
				Mass = nextState.IceMass,
				FluxEnergyIce = fluxEnergyIce,
				IceHeatingDepth = worldData.IceHeatingDepth,
			}, stateChangeJobHandle);
			energyJobHandles.Add(energyIceHandle);

			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layerIndex = _airLayer0 + j;

				var airDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					pgfJobHandles[j],
					windFrictionJobHandle,
					stateChangeJobHandle,
				};
				jobHandleDependencies.Add(airDependencies);
				energyJobHandles.Add(EnergyAirJob.Run(new EnergyAirJob()
				{
					Temperature = nextState.AirTemperature[j],
					Vapor = nextState.AirVapor[j],
					Wind = nextState.Wind[j],
					LastWind = lastState.Wind[j],
					AirMass = dependent.AirMass[j],
					Energy = fluxEnergyAir[j],
					PressureGradientForce = pressureGradientForce[j],
					Position = staticState.SphericalPosition,
					CoriolisMultiplier = staticState.CoriolisMultiplier,
					CoriolisTerm = coriolisTerm,
					SecondsPerTick = worldData.SecondsPerTick,
					WindFriction = windFriction,
					WindFrictionMultiplier = j==1 ? 1 : 0,
				}, JobHandle.CombineDependencies(airDependencies)));
			}

			JobHandle surfaceWaterJobHandle = default(JobHandle);
			for (int j = 1; j < _waterLayers - 1; j++)
			{
				int layerIndex = _waterLayer0 + j;

				var waterDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
					conductionWaterTerrainJobHandle,
					waterDensityForceJobHandles[j]
				};
				if (j == _surfaceWaterLayer)
				{
					waterDependencies.Add(stateChangeJobHandle);
				}
				jobHandleDependencies.Add(waterDependencies);
				var energyWaterJobHandle = EnergyWaterJob.Run(new EnergyWaterJob()
				{
					Temperature = nextState.WaterTemperature[j],
					SaltMass = nextState.SaltMass[j],
					Velocity = nextState.WaterVelocity[j],
					Mass = nextState.WaterMass[j],
					LastTemperature = lastState.WaterTemperature[j],
					LastSaltMass = lastState.SaltMass[j],
					LastVelocity = lastState.WaterVelocity[j],
					LastWaterMass = lastState.WaterMass[j],
					EvaporatedMass = evaporationMass[j],
					FrozenMass = frozenMass[j],
					Energy = fluxEnergyWater[j],
					Force = waterForce[j],
					SecondsPerTick = worldData.SecondsPerTick,
				}, JobHandle.CombineDependencies(waterDependencies));
				energyJobHandles.Add(energyWaterJobHandle);
				if (j == _surfaceWaterLayer)
				{
					surfaceWaterJobHandle = energyWaterJobHandle;
				}
			}


			energyJobHandles.Add(updateTerrainJobHandle);
			var energyJobHandle = JobHandle.CombineDependencies(energyJobHandles);

			#endregion


			// Diffuse from last time step
			// Air, Water, Cloud
			#region Diffusion

			JobHandle[] diffusionJobHandles = new JobHandle[_layerCount];
			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layer = _airLayer0 + j;
				// TODO: is it a problem that we are using the dependent variables from last frame while referencing our newly calculated next frame values for temperature and such?
				diffusionJobHandles[layer] = DiffusionAirJob.Run(new DiffusionAirJob()
				{
					Delta = diffusionAir[j],

					LastTemperature = nextState.AirTemperature[j],
					LastVapor = nextState.AirVapor[j],
					LastWind = nextState.Wind[j],
					Neighbors = staticState.Neighbors,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					AirMass = dependent.AirMass[j],
					UpTemperature = nextState.AirTemperature[j + 1],
					UpHumidity = nextState.AirVapor[j + 1],
					UpWind = nextState.Wind[j + 1],
					UpAirMass = dependent.AirMass[j + 1],
					UpLayerElevation = dependent.LayerElevation[j + 1],
					UpLayerHeight = dependent.LayerHeight[j + 1],
					DownTemperature = nextState.AirTemperature[j - 1],
					DownHumidity = nextState.AirVapor[j - 1],
					DownWind = nextState.Wind[j - 1],
					DownAirMass = dependent.AirMass[j - 1],
					DownLayerElevation = dependent.LayerElevation[j - 1],
					DownLayerHeight = dependent.LayerHeight[j - 1],
					IsTop = j == _airLayers - 2,
					IsBottom = j == 1,
					DiffusionCoefficientHoriztonal = worldData.AirDiffusionCoefficientHorizontal,
					DiffusionCoefficientVertical = worldData.AirDiffusionCoefficientVertical,
				}, energyJobHandle);
			}
			for (int j = 1; j < _waterLayers - 1; j++)
			{
				int layer = _waterLayer0 + j;
				diffusionJobHandles[layer] = DiffusionWaterJob.Run(new DiffusionWaterJob()
				{
					Delta = diffusionWater[j],

					LastTemperature = nextState.WaterTemperature[j],
					LastSalt = nextState.SaltMass[j],
					LastCurrent = nextState.WaterVelocity[j],
					LastMass = nextState.WaterMass[j],
					UpTemperature = nextState.WaterTemperature[j + 1],
					UpSalt = nextState.SaltMass[j + 1],
					UpCurrent = nextState.WaterVelocity[j + 1],
					UpMass = nextState.WaterMass[j + 1],
					DownTemperature = nextState.WaterTemperature[j - 1],
					DownSalt = nextState.SaltMass[j - 1],
					DownCurrent = nextState.WaterVelocity[j - 1],
					DownMass = nextState.WaterMass[j - 1],
					Neighbors = staticState.Neighbors,
					DiffusionCoefficientHoriztonal = worldData.WaterDiffusionCoefficientHorizontal,
					DiffusionCoefficientVertical = worldData.WaterDiffusionCoefficientVertical,
				}, energyJobHandle);
			}

			var diffusionCloudHandle = DiffusionCloudJob.Run(new DiffusionCloudJob()
			{
				Delta = diffusionCloud,

				LastMass = nextState.CloudMass,
				LastDropletMass = nextState.CloudDropletMass,
				LastVelocity = nextState.CloudVelocity,
				Neighbors = staticState.Neighbors,
				DiffusionCoefficient = worldData.CloudDiffusionCoefficient,
			}, energyJobHandle);

			#endregion

			#region Apply Diffusion

			JobHandle diffusionJobHandle = default(JobHandle);
			diffusionJobHandle = JobHandle.CombineDependencies(diffusionJobHandle, ApplyAdvectionCloudJob.Run(new ApplyAdvectionCloudJob()
			{
				Advection = diffusionCloud,
				CloudMass = nextState.CloudMass,
				DropletMass = nextState.CloudDropletMass,
				Velocity = nextState.CloudVelocity
			}, diffusionCloudHandle));

			for (int i = 1; i < _waterLayers - 1; i++)
			{
				diffusionJobHandle = JobHandle.CombineDependencies(diffusionJobHandle, ApplyAdvectionWaterJob.Run(new ApplyAdvectionWaterJob()
				{
					Advection = diffusionWater[i],
					SaltMass = nextState.SaltMass[i],
					Temperature = nextState.WaterTemperature[i],
					Velocity = nextState.WaterVelocity[i],
					Mass = nextState.WaterMass[i]
				}, JobHandle.CombineDependencies(diffusionJobHandles[i + _waterLayer0], diffusionJobHandles[i + _waterLayer0 + 1], diffusionJobHandles[i + _waterLayer0 - 1])));
			}

			for (int i = 1; i < _airLayers - 1; i++)
			{
				diffusionJobHandle = JobHandle.CombineDependencies(diffusionJobHandle, ApplyAdvectionAirJob.Run(new ApplyAdvectionAirJob()
				{
					Advection = diffusionAir[i],
					Vapor = nextState.AirVapor[i],
					Temperature = nextState.AirTemperature[i],
					Wind = nextState.Wind[i],
				}, JobHandle.CombineDependencies(diffusionJobHandles[i + _airLayer0], diffusionJobHandles[i + _airLayer0 - 1], diffusionJobHandles[i + _airLayer0 + 1])));
			}

			#endregion




			// Wind and currents move temperature and trace elements horizontally
			// Air, Water, Cloud
			// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
			#region Advection

			JobHandle[] advectionJobHandles = new JobHandle[_layerCount];
			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layer = _airLayer0 + j;
				var airDestJob = GetVectorDestCoordsJob.Run(new GetVectorDestCoordsJob()
				{
					Destination = destinationAir[j],
					Neighbors = staticState.Neighbors,
					Position = staticState.SphericalPosition,
					Velocity = nextState.Wind[j],
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick
				}, diffusionJobHandle);
				advectionJobHandles[layer] = AdvectionAirJob.Run(new AdvectionAirJob()
				{
					Delta = advectionAir[j],
					Temperature = nextState.AirTemperature[j],
					Vapor = nextState.AirVapor[j],
					Velocity = nextState.Wind[j],
					Neighbors = staticState.Neighbors,
					Destination = destinationAir[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					UpTemperature = nextState.AirTemperature[j + 1],
					UpHumidity = nextState.AirVapor[j + 1],
					UpAirMass = dependent.AirMass[j + 1],
					UpLayerElevation = dependent.LayerElevation[j + 1],
					UpLayerHeight = dependent.LayerHeight[j + 1],
					DownTemperature = nextState.AirTemperature[j - 1],
					DownHumidity = nextState.AirVapor[j - 1],
					DownAirMass = dependent.AirMass[j - 1],
					DownLayerElevation = dependent.LayerElevation[j - 1],
					DownLayerHeight = dependent.LayerHeight[j - 1],
					IsTop = j == _airLayers - 2,
					IsBottom = j == 1,
				}, airDestJob);
			}
			for (int j = 1; j < _waterLayers - 1; j++)
			{
				int layer = _waterLayer0 + j;
				var waterDestJob = GetVectorDestCoordsJob.Run(new GetVectorDestCoordsJob()
				{
					Destination = destinationWater[j],
					Neighbors = staticState.Neighbors,
					Position = staticState.SphericalPosition,
					Velocity = nextState.WaterVelocity[j],
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick
				}, diffusionJobHandle);
				advectionJobHandles[layer] = AdvectionWaterJob.Run(new AdvectionWaterJob()
				{
					Delta = advectionWater[j],
					Destination = destinationWater[j],
					Velocity = nextState.WaterVelocity[j],
					Temperature = nextState.WaterTemperature[j],
					Mass = nextState.WaterMass[j],
					Salt = nextState.SaltMass[j],
					Neighbors = staticState.Neighbors,
				}, waterDestJob);
			}

			var cloudDestJob = GetVectorDestCoordsJob.Run(new GetVectorDestCoordsJob()
			{
				Destination = destinationCloud,
				Neighbors = staticState.Neighbors,
				Position = staticState.SphericalPosition,
				Velocity = nextState.CloudVelocity,
				PlanetRadius = staticState.PlanetRadius,
				SecondsPerTick = worldData.SecondsPerTick
			}, diffusionJobHandle);
			var advectionJobHandleCloud = AdvectionCloudJob.Run(new AdvectionCloudJob()
			{
				Delta = advectionCloud,
				Destination = destinationCloud,
				Mass = nextState.CloudMass,
				DropletMass = nextState.CloudDropletMass,
				Velocity = nextState.CloudVelocity,
				Neighbors = staticState.Neighbors,
			}, cloudDestJob);


			#endregion

			#region Apply Advection

			NativeList<JobHandle> applyAdvectionJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
			applyAdvectionJobHandles.Add(ApplyAdvectionCloudJob.Run(new ApplyAdvectionCloudJob()
			{
				Advection = advectionCloud,
				CloudMass = nextState.CloudMass,
				DropletMass = nextState.CloudDropletMass,
				Velocity = nextState.CloudVelocity
			}, advectionJobHandleCloud));

			for (int i = 1; i < _waterLayers - 1; i++)
			{
				applyAdvectionJobHandles.Add(ApplyAdvectionWaterJob.Run(new ApplyAdvectionWaterJob()
				{
					Advection = advectionWater[i],
					SaltMass = nextState.SaltMass[i],
					Temperature = nextState.WaterTemperature[i],
					Velocity = nextState.WaterVelocity[i],
					Mass = nextState.WaterMass[i]
				}, advectionJobHandles[i+_waterLayer0]));
			}

			for (int i = 1; i < _airLayers - 1; i++)
			{
				applyAdvectionJobHandles.Add(ApplyAdvectionAirJob.Run(new ApplyAdvectionAirJob()
				{
					Advection = advectionAir[i],
					Vapor = nextState.AirVapor[i],
					Temperature = nextState.AirTemperature[i],
					Wind = nextState.Wind[i],
				}, JobHandle.CombineDependencies( advectionJobHandles[i+_airLayer0], advectionJobHandles[i + _airLayer0 - 1], advectionJobHandles[i + _airLayer0 + 1])));
			}
			lastJobHandle = JobHandle.CombineDependencies(applyAdvectionJobHandles);

			#endregion

			lastJobHandle.Complete();

			#region Update dependent variables

			var waterMassTotal = new NativeArray<float>(nextState.IceMass, Allocator.TempJob);

			NativeList<JobHandle> updateDependenciesJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
			var updateDependentWaterLayerJobHandle = lastJobHandle;
			for (int j = _waterLayers - 2; j > 0; j--)
			{
				updateDependentWaterLayerJobHandle = UpdateDependentWaterLayerJob.Run(new UpdateDependentWaterLayerJob()
				{
					Salinity = dependent.WaterSalinity[j],
					WaterCoverage = dependent.WaterCoverage[j],
					PotentialEnergy = dependent.WaterPotentialEnergy[j],
					Density = dependent.WaterDensity[j],
					Pressure = dependent.WaterPressure[j],
					WaterDepthTotal = waterDepthTotal,
					WaterMassTotal = waterMassTotal,
					LayerDepth = dependent.WaterLayerDepth[j],
					LayerHeight = dependent.WaterLayerHeight[j],

					UpLayerDepth = dependent.WaterLayerDepth[j + 1],
					UpLayerHeight = dependent.WaterLayerHeight[j + 1],
					Temperature = nextState.WaterTemperature[j],
					SaltMass = nextState.SaltMass[j],
					WaterMass = nextState.WaterMass[j],
					Terrain = nextState.Terrain,
					WaterDensityPerDegree = worldData.WaterDensityPerDegree,
					WaterDensityPerSalinity = worldData.WaterDensityPerSalinity,
					Gravity = nextState.PlanetState.Gravity,
				}, updateDependentWaterLayerJobHandle);
				updateDependenciesJobHandles.Add(updateDependentWaterLayerJobHandle);
			}

			var updateDependenciesJobHandle = UpdateDependentJob.Run(new UpdateDependentStateJob()
			{
				CloudCoverage = dependent.CloudCoverage,
				IceCoverage = dependent.IceCoverage,
				IceEnergy = dependent.IceEnergy,
				SurfaceElevation = dependent.SurfaceElevation,
				VegetationCoverage = dependent.VegetationCoverage,
				SurfaceAirTemperature = dependent.SurfaceAirTemperature,
				WaterDepth = dependent.WaterDepth,

				WaterDepthTotal = waterDepthTotal,
				CloudMass = nextState.CloudMass,
				IceMass = nextState.IceMass,
				IceTemperature = nextState.IceTemperature,
				Terrain = nextState.Terrain,
				worldData = worldData,
			}, JobHandle.CombineDependencies(updateDependentWaterLayerJobHandle, lastJobHandle));
			updateDependenciesJobHandles.Add(updateDependenciesJobHandle);


			JobHandle updateDependentAirLayerJobHandle = updateDependenciesJobHandle;
			for (int j = _airLayers - 2; j > 0; j--)
			{
				updateDependentAirLayerJobHandle = UpdateDependentAirLayerJob.Run(new UpdateDependentAirLayerJob()
				{
					Pressure = dependent.AirPressure[j],
					RelativeHumidity = dependent.AirHumidityRelative[j],
					AbsoluteHumidity = dependent.AirHumidityAbsolute[j],
					AirMass = dependent.AirMass[j],
					PotentialEnergy = dependent.AirPotentialEnergy[j],
					AirMassCloud = dependent.AirMassCloud,
					AirVaporCloud = dependent.AirVaporCloud,
					AirPressureCloud = dependent.AirPressureCloud,
					AirHumidityRelativeCloud = dependent.AirHumidityRelativeCloud,
					AirLayerCloud = dependent.AirLayerCloud,
					CloudElevation = dependent.CloudElevation,
					DewPoint = dependent.DewPoint,
					AirMassTotal = airMassTotal,

					AirTemperature = lastState.AirTemperature[j],
					CloudDropletMass = nextState.CloudDropletMass,
					CloudMass = nextState.CloudMass,
					VaporMass = nextState.AirVapor[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					SurfaceElevation = dependent.SurfaceElevation,
					IceMass = nextState.IceMass,
					Gravity = nextState.PlanetState.Gravity,
					DewPointZero = worldData.DewPointZero,
					InverseDewPointTemperatureRange = worldData.inverseDewPointTemperatureRange,
					WaterVaporMassToAirMassAtDewPoint = worldData.WaterVaporMassToAirMassAtDewPoint,
					LayerIndex = j,
				}, updateDependentAirLayerJobHandle);
				updateDependenciesJobHandles.Add(updateDependentAirLayerJobHandle);
			}

			var updateSurfaceDependenciesJobHandle = UpdateDependentJob.Run(new UpdateSurfaceDependentStateJob()
			{
				SurfaceAirTemperature = dependent.SurfaceAirTemperature,

				LowerAirTemperature = nextState.AirTemperature[1],
				lowerAirHeight = dependent.LayerHeight[1]
			}, updateDependentAirLayerJobHandle);
			updateDependenciesJobHandles.Add(updateSurfaceDependenciesJobHandle);

			lastJobHandle = JobHandle.CombineDependencies(updateDependenciesJobHandles);

			#endregion

			lastJobHandle.Complete();

			#region Debug
			if (checkForDegeneracy)
			{
				bool degen = false;
				SortedSet<int> degenIndices = new SortedSet<int>();
				List<string> degenVarNames = new List<string>();
				degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "TerrainTemperature", nextState.TerrainTemperature, 0, 1000, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "CloudMass", nextState.CloudMass, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "CloudDropletMass", nextState.CloudDropletMass, degenVarNames);
				degen |= CheckDegen(_cellCount, degenIndices, "CloudVelocity", nextState.CloudVelocity, degenVarNames);
				degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "CloudElevation", dependent.CloudElevation, -100000, 100000, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "IceMass", nextState.IceMass, degenVarNames);
				degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "IceTemperature", nextState.IceTemperature, 0, 300, degenVarNames);
				for (int i = 1; i < _airLayers - 1; i++) {
					degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "AirTemperature" + i, nextState.AirTemperature[i], 0, 1000, degenVarNames);
					degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "AirVapor" + i, nextState.AirVapor[i], 0, 1000, degenVarNames);
					degen |= CheckDegen(_cellCount, degenIndices, "Wind" + i, nextState.Wind[i], degenVarNames);
				}
				for (int i=1;i<_waterLayers - 1;i++)
				{
					degen |= CheckDegenPosValues(_cellCount, degenIndices, "WaterMass" + i, nextState.WaterMass[i], degenVarNames);
					degen |= CheckDegenPosValues(_cellCount, degenIndices, "SaltMass" + i, nextState.SaltMass[i], degenVarNames);
					degen |= CheckDegenPosValues(_cellCount, degenIndices, "WaterTemperature" + i, nextState.WaterTemperature[i], degenVarNames);
					degen |= CheckDegen(_cellCount, degenIndices, "Current" + i, nextState.WaterVelocity[i], degenVarNames);
				}
				if (degen)
				{
					foreach (var i in degenIndices)
					{
						PrintState("Degenerate", i, staticState, nextState, degenVarNames);
						PrintDependentState("Dependent Vars", i, dependent);
						PrintState("Last State", i, staticState, lastState, new List<string>());
					}
					Debug.Break();
				}
			}

			if (logState)
			{
				PrintState("State", logStateIndex, staticState, nextState, new List<string>());
				PrintDependentState("Dependent Vars", logStateIndex, dependent);
			}
			#endregion

			#region Update Display
			if (tick == ticksToAdvance-1)
			{
				var curState = states[curStateIndex];

				display.Dispose();
				display = new DisplayState();
				display.Init(_cellCount, _airLayers, _waterLayers);
				JobHandle initDisplayHandle = default(JobHandle);
				for (int i = 1; i < _airLayers - 1; i++)
				{
					initDisplayHandle = JobHandle.CombineDependencies(initDisplayHandle, (new InitDisplayAirLayerJob()
					{
						DisplayPotentialTemperature = display.PotentialTemperature[i],
						DisplayPressure = display.Pressure[i],
						DisplayPressureGradientForce = display.PressureGradientForce[i],

						Gravity = curState.PlanetState.Gravity,
						AirTemperature = curState.AirTemperature[i],
						AirLayerElevation = dependent.LayerElevation[i],
						AirLayerHeight = dependent.LayerHeight[i],
						AirPressure = dependent.AirPressure[i],
						PressureGradientForce = pressureGradientForce[i],
					}).Schedule(_cellCount, _batchCount));
				}

				var updateDisplayJob = new UpdateDisplayJob()
				{
					SolarRadiationAbsorbedSurface = display.SolarRadiationAbsorbedSurface,
					DisplayEvaporation = display.Evaporation,
					DisplayRainfall = display.Rainfall,

					SolarRadiationInTerrain = solarRadiationIn[_terrainLayer],
					SolarRadiationInIce = solarRadiationIn[_iceLayer],
					SolarRadiationInWaterSurface = solarRadiationIn[_waterLayer0 + _waterLayers - 2],
					Evaporation = evaporationMass[_surfaceWaterLayer],
					RainfallWater = precipitationMass,
				};
				var updateDisplayJobHandle = updateDisplayJob.Schedule(_cellCount, _batchCount);
				var displayHandles = JobHandle.CombineDependencies(initDisplayHandle, updateDisplayJobHandle);
				displayHandles.Complete();

				for (int i = 0; i < _cellCount; i++)
				{
					display.SolarRadiation += displaySolarRadiation[i];
					display.GlobalCloudCoverage += dependent.CloudCoverage[i];
					display.GlobalCloudMass += curState.CloudMass[i];
					display.GlobalIceMass += curState.IceMass[i];
					display.GlobalOceanCoverage += dependent.WaterCoverage[_waterLayers-2][i];
					display.GlobalTemperature += curState.AirTemperature[1][i];
					display.GlobalWaterVapor += curState.AirVapor[1][i];
					display.GlobalOceanVolume += dependent.WaterDepth[i];
					display.GlobalSeaLevel += dependent.SurfaceElevation[i];
					display.GlobalEvaporation += display.Evaporation[i];
					display.GlobalRainfall += display.Rainfall[i];
					for (int j = 1; j < _airLayers - 1; j++)
					{
						display.EnergySolarReflectedAtmosphere += solarReflected[j + _airLayer0][i];
						display.EnergySolarAbsorbedAtmosphere += solarRadiationIn[j + _airLayer0][i];
					}
					for (int j = 1; j < _waterLayers - 1; j++)
					{
						float absorbed = solarRadiationIn[j + _waterLayer0][i];
						display.EnergySolarAbsorbedOcean += absorbed;
						display.EnergySolarAbsorbedSurface += absorbed;
						display.EnergySolarReflectedSurface += solarReflected[j + _waterLayer0][i];
					}
					display.EnergySolarAbsorbedSurface += solarRadiationIn[_terrainLayer][i] + solarRadiationIn[_iceLayer][i];
					display.EnergySolarReflectedSurface += solarReflected[_terrainLayer][i] + solarReflected[_iceLayer][i];
					display.EnergySurfaceConduction += conductionSurfaceAirIce[i] + conductionSurfaceAirTerrain[i] + conductionSurfaceAirWater[i];
					display.EnergyOceanConduction += conductionSurfaceAirWater[i];
					display.EnergyEvapotranspiration += evaporationMass[_surfaceWaterLayer][i] * WorldData.LatentHeatWaterVapor;
					//display.EnergyThermalAbsorbedAtmosphere += ;
					display.EnergyThermalBackRadiation += windowRadiationTransmittedDown[_airLayer0 + 1][i] + thermalRadiationTransmittedDown[_airLayer0 + 1][i];
					display.EnergyThermalOceanRadiation += (windowRadiationTransmittedUp[_waterLayer0 + _waterLayers - 2][i] + thermalRadiationTransmittedUp[_waterLayer0 + _waterLayers - 2][i]) * dependent.WaterCoverage[_waterLayers - 2][i];
					display.EnergyThermalOutAtmosphere += thermalRadiationTransmittedUp[_airLayer0 + _airLayers - 2][i];
					display.EnergyThermalSurfaceOutAtmosphericWindow += windowRadiationTransmittedUp[_iceLayer][i];
					display.EnergyThermalSurfaceRadiation += windowRadiationTransmittedUp[_iceLayer][i] + thermalRadiationTransmittedUp[_iceLayer][i];
			
				}

			}

			#endregion

			#region Dispose Temporary Arrays
			updateDependenciesJobHandles.Dispose();
			energyJobHandles.Dispose();
			atmosphericWindowUp.Dispose();
			atmosphericWindowDown.Dispose();
			dropletDelta.Dispose();
			precipitationMass.Dispose();
			precipitationTemperature.Dispose();
			cloudEvaporationMass.Dispose();
			iceMeltedMass.Dispose();
			airMassTotal.Dispose();
			waterMassTotal.Dispose();
			waterDepthTotal.Dispose();
			fluxWaterDependencies.Dispose();
			applyAdvectionJobHandles.Dispose();
			foreach (var d in jobHandleDependencies)
			{
				d.Dispose();
			}
			thermalRadiationDeltaIceTop.Dispose();
			thermalRadiationDeltaIceBottom.Dispose();
			thermalRadiationDeltaSurfaceWater.Dispose();
			conductionWaterTerrainTotal.Dispose();
			for (int i = 0; i < _layerCount; i++)
			{
				thermalRadiationDelta[i].Dispose();
				thermalRadiationTransmittedUp[i].Dispose();
				thermalRadiationTransmittedDown[i].Dispose();
				windowRadiationTransmittedUp[i].Dispose();
				windowRadiationTransmittedDown[i].Dispose();
				solarReflected[i].Dispose();
				condensationGroundMass[i].Dispose();
				condensationCloudMass[i].Dispose();
			}
			#endregion


		}


	}

	private static bool CheckDegen(int count, SortedSet<int> degenIndices, string name, NativeArray<float> values, List<string> degenVarNames)
	{
		for (int i = 0; i < count; i++)
		{
			float v = values[i];
			if (!math.isfinite(v))
			{
				degenVarNames.Add(name);
				degenIndices.Add(i);
				return true;
			}
		}
		return false;
	}
	private static bool CheckDegen(int count, SortedSet<int> degenIndices, string name, NativeArray<float3> values, List<string> degenVarNames)
	{
		for (int i = 0; i < count; i++)
		{
			if (!math.isfinite(values[i].x) || !math.isfinite(values[i].y) || !math.isfinite(values[i].z))
			{
				degenVarNames.Add(name);
				degenIndices.Add(i);
				return true;
			}
		}
		return false;
	}

	private static bool CheckDegenPosValues(int count, SortedSet<int> degenIndices, string name, NativeArray<float> values, List<string> degenVarNames)
	{
		for (int i = 0; i < count; i++)
		{
			if (!math.isfinite(values[i]) || values[i] < 0)
			{
				degenVarNames.Add(name);
				degenIndices.Add(i);
				return true;
			}
		}
		return false;
	}

	private static bool CheckDegenMinMaxValues(int count, SortedSet<int> degenIndices, string name, NativeArray<float> values, float min, float max, List<string> degenVarNames)
	{
		for (int i = 0; i < count; i++)
		{
			if (!math.isfinite(values[i]) || values[i] < min || values[i] > max)
			{
				degenVarNames.Add(name);
				degenIndices.Add(i);
				return true;
			}
		}
		return false;
	}

	public void PrintState(string title, int i, StaticState staticState, SimState state, List<string> degenVarNames)
	{
		StringBuilder s = new StringBuilder();
		s.AppendFormat("{0} Index: {1} Time: {2}", title, i, state.PlanetState.Ticks);
		foreach (var n in degenVarNames)
		{
			s.AppendFormat(" | {0}", n);
		}
		s.AppendLine("");
		s.AppendFormat("X: {0} Y: {1}\n", staticState.Coordinate[i].x, staticState.Coordinate[i].y);
		s.AppendFormat("Elevation: {0}\n", state.Terrain[i].Elevation);
		s.AppendFormat("Roughness: {0}\n", state.Terrain[i].Roughness);
		s.AppendFormat("SoilFertility: {0}\n", state.Terrain[i].SoilFertility);
		s.AppendFormat("Vegetation: {0}\n", state.Terrain[i].Vegetation);
		s.AppendFormat("TerrainTemperature: {0}\n", state.TerrainTemperature[i]);
		s.AppendFormat("CloudMass: {0}\n", state.CloudMass[i]);
		s.AppendFormat("CloudDropletMass: {0}\n", state.CloudDropletMass[i]);
		s.AppendFormat("CloudVelocity: {0}\n", state.CloudVelocity[i]);
		s.AppendFormat("IceMass: {0}\n", state.IceMass[i]);
		s.AppendFormat("IceTemperature: {0}\n", state.IceTemperature[i]);
		for (int j = 1; j < _waterLayers - 1; j++)
		{
			s.AppendFormat("WaterMass{0}: {1}\n", j, state.WaterMass[j][i]);
			s.AppendFormat("SaltMass{0}: {1}\n", j, state.SaltMass[j][i]);
			s.AppendFormat("WaterTemperature{0}: {1}\n", j, state.WaterTemperature[j][i]);
			s.AppendFormat("WaterVelocity{0}: {1}\n", j, state.WaterVelocity[j][i]);
		}
		for (int j = 1; j < _airLayers - 1; j++)
		{
			s.AppendFormat("AirTemperature{0}: {1}\n", j, state.AirTemperature[j][i]);
			s.AppendFormat("AirVapor{0}: {1}\n", j, state.AirVapor[j][i]);
			s.AppendFormat("Wind{0}: {1}\n", j, state.Wind[j][i]);
		}
		Debug.Log(s);
	}
	public void PrintDependentState(string title, int i, DependentState dependent)
	{
		StringBuilder s = new StringBuilder();
		s.AppendFormat("{0} Index: {1}", title, i);
		s.AppendFormat("Surface Elevation: {0}\n", dependent.SurfaceElevation[i]);
		s.AppendFormat("Water Depth: {0}\n", dependent.WaterDepth[i]);
		s.AppendFormat("Ice Coverage: {0}\n", dependent.IceCoverage[i]);
		s.AppendFormat("Cloud Coverage: {0}\n", dependent.CloudCoverage[i]);
		s.AppendFormat("Cloud Elevation: {0}\n", dependent.CloudElevation[i]);
		for (int j = 1; j < _waterLayers - 1; j++)
		{
			s.AppendFormat("Water Coverage: {0}\n", dependent.WaterCoverage[j][i]);
		}
		for (int j = 1; j < _airLayers - 1; j++)
		{
		}
		Debug.Log(s);
	}
}
