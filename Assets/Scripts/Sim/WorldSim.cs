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

public class WorldSim {

	private int _cellCount;
	private int _batchCount = 100;
	private const int _terrainLayers = 3;

	private int _airLayers;
	private int _waterLayers;
	private int _layerCount;
	private int _terrainLayer;
	private int _waterLayer0;
	private int _iceLayer;
	private int _cloudLayer;
	private int _airLayer0;


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
	private NativeArray<float>[] latentHeat;
	private NativeArray<float> conductionCloudAir;
	private NativeArray<float> conductionAirIce;
	private NativeArray<float> conductionAirWater;
	private NativeArray<float> conductionAirTerrain;
	private NativeArray<float> conductionIceWater;
	private NativeArray<float> conductionIceTerrain;
	private NativeArray<float>[] conductionWaterTerrain;
	private NativeArray<float> cloudEvaporation;
	private NativeArray<float> surfaceEvaporation;
	private NativeArray<float> cloudCondensation;
	private NativeArray<float> surfaceCondensation;
	private NativeArray<float2>[] windHorizontal;
	private NativeArray<float2>[] currentHorizontal;

	private NativeArray<float> displaySolarRadiation;

	public WorldSim(int cellCount, int atmosphericLayers, int waterLayers)
	{
		_cellCount = cellCount;
		_airLayers = atmosphericLayers;
		_waterLayers = waterLayers;
		_layerCount = _airLayers + _waterLayers + _terrainLayers;
		_terrainLayer = 0;
		_waterLayer0 = _terrainLayer + 1;
		_iceLayer = _waterLayer0 + _waterLayers;
		_airLayer0 = _iceLayer + 1;
		_cloudLayer = _layerCount - 1;

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
		windHorizontal = new NativeArray<float2>[_airLayers];
		for (int i = 0; i < _airLayers; i++)
		{
			diffusionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
			advectionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
			windHorizontal[i] = new NativeArray<float2>(_cellCount, Allocator.Persistent);
		}
		diffusionWater = new NativeArray<DiffusionWater>[_waterLayers];
		advectionWater = new NativeArray<DiffusionWater>[_waterLayers];
		conductionWaterTerrain = new NativeArray<float>[_waterLayers];
		currentHorizontal = new NativeArray<float2>[_waterLayers];
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
			advectionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
			conductionWaterTerrain[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			currentHorizontal[i] = new NativeArray<float2>(_cellCount, Allocator.Persistent);
		}
		diffusionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		advectionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		conductionCloudAir = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirIce = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionIceWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionIceTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		displaySolarRadiation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		cloudEvaporation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		surfaceEvaporation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		cloudCondensation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		surfaceCondensation = new NativeArray<float>(_cellCount, Allocator.Persistent);

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
			diffusionAir[i].Dispose();
			advectionAir[i].Dispose();
			windHorizontal[i].Dispose();
		}
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i].Dispose();
			advectionWater[i].Dispose();
			conductionWaterTerrain[i].Dispose();
			currentHorizontal[i].Dispose();
		}
		diffusionCloud.Dispose();
		advectionCloud.Dispose();
		conductionCloudAir.Dispose();
		conductionAirIce.Dispose();
		conductionAirWater.Dispose();
		conductionAirTerrain.Dispose();
		conductionIceWater.Dispose();
		conductionIceTerrain.Dispose();
		cloudEvaporation.Dispose();
		surfaceEvaporation.Dispose();
		cloudCondensation.Dispose();
		surfaceCondensation.Dispose();

		displaySolarRadiation.Dispose();

	}

	public void Tick(SimState[] states, int stateCount, int ticksToAdvance, ref DependentState dependent, ref DisplayState display, ref StaticState staticState, ref WorldData worldData, ref int curStateIndex)
	{
		JobHandle lastJobHandle = default(JobHandle);
		for (int tick = 0; tick < ticksToAdvance; tick++)
		{

			var thermalRadiationDeltaIceBottom = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var thermalRadiationDeltaIceTop = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var thermalRadiationDeltaSurfaceWater = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var thermalRadiationDelta = new NativeArray<float>[_layerCount];
			var thermalRadiationTransmittedUp = new NativeArray<float>[_layerCount];
			var thermalRadiationTransmittedDown = new NativeArray<float>[_layerCount];
			var windowRadiationTransmittedUp = new NativeArray<float>[_layerCount];
			var windowRadiationTransmittedDown = new NativeArray<float>[_layerCount];
			var condensationGroundMass = new NativeArray<float>[_layerCount];
			var condensationGroundEnergy = new NativeArray<float>[_layerCount];
			var condensationCloudMass = new NativeArray<float>[_layerCount];
			var condensationCloudEnergy = new NativeArray<float>[_layerCount];
			var solarReflected = new NativeArray<float>[_layerCount];
			var conductionWaterTerrainTotal = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var evaporationMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var evaporationLatentHeat = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var frozenTopMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var frozenBottomMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var rainfallWaterMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var rainfallEvaporationMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var iceMeltedTopMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var iceMeltedBottomMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			for (int i = 0; i < _layerCount; i++)
			{
				solarReflected[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				thermalRadiationDelta[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				thermalRadiationTransmittedUp[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				thermalRadiationTransmittedDown[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				windowRadiationTransmittedUp[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				windowRadiationTransmittedDown[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				condensationGroundMass[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				condensationGroundEnergy[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				condensationCloudMass[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
				condensationCloudEnergy[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
			}

			bool updateDisplay = tick == ticksToAdvance - 1;

			ref var lastState = ref states[curStateIndex];
			curStateIndex = (curStateIndex + 1) % stateCount;
			ref var nextState = ref states[curStateIndex];

			nextState.PlanetState = lastState.PlanetState;
			nextState.PlanetState.Ticks++;

			// TODO: update
			float distanceToSun = math.length(lastState.PlanetState.Position);
			float angleToSun = math.atan2(lastState.PlanetState.Position.z, lastState.PlanetState.Position.x);
			angleToSun += lastState.PlanetState.OrbitSpeed;
			nextState.PlanetState.Position = new float3(math.cos(angleToSun), 0, math.sin(angleToSun)) * distanceToSun;
			nextState.PlanetState.Rotation = new float3(lastState.PlanetState.Rotation.x, Mathf.Repeat(lastState.PlanetState.Rotation.y + lastState.PlanetState.SpinSpeed * worldData.SecondsPerTick, math.PI * 2), 0);

			// SOLAR RADIATION INITIALIZATION

			var solarRadiationJob = new SolarRadiationJob()
			{
				SolarRadiation = solarRadiation,
				DisplaySolarRadiation = displaySolarRadiation,
				WaterSlopeAlbedo = waterSlopeAlbedo,

				SphericalPosition = staticState.SphericalPosition,
				IncomingSolarRadiation = lastState.PlanetState.SolarRadiation * worldData.SecondsPerTick,
				PlanetRotation = quaternion.Euler(lastState.PlanetState.Rotation),
				SunToPlanetDir = math.normalize(lastState.PlanetState.Position),
			};
			var solarInJobHandle = solarRadiationJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			var updateTerrainJob = new UpdateTerrainJob()
			{
				Terrain = nextState.Terrain,

				LastTerrain = lastState.Terrain,
			};
			var updateTerrainJobHandle = updateTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			
			// EMISSIVITY

			JobHandle[] emissivityJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				int layerIndex = _airLayer0 + j;
				var emissivityAirJob = new EmissivityAirJob()
				{
					Emissivity = emissivity[layerIndex],
					AirMass = dependent.AirMass[j],
					VaporMass = lastState.AirVapor[j],
					AbsorptivityAir = worldData.AbsorptivityAir,
					AbsorptivityCarbonDioxide = worldData.AbsorptivityCarbonDioxide,
					AbsorptivityWaterVapor = worldData.AbsorptivityWaterVapor,
					GreenhouseGasConcentration = lastState.PlanetState.CarbonDioxide,
				};
				emissivityJobHandles[layerIndex] = emissivityAirJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}

			for (int j = 0; j < _waterLayers; j++)
			{
				int layerIndex = _waterLayer0 + j;
				var emissivityWaterJob = new EmissivityWaterJob()
				{
					Emissivity = emissivity[layerIndex], 
					WaterMass = lastState.WaterMass[j],
					WaterSaltMass = lastState.WaterSaltMass[j]
				};
				emissivityJobHandles[j + _waterLayer0] = emissivityWaterJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			var emissivityTerrainJob = new EmissivityTerrainJob()
			{
				Emissivity = emissivity[_terrainLayer],
				Terrain = lastState.Terrain,
				VegetationCoverage = dependent.VegetationCoverage
			};
			emissivityJobHandles[_terrainLayer] = emissivityTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle);
		
			
			// SOLAR RADIATION ABSORPTION
			// process each vertical layer in order

			// atmosphere
			JobHandle[] solarInJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				int layerIndex = _airLayer0 + _airLayers - 1 - j;
				int airLayerIndex = layerIndex - _airLayer0;
				var solarInAtmosphereJob = new SolarRadiationAbsorbedAirJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[layerIndex],
					SolarRadiationIncoming = solarRadiation,
					SolarRadiationAbsorbedCloud = solarRadiationIn[_cloudLayer],
					SolarRadiationReflected = solarReflected[layerIndex],
					SolarRadiationReflectedCloud = solarReflected[_cloudLayer],
					AirMass = dependent.AirMass[airLayerIndex],
					VaporMass = lastState.AirVapor[airLayerIndex],
					CloudMass = lastState.CloudMass,
					CloudTemperature = lastState.CloudTemperature,
					CloudDropletMass = lastState.CloudDropletMass,
					CloudElevation = lastState.CloudElevation,
					CloudCoverage = dependent.CloudCoverage,
					WaterSlopeAlbedo = solarRadiationJob.WaterSlopeAlbedo,
					SolarReflectivityAir = worldData.SolarReflectivityAir,
					SolarAbsorptivityAir = worldData.SolarAbsorptivityAir,
					SolarAbsorptivityWaterVapor = worldData.SolarAbsorptivityWaterVapor,
					SolarAbsorptivityCloud = worldData.SolarAbsorptivityCloud,
					LayerElevation = dependent.LayerElevation[airLayerIndex],
					LayerHeight = dependent.LayerHeight[airLayerIndex],
					LayerIndex = airLayerIndex,
					worldData = worldData
				};
				solarInJobHandles[layerIndex] = solarInAtmosphereJob.Schedule(_cellCount, _batchCount, airLayerIndex == _airLayers - 1 ? solarInJobHandle : solarInJobHandles[layerIndex + 1]);
			}

			
			// ice
			var solarInIceJob = new SolarRadiationAbsorbedIceJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_iceLayer],
				SolarRadiationIncoming = solarRadiation,
				SolarRadiationReflected = solarReflected[_iceLayer],
				AlbedoIce = WorldData.AlbedoIce,
				IceCoverage = dependent.IceCoverage
			};
			solarInJobHandles[_iceLayer] = solarInIceJob.Schedule(_cellCount, _batchCount, solarInJobHandles[_airLayer0]);
			
			for (int j = 0; j < _waterLayers; j++)
			{
				int layerIndex = _waterLayer0 + _waterLayers - 1 - j;
				var solarRadiationAbsorbedWaterJob = new SolarRadiationAbsorbedWaterJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[layerIndex],
					SolarRadiationIncoming = solarRadiation,
					SolarRadiationReflected = solarReflected[layerIndex],
					WaterCoverage = dependent.WaterCoverage[layerIndex - _waterLayer0],
					WaterSlopeAlbedo = waterSlopeAlbedo,
				};
				solarInJobHandles[layerIndex] = solarRadiationAbsorbedWaterJob.Schedule(_cellCount, _batchCount, solarInJobHandles[layerIndex + 1]);
			}

			var solarRadiationAbsorbedTerrainJob = new SolarRadiationAbsorbedTerrainJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_terrainLayer],
				SolarRadiationIncoming = solarRadiation,
				SolarRadiationReflected = solarReflected[_terrainLayer],
				VegetationCoverage = dependent.VegetationCoverage,
				worldData = worldData,
				LastTerrain = lastState.Terrain,
			};
			solarInJobHandles[_terrainLayer] = solarRadiationAbsorbedTerrainJob.Schedule(_cellCount, _batchCount, solarInJobHandles[_terrainLayer + 1]);
			
			// THERMAL RADIATION
			JobHandle[] thermalOutJobHandles = new JobHandle[_layerCount];

			// ICE
			var thermalOutJob = new ThermalEnergyRadiatedConstantEmissivityJob()
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
				SecondsPerTick = worldData.SecondsPerTick
			};
			thermalOutJobHandles[_iceLayer] = thermalOutJob.Schedule(_cellCount, 100, lastJobHandle);

			// CLOUD
			thermalOutJob = new ThermalEnergyRadiatedConstantEmissivityJob()
			{
				ThermalRadiationDelta = thermalRadiationDelta[_cloudLayer],
				ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[_cloudLayer],
				ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[_cloudLayer],
				WindowRadiationTransmittedUp = windowRadiationTransmittedUp[_cloudLayer],
				WindowRadiationTransmittedDown = windowRadiationTransmittedDown[_cloudLayer],

				PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
				Emissivity = WorldData.EmissivityWater,
				Energy = dependent.CloudEnergy,
				Temperature = lastState.CloudTemperature,
				SecondsPerTick = worldData.SecondsPerTick
			};
			thermalOutJobHandles[_cloudLayer] = thermalOutJob.Schedule(_cellCount, 100, lastJobHandle);


			// TERRAIN
			var thermalOutTerrainJob = new ThermalEnergyRadiatedTerrainJob()
			{
				ThermalRadiationDelta = thermalRadiationDelta[_terrainLayer],
				ThermalRadiationTransmitted = thermalRadiationTransmittedUp[_terrainLayer],
				WindowRadiationTransmitted = windowRadiationTransmittedUp[_terrainLayer],

				PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
				Emissivity = emissivity[_terrainLayer],
				Temperature = lastState.TerrainTemperature,
				SecondsPerTick = worldData.SecondsPerTick
			};
			thermalOutJobHandles[_terrainLayer] = thermalOutTerrainJob.Schedule(_cellCount, 100, emissivityJobHandles[_terrainLayer]);


			// ATMOSPHERE
			for (int j = 0; j < _airLayers; j++)
			{
				int layer = _airLayer0 + j;
				var thermalOutAirJob = new ThermalEnergyRadiatedJob()
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
				};
				thermalOutJobHandles[layer] = thermalOutAirJob.Schedule(_cellCount, 100, emissivityJobHandles[layer]);
			}

			// WATER
			for (int j = 0; j < _waterLayers; j++)
			{
				int layer = _waterLayer0 + j;
				var thermalOutWaterJob = new ThermalEnergyRadiatedJob()
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
				};
				thermalOutJobHandles[layer] = thermalOutWaterJob.Schedule(_cellCount, 100, emissivityJobHandles[layer]);
			}

			
			// THERMAL RADIATION ABSORPTION
			// Start at bottom water layer and go up, then go back down
			JobHandle[] thermalInUpJobHandles = new JobHandle[_layerCount];
			JobHandle[] thermalInDownJobHandles = new JobHandle[_layerCount];
			NativeArray<float> atmosphericWindowUp = new NativeArray<float>(_cellCount, Allocator.TempJob);
			NativeArray<float> atmosphericWindowDown = new NativeArray<float>(_cellCount, Allocator.TempJob);

			// transmit up from land
			for (int j = 1; j < _cloudLayer; j++)
			{

				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[j - 1]);
				if (j >= _airLayer0 && j < _airLayer0 + _airLayers)
				{
					int airLayer = j - _airLayer0;
					thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDependenciesHandle, thermalOutJobHandles[_cloudLayer]);
					var thermalInJob = new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationDeltaCloud = thermalRadiationDelta[_cloudLayer],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[j - 1],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[j - 1],
						ThermalRadiationTransmittedCloud = thermalRadiationTransmittedUp[_cloudLayer],
						CloudElevation = lastState.CloudElevation,
						CloudMass = lastState.CloudMass,
						LayerElevation = dependent.LayerElevation[airLayer],
						LayerHeight = dependent.LayerHeight[airLayer],
						AirMass = dependent.AirMass[airLayer],
						VaporMass = lastState.AirVapor[airLayer],
						AirAbsorptivity = worldData.AbsorptivityAir,
						VaporAbsorptivity = worldData.AbsorptivityWaterVapor,
						WaterAbsorptivity = worldData.AbsorptivityWaterLiquid,
						CarbonAbsorptivity = worldData.AbsorptivityCarbonDioxide,
						LayerIndex = airLayer,
						FromTop = false
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j == _iceLayer)
				{
					var thermalInJob = new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDeltaIceBottom,
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[j - 1],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[j - 1],
						Coverage = dependent.IceCoverage,

					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j >= _waterLayer0 && j < _waterLayer0 + _waterLayers)
				{
					thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDependenciesHandle, thermalOutJobHandles[j - 1]);
					var thermalInJob = new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDeltaSurfaceWater,
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[j - 1],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[j - 1],
						Coverage = dependent.WaterCoverage[j - _waterLayer0],
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else
				{
					thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDependenciesHandle, thermalOutJobHandles[j - 1]);
					var thermalInJob = new ThermalEnergyAbsorbedJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[j - 1],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[j - 1],
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
			}

			var thermalInUpJobHandlesCombined = default(JobHandle);
			for (int j=0;j<_layerCount;j++)
			{
				thermalInUpJobHandlesCombined = JobHandle.CombineDependencies(thermalInUpJobHandlesCombined, thermalInUpJobHandles[j]);
			}
			
			// transmit down from top of atmosphere
			for (int j = _airLayer0 + _airLayers - 2; j >= 0; j--)
			{
				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[j + 1], thermalInUpJobHandlesCombined);

				if (j == _terrainLayer)
				{
					var thermalInJob = new ThermalEnergyAbsorbedTerrainJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationDelta[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[j + 1],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j == _iceLayer)
				{
					var thermalInJob = new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDeltaIceTop,
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[j + 1],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
						Coverage = dependent.IceCoverage,
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j >= _waterLayer0 && j < _waterLayer0 + _waterLayers)
				{
					var thermalInJob = new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[j + 1],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
						Coverage = dependent.WaterCoverage[j-_waterLayer0],
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j >= _airLayer0 && j < _airLayer0 + _airLayers)
				{
					int airLayer = j - _airLayer0;
					thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDependenciesHandle, thermalOutJobHandles[_cloudLayer], thermalInUpJobHandles[_cloudLayer]);
					var thermalInJob = new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationDeltaCloud = thermalRadiationDelta[_cloudLayer],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[j + 1],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
						ThermalRadiationTransmittedCloud = thermalRadiationTransmittedDown[_cloudLayer],
						CloudElevation = lastState.CloudElevation,
						LayerElevation = dependent.LayerElevation[airLayer],
						LayerHeight = dependent.LayerHeight[airLayer],
						AirMass = dependent.AirMass[airLayer],
						CloudMass = lastState.CloudMass,
						VaporMass = lastState.AirVapor[airLayer],
						AirAbsorptivity = worldData.AbsorptivityAir,
						VaporAbsorptivity = worldData.AbsorptivityWaterVapor,
						WaterAbsorptivity = worldData.AbsorptivityWaterLiquid,
						CarbonAbsorptivity = worldData.AbsorptivityCarbonDioxide,
						LayerIndex = airLayer,
						FromTop = true
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else
				{
					var thermalInJob = new ThermalEnergyAbsorbedJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[j + 1],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
			}




			// DIFFUSION

			JobHandle[] diffusionJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				var diffusionJob = new DiffusionAirJob()
				{
					Delta = diffusionAir[j],
					Temperature = lastState.AirTemperature[j],
					Humidity = lastState.AirVapor[j],
					Velocity = lastState.AirVelocity[j],
					Neighbors = staticState.Neighbors,
					DiffusionCoefficient = worldData.AirMassDiffusionSpeedHorizontal,
				};
				diffusionJobHandles[_airLayer0 + j] = diffusionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			for (int j = 0; j < _waterLayers; j++)
			{
				var diffusionJob = new DiffusionWaterJob()
				{
					Delta = diffusionWater[j],
					Temperature = lastState.WaterTemperature[j],
					Salt = lastState.WaterSaltMass[j],
					Velocity = lastState.WaterVelocity[j],
					WaterMass = lastState.WaterMass[j],
					Neighbors = staticState.Neighbors,
					DiffusionCoefficient = worldData.WaterDiffuseSpeed,
				};
				diffusionJobHandles[_waterLayer0 + j] = diffusionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			var diffusionCloudJob = new DiffusionCloudJob()
			{
				Delta = diffusionCloud,
				Mass = lastState.CloudMass,
				Temperature = lastState.CloudTemperature,
				Elevation = lastState.CloudElevation,
				DropletSize = lastState.CloudDropletMass,
				Velocity = lastState.CloudVelocity,
				Neighbors = staticState.Neighbors,
				DiffusionCoefficient = worldData.CloudDiffusionCoefficient,
			};
			diffusionJobHandles[_cloudLayer] = diffusionCloudJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// ADVECTION

			JobHandle[] advectionJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				int layer = _airLayer0 + j;
				var advectionJob = new AdvectionAirJob()
				{
					Delta = advectionAir[j],
					Temperature = lastState.AirTemperature[j],
					Vapor = lastState.AirVapor[j],
					Velocity = lastState.AirVelocity[j],
					Neighbors = staticState.Neighbors,
					Coords = staticState.Coordinate,
					InverseCellDiameter = staticState.InverseCellDiameter,
					InverseCoordDiff = staticState.InverseCoordDiameter,
					SecondsPerTick = worldData.SecondsPerTick
				};
				advectionJobHandles[layer] = advectionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			for (int j = 0; j < _waterLayers; j++)
			{
				int layer = _waterLayer0 + j;
				var advectionJob = new AdvectionWaterJob()
				{
					Delta = advectionWater[j],
					Temperature = lastState.WaterTemperature[j],
					Salt = lastState.WaterSaltMass[j],
					Velocity = lastState.WaterVelocity[j],
					Neighbors = staticState.Neighbors,
				};
				advectionJobHandles[layer] = advectionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			
			// CONDUCTION

			// cloud to air
			//var conductionCloudAirJob = new ConductionJob()
			//{
			//	EnergyDelta = conductionCloudAir,
			//	TemperatureA = dependent.CloudTemperature,
			//	TemperatureB = dependent.AirTemperature[cloudAirLayer],
			//	EnergyA = lastState.CloudEnergy,
			//	EnergyB = lastState.AirEnergy[cloudAirLayer],
			//	ConductionCoefficient = worldData.AirWaterConductionPositive
			//};
			//var conductionCloudAirJobHandle = conductionCloudAirJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// air to ice
			var conductionAirIceJob = new ConductionAirIceJob()
			{
				EnergyDelta = conductionAirIce,
				TemperatureA = dependent.SurfaceAirTemperature,
				TemperatureB = lastState.IceTemperature,
				EnergyB = dependent.IceEnergy,
				ConductionCoefficient = WorldData.ConductivityAirIce,
				Coverage = dependent.IceCoverage,
				SecondsPerTick = worldData.SecondsPerTick
			};
			var conductionAirIceJobHandle = conductionAirIceJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// air to water
			int surfaceWaterLayer = _waterLayers - 1;
			var conductionAirWaterJob = new ConductionAirWaterJob()
			{
				EnergyDelta = conductionAirWater,
				TemperatureA = dependent.SurfaceAirTemperature,
				TemperatureB = lastState.WaterTemperature[surfaceWaterLayer],
				EnergyA = dependent.IceEnergy,
				EnergyB = dependent.WaterPotentialEnergy[surfaceWaterLayer],
				ConductionCoefficient = WorldData.ConductivityAirWater,
				CoverageIce = dependent.IceCoverage,
				CoverageWater = dependent.WaterCoverage[surfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			};
			var conductionAirWaterJobHandle = conductionAirWaterJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// air to terrain
			var conductionAirTerrainJob = new ConductionAirTerrainJob()
			{
				EnergyDelta = conductionAirTerrain,
				TemperatureA = dependent.SurfaceAirTemperature,
				TemperatureB = lastState.TerrainTemperature,
				ConductionCoefficient = WorldData.ConductivityAirTerrain,
				CoverageIce = dependent.IceCoverage,
				CoverageWater = dependent.WaterCoverage[surfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			};
			var conductionAirTerrainJobHandle = conductionAirTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// ice to water
			var conductionIceWaterJob = new ConductionIceWaterJob()
			{
				EnergyDelta = conductionIceWater,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.WaterTemperature[surfaceWaterLayer],
				EnergyA = dependent.IceEnergy,
				EnergyB = dependent.WaterPotentialEnergy[surfaceWaterLayer],
				ConductionCoefficient = WorldData.ConductivityIceWater,
				CoverageA = dependent.IceCoverage,
				CoverageB = dependent.WaterCoverage[surfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			};
			var conductionIceWaterJobHandle = conductionIceWaterJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// ice to terrain
			var conductionIceTerrainJob = new ConductionIceTerrainJob()
			{
				EnergyDelta = conductionIceTerrain,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.TerrainTemperature,
				EnergyA = dependent.IceEnergy,
				ConductionCoefficient = WorldData.ConductivityIceTerrain,
				CoverageIce = dependent.IceCoverage,
				CoverageWater = dependent.WaterCoverage[surfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			};
			var conductionIceTerrainJobHandle = conductionIceTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// water to terrain
			JobHandle conductionWaterTerrainJobHandle = default(JobHandle);
			{
				var conductionWaterTerrainJob = new ConductionWaterBottomJob()
				{
					EnergyDelta = conductionWaterTerrain[0],
					EnergyDeltaWaterTotal = conductionWaterTerrainTotal,
					TemperatureA = lastState.WaterTemperature[0],
					TemperatureB = lastState.TerrainTemperature,
					EnergyA = dependent.WaterPotentialEnergy[0],
					ConductionCoefficient = WorldData.ConductivityWaterTerrain,
					Coverage = dependent.WaterCoverage[0],
					SecondsPerTick = worldData.SecondsPerTick
				};
				conductionWaterTerrainJobHandle = JobHandle.CombineDependencies(conductionWaterTerrainJobHandle, conductionWaterTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle));
			}
			for (int i = 1; i < _waterLayers; i++) {
				var conductionWaterTerrainJob = new ConductionWaterTerrainJob()
				{
					EnergyDelta = conductionWaterTerrain[i],
					EnergyDeltaWaterTotal = conductionWaterTerrainTotal,
					TemperatureA = lastState.WaterTemperature[i],
					TemperatureB = lastState.TerrainTemperature,
					EnergyA = dependent.WaterPotentialEnergy[i],
					ConductionCoefficient = WorldData.ConductivityWaterTerrain,
					Coverage = dependent.WaterCoverage[i],
					CoverageBelow = dependent.WaterCoverage[i-1],
					SecondsPerTick = worldData.SecondsPerTick
				};
				conductionWaterTerrainJobHandle = JobHandle.CombineDependencies(conductionWaterTerrainJobHandle, conductionWaterTerrainJob.Schedule(_cellCount, _batchCount, conductionWaterTerrainJobHandle));
			}

			// COMBINE ADVECTION, DIFFUSION, SOLAR, THERMAL DELTA

			var energyJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
			var energyJobHandleDependencies = new List<NativeList<JobHandle>>();
			var energyTerrainJob = new EnergyTerrainJob()
			{
				Temperature = nextState.TerrainTemperature,
				LastTemperature = lastState.TerrainTemperature,
				Terrain = lastState.Terrain,
				SolarRadiationIn = solarRadiationIn[_terrainLayer],
				ThermalRadiationDelta = thermalRadiationDelta[_terrainLayer],
				ConductionEnergyAir = conductionAirTerrain,
				ConductionEnergyIce = conductionIceTerrain,
				ConductionEnergyWater = conductionWaterTerrainTotal,
				VegetationCoverage = dependent.VegetationCoverage,
				GeothermalEnergy = nextState.PlanetState.GeothermalHeat * worldData.SecondsPerTick,
				HeatingDepth = worldData.SoilHeatDepth
			};
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
			energyJobHandleDependencies.Add(terrainEnergyJobHandleDependencies);
			energyJobHandles.Add(energyTerrainJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(terrainEnergyJobHandleDependencies)));

			var energyIceJob = new EnergyIceJob()
			{
				Temperature = nextState.IceTemperature,
				Mass = nextState.IceMass,
				MeltedTopMass = iceMeltedTopMass,
				MeltedBottomMass = iceMeltedBottomMass,
				LastMass = lastState.IceMass,
				SolarRadiationIn = solarRadiationIn[_iceLayer],
				ThermalRadiationDeltaBottom = thermalRadiationDeltaIceBottom,
				ThermalRadiationDeltaTop = thermalRadiationDeltaIceTop,
				ConductionEnergyAir = conductionAirIce,
				ConductionEnergyTerrain = conductionIceTerrain,
				ConductionEnergyWater = conductionIceWater,
				LastTemperature = lastState.IceTemperature,
			};
			var iceEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_iceLayer],
				thermalOutJobHandles[_iceLayer],
				thermalInDownJobHandles[_iceLayer],
				thermalInUpJobHandles[_iceLayer],
				conductionAirIceJobHandle,
				conductionIceWaterJobHandle,
				conductionIceTerrainJobHandle,
			};
			energyJobHandleDependencies.Add(iceEnergyJobHandleDependencies);
			energyJobHandles.Add(energyIceJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(iceEnergyJobHandleDependencies)));

			var energyCloudJob = new EnergyCloudJob()
			{
				CloudTemperature = nextState.CloudTemperature,
				CloudMass = nextState.CloudMass,
				DropletMass = nextState.CloudDropletMass,
				Elevation = nextState.CloudElevation,
				Velocity = nextState.CloudVelocity,
				EvaporationMass = rainfallEvaporationMass,
				RainfallWaterMass = rainfallWaterMass,
				SolarRadiationIn = solarRadiationIn[_cloudLayer],
				ThermalRadiationDelta = thermalRadiationDelta[_cloudLayer],
				ConductionEnergyAir = conductionCloudAir,
				LastCloudMass = lastState.CloudMass,
				LastCloudTemperature = lastState.CloudTemperature,
				LastVelocity = lastState.CloudVelocity,
				Advection = advectionCloud,
				Diffusion = diffusionCloud,
				LastDropletMass = lastState.CloudDropletMass,
				LastElevation = lastState.CloudElevation,
				Gravity =lastState.PlanetState.Gravity,
				RainDropDragCoefficient = worldData.rainDropDragCoefficient,
				RainDropMaxSize =worldData.rainDropMaxSize,
				RainDropMinSize = worldData.rainDropMinSize,
				RainfallRate = worldData.RainfallRate,
				WindVertical = dependent.WindVertical[1], // TODO: make this be "wind at cloud elevation"
				LastAirMass = dependent.AirMass[1],
				LastWaterVapor = lastState.AirVapor[1],
				LastAirTemperature = lastState.AirTemperature[1],
				LastAirPressure = dependent.AirPressure[1],
				SecondsPerTick = worldData.SecondsPerTick,
				CloudDissapationRateDryAir = worldData.CloudDissapationRateDryAir,
				CloudDissapationRateWind = worldData.CloudDissapationRateWind,
				InverseCellDiameter = staticState.InverseCellDiameter,
				RelativeHumidityAtCloudElevation = dependent.AirHumidityRelative[1], // TODO: make a cloud elevation version of this
				TicksPerSecond = worldData.TicksPerSecond,
			};
			var cloudEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_cloudLayer],
				thermalOutJobHandles[_cloudLayer],
				thermalInDownJobHandles[_cloudLayer],
				thermalInUpJobHandles[_cloudLayer],
				diffusionJobHandles[_cloudLayer],
				//conductionCloudAirJobHandle,
			};
			for (int j=_airLayer0;j<_airLayer0+_airLayers;j++)
			{
				cloudEnergyJobHandleDependencies.Add(solarInJobHandles[j]);
				cloudEnergyJobHandleDependencies.Add(thermalInUpJobHandles[j]);
				cloudEnergyJobHandleDependencies.Add(thermalInDownJobHandles[j]);
			}
			energyJobHandleDependencies.Add(cloudEnergyJobHandleDependencies);
			energyJobHandles.Add(energyCloudJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(cloudEnergyJobHandleDependencies)));

			{
				var energyJob = new EnergySurfaceAirJob()
				{
					Temperature = nextState.AirTemperature[0],
					Vapor = nextState.AirVapor[0],
					Velocity = nextState.AirVelocity[0],
					CondensationCloudEnergy = condensationCloudEnergy[0],
					CondensationCloudMass = condensationCloudMass[0],
					CondensationGroundEnergy = condensationGroundEnergy[0],
					CondensationGroundMass = condensationGroundMass[0],
					LastTemperature = lastState.AirTemperature[0],
					LastVapor = lastState.AirVapor[0],
					LastVelocity = lastState.AirVelocity[0],
					AirMass = dependent.AirMass[0],
					Advection = advectionAir[0],
					Diffusion = diffusionAir[0],
					SolarRadiationIn = solarRadiationIn[_airLayer0],
					ThermalRadiationDelta = thermalRadiationDelta[_airLayer0],
					ConductionEnergyCloud = conductionCloudAir,
					ConductionEnergyWater = conductionAirWater,
					ConductionEnergyIce = conductionAirIce,
					ConductionEnergyTerrain = conductionAirTerrain,
					SecondsPerTick = worldData.SecondsPerTick,
				};

				var airDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					advectionJobHandles[_airLayer0],
					diffusionJobHandles[_airLayer0],
					solarInJobHandles[_airLayer0],
					thermalOutJobHandles[_airLayer0],
					thermalInDownJobHandles[_airLayer0],
					thermalInUpJobHandles[_airLayer0],
					//conductionCloudAirJobHandle,
					conductionAirWaterJobHandle,
					conductionAirIceJobHandle,
					conductionAirTerrainJobHandle
				};
				energyJobHandleDependencies.Add(airDependencies);
				energyJobHandles.Add(energyJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(airDependencies)));
			}
			for (int j = 1; j < _airLayers; j++)
			{
				int layerIndex = _airLayer0 + j;
				var energyJob = new EnergyAirJob()
				{
					Temperature = nextState.AirTemperature[j],
					Vapor = nextState.AirVapor[j],
					Velocity = nextState.AirVelocity[j],	
					CondensationCloudEnergy = condensationCloudEnergy[j],
					CondensationCloudMass = condensationCloudMass[j],
					CondensationGroundEnergy = condensationGroundEnergy[j],
					CondensationGroundMass = condensationGroundMass[j],
					LastTemperature = lastState.AirTemperature[j],
					LastVapor = lastState.AirVapor[j],
					LastVelocity = lastState.AirVelocity[j],
					AirMass = dependent.AirMass[j],
					Advection = advectionAir[j],
					Diffusion = diffusionAir[j],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationDelta = thermalRadiationDelta[layerIndex],
					ConductionEnergyCloud = conductionCloudAir,
					SecondsPerTick = worldData.SecondsPerTick,
				};

				var airDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					advectionJobHandles[layerIndex],
					diffusionJobHandles[layerIndex],
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
					//conductionCloudAirJobHandle,
				};
				energyJobHandleDependencies.Add(airDependencies);
				energyJobHandles.Add(energyJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(airDependencies)));
			}

			for (int j = 0; j < _waterLayers-1; j++)
			{
				int layerIndex = _waterLayer0 + j;
				var energyJob = new EnergyWaterJob()
				{
					Temperature = nextState.WaterTemperature[j],
					SaltMass = nextState.WaterSaltMass[j],
					Velocity = nextState.WaterVelocity[j],
					Mass = nextState.WaterMass[j],
					LastTemperature = lastState.WaterTemperature[j],
					LastSaltMass = lastState.WaterSaltMass[j],
					LastVelocity = lastState.WaterVelocity[j],
					LastMass = lastState.WaterMass[j],
					Advection = advectionWater[j],
					Diffusion = diffusionWater[j],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationDelta = thermalRadiationDelta[layerIndex],
					ConductionEnergyTerrain = conductionWaterTerrain[j],
					SecondsPerTick = worldData.SecondsPerTick,
				};

				var waterDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					advectionJobHandles[layerIndex],
					diffusionJobHandles[layerIndex],
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
					conductionAirWaterJobHandle,
					conductionIceWaterJobHandle,
					conductionWaterTerrainJobHandle,
				};
				energyJobHandleDependencies.Add(waterDependencies);
				energyJobHandles.Add(energyJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(waterDependencies)));
			}

			// surface water
			{
				int waterLayer = _waterLayers - 1;
				int layerIndex = _waterLayer0 + waterLayer;
				var energyJob = new EnergyWaterJobSurface()
				{
					Temperature = nextState.WaterTemperature[waterLayer],
					SaltMass = nextState.WaterSaltMass[waterLayer],
					Velocity = nextState.WaterVelocity[waterLayer],
					WaterMass = nextState.WaterMass[waterLayer],
					EvaporatedWaterMass = evaporationMass,
					FrozenBottomMass = frozenBottomMass,
					FrozenTopMass = frozenTopMass,

					LastMass = lastState.WaterMass[waterLayer],
					LastSaltMass = lastState.WaterSaltMass[waterLayer],
					LastVelocity = lastState.WaterVelocity[waterLayer],
					LastTemperature = lastState.WaterTemperature[waterLayer],
					LastSurfaceWind = lastState.AirVelocity[0],
					Advection = advectionWater[waterLayer],
					Diffusion = diffusionWater[waterLayer],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationDeltaTop = thermalRadiationDelta[layerIndex],
					ThermalRadiationDeltaBottom = thermalRadiationDeltaSurfaceWater,
					ConductionEnergyAir = conductionAirWater,
					ConductionEnergyIce = conductionIceWater,
					ConductionEnergyTerrain = conductionWaterTerrain[waterLayer],
					RelativeHumidity = dependent.AirHumidityRelative[0],
					IceCoverage = dependent.IceCoverage,
					WaterCoverage = dependent.WaterCoverage[waterLayer],
					EvaporationRate = worldData.EvaporationRate,
					EvapTemperatureMax = worldData.EvapMaxTemperature,
					EvapTemperatureMin = worldData.EvapMinTemperature,
					SecondsPerTick = worldData.SecondsPerTick,
				};

				var waterDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					advectionJobHandles[layerIndex],
					diffusionJobHandles[layerIndex],
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
					conductionAirWaterJobHandle,
					conductionIceWaterJobHandle,
					conductionWaterTerrainJobHandle,
				};
				energyJobHandleDependencies.Add(waterDependencies);
				energyJobHandles.Add(energyJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(waterDependencies)));
			}

			energyJobHandles.Add(updateTerrainJobHandle);
			var energyJobHandle = JobHandle.CombineDependencies(energyJobHandles);


			var stateChangeJob = new StateChangeJob()
			{
				IceTemperature = nextState.IceTemperature,
				IceMass = nextState.IceMass,
				SurfaceWaterTemperature = nextState.WaterTemperature[surfaceWaterLayer],
				SurfaceWaterMass = nextState.WaterMass[surfaceWaterLayer],
				SurfaceAirTemperature = nextState.AirTemperature[0],
				SurfaceAirVapor = nextState.AirVapor[0],

				RainfallEvaporatedMass = rainfallEvaporationMass,
				SurfaceAirMass = dependent.AirMass[0],
				SurfaceSaltMass = lastState.WaterSaltMass[surfaceWaterLayer],
				WaterEvaporatedMass = evaporationMass,
				WaterFrozenTopMass = frozenTopMass,
				WaterFrozenBottomMass = frozenBottomMass,
				IceMeltedTopMass = iceMeltedTopMass,
				IceMeltedBottomMass = iceMeltedBottomMass,
				RainfallTemperature = lastState.CloudTemperature,
				RainfallWaterMass = rainfallWaterMass,
			};
			var stateChangeJobHandle = stateChangeJob.Schedule(_cellCount, _batchCount, energyJobHandle);

			for (int j = 0; j < _airLayers; j++)
			{
				int layerIndex = _airLayer0 + j;
				var stateChangeAirLayerJob = new StateChangeAirLayerJob()
				{
					AirTemperature = nextState.AirTemperature[j],
					VaporMass = nextState.AirVapor[j],
					CloudDropletMass = nextState.CloudDropletMass,
					CloudElevation = nextState.CloudElevation,
					CloudEvaporationMass = cloudEvaporation,
					CloudMass = nextState.CloudMass,
					CloudTemperature = nextState.CloudTemperature,
					SurfaceWaterMass = nextState.WaterMass[surfaceWaterLayer],
					SurfaceWaterTemperature = nextState.WaterTemperature[surfaceWaterLayer],

					SurfaceSaltMass = nextState.WaterSaltMass[surfaceWaterLayer],
					CloudCondensationMass = cloudCondensation,
					GroundCondensationMass = condensationGroundMass[j],
					AirMass = dependent.AirMass[j],
					LayerIndex = layerIndex,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],

				};

				stateChangeJobHandle = stateChangeAirLayerJob.Schedule(_cellCount, _batchCount, stateChangeJobHandle);
			}


			NativeList<JobHandle> updateDependenciesJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
			for (int j = 0; j < _waterLayers; j++)
			{
				var updateDependentWaterLayerJob = new UpdateDependentWaterLayerJob()
				{
					Salinity = dependent.WaterSalinity[j],
					WaterCoverage = dependent.WaterCoverage[j],
					PotentialEnergy = dependent.WaterPotentialEnergy[j],

					Temperature = nextState.WaterTemperature[j],
					SaltMass = nextState.WaterSaltMass[j],
					WaterMass = nextState.WaterMass[j],
					Terrain = nextState.Terrain,
					worldData = worldData,
				};
				updateDependenciesJobHandles.Add(updateDependentWaterLayerJob.Schedule(_cellCount, _batchCount, stateChangeJobHandle));
			}
			JobHandle surfaceAirJobHandle = default(JobHandle);
			for (int j = 0; j < _airLayers; j++)
			{
				var updateDependentAirLayerJob = new UpdateDependentAirLayerJob()
				{
					Temperature = lastState.AirTemperature[j],
					Pressure = dependent.AirPressure[j],
					RelativeHumidity = dependent.AirHumidityRelative[j],
					AbsoluteHumidity = dependent.AirHumidityAbsolute[j],
					AirMass = dependent.AirMass[j],
					PotentialEnergy = dependent.AirPotentialEnergy[j],
					WindVertical = dependent.WindVertical[j],					

					VaporMass = nextState.AirVapor[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					IceMass = nextState.IceMass,
					Gravity = nextState.PlanetState.Gravity,
					DewPointZero = worldData.DewPointZero,
					InverseDewPointTemperatureRange = worldData.inverseDewPointTemperatureRange,
					WaterVaporMassToAirMassAtDewPoint = worldData.WaterVaporMassToAirMassAtDewPoint
				};
				var h = updateDependentAirLayerJob.Schedule(_cellCount, _batchCount, stateChangeJobHandle);
				updateDependenciesJobHandles.Add(h);
				if (j == 0)
				{
					surfaceAirJobHandle = h;
				}
			}


			JobHandle summationHandle = stateChangeJobHandle;
			var waterSaltMass = new NativeArray<WaterSaltMass>(_cellCount, Allocator.TempJob);
			for (int j = 0; j < _waterLayers; j++)
			{
				var updateWaterSaltMassJob = new UpdateWaterSaltMassJob()
				{
					WaterSaltMass = waterSaltMass,
					WaterLayerMass = nextState.WaterMass[j],
					SaltLayerMass = nextState.WaterSaltMass[j]
				};
				summationHandle = updateWaterSaltMassJob.Schedule(_cellCount, _batchCount, summationHandle);
			}
			var updateDependentStateJob = new UpdateDependentStateJob()
			{
				CloudCoverage = dependent.CloudCoverage,
				IceCoverage = dependent.IceCoverage,
				SurfaceElevation = dependent.SurfaceElevation,
				VegetationCoverage = dependent.VegetationCoverage,
				WaterDepth = dependent.WaterDepth,
				SurfaceAirTemperature = dependent.SurfaceAirTemperature,

				CloudMass = nextState.CloudMass,
				IceMass = nextState.IceMass,
				Terrain = nextState.Terrain,
				WaterSaltMass = waterSaltMass,
				worldData = worldData,
				LowerAirTemperature = lastState.AirTemperature[0],
				lowerAirHeight = dependent.LayerHeight[0],
			};

			var updateDependenciesJobHandle = updateDependentStateJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(surfaceAirJobHandle, summationHandle));
			updateDependenciesJobHandles.Add(updateDependenciesJobHandle);

			lastJobHandle = JobHandle.CombineDependencies(updateDependenciesJobHandles);
			lastJobHandle.Complete();


			if (tick == ticksToAdvance-1)
			{
				display.Dispose();
				display = new DisplayState();
				display.Init(_cellCount, _airLayers, _waterLayers);

				var updateDisplayJob = new UpdateDisplayJob()
				{
					SolarRadiationAbsorbedSurface = display.SolarRadiationAbsorbedSurface,
					DisplayEvaporation = display.Evaporation,
					DisplayRainfall = display.Rainfall,
					DisplayPressure = display.Pressure,

					SolarRadiationInTerrain = solarRadiationIn[_terrainLayer],
					SolarRadiationInIce = solarRadiationIn[_iceLayer],
					SolarRadiationInWaterSurface = solarRadiationIn[_waterLayer0 + _waterLayers - 1],
					Evaporation = evaporationMass,
					RainfallWater = rainfallWaterMass,
					AirMass = dependent.AirMass[0],
					Gravity = nextState.PlanetState.Gravity,
					AirLayerElevation = dependent.LayerElevation[0],
					AirLayerHeight = dependent.LayerHeight[0],
					AirPressure = dependent.AirPressure[0],
					AirTemperature = nextState.AirTemperature[0],
					WaterVaporMass = nextState.AirVapor[0],

				};
				var updateDisplayJobHandle = updateDisplayJob.Schedule(_cellCount, _batchCount);
				updateDisplayJobHandle.Complete();

				var curState = states[curStateIndex];
				for (int i = 0; i < _cellCount; i++)
				{
					display.SolarRadiation += displaySolarRadiation[i];
					display.GlobalCloudCoverage += dependent.CloudCoverage[i];
					display.GlobalCloudMass += curState.CloudMass[i];
					display.GlobalIceMass += curState.IceMass[i];
					display.GlobalOceanCoverage += dependent.WaterCoverage[_waterLayers-1][i];
					display.GlobalTemperature += curState.AirTemperature[0][i];
					display.GlobalWaterVapor += curState.AirVapor[0][i];
					display.GlobalOceanVolume += dependent.WaterDepth[i];
					display.GlobalSeaLevel += dependent.SurfaceElevation[i];
					display.GlobalEvaporation += display.Evaporation[i];
					display.GlobalRainfall += display.Rainfall[i];
					for (int j = 0; j < _airLayers; j++)
					{
						display.EnergySolarReflectedAtmosphere += solarReflected[j + _airLayer0][i];
						display.EnergySolarAbsorbedAtmosphere += solarRadiationIn[j + _airLayer0][i];
					}
					for (int j = 0; j < _waterLayers; j++)
					{
						float absorbed = solarRadiationIn[j + _waterLayer0][i];
						display.EnergySolarAbsorbedOcean += absorbed;
						display.EnergySolarAbsorbedSurface += absorbed;
						display.EnergySolarReflectedSurface += solarReflected[j + _waterLayer0][i];
					}
					display.EnergySolarAbsorbedCloud += solarRadiationIn[_cloudLayer][i];
					display.EnergySolarAbsorbedSurface += solarRadiationIn[_terrainLayer][i] + solarRadiationIn[_iceLayer][i];
					display.EnergySolarReflectedCloud += solarReflected[_cloudLayer][i];
					display.EnergySolarReflectedSurface += solarReflected[_terrainLayer][i] + solarReflected[_iceLayer][i];
					display.EnergySurfaceConduction += conductionAirIce[i] + conductionAirTerrain[i] + conductionAirWater[i];
					display.EnergyOceanConduction += conductionAirWater[i];
					display.EnergyEvapotranspiration += evaporationMass[i] * WorldData.LatentHeatWaterVapor;
					//display.EnergyThermalAbsorbedAtmosphere += ;
					display.EnergyThermalBackRadiation += windowRadiationTransmittedDown[_airLayer0][i] + thermalRadiationTransmittedDown[_airLayer0][i];
					display.EnergyThermalOceanRadiation += (windowRadiationTransmittedUp[_waterLayer0 + _waterLayers - 1][i] + thermalRadiationTransmittedUp[_waterLayer0 + _waterLayers - 1][i]) * dependent.WaterCoverage[_waterLayers - 1][i];
					display.EnergyThermalOutAtmosphere += thermalRadiationTransmittedUp[_airLayer0 + _airLayers - 1][i];
					display.EnergyThermalSurfaceOutAtmosphericWindow += windowRadiationTransmittedUp[_iceLayer][i];
					display.EnergyThermalSurfaceRadiation += windowRadiationTransmittedUp[_iceLayer][i] + thermalRadiationTransmittedUp[_iceLayer][i];
			
				}

			}

			updateDependenciesJobHandles.Dispose();
			waterSaltMass.Dispose();
			energyJobHandles.Dispose();
			atmosphericWindowUp.Dispose();
			atmosphericWindowDown.Dispose();
			evaporationMass.Dispose();
			evaporationLatentHeat.Dispose();
			frozenTopMass.Dispose();
			frozenBottomMass.Dispose();
			rainfallWaterMass.Dispose();
			rainfallEvaporationMass.Dispose();
			iceMeltedTopMass.Dispose();
			iceMeltedBottomMass.Dispose();
			foreach (var d in energyJobHandleDependencies)
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
				condensationGroundEnergy[i].Dispose();
				condensationCloudMass[i].Dispose();
				condensationCloudEnergy[i].Dispose();
			}



		}


	}
}
