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
	private NativeArray<float> thermalRadiationDeltaIceBottom;
	private NativeArray<float> thermalRadiationDeltaIceTop;
	private NativeArray<float> thermalRadiationDeltaSurfaceWater;
	private NativeArray<float>[] thermalRadiationDelta;
	private NativeArray<float>[] thermalRadiationTransmittedUp;
	private NativeArray<float>[] thermalRadiationTransmittedDown;
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
	private NativeArray<float> conductionWaterTerrain;
	private NativeArray<float> rainfall;
	private NativeArray<float> iceMelt;
	private NativeArray<float> waterFrozen;
	private NativeArray<float> cloudEvaporation;
	private NativeArray<float> surfaceEvaporation;
	private NativeArray<float> cloudCondensation;
	private NativeArray<float> surfaceCondensation;


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
		thermalRadiationDeltaIceBottom = new NativeArray<float>(_cellCount, Allocator.Persistent);
		thermalRadiationDeltaIceTop = new NativeArray<float>(_cellCount, Allocator.Persistent);
		thermalRadiationDeltaSurfaceWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		emissivity = new NativeArray<float>[_layerCount];
		solarRadiationIn = new NativeArray<float>[_layerCount];
		thermalRadiationDelta = new NativeArray<float>[_layerCount];
		thermalRadiationTransmittedUp = new NativeArray<float>[_layerCount];
		thermalRadiationTransmittedDown = new NativeArray<float>[_layerCount];
		latentHeat = new NativeArray<float>[_layerCount];
		for (int i=0;i<_layerCount;i++)
		{
			latentHeat[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			emissivity[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			solarRadiationIn[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			thermalRadiationDelta[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			thermalRadiationTransmittedUp[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			thermalRadiationTransmittedDown[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
		}
		diffusionAir = new NativeArray<DiffusionAir>[_airLayers];
		advectionAir = new NativeArray<DiffusionAir>[_airLayers];
		for (int i = 0; i < _airLayers; i++)
		{
			diffusionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
			advectionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
		}
		diffusionWater = new NativeArray<DiffusionWater>[_waterLayers];
		advectionWater = new NativeArray<DiffusionWater>[_waterLayers];
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
			advectionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
		}
		diffusionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		advectionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		conductionCloudAir = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirIce = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionIceWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionIceTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionWaterTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);

	}

	public void Dispose()
	{
		solarRadiation.Dispose();
		waterSlopeAlbedo.Dispose();
		thermalRadiationDeltaIceTop.Dispose();
		thermalRadiationDeltaIceBottom.Dispose();
		thermalRadiationDeltaSurfaceWater.Dispose();
		for (int i=0;i<_layerCount;i++)
		{
			latentHeat[i].Dispose();
			emissivity[i].Dispose();
			solarRadiationIn[i].Dispose();
			thermalRadiationDelta[i].Dispose();
			thermalRadiationTransmittedUp[i].Dispose();
			thermalRadiationTransmittedDown[i].Dispose();
		}
		for (int i = 0; i < _airLayers; i++)
		{
			diffusionAir[i].Dispose();
			advectionAir[i].Dispose();
		}
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i].Dispose();
			advectionWater[i].Dispose();
		}
		diffusionCloud.Dispose();
		advectionCloud.Dispose();
		conductionCloudAir.Dispose();
		conductionAirIce.Dispose();
		conductionAirWater.Dispose();
		conductionAirTerrain.Dispose();
		conductionIceWater.Dispose();
		conductionIceTerrain.Dispose();
		conductionWaterTerrain.Dispose();
	}

	public void Tick(SimState[] states, int stateCount, int ticksToAdvance, ref DependentState dependent, ref StaticState staticState, ref WorldData worldData, ref int curStateIndex)
	{
		JobHandle lastJobHandle = default(JobHandle);

		for (int tick = 0; tick < ticksToAdvance; tick++)
		{
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
				WaterSlopeAlbedo = waterSlopeAlbedo,

				SphericalPosition = staticState.SphericalPosition,
				IncomingSolarRadiation = lastState.PlanetState.SolarRadiation,
				PlanetRotation = quaternion.Euler(lastState.PlanetState.Rotation),
				SecondsPerTick = worldData.SecondsPerTick,
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
				var solarInAtmosphereJob = new SolarRadiationAbsorbedAirJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[layerIndex],
					SolarRadiationIncoming = solarRadiation,
					SolarRadiationAbsorbedCloud = solarRadiationIn[_cloudLayer],
					AirMass = dependent.AirMass[j],
					VaporMass = lastState.AirVapor[j],
					CloudMass = lastState.CloudMass,
					CloudTemperature = lastState.CloudEnergy,
					CloudDropletMass = lastState.CloudDropletMass,
					CloudElevation = lastState.CloudElevation,
					CloudCoverage = dependent.CloudCoverage,
					WaterSlopeAlbedo = solarRadiationJob.WaterSlopeAlbedo,
					SolarReflectivityAir = worldData.SolarReflectivityAir,
					SolarAbsorptivityAir = worldData.SolarAbsorptivityAir,
					SolarAbsorptivityWaterVapor = worldData.SolarAbsorptivityWaterVapor,
					SolarAbsorptivityCloud = worldData.SolarAbsorptivityCloud,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					worldData = worldData
				};
				solarInJobHandles[layerIndex] = solarInAtmosphereJob.Schedule(_cellCount, _batchCount, j == 0 ? solarInJobHandle : solarInJobHandles[layerIndex + 1]);
			}

			
			// ice
			var solarInIceJob = new SolarRadiationAbsorbedIceJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_iceLayer],
				SolarRadiationIncoming = solarRadiation,
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
					WaterCoverage = dependent.WaterCoverage,
					WaterSlopeAlbedo = waterSlopeAlbedo,
				};
				solarInJobHandles[layerIndex] = solarRadiationAbsorbedWaterJob.Schedule(_cellCount, _batchCount, solarInJobHandles[layerIndex + 1]);
			}

			var solarRadiationAbsorbedTerrainJob = new SolarRadiationAbsorbedTerrainJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_terrainLayer],
				SolarRadiationIncoming = solarRadiation,
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
				Emissivity = WorldData.EmissivityIce,
				Temperature = lastState.IceEnergy, 
			};
			thermalOutJobHandles[_iceLayer] = thermalOutJob.Schedule(_cellCount, 100, lastJobHandle);

			// CLOUD
			thermalOutJob = new ThermalEnergyRadiatedConstantEmissivityJob()
			{
				ThermalRadiationDelta = thermalRadiationDelta[_cloudLayer],
				ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[_cloudLayer],
				ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[_cloudLayer],
				Emissivity = WorldData.EmissivityWater,
				Temperature = lastState.CloudEnergy,
			};
			thermalOutJobHandles[_cloudLayer] = thermalOutJob.Schedule(_cellCount, 100, lastJobHandle);


			// TERRAIN
			var thermalOutTerrainJob = new ThermalEnergyRadiatedTerrainJob()
			{
				ThermalRadiationDelta = thermalRadiationDelta[_terrainLayer],
				ThermalRadiationTransmitted = thermalRadiationTransmittedUp[_terrainLayer],
				Emissivity = emissivity[_terrainLayer],
				Temperature = lastState.TerrainEnergy,
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
					Emissivity = emissivity[layer],
					Temperature = lastState.AirEnergy[j]			
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
					Emissivity = emissivity[layer],
					Temperature = lastState.WaterEnergy[j],
				};
				thermalOutJobHandles[layer] = thermalOutWaterJob.Schedule(_cellCount, 100, emissivityJobHandles[layer]);
			}

			
			// THERMAL RADIATION ABSORPTION
			// Start at bottom water layer and go up, then go back down
			JobHandle[] thermalInUpJobHandles = new JobHandle[_layerCount];
			JobHandle[] thermalInDownJobHandles = new JobHandle[_layerCount];

			// transmit up from land
			for (int j = 1; j < _cloudLayer; j++)
			{

				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[j - 1]);
				if (j >= _airLayer0 && j < _airLayer0 + _airLayers)
				{
					thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDependenciesHandle, thermalOutJobHandles[_cloudLayer]);
					var thermalInJob = new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationDeltaCloud = thermalRadiationDelta[_cloudLayer],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[j - 1],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						ThermalRadiationEmittedCloudUp = thermalRadiationTransmittedUp[_cloudLayer],
						ThermalRadiationEmittedCloudDown = thermalRadiationTransmittedDown[_cloudLayer],
						Emissivity = emissivity[j],
						CloudElevation = lastState.CloudElevation,
						CloudCoverage = dependent.CloudCoverage,
						LayerElevation = dependent.LayerElevation[j - _airLayer0],
						LayerHeight = dependent.LayerHeight[j - _airLayer0],
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j == _iceLayer)
				{
					var thermalInJob = new ThermalEnergyAbsorbedPartialCoverageConstantEmissivityJob()
					{
						ThermalRadiationDelta = thermalRadiationDeltaIceBottom,
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[j - 1],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						Emissivity = WorldData.EmissivityIce,
						Coverage = dependent.IceCoverage
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j == _waterLayer0 + _waterLayers - 1)
				{
					var thermalInJob = new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[j - 1],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						Emissivity = emissivity[j],
						Coverage = dependent.WaterCoverage
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else
				{
					thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDependenciesHandle, thermalOutJobHandles[j - 1]);
					var thermalInJob = new ThermalEnergyAbsorbedJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[j - 1],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						Emissivity = emissivity[j]
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
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j == _iceLayer)
				{
					var thermalInJob = new ThermalEnergyAbsorbedPartialCoverageConstantEmissivityJob()
					{
						ThermalRadiationDelta = thermalRadiationDeltaIceTop,
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
						Emissivity = WorldData.EmissivityIce,
						Coverage = dependent.IceCoverage
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j >= _airLayer0 && j < _airLayer0 + _airLayers)
				{
					thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDependenciesHandle, thermalOutJobHandles[_cloudLayer], thermalInUpJobHandles[_cloudLayer]);
					var thermalInJob = new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationDeltaCloud = thermalRadiationDelta[_cloudLayer],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						ThermalRadiationEmittedCloudUp = thermalRadiationTransmittedUp[_cloudLayer],
						ThermalRadiationEmittedCloudDown = thermalRadiationTransmittedDown[_cloudLayer],
						Emissivity = emissivity[j],
						CloudElevation = lastState.CloudElevation,
						CloudCoverage = dependent.CloudCoverage,
						LayerElevation = dependent.LayerElevation[j-_airLayer0],
						LayerHeight = dependent.LayerHeight[j-_airLayer0],
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else
				{
					var thermalInJob = new ThermalEnergyAbsorbedJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						Emissivity = emissivity[j]
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
					Temperature = dependent.AirTemperature[j],
					Humidity = lastState.AirVapor[j],
					Velocity = lastState.AirVelocity[j],
					Neighbors = staticState.Neighbors,
					DiffusionCoefficientHumidity = worldData.AirMassDiffusionSpeedHorizontal,
					DiffusionCoefficientTemperature = worldData.AirMassDiffusionSpeedHorizontal,
					DiffusionCoefficientVelocity = worldData.AirMassDiffusionSpeedHorizontal,
				};
				diffusionJobHandles[_airLayer0 + j] = diffusionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			for (int j = 0; j < _waterLayer0; j++)
			{
				var diffusionJob = new DiffusionWaterJob()
				{
					Delta = diffusionWater[j],
					Temperature = dependent.WaterTemperature[j],
					Salt = lastState.WaterSaltMass[j],
					Velocity = lastState.WaterVelocity[j],
					Neighbors = staticState.Neighbors,
					DiffusionCoefficientSalinity = worldData.WaterDiffuseSpeed,
					DiffusionCoefficientTemperature = worldData.WaterDiffuseSpeed,
					DiffusionCoefficientVelocity = worldData.WaterDiffuseSpeed
				};
				diffusionJobHandles[_waterLayer0 + j] = diffusionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			
			// ADVECTION

			JobHandle[] advectionJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				int layer = _airLayer0 + j;
				var advectionJob = new AdvectionAirJob()
				{
					Delta = advectionAir[j],
					Temperature = dependent.AirTemperature[j],
					Vapor = lastState.AirVapor[j],
					Velocity = lastState.AirVelocity[j],
					Neighbors = staticState.Neighbors,
				};
				advectionJobHandles[layer] = advectionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			for (int j = 0; j < _waterLayers; j++)
			{
				int layer = _waterLayer0 + j;
				var advectionJob = new AdvectionWaterJob()
				{
					Delta = advectionWater[j],
					Energy = lastState.WaterEnergy[j],
					Salt = lastState.WaterSaltMass[j],
					Velocity = lastState.WaterVelocity[j],
					Neighbors = staticState.Neighbors,
				};
				advectionJobHandles[layer] = advectionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			
			// CONDUCTION

			// cloud to air
			int cloudAirLayer = 1;
			var conductionCloudAirJob = new ConductionJob()
			{
				EnergyDelta = conductionCloudAir,
				MassA = lastState.CloudMass,
				MassB = dependent.AirMass[cloudAirLayer],
				EnergyA = lastState.CloudEnergy,
				EnergyB = lastState.AirEnergy[cloudAirLayer],
				ConductionCoefficient = worldData.AirWaterConductionPositive
			};
			var conductionCloudAirJobHandle = conductionCloudAirJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// air to ice
			var conductionAirIceJob = new ConductionAirIceJob()
			{
				EnergyDelta = conductionAirIce,
				MassA = dependent.AirMass[0],
				MassB = lastState.IceMass,
				EnergyA = lastState.AirEnergy[0],
				EnergyB = lastState.IceEnergy,
				ConductionCoefficient = worldData.AirIceConduction,
				IceCoverage = dependent.IceCoverage
			};
			var conductionAirIceJobHandle = conductionAirIceJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// air to water
			int surfaceWaterLayer = _waterLayers - 1;
			var conductionAirWaterJob = new ConductionAirWaterJob()
			{
				EnergyDelta = conductionAirWater,
				MassA = dependent.AirMass[0],
				MassB = lastState.WaterMass[surfaceWaterLayer],
				EnergyA = lastState.IceEnergy,
				EnergyB = lastState.WaterEnergy[surfaceWaterLayer],
				ConductionCoefficientPositive = worldData.AirWaterConductionPositive,
				ConductionCoefficientNegative = worldData.AirWaterConductionNegative,
				IceCoverage = dependent.IceCoverage,
				WaterCoverage = dependent.WaterCoverage
				
			};
			var conductionAirWaterJobHandle = conductionAirWaterJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			//// air to terrain
			//var conductionAirTerrainJob = new ConductionAirTerrainJob()
			//{
			//	EnergyDelta = airTerrainConduction,
			//	MassA = dependent.AirMass[0],
			//	MassB = ,
			//	TemperatureA = lastState.AirTemperature[0],
			//	TemperatureB = lastState.TerrainTemperature,
			//	ConductionCoefficient = worldData.AirTerrainConduction,
			//	IceCoverage = dependent.IceCoverage,
			//	WaterCoverage = dependent.WaterCoverage
			//};
			//var conductionAirTerrainJobHandle = conductionAirTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// ice to water
			var conductionIceWaterJob = new ConductionIceWaterJob()
			{
				EnergyDelta = conductionIceWater,
				MassA = lastState.IceMass,
				MassB = lastState.WaterMass[surfaceWaterLayer],
				EnergyA = lastState.IceEnergy,
				EnergyB = lastState.WaterEnergy[surfaceWaterLayer],
				ConductionCoefficient = worldData.IceWaterConduction,
				IceCoverage = dependent.IceCoverage,
				WaterCoverage = dependent.WaterCoverage
			};
			var conductionIceWaterJobHandle = conductionIceWaterJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			//// ice to terrain
			//var conductionIceTerrainJob = new ConductionIceTerrainJob()
			//{
			//	EnergyDelta = iceTerrainConduction,
			//	MassA = lastState.IceMass,
			//	MassB = ,
			//	TemperatureA = lastState.IceTemperature,
			//	TemperatureB = lastState.TerrainTemperature,
			//	ConductionCoefficient = worldData.IceTerrainConduction,
			//	IceCoverage = dependent.IceCoverage,
			//	WaterCoverage = dependent.WaterCoverage
			//};
			//var conductionIceTerrainJobHandle = conductionIceToTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			//// water to terrain
			//var conductionWaterTerrainJob = new ConductionWaterTerrainJob()
			//{
			//	EnergyDelta = waterTerrainConduction,
			//	MassA = lastState.WaterMass[0],
			//	MassB = ,
			//	TemperatureA = lastState.WaterTemperature[0],
			//	TemperatureB = lastState.TerrainTemperature,
			//	ConductionCoefficient = worldData.WaterTerrainConduction,
			//	WaterCoverage = dependent.WaterCoverage
			//};
			//var conductionWaterTerrainJobHandle = conductionWaterToTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle);


			// COMBINE ADVECTION, DIFFUSION, SOLAR, THERMAL DELTA

			var energyJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
			var energyJobHandleDependencies = new List<NativeList<JobHandle>>();
			var energyTerrainJob = new EnergyTerrainJob()
			{
				Energy = nextState.TerrainEnergy,
				Terrain = lastState.Terrain,
				SolarRadiationIn = solarRadiationIn[_terrainLayer],
				ThermalRadiationDelta = thermalRadiationDelta[_terrainLayer],
				ConductionEnergyAir = conductionAirTerrain,
				ConductionEnergyIce = conductionIceTerrain,
				ConductionEnergyWater = conductionWaterTerrain
			};
			var terrainEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_terrainLayer],
				thermalOutJobHandles[_terrainLayer],
				thermalInDownJobHandles[_terrainLayer],
				thermalInUpJobHandles[_terrainLayer],
				//conductionAirTerrainJobHandle,
				//conductionIceTerrainJobHandle,
				//conductionWaterTerrainHandle
			};
			energyJobHandleDependencies.Add(terrainEnergyJobHandleDependencies);
			energyJobHandles.Add(energyTerrainJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(terrainEnergyJobHandleDependencies)));

			var energyIceJob = new EnergyIceJob()
			{
				Energy = nextState.IceEnergy,
				SolarRadiationIn = solarRadiationIn[_iceLayer],
				ThermalRadiationDeltaBottom = thermalRadiationDeltaIceBottom,
				ThermalRadiationDeltaTop = thermalRadiationDeltaIceTop,
				ConductionEnergyAir = conductionAirIce,
				ConductionEnergyTerrain = conductionIceTerrain,
				ConductionEnergyWater = conductionIceWater,
				LastEnergy = lastState.IceEnergy,
			};
			var iceEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_iceLayer],
				thermalOutJobHandles[_iceLayer],
				thermalInDownJobHandles[_iceLayer],
				thermalInUpJobHandles[_iceLayer],
				conductionAirIceJobHandle,
				conductionIceWaterJobHandle,
				//conductionIceTerrainJobHandle,
			};
			energyJobHandleDependencies.Add(iceEnergyJobHandleDependencies);
			energyJobHandles.Add(energyIceJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(iceEnergyJobHandleDependencies)));

			var energyCloudJob = new EnergyCloudJob()
			{
				Energy = nextState.CloudEnergy,
				Mass = nextState.CloudMass,
				DropletMass = nextState.CloudDropletMass,
				Elevation = nextState.CloudElevation,
				Velocity = nextState.CloudVelocity,
				SolarRadiationIn = solarRadiationIn[_cloudLayer],
				ThermalRadiationDelta = thermalRadiationDelta[_cloudLayer],
				ConductionEnergyAir = conductionCloudAir,
				LastMass = lastState.CloudMass,
				LastEnergy = lastState.CloudEnergy,
				LastVelocity = lastState.CloudVelocity,
				Advection = advectionCloud,
				Diffusion = diffusionCloud,
				LastDropletMass = lastState.CloudDropletMass,
				LastElevation = lastState.CloudElevation,
			};
			var cloudEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_cloudLayer],
				thermalOutJobHandles[_cloudLayer],
				thermalInDownJobHandles[_cloudLayer],
				thermalInUpJobHandles[_cloudLayer],
				conductionCloudAirJobHandle,
			};
			for (int j=_airLayer0;j<_airLayer0+_airLayers;j++)
			{
				cloudEnergyJobHandleDependencies.Add(solarInJobHandles[j]);
				cloudEnergyJobHandleDependencies.Add(thermalInUpJobHandles[j]);
				cloudEnergyJobHandleDependencies.Add(thermalInDownJobHandles[j]);
			}
			energyJobHandleDependencies.Add(cloudEnergyJobHandleDependencies);
			energyJobHandles.Add(energyCloudJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(cloudEnergyJobHandleDependencies)));

			for (int j = 0; j < _airLayers; j++)
			{
				int layerIndex = _airLayer0 + j;
				var energyJob = new EnergyAirJob()
				{
					Energy = nextState.AirEnergy[j],
					Vapor = nextState.AirVapor[j],
					Velocity = nextState.AirVelocity[j],
					LastEnergy = lastState.AirEnergy[j],
					LastVapor = lastState.AirVapor[j],
					LastVelocity = lastState.AirVelocity[j],
					AirMass = dependent.AirMass[j],
					Advection = advectionAir[j],
					Diffusion = diffusionAir[j],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationDelta = thermalRadiationDelta[layerIndex],
					ConductionEnergyCloud = conductionCloudAir,
					ConductionEnergyWater = conductionAirWater,
					ConductionEnergyIce = conductionAirIce,
					ConductionEnergyTerrain = conductionAirTerrain
				};

				var airDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					advectionJobHandles[layerIndex],
					diffusionJobHandles[layerIndex],
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
					conductionCloudAirJobHandle,
					conductionAirWaterJobHandle,
					conductionAirIceJobHandle, 
					//conductionAirTerrainJobHandle;
				};
				energyJobHandleDependencies.Add(airDependencies);
				energyJobHandles.Add(energyJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(airDependencies)));
			}

			for (int j = 0; j < _waterLayers-1; j++)
			{
				int layerIndex = _waterLayer0 + j;
				var energyJob = new EnergyWaterJob()
				{
					Energy = nextState.WaterEnergy[j],
					Salinity = nextState.WaterSaltMass[j],
					Velocity = nextState.WaterVelocity[j],
					Mass = nextState.WaterMass[j],
					LastTemperature = lastState.WaterEnergy[j],
					LastSalinity = lastState.WaterSaltMass[j],
					LastVelocity = lastState.WaterVelocity[j],
					LastMass = lastState.WaterMass[j],
					Advection = advectionWater[j],
					Diffusion = diffusionWater[j],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationDelta = thermalRadiationDelta[layerIndex],
					ConductionEnergyTerrain = conductionAirTerrain
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
					//conductionWaterTerrainJobHandle,
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
					Energy = nextState.WaterEnergy[waterLayer],
					Salinity = nextState.WaterSaltMass[waterLayer],
					Velocity = nextState.WaterVelocity[waterLayer],
					Mass = nextState.WaterMass[waterLayer],
					LastEnergy = lastState.WaterEnergy[waterLayer],
					LastSalinity = lastState.WaterSaltMass[waterLayer],
					LastVelocity = lastState.WaterVelocity[waterLayer],
					Advection = advectionWater[waterLayer],
					Diffusion = diffusionWater[waterLayer],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationDeltaTop = thermalRadiationDelta[layerIndex],
					ThermalRadiationDeltaBottom = thermalRadiationDeltaSurfaceWater,
					ConductionEnergyAir = conductionAirWater,
					ConductionEnergyIce = conductionIceWater,
					ConductionEnergyTerrain = conductionAirTerrain
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
					//conductionWaterTerrainJobHandle,
				};
				energyJobHandleDependencies.Add(waterDependencies);
				energyJobHandles.Add(energyJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(waterDependencies)));
			}

			energyJobHandles.Add(updateTerrainJobHandle);
			var energyJobHandle = JobHandle.CombineDependencies(energyJobHandles);

			NativeList<JobHandle> updateDependenciesJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
			for (int j = 0; j < _waterLayers; j++)
			{
				var updateDependentWaterLayerJob = new UpdateDependentWaterLayerJob()
				{
					Temperature = dependent.WaterTemperature[j],
					Salinity = dependent.WaterSalinity[j],
					Energy = nextState.WaterEnergy[j],
					SaltMass = nextState.WaterSaltMass[j],
					WaterMass = nextState.WaterMass[j],
					worldData = worldData
				};
				updateDependenciesJobHandles.Add(updateDependentWaterLayerJob.Schedule(_cellCount, _batchCount, energyJobHandle));
			}
			for (int j = 0; j < _airLayers; j++)
			{
				var updateDependentAirLayerJob = new UpdateDependentAirLayerJob()
				{
					Temperature = dependent.AirTemperature[j],
					Pressure = dependent.AirPressure[j],
					RelativeHumidity = dependent.AirHumidityRelative[j],
					AbsoluteHumidity = dependent.AirHumidityAbsolute[j],
					AirMass = dependent.AirMass[j],
					AirEnergy = nextState.AirEnergy[j],
					AirVapor = nextState.AirVapor[j],
					VaporMass = nextState.AirVapor[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					IceMass = nextState.IceMass,
					Gravity = nextState.PlanetState.Gravity,
					WorldData = worldData,
				};
				updateDependenciesJobHandles.Add(updateDependentAirLayerJob.Schedule(_cellCount, _batchCount, energyJobHandle));
			}


			JobHandle summationHandle = energyJobHandle;
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
				WaterCoverage = dependent.WaterCoverage,
				WaterDepth = dependent.WaterDepth,

				CloudMass = nextState.CloudMass,
				IceMass = nextState.IceMass,
				Terrain = nextState.Terrain,
				WaterSaltMass = waterSaltMass,
				worldData = worldData
			};
			updateDependenciesJobHandles.Add(updateDependentStateJob.Schedule(_cellCount, _batchCount, summationHandle));
			lastJobHandle = JobHandle.CombineDependencies(updateDependenciesJobHandles);
			lastJobHandle.Complete();

			updateDependenciesJobHandles.Dispose();
			waterSaltMass.Dispose();
			energyJobHandles.Dispose();
			foreach (var d in energyJobHandleDependencies)
			{
				d.Dispose();
			}

		}

		int outputBuffer = ticksToAdvance % 2;
		dependent.DisplayPlanet = new PlanetDisplay();
		var curState = states[curStateIndex];
		for (int i = 0; i < _cellCount; i++)
		{

			dependent.DisplayPlanet.CloudCoverage += dependent.CloudCoverage[i];
			dependent.DisplayPlanet.CloudMass += curState.CloudMass[i];
			dependent.DisplayPlanet.IceMass += curState.IceMass[i];
			dependent.DisplayPlanet.OceanCoverage += dependent.WaterCoverage[i];
			dependent.DisplayPlanet.Temperature += curState.AirEnergy[0][i];
			dependent.DisplayPlanet.WaterVapor += curState.AirVapor[0][i];
			//dependent.DisplayPlanet.EnergyDelta += nextState.CellDisplays[i].EnergyDelta;
			//dependent.DisplayPlanet.EnergyEvapotranspiration += nextState.CellDisplays[i].EnergyEvapotranspiration;
			//dependent.DisplayPlanet.EnergyIncoming += nextState.CellDisplays[i].EnergyIncoming;
			//dependent.DisplayPlanet.EnergyLand += nextState.CellStates[i].GroundEnergy;
			//dependent.DisplayPlanet.EnergyLowerAir += nextState.CellStates[i].AirEnergy;
			//dependent.DisplayPlanet.EnergyOceanConduction += nextState.CellDisplays[i].EnergyOceanConduction;
			//dependent.DisplayPlanet.EnergyShallowWater += nextState.CellStates[i].WaterEnergy;
			//dependent.DisplayPlanet.EnergySolarAbsorbedAtmosphere += nextState.CellDisplays[i].EnergySolarAbsorbedAtmosphere;
			//dependent.DisplayPlanet.EnergySolarAbsorbedCloud += nextState.CellDisplays[i].EnergySolarAbsorbedCloud;
			//dependent.DisplayPlanet.EnergySolarAbsorbedOcean += nextState.CellDisplays[i].EnergySolarAbsorbedOcean;
			//dependent.DisplayPlanet.EnergySolarAbsorbedSurface += nextState.CellDisplays[i].EnergySolarAbsorbedSurface;
			//dependent.DisplayPlanet.EnergySolarReflectedAtmosphere += nextState.CellDisplays[i].EnergySolarReflectedAtmosphere;
			//dependent.DisplayPlanet.EnergySolarReflectedCloud += nextState.CellDisplays[i].EnergySolarReflectedCloud;
			//dependent.DisplayPlanet.EnergySolarReflectedSurface += nextState.CellDisplays[i].EnergySolarReflectedSurface;
			//dependent.DisplayPlanet.EnergySurfaceConduction += nextState.CellDisplays[i].EnergySurfaceConduction;
			//dependent.DisplayPlanet.EnergyThermalAbsorbedAtmosphere += nextState.CellDisplays[i].EnergySolarAbsorbedAtmosphere;
			//dependent.DisplayPlanet.EnergyThermalBackRadiation += nextState.CellDisplays[i].EnergyThermalBackRadiation;
			//dependent.DisplayPlanet.EnergyThermalOceanRadiation += nextState.CellDisplays[i].EnergyThermalOceanRadiation;
			//dependent.DisplayPlanet.EnergyThermalOutAtmosphere += nextState.CellDisplays[i].EnergyThermalOutAtmosphere;
			//dependent.DisplayPlanet.EnergyThermalSurfaceOutAtmosphericWindow += nextState.CellDisplays[i].EnergyThermalSurfaceOutAtmosphericWindow;
			//dependent.DisplayPlanet.EnergyThermalSurfaceRadiation += nextState.CellDisplays[i].EnergyThermalSurfaceRadiation;
			//dependent.DisplayPlanet.Evaporation += nextState.CellDisplays[i].Evaporation;
			//dependent.DisplayPlanet.OceanVolume += nextState.CellDependents[i].WaterAndIceDepth;
			//dependent.DisplayPlanet.Rainfall += nextState.CellDisplays[i].Rainfall;
			//dependent.DisplayPlanet.SeaLevel += nextState.CellTerrains[i].Elevation + nextState.CellDependents[i].WaterAndIceDepth;

		}

	}
}
