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
	private int _atmosphereLayer0;

	private NativeArray<CellTerrain>[] terrain;
	private NativeArray<float>[] terrainTemperature;
	private NativeArray<float>[] iceMass;
	private NativeArray<float>[] iceTemperature;
	private NativeArray<float>[] cloudMass;
	private NativeArray<float>[] cloudTemperature;
	private NativeArray<float>[] cloudDropletMass;
	private NativeArray<float>[] cloudElevation;
	private NativeArray<float2>[] cloudVelocity;
	private NativeArray<float>[][] airTemperature;
	private NativeArray<float>[][] airHumidity;
	private NativeArray<float2>[][] airVelocity;
	private NativeArray<float>[][] waterTemperature;
	private NativeArray<float>[][] waterSalinity;
	private NativeArray<float2>[][] waterVelocity;

	private NativeArray<float> surfaceElevation;
	private NativeArray<float> solarRadiation;
	private NativeArray<float> iceCoverage;
	private NativeArray<float> waterCoverage;
	private NativeArray<float> cloudCoverage;
	private NativeArray<float> vegetationCoverage;
	private NativeArray<float> dewPoint;
	private NativeArray<float> waterSlopeAlbedo;

	private NativeArray<float>[] emissivity;
	private NativeArray<float>[] solarRadiationIn;
	private NativeArray<float>[] thermalRadiationIn;
	private NativeArray<float>[] thermalRadiationOut;
	private NativeArray<float>[] thermalRadiationTransmittedUp;
	private NativeArray<float>[] thermalRadiationTransmittedDown;

	private NativeArray<float> cloudAirConduction;
	private NativeArray<float> airIceConduction;
	private NativeArray<float> airWaterConduction;
	private NativeArray<float> airTerrainConduction;
	private NativeArray<float> iceWaterConduction;
	private NativeArray<float> iceTerrainConduction;
	private NativeArray<float> waterTerrainConduction;


	private NativeArray<DiffusionAir>[] diffusionAir;
	private NativeArray<DiffusionAir>[] advectionAir;
	private NativeArray<DiffusionWater>[] diffusionWater;
	private NativeArray<DiffusionWater>[] advectionWater;

	private NativeArray<JobHandle> terrainEnergyJobHandleDependencies;
	private NativeArray<JobHandle> energyJobHandles;
	private NativeArray<JobHandle>[] energyJobDependencies;


	public WorldSim(int cellCount, int atmosphericLayers, int waterLayers)
	{
		_cellCount = cellCount;
		_airLayers = atmosphericLayers;
		_waterLayers = waterLayers;
		_layerCount = _airLayers + _waterLayers + _terrainLayers;
		_terrainLayer = 0;
		_waterLayer0 = _terrainLayer + 1;
		_iceLayer = _waterLayer0 + _waterLayers;
		_atmosphereLayer0 = _iceLayer + 1;
		_cloudLayer = _layerCount - 1;

		terrain = new NativeArray<CellTerrain>[2];
		airTemperature = new NativeArray<float>[2][];
		airHumidity = new NativeArray<float>[2][];
		airVelocity = new NativeArray<float2>[2][];
		waterTemperature = new NativeArray<float>[2][];
		waterSalinity = new NativeArray<float>[2][];
		waterVelocity = new NativeArray<float2>[2][];
		for (int i = 0; i < 2; i++)
		{
			terrain[i] = new NativeArray<CellTerrain>(_cellCount, Allocator.Persistent);
			terrainTemperature[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			iceMass[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			iceTemperature[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			cloudMass[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			cloudTemperature[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			cloudDropletMass[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			cloudElevation[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			cloudVelocity[i] = new NativeArray<float2>(_cellCount, Allocator.Persistent);
			airTemperature[i] = new NativeArray<float>[_waterLayers];
			airHumidity[i] = new NativeArray<float>[_waterLayers];
			airVelocity[i] = new NativeArray<float2>[_waterLayers];
			waterTemperature[i] = new NativeArray<float>[_waterLayers];
			waterSalinity[i] = new NativeArray<float>[_waterLayers];
			waterVelocity[i] = new NativeArray<float2>[_waterLayers];
			for (int j = 0; j < _airLayers; j++)
			{
				airTemperature[i][j] = new NativeArray<float>(_cellCount * _airLayers, Allocator.Persistent);
				airHumidity[i][j] = new NativeArray<float>(_cellCount * _airLayers, Allocator.Persistent);
				airVelocity[i][j] = new NativeArray<float2>(_cellCount * _airLayers, Allocator.Persistent);
			}
			for (int j = 0; j < _waterLayers; j++)
			{
				waterTemperature[i][j] = new NativeArray<float>(_cellCount * _waterLayers, Allocator.Persistent);
				waterSalinity[i][j] = new NativeArray<float>(_cellCount * _waterLayers, Allocator.Persistent);
				waterVelocity[i][j] = new NativeArray<float2>(_cellCount * _waterLayers, Allocator.Persistent);
			}
		}
		surfaceElevation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		solarRadiation = new NativeArray<float>(_cellCount, Allocator.Persistent);

		cloudAirConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		airIceConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		airWaterConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		airTerrainConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		iceWaterConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		iceTerrainConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		waterTerrainConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);

		iceCoverage = new NativeArray<float>(_cellCount, Allocator.Persistent);
		waterCoverage = new NativeArray<float>(_cellCount, Allocator.Persistent);
		cloudCoverage = new NativeArray<float>(_cellCount, Allocator.Persistent);
		vegetationCoverage = new NativeArray<float>(_cellCount, Allocator.Persistent);
		dewPoint = new NativeArray<float>(_cellCount, Allocator.Persistent);
		waterSlopeAlbedo = new NativeArray<float>(_cellCount, Allocator.Persistent);

		emissivity = new NativeArray<float>[_layerCount];
		solarRadiationIn = new NativeArray<float>[_layerCount];
		thermalRadiationIn = new NativeArray<float>[_layerCount];
		thermalRadiationOut = new NativeArray<float>[_layerCount];
		thermalRadiationTransmittedUp = new NativeArray<float>[_layerCount];
		thermalRadiationTransmittedDown = new NativeArray<float>[_layerCount];

		diffusionAir = new NativeArray<DiffusionAir>[_airLayers];
		advectionAir = new NativeArray<DiffusionAir>[_airLayers];
		diffusionWater = new NativeArray<DiffusionWater>[_waterLayers];
		advectionWater = new NativeArray<DiffusionWater>[_waterLayers];

		energyJobDependencies = new NativeArray<JobHandle>[_layerCount];
		for (int i = 0; i < _layerCount; i++)
		{
			solarRadiationIn[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			emissivity[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			thermalRadiationIn[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			thermalRadiationOut[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			thermalRadiationTransmittedUp[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			thermalRadiationTransmittedDown[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			energyJobDependencies[i] = new NativeArray<JobHandle>(6, Allocator.Persistent);
		}

		terrainEnergyJobHandleDependencies = new NativeArray<JobHandle>(4, Allocator.Persistent);
		energyJobHandles = new NativeArray<JobHandle>(_layerCount, Allocator.Persistent);

	}

	public void Dispose()
	{
		for (int i = 0; i < 2; i++)
		{
			terrain[i].Dispose();
			terrainTemperature[i].Dispose();
			iceMass[i].Dispose();
			iceTemperature[i].Dispose();
			cloudMass[i].Dispose();
			cloudTemperature[i].Dispose();
			cloudDropletMass[i].Dispose();
			cloudVelocity[i].Dispose();
			for (int j = 0; j < _layerCount; j++)
			{
				airTemperature[i][j].Dispose();
				airHumidity[i][j].Dispose();
				airVelocity[i][j].Dispose();
				waterTemperature[i][j].Dispose();
				waterSalinity[i][j].Dispose();
				waterVelocity[i][j].Dispose();
			}
		}
		surfaceElevation.Dispose();
		solarRadiation.Dispose();
		iceCoverage.Dispose();
		waterCoverage.Dispose();
		cloudCoverage.Dispose();
		vegetationCoverage.Dispose();
		dewPoint.Dispose();
		waterSlopeAlbedo.Dispose();

		cloudAirConduction.Dispose();
		airIceConduction.Dispose();
		airWaterConduction.Dispose();
		airTerrainConduction.Dispose();
		iceWaterConduction.Dispose();
		iceTerrainConduction.Dispose();
		waterTerrainConduction.Dispose();


		for (int i = 0; i < _layerCount; i++)
		{
			emissivity[i].Dispose();
			solarRadiationIn[i].Dispose();
			thermalRadiationIn[i].Dispose();
			thermalRadiationOut[i].Dispose();
			thermalRadiationTransmittedUp[i].Dispose();
			thermalRadiationTransmittedDown[i].Dispose();
			energyJobDependencies[i].Dispose();
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

		terrainEnergyJobHandleDependencies.Dispose();
		energyJobHandles.Dispose();

	}

	public void Tick(ref SimState state, ref SimState nextState, int ticksToAdvance, ref StaticState staticState, ref WorldData worldData)
	{
		JobHandle lastJobHandle = default(JobHandle);

		terrain[0].CopyFrom(state.CellTerrains);
		terrainTemperature[0].CopyFrom(state.TerrainTemperature);
		iceMass[0].CopyFrom(state.IceMass);
		iceTemperature[0].CopyFrom(state.IceTemperature);
		cloudMass[0].CopyFrom(state.CloudMass);
		cloudTemperature[0].CopyFrom(state.CloudTemperature);
		cloudElevation[0].CopyFrom(state.CloudElevation);
		cloudDropletMass[0].CopyFrom(state.CloudDropletMass);
		cloudVelocity[0].CopyFrom(state.CloudVelocity);
		for (int i = 0; i < _airLayers; i++) {
			airTemperature[0][i].CopyFrom(state.AirTemperature[i]);
			airHumidity[0][i].CopyFrom(state.AirHumidity[i]);
			airVelocity[0][i].CopyFrom(state.AirVelocity[i]);
		}
		for (int i=0;i<_waterLayers;i++)
		{
			waterTemperature[0][i].CopyFrom(state.WaterTemperature[i]);
			waterSalinity[0][i].CopyFrom(state.WaterSaltMass[i]);
			waterVelocity[0][i].CopyFrom(state.WaterVelocity[i]);
		}

		for (int i = 0; i < ticksToAdvance; i++)
		{
			state.PlanetState.Ticks++;

			// TODO: update
			float distanceToSun = math.length(state.PlanetState.Position);
			float angleToSun = math.atan2(state.PlanetState.Position.z, state.PlanetState.Position.x);
			angleToSun += state.PlanetState.OrbitSpeed;
			state.PlanetState.Position = new float3(math.cos(angleToSun), 0, math.sin(angleToSun)) * distanceToSun;
			state.PlanetState.Rotation = new float3(state.PlanetState.Rotation.x, Mathf.Repeat(state.PlanetState.Rotation.y + state.PlanetState.SpinSpeed * worldData.SecondsPerTick, math.PI * 2), 0);

			int lastStateIndex = i % 2;

			// SOLAR RADIATION INITIALIZATION

			var solarRadiationJob = new SolarRadiationJob()
			{
				SurfaceElevation = surfaceElevation,
				SolarRadiation = solarRadiation,
				CloudCoverage = cloudCoverage,
				IceCoverage = iceCoverage,
				VegetationCoverage = vegetationCoverage,
				WaterCoverage = waterCoverage,
				WaterSlopeAlbedo = waterSlopeAlbedo,

				IncomingSolarRadiation = state.PlanetState.SolarRadiation,
				PlanetRotation = quaternion.Euler(state.PlanetState.Rotation),
				SecondsPerTick = worldData.SecondsPerTick,
				SphericalPosition = staticState.SphericalPosition,
				SunToPlanetDir = math.normalize(state.PlanetState.Position),
				Terrain = terrain[lastStateIndex],
				worldData = worldData
			};
			var solarInJobHandle = solarRadiationJob.Schedule(_cellCount, _batchCount, lastJobHandle);



			// EMISSIVITY

			JobHandle[] emissivityJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				var emissivityAirJob = new EmissivityAirJob()
				{
					Emissivity = emissivity[_atmosphereLayer0 + j],
					AbsorptivityAir = worldData.AbsorptivityAir,
					AbsorptivityCarbonDioxide = worldData.AbsorptivityCarbonDioxide,
					AbsorptivityWaterVapor = worldData.AbsorptivityWaterVapor,
					AirMass = airMass[j],
					GreenhouseGasConcentration = state.PlanetState.CarbonDioxide,
					Humidity = airHumidity[lastStateIndex][j]
				};
				emissivityJobHandles[j] = emissivityAirJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}

			for (int j = 0; j < _waterLayers; j++)
			{
				var emissivityWaterJob = new EmissivityWaterJob()
				{
					Emissivity = emissivity[_waterLayer0 + j], 
					WaterMass = waterMass[lastStateIndex][j],
					WaterSaltMass = waterSalinity[lastStateIndex][j]
				};
				emissivityJobHandles[j + _waterLayer0] = emissivityWaterJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			var emissivityTerrainJob = new EmissivityTerrainJob()
			{
				Emissivity = emissivity[_terrainLayer],
				Terrain = terrain[lastStateIndex],
				VegetationCoverage = vegetationCoverage,
			};
			emissivityJobHandles[_terrainLayer] = emissivityTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle);



			// SOLAR RADIATION ABSORPTION
			// process each vertical layer in order

			// atmosphere
			JobHandle[] solarInJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				int layerIndex = _atmosphereLayer0 + _airLayers - j - 1;
				var solarInAtmosphereJob = new SolarRadiationAbsorbedAirJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[layerIndex],
					SolarRadiationIncoming = solarRadiation,
					SolarRadiationAbsorbedCloud = solarRadiationInCloud,
					AirMass = airMass[j],
					Humidity = airHumidity[lastStateIndex][j],
					CloudMass = cloudMass[lastStateIndex],
					CloudTemperature = cloudTemperature[lastStateIndex],
					CloudDropletMass = cloudDropletMass[lastStateIndex],
					CloudElevation = cloudElevation[lastStateIndex],
					CloudCoverage = cloudCoverage,
					WaterSlopeAlbedo = waterSlopeAlbedo,
					SolarReflectivityAir = worldData.SolarReflectivityAir,
					SolarAbsorptivityAir = worldData.SolarAbsorptivityAir,
					SolarAbsorptivityWaterVapor = worldData.SolarAbsorptivityWaterVapor,
					SolarAbsorptivityCloud = worldData.SolarAbsorptivityCloud,
					LayerElevation = ,
					LayerHeight = ,
					worldData = worldData
				};
				solarInJobHandles[_atmosphereLayer0 + j] = solarInAtmosphereJob.Schedule(_cellCount, _batchCount, j == _airLayers - 1 ? solarInJobHandle : solarInJobHandles[_atmosphereLayer0 + j + 1]);
			}

			// ice
			var solarInIceJob = new SolarRadiationAbsorbedIceJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_iceLayer],
				SolarRadiationIncoming = solarRadiation,
				AlbedoIce = WorldData.AlbedoIce,
				IceCoverage = iceCoverage
			};
			solarInJobHandles[_iceLayer] = solarInIceJob.Schedule(_cellCount, _batchCount, solarInJobHandles[_atmosphereLayer0]);

			for (int j = 0; j < _waterLayers; j++)
			{
				var solarRadiationAbsorbedWaterJob = new SolarRadiationAbsorbedWaterJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[_waterLayer0 + _waterLayers - j - 1],
					SolarRadiationIncoming = solarRadiation,
					WaterCoverage = waterCoverage,
					WaterSlopeAlbedo = waterSlopeAlbedo,
				};
				solarInJobHandles[_waterLayer0 + j] = solarRadiationAbsorbedWaterJob.Schedule(_cellCount, _batchCount, solarInJobHandles[_waterLayer0 + j + 1]);
			}

			var solarRadiationAbsorbedTerrainJob = new SolarRadiationAbsorbedTerrainJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_terrainLayer],
				SolarRadiationIncoming = solarRadiation,
				VegetationCoverage = vegetationCoverage,
				worldData = worldData,
				LastTerrain = terrain[lastStateIndex],
			};
			solarInJobHandles[_terrainLayer] = solarRadiationAbsorbedTerrainJob.Schedule(_cellCount, _batchCount, solarInJobHandles[_terrainLayer + 1]);

			// THERMAL RADIATION
			JobHandle[] thermalOutJobHandles = new JobHandle[_layerCount];

			// ICE
			var thermalOutJob = new ThermalEnergyRadiatedConstantEmissivityJob()
			{
				ThermalRadiationEmitted = thermalRadiationOut[_iceLayer],
				ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[_iceLayer],
				ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[_iceLayer],
				Emissivity = WorldData.EmissivityIce,
				Temperature = iceTemperature[lastStateIndex], 
			};
			thermalOutJobHandles[_iceLayer] = thermalOutJob.Schedule(_cellCount, 100, lastJobHandle);

			// CLOUD
			thermalOutJob = new ThermalEnergyRadiatedConstantEmissivityJob()
			{
				ThermalRadiationEmitted = thermalRadiationOut[_cloudLayer],
				ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[_cloudLayer],
				ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[_cloudLayer],
				Emissivity = WorldData.EmissivityWater,
				Temperature = cloudTemperature[lastStateIndex],
			};
			thermalOutJobHandles[_cloudLayer] = thermalOutJob.Schedule(_cellCount, 100, lastJobHandle);


			// TERRAIN
			var thermalOutTerrainJob = new ThermalEnergyRadiatedTerrainJob()
			{
				ThermalRadiationEmitted = thermalRadiationOut[_terrainLayer],
				ThermalRadiationTransmitted = thermalRadiationTransmittedUp[_terrainLayer],
				Emissivity = emissivity[_terrainLayer],
				Temperature = terrainTemperature[lastStateIndex],
			};
			thermalOutJobHandles[_terrainLayer] = thermalOutJob.Schedule(_cellCount, 100, emissivityJobHandles[_terrainLayer]);


			// ATMOSPHERE
			for (int j = 0; j < _airLayers; j++)
			{
				int layer = _waterLayer0 + j;
				var thermalOutAirJob = new ThermalEnergyRadiatedJob()
				{
					ThermalRadiationEmitted = thermalRadiationOut[layer],
					ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[layer],
					ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[layer],
					Emissivity = emissivity[_atmosphereLayer0 + j],
					Temperature = airTemperature[lastStateIndex][j]			
				};
				thermalOutJobHandles[_atmosphereLayer0 + j] = thermalOutAirJob.Schedule(_cellCount, 100, emissivityJobHandles[_atmosphereLayer0 + j]);
			}

			// WATER
			for (int j = 0; j < _waterLayers; j++)
			{
				int layer = _waterLayer0 + j;
				var thermalOutWaterJob = new ThermalEnergyRadiatedJob()
				{
					ThermalRadiationEmitted = thermalRadiationOut[layer],
					ThermalRadiationTransmittedDown = thermalRadiationTransmittedUp[layer],
					ThermalRadiationTransmittedUp =thermalRadiationTransmittedUp[layer],
					Emissivity = emissivity[layer],
					Temperature = waterTemperature[lastStateIndex][j],
				};
				thermalOutJobHandles[_waterLayer0 + j] = thermalOutWaterJob.Schedule(_cellCount, 100, emissivityJobHandles[_waterLayer0 + j]);
			}


			// THERMAL RADIATION ABSORPTION
			// Start at bottom water layer and go up, then go back down
			JobHandle[] thermalInUpJobHandles = new JobHandle[_layerCount];
			JobHandle[] thermalInDownJobHandles = new JobHandle[_layerCount];

			// transmit up from land
			for (int j = 1; j < _layerCount; j++)
			{
				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInUpJobHandles[j-1], thermalOutJobHandles[j]);

				if (j >= _atmosphereLayer0 && j < _atmosphereLayer0 + _airLayers)
				{
					var thermalInJob = new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationIn[j],
						ThermalRadiationAbsorbedCloud = thermalRadiationIn[_cloudLayer],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j - 1],
						ThermalRadiationEmitted = thermalRadiationTransmittedUp[j],
						ThermalRadiationEmittedCloudUp = thermalRadiationTransmittedUp[_cloudLayer],
						ThermalRadiationEmittedCloudDown = thermalRadiationTransmittedDown[_cloudLayer],
						Emissivity = emissivity[j],
						CloudElevation = cloudElevation[lastStateIndex],
						CloudCoverage = cloudCoverage,
						LayerElevation = ,
						LayerHeight = ,
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j == _iceLayer)
				{
					var thermalInJob = new ThermalEnergyAbsorbedPartialCoverageConstantEmissivityJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationIn[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j - 1],
						ThermalRadiationEmitted = thermalRadiationTransmittedUp[j],
						Emissivity = WorldData.EmissivityIce,
						Coverage = iceCoverage
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j == _waterLayer0 + _waterLayers - 1)
				{
					var thermalInJob = new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationIn[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j - 1],
						ThermalRadiationEmitted = thermalRadiationTransmittedUp[j],
						Emissivity = emissivity[j],
						Coverage = waterCoverage
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else
				{
					var thermalInJob = new ThermalEnergyAbsorbedJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationIn[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j - 1],
						ThermalRadiationEmitted = thermalRadiationTransmittedUp[j],
						Emissivity = emissivity[j]
					};
					thermalInUpJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
			}

			// transmit down from top of atmosphere
			for (int j = _layerCount-2; j >= 0; j--)
			{
				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[j + 1], thermalOutJobHandles[j]);

				if (j == _terrainLayer) {
					var thermalInJob = new ThermalEnergyAbsorbedTerrainJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationIn[j],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[j + 1],
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else if (j >= _atmosphereLayer0 && j < _atmosphereLayer0 + _airLayers)
				{
					var thermalInJob = new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationIn[j],
						ThermalRadiationAbsorbedCloud = thermalRadiationIn[_cloudLayer],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j + 1],
						ThermalRadiationEmitted = thermalRadiationTransmittedDown[j],
						ThermalRadiationEmittedCloudUp = thermalRadiationTransmittedUp[_cloudLayer],
						ThermalRadiationEmittedCloudDown = thermalRadiationTransmittedDown[_cloudLayer],
						Emissivity = emissivity[j],
						CloudElevation = cloudElevation[lastStateIndex],
						CloudCoverage = cloudCoverage,
						LayerElevation = ,
						LayerHeight = ,
					};
					thermalInDownJobHandles[j] = thermalInJob.Schedule(_cellCount, _batchCount, thermalInDependenciesHandle);
				}
				else
				{
					var thermalInJob = new ThermalEnergyAbsorbedJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationIn[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j + 1],
						ThermalRadiationEmitted = thermalRadiationTransmittedDown[j],
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
					Temperature = airTemperature[lastStateIndex][j],
					Humidity = airHumidity[lastStateIndex][j],
					Velocity = airVelocity[lastStateIndex][j],
					Elevation = , 
					Neighbors = staticState.Neighbors,
					DiffusionCoefficientHumidity = worldData.AirMassDiffusionSpeedHorizontal,
					DiffusionCoefficientTemperature = worldData.AirMassDiffusionSpeedHorizontal,
					DiffusionCoefficientVelocity = worldData.AirMassDiffusionSpeedHorizontal,
				};
				diffusionJobHandles[_atmosphereLayer0 + j] = diffusionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			for (int j = 0; j < _waterLayer0; j++)
			{
				var diffusionJob = new DiffusionWaterJob()
				{
					Delta = diffusionWater[j],
					Temperature = waterTemperature[lastStateIndex][j],
					Salinity = waterSalinity[lastStateIndex][j],
					Velocity = waterVelocity[lastStateIndex][j],
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
				var advectionJob = new AdvectionAirJob()
				{
					Delta = advectionAir[j],
					Temperature = airTemperature[lastStateIndex][j],
					Humidity = airHumidity[lastStateIndex][j],
					Velocity = airVelocity[lastStateIndex][j],
					Neighbors = staticState.Neighbors,
				};
				advectionJobHandles[j] = advectionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			for (int j = 0; j < _waterLayers; j++)
			{
				var advectionJob = new AdvectionWaterJob()
				{
					Delta = advectionAir[j],
					Temperature = waterTemperature[lastStateIndex][j],
					Salinity = waterSalinity[lastStateIndex][j],
					Velocity = waterVelocity[lastStateIndex][j],
					Neighbors = staticState.Neighbors,
				};
				advectionJobHandles[j] = advectionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}

			// CONDUCTION
			List<JobHandle> conductionJobHandles = new List<JobHandle>();

			// cloud to air
			var conductionCloudToAirJob = new ConductionJob()
			{
				EnergyDelta = cloudAirConduction
			};
			conductionJobHandles.Add(conductionCloudToAirJob.Schedule(_cellCount, _batchCount, lastJobHandle));

			// air to ice
			var conductionAirToIceJob = new ConductionJob()
			{
				EnergyDelta = airIceConduction
			};
			conductionJobHandles.Add(conductionAirToIceJob.Schedule(_cellCount, _batchCount, lastJobHandle));

			// air to water
			var conductionAirToWaterJob = new ConductionAirWaterJob()
			{
				EnergyDelta = airWaterConduction
			};
			conductionJobHandles.Add(conductionAirToWaterJob.Schedule(_cellCount, _batchCount, lastJobHandle));

			// air to terrain
			var conductionAirToTerrainJob = new ConductionAirTerrainJob()
			{
				EnergyDelta = airTerrainConduction
			};
			conductionJobHandles.Add(conductionAirToTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle));

			// ice to water
			var conductionIceToWaterJob = new ConductionIceWaterJob()
			{
				EnergyDelta = iceWaterConduction
			};
			conductionJobHandles.Add(conductionIceToWaterJob.Schedule(_cellCount, _batchCount, lastJobHandle));

			// ice to terrain
			var conductionIceToTerrainJob = new ConductionIceTerrainJob()
			{
				EnergyDelta = iceTerrainConduction
			};
			conductionJobHandles.Add(conductionIceToTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle));

			// water to terrain
			var conductionWaterToTerrainJob = new ConductionWaterTerrainJob()
			{
				EnergyDelta = waterTerrainConduction
			};
			conductionJobHandles.Add(conductionWaterToTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle));


			// COMBINE ADVECTION, DIFFUSION, SOLAR, THERMAL DELTA

			var energyTerrainJob = new EnergyTerrainJob()
			{
				Temperature = terrainTemperature[lastStateIndex],
				Terrain = terrain[i],
				SolarRadiationInTerrain = solarRadiationIn[_terrainLayer],
				ThermalRadiationInTerrain = thermalRadiationIn[_terrainLayer],
				ThermalRadiationOutTerrain = thermalRadiationOut[_terrainLayer],
			};


			for (int j = 0; j < _airLayers; j++)
			{
				int layerIndex = _atmosphereLayer0 + j;
				var energyJob = new EnergyAirJob()
				{
					Temperature = airTemperature[i][j],
					Humidity = airHumidity[i][j],
					Velocity = airVelocity[i][j],
					LastTemperature = airTemperature[lastStateIndex][j],
					LastHumidity = airHumidity[lastStateIndex][j],
					LastVelocity = airVelocity[lastStateIndex][j],
					Advection = advectionAir[j],
					Diffusion = diffusionAir[j],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationIn = thermalRadiationIn[layerIndex],
					ThermalRadiationOut = thermalRadiationOut[layerIndex],
				};

				energyJobDependencies[layerIndex][0] = advectionJobHandles[j];
				energyJobDependencies[layerIndex][1] = diffusionJobHandles[j];
				energyJobDependencies[layerIndex][2] = solarInJobHandles[layerIndex];
				energyJobDependencies[layerIndex][3] = thermalOutJobHandles[layerIndex];
				energyJobDependencies[layerIndex][4] = thermalInDownJobHandles[layerIndex];
				energyJobDependencies[layerIndex][5] = thermalInUpJobHandles[layerIndex];
				var energyJobDependencyHandles = JobHandle.CombineDependencies(energyJobDependencies[j]);
				energyJobHandles[layerIndex] = energyJob.Schedule(_cellCount, _batchCount, energyJobDependencyHandles);
			}

			for (int j = 0; j < _waterLayers; j++)
			{
				int layerIndex = _waterLayer0 + j;
				var energyJob = new EnergyWaterJob()
				{
					Temperature = waterTemperature[i][j],
					Salinity = waterSalinity[i][j],
					Velocity = waterVelocity[i][j],
					LastTemperature = waterTemperature[lastStateIndex][j],
					LastSalinity = waterSalinity[lastStateIndex][j],
					LastVelocity = waterVelocity[lastStateIndex][j],
					Advection = advectionWater[j],
					Diffusion = diffusionWater[j],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationIn = thermalRadiationIn[layerIndex],
					ThermalRadiationOut = thermalRadiationOut[layerIndex],
				};

				energyJobDependencies[layerIndex][0] = advectionJobHandles[j];
				energyJobDependencies[layerIndex][1] = diffusionJobHandles[j];
				energyJobDependencies[layerIndex][2] = solarInJobHandles[layerIndex];
				energyJobDependencies[layerIndex][3] = thermalOutJobHandles[layerIndex];
				energyJobDependencies[layerIndex][4] = thermalInDownJobHandles[layerIndex];
				energyJobDependencies[layerIndex][5] = thermalInUpJobHandles[layerIndex];
				var energyJobDependencyHandles = JobHandle.CombineDependencies(energyJobDependencies[j]);
				energyJobHandles[layerIndex] = energyJob.Schedule(_cellCount, _batchCount, energyJobDependencyHandles);
			}

			terrainEnergyJobHandleDependencies[0] = solarInJobHandles[_terrainLayer];
			terrainEnergyJobHandleDependencies[1] = thermalOutJobHandles[_terrainLayer];
			terrainEnergyJobHandleDependencies[2] = thermalInDownJobHandles[_terrainLayer];
			terrainEnergyJobHandleDependencies[3] = thermalInUpJobHandles[_terrainLayer];
			var terrainEnergyDependenciesHandle = JobHandle.CombineDependencies(terrainEnergyJobHandleDependencies);
			energyJobHandles[_terrainLayer] = energyTerrainJob.Schedule(_cellCount, _batchCount, terrainEnergyDependenciesHandle);

			lastJobHandle = JobHandle.CombineDependencies(energyJobHandles);

		}
		lastJobHandle.Complete();

		int outputBuffer = ticksToAdvance % 2;

		nextState.PlanetState = state.PlanetState;

		nextState.DisplayPlanet = new PlanetDisplay();
		terrain[outputBuffer].CopyTo(nextState.CellTerrains);
		for (int i = 0; i < _airLayers; i++) {
			airTemperature[outputBuffer][i].CopyTo(nextState.AirTemperature[i]);
			airHumidity[outputBuffer][i].CopyTo(nextState.AirHumidity[i]);
			airVelocity[outputBuffer][i].CopyTo(nextState.AirVelocity[i]);
		}
		for (int i = 0; i < _waterLayers; i++) {
			waterTemperature[outputBuffer][i].CopyTo(nextState.WaterTemperature[i]);
			waterSalinity[outputBuffer][i].CopyTo(nextState.WaterSaltMass[i]);
			waterVelocity[outputBuffer][i].CopyTo(nextState.WaterVelocity[i]);
		}

		//for (int i = 0; i < _cellCount; i++)
		//{

		//nextState.DisplayPlanet.CloudCoverage += nextState.CellDependents[i].CloudCoverage;
		//nextState.DisplayPlanet.CloudMass += nextState.CellStates[i].CloudMass;
		//nextState.DisplayPlanet.EnergyDelta += nextState.CellDisplays[i].EnergyDelta;
		//nextState.DisplayPlanet.EnergyEvapotranspiration += nextState.CellDisplays[i].EnergyEvapotranspiration;
		//nextState.DisplayPlanet.EnergyIncoming += nextState.CellDisplays[i].EnergyIncoming;
		//nextState.DisplayPlanet.EnergyLand += nextState.CellStates[i].GroundEnergy;
		//nextState.DisplayPlanet.EnergyLowerAir += nextState.CellStates[i].AirEnergy;
		//nextState.DisplayPlanet.EnergyOceanConduction += nextState.CellDisplays[i].EnergyOceanConduction;
		//nextState.DisplayPlanet.EnergyShallowWater += nextState.CellStates[i].WaterEnergy;
		//nextState.DisplayPlanet.EnergySolarAbsorbedAtmosphere += nextState.CellDisplays[i].EnergySolarAbsorbedAtmosphere;
		//nextState.DisplayPlanet.EnergySolarAbsorbedCloud += nextState.CellDisplays[i].EnergySolarAbsorbedCloud;
		//nextState.DisplayPlanet.EnergySolarAbsorbedOcean += nextState.CellDisplays[i].EnergySolarAbsorbedOcean;
		//nextState.DisplayPlanet.EnergySolarAbsorbedSurface += nextState.CellDisplays[i].EnergySolarAbsorbedSurface;
		//nextState.DisplayPlanet.EnergySolarReflectedAtmosphere += nextState.CellDisplays[i].EnergySolarReflectedAtmosphere;
		//nextState.DisplayPlanet.EnergySolarReflectedCloud += nextState.CellDisplays[i].EnergySolarReflectedCloud;
		//nextState.DisplayPlanet.EnergySolarReflectedSurface += nextState.CellDisplays[i].EnergySolarReflectedSurface;
		//nextState.DisplayPlanet.EnergySurfaceConduction += nextState.CellDisplays[i].EnergySurfaceConduction;
		//nextState.DisplayPlanet.EnergyThermalAbsorbedAtmosphere += nextState.CellDisplays[i].EnergySolarAbsorbedAtmosphere;
		//nextState.DisplayPlanet.EnergyThermalBackRadiation += nextState.CellDisplays[i].EnergyThermalBackRadiation;
		//nextState.DisplayPlanet.EnergyThermalOceanRadiation += nextState.CellDisplays[i].EnergyThermalOceanRadiation;
		//nextState.DisplayPlanet.EnergyThermalOutAtmosphere += nextState.CellDisplays[i].EnergyThermalOutAtmosphere;
		//nextState.DisplayPlanet.EnergyThermalSurfaceOutAtmosphericWindow += nextState.CellDisplays[i].EnergyThermalSurfaceOutAtmosphericWindow;
		//nextState.DisplayPlanet.EnergyThermalSurfaceRadiation += nextState.CellDisplays[i].EnergyThermalSurfaceRadiation;
		//nextState.DisplayPlanet.Evaporation += nextState.CellDisplays[i].Evaporation;
		//nextState.DisplayPlanet.IceMass += nextState.CellStates[i].IceMass;
		//nextState.DisplayPlanet.OceanCoverage += math.saturate(nextState.CellDependents[i].WaterDepth / nextState.CellTerrains[i].Roughness);
		//nextState.DisplayPlanet.OceanVolume += nextState.CellDependents[i].WaterAndIceDepth;
		//nextState.DisplayPlanet.Rainfall += nextState.CellDisplays[i].Rainfall;
		//nextState.DisplayPlanet.SeaLevel += nextState.CellTerrains[i].Elevation + nextState.CellDependents[i].WaterAndIceDepth;
		//nextState.DisplayPlanet.Temperature += nextState.CellDependents[i].AirTemperature;
		//nextState.DisplayPlanet.WaterVapor += nextState.CellStates[i].AirWaterMass;

		//}

	}
}
