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

	private NativeArray<float> solarRadiation;
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

		solarRadiation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		dewPoint = new NativeArray<float>(_cellCount, Allocator.Persistent);
		waterSlopeAlbedo = new NativeArray<float>(_cellCount, Allocator.Persistent);

		cloudAirConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		airIceConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		airWaterConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		airTerrainConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		iceWaterConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		iceTerrainConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		waterTerrainConduction = new NativeArray<float>(_cellCount, Allocator.Persistent);

		emissivity = new NativeArray<float>[_layerCount];
		solarRadiationIn = new NativeArray<float>[_layerCount];
		thermalRadiationIn = new NativeArray<float>[_layerCount];
		thermalRadiationOut = new NativeArray<float>[_layerCount];
		thermalRadiationTransmittedUp = new NativeArray<float>[_layerCount];
		thermalRadiationTransmittedDown = new NativeArray<float>[_layerCount];
		energyJobDependencies = new NativeArray<JobHandle>[_layerCount];
		energyJobHandles = new NativeArray<JobHandle>(_layerCount, Allocator.Persistent);

		diffusionAir = new NativeArray<DiffusionAir>[_airLayers];
		advectionAir = new NativeArray<DiffusionAir>[_airLayers];
		diffusionWater = new NativeArray<DiffusionWater>[_waterLayers];
		advectionWater = new NativeArray<DiffusionWater>[_waterLayers];

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

	}

	public void Dispose()
	{
		solarRadiation.Dispose();
		dewPoint.Dispose();
		waterSlopeAlbedo.Dispose();

		cloudAirConduction.Dispose();
		airIceConduction.Dispose();
		airWaterConduction.Dispose();
		airTerrainConduction.Dispose();
		iceWaterConduction.Dispose();
		iceTerrainConduction.Dispose();
		waterTerrainConduction.Dispose();
		terrainEnergyJobHandleDependencies.Dispose();
		energyJobHandles.Dispose();

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
	}

	public void Tick(SimState[] states, int stateCount, int ticksToAdvance, ref DependentState dependent, ref StaticState staticState, ref WorldData worldData, ref int curStateIndex)
	{
		JobHandle lastJobHandle = default(JobHandle);

		for (int i = 0; i < ticksToAdvance; i++)
		{
			ref var lastState = ref states[curStateIndex];
			curStateIndex = (curStateIndex + 1) % stateCount;
			ref var nextState = ref states[curStateIndex];

			nextState.PlanetState.Ticks = lastState.PlanetState.Ticks + 1;

			// TODO: update
			float distanceToSun = math.length(lastState.PlanetState.Position);
			float angleToSun = math.atan2(lastState.PlanetState.Position.z, lastState.PlanetState.Position.x);
			angleToSun += lastState.PlanetState.OrbitSpeed;
			nextState.PlanetState.Position = new float3(math.cos(angleToSun), 0, math.sin(angleToSun)) * distanceToSun;
			nextState.PlanetState.Rotation = new float3(lastState.PlanetState.Rotation.x, Mathf.Repeat(lastState.PlanetState.Rotation.y + lastState.PlanetState.SpinSpeed * worldData.SecondsPerTick, math.PI * 2), 0);

			int lastStateIndex = i % 2;

			/*

			// SOLAR RADIATION INITIALIZATION

			var solarRadiationJob = new SolarRadiationJob()
			{
				SurfaceElevation = dependent.SurfaceElevation,
				SolarRadiation = solarRadiation,
				CloudCoverage = dependent.CloudCoverage,
				IceCoverage = dependent.IceCoverage,
				VegetationCoverage = dependent.VegetationCoverage,
				WaterCoverage = dependent.WaterCoverage,
				WaterSlopeAlbedo = waterSlopeAlbedo,

				IncomingSolarRadiation = lastState.PlanetState.SolarRadiation,
				PlanetRotation = quaternion.Euler(lastState.PlanetState.Rotation),
				SecondsPerTick = worldData.SecondsPerTick,
				SphericalPosition = staticState.SphericalPosition,
				SunToPlanetDir = math.normalize(lastState.PlanetState.Position),
				Terrain = lastState.Terrain,
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
					AirMass = dependent.AirMass[j],
					GreenhouseGasConcentration = lastState.PlanetState.CarbonDioxide,
					Humidity = lastState.AirHumidity[j]
				};
				emissivityJobHandles[j] = emissivityAirJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}

			for (int j = 0; j < _waterLayers; j++)
			{
				var emissivityWaterJob = new EmissivityWaterJob()
				{
					Emissivity = emissivity[_waterLayer0 + j], 
					WaterMass = lastState.WaterMass[j],
					WaterSaltMass = lastState.WaterSaltMass[j]
				};
				emissivityJobHandles[j + _waterLayer0] = emissivityWaterJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			var emissivityTerrainJob = new EmissivityTerrainJob()
			{
				Emissivity = emissivity[_terrainLayer],
				Terrain = lastState.Terrain,
				VegetationCoverage = dependent.VegetationCoverage,
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
					SolarRadiationAbsorbedCloud = solarRadiationIn[_cloudLayer],
					AirMass = dependent.AirMass[j],
					Humidity = lastState.AirHumidity[j],
					CloudMass = lastState.CloudMass,
					CloudTemperature = lastState.CloudTemperature,
					CloudDropletMass = lastState.CloudDropletMass,
					CloudElevation = lastState.CloudElevation,
					CloudCoverage = dependent.CloudCoverage,
					WaterSlopeAlbedo = waterSlopeAlbedo,
					SolarReflectivityAir = worldData.SolarReflectivityAir,
					SolarAbsorptivityAir = worldData.SolarAbsorptivityAir,
					SolarAbsorptivityWaterVapor = worldData.SolarAbsorptivityWaterVapor,
					SolarAbsorptivityCloud = worldData.SolarAbsorptivityCloud,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
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
				IceCoverage = dependent.IceCoverage
			};
			solarInJobHandles[_iceLayer] = solarInIceJob.Schedule(_cellCount, _batchCount, solarInJobHandles[_atmosphereLayer0]);

			for (int j = 0; j < _waterLayers; j++)
			{
				var solarRadiationAbsorbedWaterJob = new SolarRadiationAbsorbedWaterJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[_waterLayer0 + _waterLayers - j - 1],
					SolarRadiationIncoming = solarRadiation,
					WaterCoverage = dependent.WaterCoverage,
					WaterSlopeAlbedo = waterSlopeAlbedo,
				};
				solarInJobHandles[_waterLayer0 + j] = solarRadiationAbsorbedWaterJob.Schedule(_cellCount, _batchCount, solarInJobHandles[_waterLayer0 + j + 1]);
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
				ThermalRadiationEmitted = thermalRadiationOut[_iceLayer],
				ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[_iceLayer],
				ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[_iceLayer],
				Emissivity = WorldData.EmissivityIce,
				Temperature = lastState.IceTemperature, 
			};
			thermalOutJobHandles[_iceLayer] = thermalOutJob.Schedule(_cellCount, 100, lastJobHandle);

			// CLOUD
			thermalOutJob = new ThermalEnergyRadiatedConstantEmissivityJob()
			{
				ThermalRadiationEmitted = thermalRadiationOut[_cloudLayer],
				ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[_cloudLayer],
				ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[_cloudLayer],
				Emissivity = WorldData.EmissivityWater,
				Temperature = lastState.CloudTemperature,
			};
			thermalOutJobHandles[_cloudLayer] = thermalOutJob.Schedule(_cellCount, 100, lastJobHandle);


			// TERRAIN
			var thermalOutTerrainJob = new ThermalEnergyRadiatedTerrainJob()
			{
				ThermalRadiationEmitted = thermalRadiationOut[_terrainLayer],
				ThermalRadiationTransmitted = thermalRadiationTransmittedUp[_terrainLayer],
				Emissivity = emissivity[_terrainLayer],
				Temperature = lastState.TerrainTemperature,
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
					Temperature = lastState.AirTemperature[j]			
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
					Temperature = lastState.WaterTemperature[j],
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
						CloudElevation = lastState.CloudElevation,
						CloudCoverage = dependent.CloudCoverage,
						LayerElevation = dependent.LayerElevation[j],
						LayerHeight = dependent.LayerHeight[j],
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
						Coverage = dependent.IceCoverage
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
						Coverage = dependent.WaterCoverage
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
						CloudElevation = lastState.CloudElevation,
						CloudCoverage = dependent.CloudCoverage,
						LayerElevation = dependent.LayerElevation[j],
						LayerHeight = dependent.LayerHeight[j],
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
					Temperature = dependent.AirTemperaturePotential[j],
					Humidity = lastState.AirHumidity[j],
					Velocity = lastState.AirVelocity[j],
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
					Temperature = lastState.WaterTemperature[j],
					Salinity = lastState.WaterSaltMass[j],
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
				var advectionJob = new AdvectionAirJob()
				{
					Delta = advectionAir[j],
					Temperature = dependent.AirTemperaturePotential[j],
					Humidity = lastState.AirHumidity[j],
					Velocity = lastState.AirVelocity[j],
					Neighbors = staticState.Neighbors,
				};
				advectionJobHandles[j] = advectionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			for (int j = 0; j < _waterLayers; j++)
			{
				var advectionJob = new AdvectionWaterJob()
				{
					Delta = advectionAir[j],
					Temperature = lastState.WaterTemperature[j],
					Salinity = lastState.WaterSaltMass[j],
					Velocity = lastState.WaterVelocity[j],
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
				Temperature = nextState.TerrainTemperature,
				Terrain = lastState.Terrain,
				SolarRadiationInTerrain = solarRadiationIn[_terrainLayer],
				ThermalRadiationInTerrain = thermalRadiationIn[_terrainLayer],
				ThermalRadiationOutTerrain = thermalRadiationOut[_terrainLayer],
			};


			for (int j = 0; j < _airLayers; j++)
			{
				int layerIndex = _atmosphereLayer0 + j;
				var energyJob = new EnergyAirJob()
				{
					Temperature = nextState.AirTemperature[j],
					Humidity = nextState.AirHumidity[j],
					Velocity = nextState.AirVelocity[j],
					LastTemperature = lastState.AirTemperature[j],
					LastHumidity = lastState.AirHumidity[j],
					LastVelocity = lastState.AirVelocity[j],
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
					Temperature = nextState.WaterTemperature[j],
					Salinity = nextState.WaterSaltMass[j],
					Velocity = nextState.WaterVelocity[j],
					LastTemperature = lastState.WaterTemperature[j],
					LastSalinity = lastState.WaterSaltMass[j],
					LastVelocity = lastState.WaterVelocity[j],
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
			
			var updateDependentsJob = new UpdateDependentStateJob()
			{
				SurfaceElevation = dependent.SurfaceElevation,
				CloudCoverage = dependent.CloudCoverage,
				IceCoverage = dependent.IceCoverage,
				VegetationCoverage = dependent.VegetationCoverage,
				WaterCoverage = dependent.WaterCoverage,

				Terrain = lastState.Terrain,
				worldData = worldData
			};
			lastJobHandle = updateDependentsJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			*/
		}
		lastJobHandle.Complete();

		int outputBuffer = ticksToAdvance % 2;

		dependent.DisplayPlanet = new PlanetDisplay();
		var curState = states[curStateIndex];
		for (int i = 0; i < _cellCount; i++)
		{

			dependent.DisplayPlanet.CloudCoverage += dependent.CloudCoverage[i];
			dependent.DisplayPlanet.CloudMass += curState.CloudMass[i];
			dependent.DisplayPlanet.IceMass += curState.IceMass[i];
			dependent.DisplayPlanet.OceanCoverage += dependent.WaterCoverage[i];
			dependent.DisplayPlanet.Temperature += curState.AirTemperature[0][i];
			dependent.DisplayPlanet.WaterVapor += curState.AirHumidity[0][i];
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
