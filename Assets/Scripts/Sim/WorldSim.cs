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
	private NativeArray<DiffusionAirVertical>[] verticalMovementAir;
	private NativeArray<float>[] latentHeat;
	private NativeArray<float2>[] pressureGradientForce;
	private NativeArray<float> windFriction;
	private NativeArray<float> conductionCloudAir;
	private NativeArray<float> conductionAirIce;
	private NativeArray<float> conductionAirWater;
	private NativeArray<float> conductionAirTerrain;
	private NativeArray<float> conductionIceWater;
	private NativeArray<float> conductionIceTerrain;
	private NativeArray<float>[] conductionWaterTerrain;

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
		_surfaceWaterLayer = _waterLayers - 1;

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
		pressureGradientForce = new NativeArray<float2>[_airLayers];
		for (int i = 0; i < _airLayers; i++)
		{
			diffusionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
			advectionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
			pressureGradientForce[i] = new NativeArray<float2>(_cellCount, Allocator.Persistent);
		}
		verticalMovementAir = new NativeArray<DiffusionAirVertical>[_airLayers];
		for (int i=0;i<_airLayers;i++)
		{
			verticalMovementAir[i] = new NativeArray<DiffusionAirVertical>(_cellCount, Allocator.Persistent);
		}
		diffusionWater = new NativeArray<DiffusionWater>[_waterLayers];
		advectionWater = new NativeArray<DiffusionWater>[_waterLayers];
		conductionWaterTerrain = new NativeArray<float>[_waterLayers];
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
			advectionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
			conductionWaterTerrain[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
		}
		windFriction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		diffusionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		advectionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		conductionCloudAir = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirIce = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
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
			diffusionAir[i].Dispose();
			advectionAir[i].Dispose();
			pressureGradientForce[i].Dispose();
		}
		for (int i=0;i<_airLayers;i++)
		{
			verticalMovementAir[i].Dispose();
		}
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i].Dispose();
			advectionWater[i].Dispose();
			conductionWaterTerrain[i].Dispose();
		}
		windFriction.Dispose();
		diffusionCloud.Dispose();
		advectionCloud.Dispose();
		conductionCloudAir.Dispose();
		conductionAirIce.Dispose();
		conductionAirWater.Dispose();
		conductionAirTerrain.Dispose();
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
			var evaporationMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var evaporationLatentHeat = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var frozenTopMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var frozenBottomMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var cloudEvaporationMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var rainfallWaterMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
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
			#endregion

			// Calculate emissivity per cell
			// TODO: combine cloud emissivity with cloud conduction
			// TODO: we can probably combine this step with the thermal radiation step
			#region Emissivity Per Cell
			JobHandle[] emissivityJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				int layerIndex = _airLayer0 + j;
				var emissivityAirJob = new EmissivityAirJob()
				{
					Emissivity = emissivity[layerIndex],
					AirMass = dependent.AirMass[j],
					VaporMass = lastState.AirVapor[j],
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
			#endregion

			// Follow the solar radiation down from the top of the atmosphere to ther terrain, and absorb some as it passes through each layer
			#region Solar Radiation Absorbed
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
			#endregion

			// Calculate how much thermal radition is being emitted out of each layer
			#region Thermal Radiation
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
				SurfaceArea = dependent.IceCoverage,
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
				Emissivity = 0, /*WorldData.EmissivityWater,*/ // TODO: limit emissiveness when mass is small to surrounding air
				Energy = dependent.CloudEnergy,
				Temperature = lastState.CloudTemperature,
				SurfaceArea = dependent.CloudCoverage,
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
			#endregion


			// Thermal radiation travels upwards, partially reflecting downwards (clouds), partially absorbed, and partially lost to space
			#region Thermal Radiation Absorbed Up
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
						FromTop = false,
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
			#endregion

			// Thermal radiation is absorbed travelling downwards, collecting and then eventually hitting the earth (back radiation)
			// TODO: we need to include the top layer of atmosphere here, since we calculate cloud absorption as part of the air layer step
			#region Thermal Radiation Absorbed Down

			// transmit down from top of atmosphere			
			for (int j = _airLayer0 + _airLayers - 2; j >= 0; j--)
			{
				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[j + 1], thermalInUpJobHandlesCombined);

				if (j == _terrainLayer)
				{
					// TERRAIN
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
					// ICE
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
					// WATER
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
						FromTop = true,
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
			#endregion

			// Bouyancy, Updrafts, and mixing occur across air layers and water layers
			// TODO: add an empty air layer on top and bottom so we can calculate up/down diffusion in a single step 
			#region Vertical mixing

			JobHandle[] airVerticalMovementJobHandles = new JobHandle[_airLayers];
			for (int j = 0; j < 1; j++)
			{
				var airVerticalMovementJob = new AirVerticalMovementUpJob()
				{
					Delta = verticalMovementAir[j],
					Temperature = lastState.AirTemperature[j],
					Humidity = lastState.AirVapor[j],
					AirMass = dependent.AirMass[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					UpTemperature = lastState.AirTemperature[j + 1],
					UpHumidity = lastState.AirTemperature[j + 1],
					UpAirMass = dependent.AirMass[j + 1],
					UpLayerElevation = dependent.LayerElevation[j + 1],
					UpLayerHeight = dependent.LayerHeight[j + 1],
					SecondsPerTick = worldData.SecondsPerTick,
					MaxVerticalMovement = worldData.MaxBouyancy,
					DiffusionCoefficient = worldData.AirMassDiffusionSpeedVertical,
					Gravity = lastState.PlanetState.Gravity
				};
				airVerticalMovementJobHandles[j] = airVerticalMovementJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			for (int j = 1; j < _airLayers; j++)
			{
				var airVerticalMovementJob = new AirVerticalMovementDownJob()
				{
					Delta = verticalMovementAir[j],
					Temperature = lastState.AirTemperature[j],
					Humidity = lastState.AirVapor[j],
					AirMass = dependent.AirMass[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					UpTemperature = lastState.AirTemperature[j - 1],
					UpHumidity = lastState.AirTemperature[j - 1],
					UpAirMass = dependent.AirMass[j - 1],
					UpLayerElevation = dependent.LayerElevation[j - 1],
					UpLayerHeight = dependent.LayerHeight[j - 1],
					SecondsPerTick = worldData.SecondsPerTick,
					MaxVerticalMovement = worldData.MaxBouyancy,
					DiffusionCoefficient = worldData.AirMassDiffusionSpeedVertical,
					Gravity = lastState.PlanetState.Gravity
				};
				airVerticalMovementJobHandles[j] = airVerticalMovementJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			#endregion

			// Temperature and trace elements diffuse into neighboring horizontal cells based on a diffusion constant
			// Air, Water, Cloud
			#region Diffusion

			JobHandle[] diffusionJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				var diffusionJob = new DiffusionAirJob()
				{
					Delta = diffusionAir[j],
					Temperature = lastState.AirTemperature[j],
					Humidity = lastState.AirVapor[j],
					Wind = lastState.Wind[j],
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
					Current = lastState.WaterVelocity[j],
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
			#endregion

			#region Pressure Gradient Force

			JobHandle[] pgfJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				var pgfJob = new PressureGradientForceAirJob()
				{
					Delta = pressureGradientForce[j],
					Pressure = dependent.AirPressure[j],
					AirMass = dependent.AirMass[j],
					Temperature = lastState.AirTemperature[j],
					VaporMass = lastState.AirVapor[j],
					Neighbors = staticState.Neighbors,
					Coords = staticState.Coordinate,
					InverseCellDiameter = staticState.InverseCellDiameter,
					InverseCoordDiff = staticState.InverseCoordDiameter,
					Gravity = lastState.PlanetState.Gravity,
				};
				pgfJobHandles[_airLayer0 + j] = pgfJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}
			var windFrictionJob = new WindFrictionJob()
			{
				Friction = windFriction,
				IceCoverage = dependent.IceCoverage,
				WaterCoverage = dependent.WaterCoverage[_surfaceWaterLayer],
				VegetationCoverage = dependent.VegetationCoverage,
				Terrain = lastState.Terrain,
				IceFriction = worldData.WindIceFriction,
				TerrainFriction = worldData.WindTerrainFriction,
				VegetationFriction = worldData.WindVegetationFriction,
				WaterFriction = worldData.WindWaterFriction,
			};
			var windFrictionJobHandle = windFrictionJob.Schedule(_cellCount, _batchCount, lastJobHandle);


			#endregion

			// Wind and currents move temperature and trace elements horizontally
			// Air, Water, Cloud
			// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
			#region Advection

			JobHandle[] advectionJobHandles = new JobHandle[_layerCount];
			for (int j = 0; j < _airLayers; j++)
			{
				int layer = _airLayer0 + j;
				var advectionJob = new AdvectionAirJob()
				{
					Delta = advectionAir[j],
					Temperature = lastState.AirTemperature[j],
					Vapor = lastState.AirVapor[j],
					Wind = lastState.Wind[j],
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
					Current = lastState.WaterVelocity[j],
					Neighbors = staticState.Neighbors,
				};
				advectionJobHandles[layer] = advectionJob.Schedule(_cellCount, _batchCount, lastJobHandle);
			}

			var advectionCloudJob = new AdvectionCloudJob()
			{
				Delta = advectionCloud,
				Temperature = lastState.CloudTemperature,
				Mass = lastState.CloudMass,
				Elevation = lastState.CloudElevation,
				DropletMass = lastState.CloudDropletMass,
				Velocity = lastState.CloudVelocity,
				Neighbors = staticState.Neighbors,
				Coords = staticState.Coordinate,
				InverseCellDiameter = staticState.InverseCellDiameter,
				InverseCoordDiff = staticState.InverseCoordDiameter,
				SecondsPerTick = worldData.SecondsPerTick
			};
			advectionJobHandles[_cloudLayer] = advectionCloudJob.Schedule(_cellCount, _batchCount, lastJobHandle);


			#endregion

			// Conduction is calculated for each Surface that might touch another surface
			// Air to Cloud, Air to Ice, Air to Water, Air to Terrain, Ice to Water, Ice to Terrain, Water to Terrain
			#region Conduction
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
			var conductionAirWaterJob = new ConductionAirWaterJob()
			{
				EnergyDelta = conductionAirWater,
				TemperatureA = dependent.SurfaceAirTemperature,
				TemperatureB = lastState.WaterTemperature[_surfaceWaterLayer],
				EnergyA = dependent.IceEnergy,
				EnergyB = dependent.WaterPotentialEnergy[_surfaceWaterLayer],
				ConductionCoefficient = WorldData.ConductivityAirWater,
				CoverageIce = dependent.IceCoverage,
				CoverageWater = dependent.WaterCoverage[_surfaceWaterLayer],
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
				CoverageWater = dependent.WaterCoverage[_surfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			};
			var conductionAirTerrainJobHandle = conductionAirTerrainJob.Schedule(_cellCount, _batchCount, lastJobHandle);

			// ice to water
			var conductionIceWaterJob = new ConductionIceWaterJob()
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
				CoverageWater = dependent.WaterCoverage[_surfaceWaterLayer],
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

			JobHandle conductionCloudAirJobHandle = default(JobHandle);
			for (int i = 0; i < _airLayers; i++)
			{
				var conductionCloudAirJob = new ConductionCloudAirJob()
				{
					EnergyDelta = conductionCloudAir,

					ConductionCoefficient = worldData.AirWaterConductionPositive,
					AirTemperature = lastState.AirTemperature[i],
					ThermalDelta = thermalRadiationDelta[_cloudLayer],
					SolarIn = solarRadiationIn[_cloudLayer],
					CloudMass = lastState.CloudMass,
					CloudDropletSize = lastState.CloudDropletMass,
					CloudTemperature = lastState.CloudTemperature,
					CloudElevation = lastState.CloudElevation,
					CloudEnergy = dependent.CloudEnergy,
					CloudCoverage = dependent.CloudCoverage,
					SecondsPerTick = worldData.SecondsPerTick,
					LayerElevation = dependent.LayerElevation[i],
					LayerHeight = dependent.LayerHeight[i],
				};
				conductionCloudAirJobHandle = conductionCloudAirJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(solarInJobHandles[i], thermalInDownJobHandles[i], conductionCloudAirJobHandle));
			}
			#endregion


			#region COMBINE ADVECTION, DIFFUSION, SOLAR, THERMAL DELTA, and State changes within each layer

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
				CloudElevation = nextState.CloudElevation,
				Velocity = nextState.CloudVelocity,
				CloudEvaporationMass = cloudEvaporationMass,
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
				LastCloudElevation = lastState.CloudElevation,
				AirMassCloud = dependent.AirMassCloud,
				WaterVaporCloud = dependent.AirVaporCloud,
				AirTemperatureCloud = dependent.AirTemperatureCloud,
				AirPressureCloud = dependent.AirPressureCloud,
				RelativeHumidityCloud = dependent.AirHumidityRelativeCloud,
				Gravity = lastState.PlanetState.Gravity,
				RainDropDragCoefficient = worldData.rainDropDragCoefficient,
				RainDropMaxSize =worldData.rainDropMaxSize,
				RainDropMinSize = worldData.rainDropMinSize,
				RainfallRate = worldData.RainfallRate,
				SecondsPerTick = worldData.SecondsPerTick,
				CloudDissapationRateDryAir = worldData.CloudDissapationRateDryAir,
				CloudDissapationRateWind = worldData.CloudDissapationRateWind,
				InverseCellDiameter = staticState.InverseCellDiameter,
				SurfaceElevation = dependent.SurfaceElevation,
				TicksPerSecond = worldData.TicksPerSecond,
				WindFriction = windFriction,
				WindFrictionMultiplier= 0,
				CoriolisMultiplier = staticState.CoriolisMultiplier,
				CoriolisTerm = coriolisTerm,
				PressureGradientForce = dependent.PressureGradientAtCloudElevation,
			};
			var cloudEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_cloudLayer],
				thermalOutJobHandles[_cloudLayer],
				thermalInDownJobHandles[_cloudLayer],
				thermalInUpJobHandles[_cloudLayer],
				diffusionJobHandles[_cloudLayer],
				conductionCloudAirJobHandle,
				advectionJobHandles[_cloudLayer]
			};
			for (int j=_airLayer0;j<_airLayer0+_airLayers;j++)
			{
				cloudEnergyJobHandleDependencies.Add(solarInJobHandles[j]);
				cloudEnergyJobHandleDependencies.Add(thermalInUpJobHandles[j]);
				cloudEnergyJobHandleDependencies.Add(thermalInDownJobHandles[j]);
				cloudEnergyJobHandleDependencies.Add(windFrictionJobHandle);				
			}
			energyJobHandleDependencies.Add(cloudEnergyJobHandleDependencies);
			energyJobHandles.Add(energyCloudJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(cloudEnergyJobHandleDependencies)));

			{
				var energyJob = new EnergySurfaceAirJob()
				{
					Temperature = nextState.AirTemperature[0],
					Vapor = nextState.AirVapor[0],
					Wind = nextState.Wind[0],
					CondensationCloudMass = condensationCloudMass[0],
					CondensationGroundMass = condensationGroundMass[0],
					LastTemperature = lastState.AirTemperature[0],
					LastVapor = lastState.AirVapor[0],
					LastWind = lastState.Wind[0],
					AirMass = dependent.AirMass[0],
					Advection = advectionAir[0],
					Diffusion = diffusionAir[0],
					VerticalAirMovement = verticalMovementAir[0],
					SolarRadiationIn = solarRadiationIn[_airLayer0],
					ThermalRadiationDelta = thermalRadiationDelta[_airLayer0],
					ConductionEnergyCloud = conductionCloudAir,
					ConductionEnergyWater = conductionAirWater,
					ConductionEnergyIce = conductionAirIce,
					ConductionEnergyTerrain = conductionAirTerrain,
					CloudElevation = lastState.CloudElevation,
					LayerElevation = dependent.LayerElevation[0],
					LayerHeight = dependent.LayerHeight[0],
					PressureGradientForce = pressureGradientForce[0],
					CoriolisMultiplier = staticState.CoriolisMultiplier,
					CoriolisTerm = coriolisTerm,
					WindFriction = windFriction,
					SecondsPerTick = worldData.SecondsPerTick,
					DewPointZero = worldData.DewPointZero,
					InverseDewPointTemperatureRange = worldData.inverseDewPointTemperatureRange,
					WaterVaporMassToAirMassAtDewPoint = worldData.WaterVaporMassToAirMassAtDewPoint,
					WindFrictionMultiplier = 1.0f,
				};

				var airDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					advectionJobHandles[_airLayer0],
					diffusionJobHandles[_airLayer0],
					solarInJobHandles[_airLayer0],
					thermalOutJobHandles[_airLayer0],
					thermalInDownJobHandles[_airLayer0],
					thermalInUpJobHandles[_airLayer0],
					conductionCloudAirJobHandle,
					conductionAirWaterJobHandle,
					conductionAirIceJobHandle,
					conductionAirTerrainJobHandle,
					airVerticalMovementJobHandles[0],
					pgfJobHandles[_airLayer0],
					windFrictionJobHandle,
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
					Wind = nextState.Wind[j],	
					CondensationCloudMass = condensationCloudMass[j],
					CondensationGroundMass = condensationGroundMass[j],
					LastTemperature = lastState.AirTemperature[j],
					LastVapor = lastState.AirVapor[j],
					LastWind = lastState.Wind[j],
					AirMass = dependent.AirMass[j],
					Advection = advectionAir[j],
					Diffusion = diffusionAir[j],
					VerticalAirMovementUp = verticalMovementAir[j],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationDelta = thermalRadiationDelta[layerIndex],
					ConductionEnergyCloud = conductionCloudAir,
					CloudElevation = lastState.CloudElevation,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					PressureGradientForce = pressureGradientForce[j],
					CoriolisMultiplier = staticState.CoriolisMultiplier,
					CoriolisTerm = coriolisTerm,
					SecondsPerTick = worldData.SecondsPerTick,
					DewPointZero = worldData.DewPointZero,
					InverseDewPointTemperatureRange = worldData.inverseDewPointTemperatureRange,
					WaterVaporMassToAirMassAtDewPoint = worldData.WaterVaporMassToAirMassAtDewPoint,
					WindFriction = windFriction,
					WindFrictionMultiplier = 0
				};

				var airDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					advectionJobHandles[layerIndex],
					diffusionJobHandles[layerIndex],
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
					airVerticalMovementJobHandles[j],
					airVerticalMovementJobHandles[j-1],
					conductionCloudAirJobHandle,
					pgfJobHandles[layerIndex],
					windFrictionJobHandle,
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
					LastSurfaceWind = lastState.Wind[0],
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

			#endregion


			#region State Changes - Evaporation, Condensation, Melting, Rainfall

			var stateChangeJob = new StateChangeJob()
			{
				IceTemperature = nextState.IceTemperature,
				IceMass = nextState.IceMass,
				SurfaceWaterTemperature = nextState.WaterTemperature[_surfaceWaterLayer],
				SurfaceWaterMass = nextState.WaterMass[_surfaceWaterLayer],
				SurfaceAirTemperature = nextState.AirTemperature[0],
				SurfaceAirVapor = nextState.AirVapor[0],

				SurfaceAirMass = dependent.AirMass[0],
				SurfaceSaltMass = lastState.WaterSaltMass[_surfaceWaterLayer],
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
					CloudEvaporationMass = cloudEvaporationMass,
					CloudMass = nextState.CloudMass,
					CloudTemperature = nextState.CloudTemperature,
					SurfaceWaterMass = nextState.WaterMass[_surfaceWaterLayer],
					SurfaceWaterTemperature = nextState.WaterTemperature[_surfaceWaterLayer],

					SurfaceSaltMass = nextState.WaterSaltMass[_surfaceWaterLayer],
					CloudCondensationMass = condensationCloudMass[j],
					GroundCondensationMass = condensationGroundMass[j],
					AirMass = dependent.AirMass[j],
					LayerIndex = layerIndex,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],

				};

				stateChangeJobHandle = stateChangeAirLayerJob.Schedule(_cellCount, _batchCount, stateChangeJobHandle);
			}
			#endregion


			#region Update dependent variables
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
			JobHandle updateDependentAirLayerJobHandle = stateChangeJobHandle;
			for (int j = 0; j < _airLayers; j++)
			{
				var updateDependentAirLayerJob = new UpdateDependentAirLayerJob()
				{
					Pressure = dependent.AirPressure[j],
					RelativeHumidity = dependent.AirHumidityRelative[j],
					AbsoluteHumidity = dependent.AirHumidityAbsolute[j],
					AirMass = dependent.AirMass[j],
					PotentialEnergy = dependent.AirPotentialEnergy[j],
					PressureGradientAtCloudElevation = dependent.PressureGradientAtCloudElevation,
					AirMassCloud = dependent.AirMassCloud,
					AirVaporCloud = dependent.AirVaporCloud,
					AirTemperatureCloud = dependent.AirTemperatureCloud,
					AirPressureCloud = dependent.AirPressureCloud,
					AirHumidityRelativeCloud = dependent.AirHumidityRelativeCloud,
					AirLayerCloud = dependent.AirLayerCloud,
					PressureGradient = pressureGradientForce[j],

					AirTemperature = lastState.AirTemperature[j],
					CloudDropletMass = nextState.CloudDropletMass,
					CloudElevation= nextState.CloudElevation,
					CloudMass = nextState.CloudMass,
					CloudTemperature = nextState.CloudTemperature,
					VaporMass = nextState.AirVapor[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					IceMass = nextState.IceMass,
					Gravity = nextState.PlanetState.Gravity,
					DewPointZero = worldData.DewPointZero,
					InverseDewPointTemperatureRange = worldData.inverseDewPointTemperatureRange,
					WaterVaporMassToAirMassAtDewPoint = worldData.WaterVaporMassToAirMassAtDewPoint,
					LayerIndex = j,
				};
				updateDependentAirLayerJobHandle = updateDependentAirLayerJob.Schedule(_cellCount, _batchCount, updateDependentAirLayerJobHandle);
				updateDependenciesJobHandles.Add(updateDependentAirLayerJobHandle);
				if (j == 0)
				{
					surfaceAirJobHandle = updateDependentAirLayerJobHandle;
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
				CloudEnergy = dependent.CloudEnergy,
				IceCoverage = dependent.IceCoverage,
				IceEnergy = dependent.IceEnergy,
				SurfaceElevation = dependent.SurfaceElevation,
				VegetationCoverage = dependent.VegetationCoverage,
				WaterDepth = dependent.WaterDepth,
				SurfaceAirTemperature = dependent.SurfaceAirTemperature,

				CloudMass = nextState.CloudMass,
				CloudTemperature = nextState.CloudTemperature,
				IceMass = nextState.IceMass,
				IceTemperature = nextState.IceTemperature,
				Terrain = nextState.Terrain,
				WaterSaltMass = waterSaltMass,
				worldData = worldData,
				LowerAirTemperature = lastState.AirTemperature[0],
				lowerAirHeight = dependent.LayerHeight[0],
			};

			var updateDependenciesJobHandle = updateDependentStateJob.Schedule(_cellCount, _batchCount, JobHandle.CombineDependencies(surfaceAirJobHandle, summationHandle));
			updateDependenciesJobHandles.Add(updateDependenciesJobHandle);
			lastJobHandle = JobHandle.CombineDependencies(updateDependenciesJobHandles);

			#endregion

			lastJobHandle.Complete();

			#region Debug
			if (checkForDegeneracy)
			{
				bool degen = false;
				SortedSet<int> degenIndices = new SortedSet<int>();
				List<string> degenVarNames = new List<string>();
				degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "TerrainTemperature", nextState.TerrainTemperature, 0, 400, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "CloudDropletMass", nextState.CloudDropletMass, degenVarNames);
				degen |= CheckDegen(_cellCount, degenIndices, "CloudElevation", nextState.CloudElevation, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "CloudMass", nextState.CloudMass, degenVarNames);
				degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "CloudTemperature", nextState.CloudTemperature, 0, 400, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "IceMass", nextState.IceMass, degenVarNames);
				degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "IceTemperature", nextState.IceTemperature, 0, 400, degenVarNames);
				for (int i = 0; i < _airLayers; i++) {
					degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "AirTemperature" + i, nextState.AirTemperature[i], 0, 400, degenVarNames);
					degen |= CheckDegenPosValues(_cellCount, degenIndices, "AirVapor" + i, nextState.AirVapor[i], degenVarNames);
				}
				for (int i=0;i<_waterLayers;i++)
				{
					degen |= CheckDegenPosValues(_cellCount, degenIndices, "WaterMass" + i, nextState.WaterMass[i], degenVarNames);
					degen |= CheckDegenPosValues(_cellCount, degenIndices, "WaterSaltMass" + i, nextState.WaterSaltMass[i], degenVarNames);
					degen |= CheckDegenPosValues(_cellCount, degenIndices, "WaterTemperature" + i, nextState.WaterTemperature[i], degenVarNames);
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

			#endregion

			#region Dispose Temporary Arrays
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
			cloudEvaporationMass.Dispose();
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

	private static bool CheckDegenPosValues(int count, SortedSet<int> degenIndices, string name, NativeArray<float> values, List<string> degenVarNames)
	{
		for (int i = 0; i < count; i++)
		{
			float v = values[i];
			if (!math.isfinite(v) || v < 0)
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
			float v = values[i];
			if (!math.isfinite(v) || v < min || v > max)
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
		s.Append(title + " Index: " + i + " Time: " + state.PlanetState.Ticks);
		foreach (var n in degenVarNames)
		{
			s.Append(" | " + n);
		}
		s.AppendLine("");
		s.AppendLine("X: " + staticState.Coordinate[i].x + " Y: " + staticState.Coordinate[i].y);
		s.AppendLine("Elevation" + ": " + state.Terrain[i].Elevation);
		s.AppendLine("Roughness" + ": " + state.Terrain[i].Roughness);
		s.AppendLine("SoilFertility" + ": " + state.Terrain[i].SoilFertility);
		s.AppendLine("Vegetation" + ": " + state.Terrain[i].Vegetation);
		s.AppendLine("TerrainTemperature" + ": " + state.TerrainTemperature[i]);
		s.AppendLine("CloudMass" + ": " + state.CloudMass[i]);
		s.AppendLine("CloudElevation" + ": " + state.CloudElevation[i]);
		s.AppendLine("CloudTemperature" + ": " + state.CloudTemperature[i]);
		s.AppendLine("CloudDropletMass" + ": " + state.CloudDropletMass[i]);
		s.AppendLine("IceMass" + ": " + state.IceMass[i]);
		s.AppendLine("IceTemperature" + ": " + state.IceTemperature[i]);
		for (int j = 0; j < _waterLayers; j++)
		{
			s.AppendLine("WaterMass" + j + ": " + state.WaterMass[j][i]);
			s.AppendLine("WaterSaltMass" + j + ": " + state.WaterSaltMass[j][i]);
			s.AppendLine("WaterTemperature" + j + ": " + state.WaterTemperature[j][i]);
			s.AppendLine("WaterVelocity" + j + ": " + state.WaterVelocity[j][i]);
		}
		for (int j = 0; j < _airLayers; j++)
		{
			s.AppendLine("AirTemperature" + j + ": " + state.AirTemperature[j][i]);
			s.AppendLine("AirVapor" + j + ": " + state.AirVapor[j][i]);
			s.AppendLine("Wind" + j + ": " + state.Wind[j][i]);
		}
		Debug.Log(s);
	}
	public void PrintDependentState(string title, int i, DependentState dependent)
	{
		StringBuilder s = new StringBuilder();
		s.AppendLine(title + " Index: " + i);
		s.AppendLine("Surface Elevation" + ": " + dependent.SurfaceElevation[i]);
		s.AppendLine("Water Depth" + ": " + dependent.WaterDepth[i]);
		s.AppendLine("Ice Coverage" + ": " + dependent.IceCoverage[i]);
		s.AppendLine("Cloud Coverage" + ": " + dependent.CloudCoverage[i]);
		for (int j = 0; j < _waterLayers; j++)
		{
			s.AppendLine("Water Coverage" + ": " + dependent.WaterCoverage[j][i]);
		}
		for (int j = 0; j < _airLayers; j++)
		{
		}
		Debug.Log(s);
	}
}
