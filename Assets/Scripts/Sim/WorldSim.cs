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

	public JobHelper SimJob;
	public JobHelper NeighborJob;


	private int _cellCount;
	private int _batchCount = 100;
	private const int _terrainLayers = 4;

	private int _airLayers;
	private int _waterLayers;
	private int _layerCount;
	private int _terrainLayer;
	private int _cloudLayer;
	private int _waterLayer0;
	private int _floraLayer;
	private int _iceLayer;
	private int _airLayer0;
	private int _surfaceWaterLayer;

	private NativeArray<float> solarRadiation;
	private NativeArray<float> waterSlopeAlbedo;
	private NativeArray<float> cloudAlbedo;
	private NativeArray<SolarAbsorptivity>[] absorptivitySolar;
	private NativeArray<ThermalAbsorptivity>[] absorptivityThermal;
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
	private NativeArray<float3>[] airAcceleration;
	private NativeArray<float> terrainGradient;
	private NativeArray<float> windFriction;
	private NativeArray<float3> waterFriction;
	private NativeArray<float> conductionAirIce;
	private NativeArray<float> conductionAirWater;
	private NativeArray<float> conductionAirFlora;
	private NativeArray<float> conductionAirTerrain;
	private NativeArray<float> conductionIceWater;
	private NativeArray<float> conductionIceFlora;
	private NativeArray<float> conductionIceTerrain;
	private NativeArray<float> conductionFloraTerrain;
	private NativeArray<float>[] conductionWaterTerrain;
	private NativeArray<float>[] conductionWaterFlora;
	private NativeArray<float> frozenTemperature;
	private NativeArray<float> frozenMass;
	private NativeArray<float> saltPlume;
	private NativeArray<float> evaporationMassWater;
	private NativeArray<float> evaporationTemperaturePotentialWater;
	private NativeArray<float> evaporationMassFlora;
	private NativeArray<float> evaporationTemperaturePotentialFlora;
	private NativeArray<float> groundWaterConsumed;
	private NativeArray<float> floraMassDelta;
	private NativeArray<float> geothermalRadiation;
	private NativeArray<float> groundWaterFlowMass;
	private NativeArray<float> groundWaterFlowTemperature;

	private NativeArray<float> displaySolarRadiation;

	public WorldSim(int cellCount, int airLayers, int waterLayers)
	{
		_cellCount = cellCount;
		_airLayers = airLayers;
		_waterLayers = waterLayers;
		_layerCount = _airLayers + _waterLayers + _terrainLayers;
		_terrainLayer = 0;
		_floraLayer = _terrainLayer + 1;
		_waterLayer0 = _floraLayer + 1;
		_iceLayer = _waterLayer0 + _waterLayers;
		_airLayer0 = _iceLayer + 1;
		_surfaceWaterLayer = _waterLayers - 2;
		_cloudLayer = _airLayer0 + _airLayers;

		SimJob = new JobHelper(_cellCount);
		NeighborJob = new JobHelper(_cellCount * 6);

		solarRadiation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		waterSlopeAlbedo = new NativeArray<float>(_cellCount, Allocator.Persistent);
		cloudAlbedo = new NativeArray<float>(_cellCount, Allocator.Persistent);
		emissivity = new NativeArray<float>[_layerCount];
		solarRadiationIn = new NativeArray<float>[_layerCount];
		for (int i=0;i<_layerCount;i++)
		{
			emissivity[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			solarRadiationIn[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
		}
		diffusionAir = new NativeArray<DiffusionAir>[_airLayers];
		advectionAir = new NativeArray<DiffusionAir>[_airLayers];
		destinationAir = new NativeArray<BarycentricValue>[_airLayers];
		airAcceleration = new NativeArray<float3>[_airLayers];
		absorptivitySolar = new NativeArray<SolarAbsorptivity>[_airLayers];
		absorptivityThermal = new NativeArray<ThermalAbsorptivity>[_airLayers];
		for (int i = 0; i < _airLayers; i++)
		{
			diffusionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
			advectionAir[i] = new NativeArray<DiffusionAir>(_cellCount, Allocator.Persistent);
			destinationAir[i] = new NativeArray<BarycentricValue>(_cellCount, Allocator.Persistent);
			airAcceleration[i] = new NativeArray<float3>(_cellCount, Allocator.Persistent);
			absorptivitySolar[i] = new NativeArray<SolarAbsorptivity>(_cellCount, Allocator.Persistent);
			absorptivityThermal[i] = new NativeArray<ThermalAbsorptivity>(_cellCount, Allocator.Persistent);
		}
		diffusionWater = new NativeArray<DiffusionWater>[_waterLayers];
		advectionWater = new NativeArray<DiffusionWater>[_waterLayers];
		destinationWater = new NativeArray<BarycentricValue>[_waterLayers];
		conductionWaterFlora = new NativeArray<float>[_waterLayers];
		conductionWaterTerrain = new NativeArray<float>[_waterLayers];
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
			advectionWater[i] = new NativeArray<DiffusionWater>(_cellCount, Allocator.Persistent);
			destinationWater[i] = new NativeArray<BarycentricValue>(_cellCount, Allocator.Persistent);
			conductionWaterFlora[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
			conductionWaterTerrain[i] = new NativeArray<float>(_cellCount, Allocator.Persistent);
		}
		frozenTemperature = new NativeArray<float>(_cellCount, Allocator.Persistent);
		frozenMass = new NativeArray<float>(_cellCount, Allocator.Persistent);
		saltPlume = new NativeArray<float>(_cellCount, Allocator.Persistent);
		evaporationMassWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		evaporationTemperaturePotentialWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		evaporationMassFlora = new NativeArray<float>(_cellCount, Allocator.Persistent);
		evaporationTemperaturePotentialFlora = new NativeArray<float>(_cellCount, Allocator.Persistent);
		groundWaterConsumed = new NativeArray<float>(_cellCount, Allocator.Persistent);
		floraMassDelta = new NativeArray<float>(_cellCount, Allocator.Persistent);
		windFriction = new NativeArray<float>(_cellCount, Allocator.Persistent);
		waterFriction = new NativeArray<float3>(_cellCount, Allocator.Persistent);
		diffusionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		advectionCloud = new NativeArray<DiffusionCloud>(_cellCount, Allocator.Persistent);
		destinationCloud = new NativeArray<BarycentricValue>(_cellCount, Allocator.Persistent);
		conductionAirIce = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirFlora = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionAirTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionIceWater = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionIceFlora = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionIceTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		conductionFloraTerrain = new NativeArray<float>(_cellCount, Allocator.Persistent);
		displaySolarRadiation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		geothermalRadiation = new NativeArray<float>(_cellCount, Allocator.Persistent);
		terrainGradient = new NativeArray<float>(_cellCount * 6, Allocator.Persistent);
		groundWaterFlowMass = new NativeArray<float>(_cellCount, Allocator.Persistent);
		groundWaterFlowTemperature = new NativeArray<float>(_cellCount, Allocator.Persistent);

}

public void Dispose()
	{
		solarRadiation.Dispose();
		waterSlopeAlbedo.Dispose();
		cloudAlbedo.Dispose();
		for (int i=0;i<_layerCount;i++)
		{
			emissivity[i].Dispose();
			solarRadiationIn[i].Dispose();
		}
		for (int i = 0; i < _airLayers; i++)
		{
			diffusionAir[i].Dispose();
			advectionAir[i].Dispose();
			destinationAir[i].Dispose();
			airAcceleration[i].Dispose();
			absorptivitySolar[i].Dispose();
			absorptivityThermal[i].Dispose();
		}
		for (int i = 0; i < _waterLayers; i++)
		{
			diffusionWater[i].Dispose();
			advectionWater[i].Dispose();
			destinationWater[i].Dispose();
			conductionWaterFlora[i].Dispose();
			conductionWaterTerrain[i].Dispose();
		}
		frozenTemperature.Dispose();
		frozenMass.Dispose();
		saltPlume.Dispose();
		evaporationMassWater.Dispose();
		evaporationTemperaturePotentialWater.Dispose();
		evaporationMassFlora.Dispose();
		evaporationTemperaturePotentialFlora.Dispose();
		groundWaterConsumed.Dispose();
		floraMassDelta.Dispose();
		windFriction.Dispose();
		waterFriction.Dispose();
		diffusionCloud.Dispose();
		advectionCloud.Dispose();
		destinationCloud.Dispose();
		conductionAirIce.Dispose();
		conductionAirWater.Dispose();
		conductionAirFlora.Dispose();
		conductionAirTerrain.Dispose();
		conductionIceWater.Dispose();
		conductionIceFlora.Dispose();
		conductionIceTerrain.Dispose();
		conductionFloraTerrain.Dispose();
		geothermalRadiation.Dispose();
		terrainGradient.Dispose();
		groundWaterFlowMass.Dispose();
		groundWaterFlowTemperature.Dispose();



		displaySolarRadiation.Dispose();

	}

	public void Tick(SimState[] states, int stateCount, int ticksToAdvance, ref DependentState dependent, ref DisplayState display, ref StaticState staticState, ref WorldData worldData, ref int curStateIndex, bool checkForDegeneracy, bool logState, int logStateIndex, bool displayGlobals)
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
			var tempArrays = new List<NativeArray<float>>();

			var thermalRadiationDelta = new NativeArray<float>[_layerCount];
			var thermalRadiationTransmittedUp = new NativeArray<float>[_layerCount];
			var thermalRadiationTransmittedDown = new NativeArray<float>[_layerCount];
			var windowRadiationTransmittedUp = new NativeArray<float>[_layerCount];
			var windowRadiationTransmittedDown = new NativeArray<float>[_layerCount];
			var condensationGroundMass = new NativeArray<float>[_layerCount];
			var condensationCloudMass = new NativeArray<float>[_layerCount];
			var solarReflected = new NativeArray<float>[_layerCount];
			var latentHeat = new NativeArray<float>[_layerCount];
			var conductionWaterTerrainTotal = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var conductionWaterFloraTotal = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var cloudEvaporationMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var dropletDelta = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var precipitationMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var precipitationTemperature = new NativeArray<float>(_cellCount, Allocator.TempJob);
			var iceMeltedMass = new NativeArray<float>(_cellCount, Allocator.TempJob);
			for (int i = 0; i < _layerCount; i++)
			{
				latentHeat[i] = new NativeArray<float>(_cellCount, Allocator.TempJob);
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
			var solarInJobHandle =SimJob.Schedule(new SolarRadiationJob()
			{
				SolarRadiation = solarRadiation,
				GeothermalRadiation = geothermalRadiation,
				DisplaySolarRadiation = displaySolarRadiation,
				WaterSlopeAlbedo = waterSlopeAlbedo,

				SphericalPosition = staticState.SphericalPosition,
				IncomingSolarRadiation = lastState.PlanetState.SolarRadiation * worldData.SecondsPerTick,
				IncomingGeothermalRadiation = lastState.PlanetState.GeothermalHeat * worldData.SecondsPerTick,
				PlanetRotation = quaternion.Euler(lastState.PlanetState.Rotation),
				SunToPlanetDir = math.normalize(lastState.PlanetState.Position),
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
				emissivityJobHandles[layerIndex] =SimJob.Schedule(new EmissivityAirJob()
				{
					Emissivity = emissivity[layerIndex],
					AirMass = dependent.AirMass[j],
					VaporMass = lastState.AirVapor[j],
					CarbonDioxide = lastState.PlanetState.CarbonDioxide,
					EmissivityAir = worldData.ThermalEmissivityAir,
					EmissivityWaterVapor = worldData.ThermalEmissivityWaterVapor,
					EmissivityCarbonDioxide = worldData.ThermalEmissivityCarbonDioxide
				});
			}

			// we only do thermal radiation upwards for the surface layer of water,
			// for the bottom we rely on conduction with the terrain for heat transfer (although this might lead to an imbalance!)
			{
				emissivityJobHandles[_surfaceWaterLayer + _waterLayer0] =SimJob.Schedule(new EmissivityWaterJob()
				{
					Emissivity = emissivity[_waterLayer0 + _surfaceWaterLayer],
					WaterMass = lastState.WaterMass[_surfaceWaterLayer],
					SaltMass = lastState.SaltMass[_surfaceWaterLayer],
					EmissivitySalt = worldData.ThermalEmissivitySalt,
					EmissivityWater = worldData.ThermalEmissivityWater
				});
			}
			emissivityJobHandles[_terrainLayer] =SimJob.Schedule(new EmissivityTerrainJob()
			{
				Emissivity = emissivity[_terrainLayer],
				SoilFertility = lastState.SoilFertility,
				EmissivityDirt = worldData.ThermalEmissivityDirt,
				EmissivitySand = worldData.ThermalEmissivitySand,
			});

			#endregion

			// Calculate how much thermal radition is being emitted out of each layer
			#region Thermal Radiation
			JobHandle[] thermalOutJobHandles = new JobHandle[_layerCount];

			// ICE
			thermalOutJobHandles[_iceLayer] = SimJob.Schedule(new ThermalEnergyRadiatedConstantEmissivityJob()
			{
				ThermalRadiationDelta = thermalRadiationDelta[_iceLayer],
				ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[_iceLayer],
				ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[_iceLayer],
				WindowRadiationTransmittedUp = windowRadiationTransmittedUp[_iceLayer],
				WindowRadiationTransmittedDown = windowRadiationTransmittedDown[_iceLayer],

				PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
				Emissivity = worldData.ThermalEmissivityIce,
				Energy = dependent.IceEnergy,
				Temperature = lastState.IceTemperature,
				SurfaceArea = dependent.IceCoverage,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);


			// FLORA
			thermalOutJobHandles[_floraLayer] = SimJob.Schedule(new ThermalEnergyRadiatedConstantEmissivityJob()
			{
				ThermalRadiationDelta = thermalRadiationDelta[_floraLayer],
				ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[_floraLayer],
				ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[_floraLayer],
				WindowRadiationTransmittedUp = windowRadiationTransmittedUp[_floraLayer],
				WindowRadiationTransmittedDown = windowRadiationTransmittedDown[_floraLayer],

				PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
				Emissivity = worldData.ThermalEmissivityFlora,
				Energy = dependent.FloraEnergy,
				Temperature = lastState.FloraTemperature,
				SurfaceArea = dependent.FloraCoverage,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);



			// TERRAIN
			thermalOutJobHandles[_terrainLayer] =SimJob.Schedule(new ThermalEnergyRadiatedTerrainJob()
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
			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layer = _airLayer0 + j;
				thermalOutJobHandles[layer] =SimJob.Schedule(new ThermalEnergyRadiatedAirJob()
				{
					ThermalRadiationDelta = thermalRadiationDelta[layer],
					ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[layer],
					ThermalRadiationTransmittedDown = thermalRadiationTransmittedDown[layer],
					WindowRadiationTransmittedUp = windowRadiationTransmittedUp[layer],
					WindowRadiationTransmittedDown = windowRadiationTransmittedDown[layer],

					PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
					Energy = dependent.AirPotentialEnergy[j],
					Emissivity = emissivity[layer],
					TemperaturePotential = lastState.AirTemperaturePotential[j],
					LayerMiddle = dependent.LayerMiddle[j],
					SecondsPerTick = worldData.SecondsPerTick
				}, emissivityJobHandles[layer]);
			}

			// WATER
			// we only do thermal radiation upwards for the surface layer of water,
			// for the bottom we rely on conduction with the terrain for heat transfer (although this might lead to an imbalance!)
			{
				int layer = _surfaceWaterLayer + _waterLayer0;
				thermalOutJobHandles[layer] =SimJob.Schedule(new ThermalEnergyRadiatedWaterJob()
				{
					ThermalRadiationDelta = thermalRadiationDelta[layer],
					ThermalRadiationTransmittedUp = thermalRadiationTransmittedUp[layer],
					WindowRadiationTransmittedUp = windowRadiationTransmittedUp[layer],

					PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
					Emissivity = emissivity[layer],
					Energy = dependent.WaterPotentialEnergy[_surfaceWaterLayer],
					TemperatureAbsolute = lastState.WaterTemperature[_surfaceWaterLayer],
					SurfaceArea = dependent.WaterCoverage[_surfaceWaterLayer],
					SecondsPerTick = worldData.SecondsPerTick
				}, emissivityJobHandles[layer]);
			}
			#endregion


			#region absorptivity

			var cloudAlbedoJobHandle =SimJob.Schedule(new CloudAlbedoJob()
			{
				CloudAlbedo = cloudAlbedo,
				CloudMass = lastState.CloudMass,
				CloudDropletMass = lastState.CloudDropletMass,
				DewPoint = dependent.DewPoint,
				CloudElevation = dependent.CloudElevation,
				WaterSlopeAlbedo = waterSlopeAlbedo,
				AlbedoCloud = worldData.AlbedoCloud,
				CloudFreezingTemperatureMax = worldData.maxCloudFreezingTemperature,
				CloudFreezingTemperatureMin = worldData.minCloudFreezingTemperature,
				CloudSlopeAlbedoMax = worldData.maxCloudSlopeAlbedo,
				RainDropSizeAlbedoMax = worldData.rainDropSizeAlbedoMax,
				RainDropSizeAlbedoMin = worldData.rainDropSizeAlbedoMin,
			}, solarInJobHandle);

			var absorptivityAirJobHandles = new JobHandle[_airLayers];
			for (int j = _airLayers - 2; j > 0; j--)
			{
				absorptivityAirJobHandles[j] =SimJob.Schedule(new AbsorptivityAirJob()
				{
					AbsorptivitySolar = absorptivitySolar[j],
					AbsorptivityThermal = absorptivityThermal[j],
					EmissivityAir = emissivity[_airLayer0 + j],
					AirMass = dependent.AirMass[j],
					VaporMass = lastState.AirVapor[j],
					CloudMass = lastState.CloudMass,
					CloudAlbedo = cloudAlbedo,
					CloudElevation = dependent.CloudElevation,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					CarbonDioxide = lastState.PlanetState.CarbonDioxide,
					AlbedoAir = worldData.AlbedoAir,
					SolarAbsorptivityAir = worldData.SolarAbsorptivityAir,
					SolarAbsorptivityCloud = worldData.SolarAbsorptivityCloud,
					SolarAbsorptivityWaterVapor = worldData.SolarAbsorptivityWaterVapor,
					ThermalAbsorptivityAir = worldData.ThermalAbsorptivityAir,
					ThermalAbsorptivityCloud = worldData.ThermalAbsorptivityCloud,
					EmissivityCloud = worldData.ThermalEmissivityWater,
				}, JobHandle.CombineDependencies(cloudAlbedoJobHandle, emissivityJobHandles[_airLayer0 + j]));
			}

			#endregion


			// Follow the solar radiation down from the top of the atmosphere to ther terrain, and absorb some as it passes through each layer
			#region Solar Radiation Absorbed

			// process each vertical layer in order

			// atmosphere
			JobHandle[] solarInJobHandles = new JobHandle[_layerCount];
			for (int j = _airLayer0 + _airLayers - 2; j > _airLayer0; j--)
			{
				int airLayerIndex = j - _airLayer0;
				solarInJobHandles[j] = solarInJobHandle =SimJob.Schedule(new SolarRadiationAbsorbedAirJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[j],
					SolarRadiationIncoming = solarRadiation,
					SolarRadiationReflected = solarReflected[j],
					SolarRadiationAbsorbedCloud = solarRadiationIn[_cloudLayer],
					SolarRadiationReflectedCloud = solarReflected[_cloudLayer],
					AbsorptivitySolar = absorptivitySolar[airLayerIndex],
				}, JobHandle.CombineDependencies(solarInJobHandle, absorptivityAirJobHandles[airLayerIndex]));
			}


			// ice
			solarInJobHandles[_iceLayer] = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedPartialCoverageConstantAlbedoJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_iceLayer],
				SolarRadiationIncoming = solarRadiation,
				SolarRadiationReflected = solarReflected[_iceLayer],
				Albedo = WorldData.AlbedoIce,
				Coverage = dependent.IceCoverage
			}, solarInJobHandle);

			for (int j = _surfaceWaterLayer; j >= 1; j--)
			{
				int layerIndex = _waterLayer0 + j;
				solarInJobHandles[layerIndex] = solarInJobHandle =SimJob.Schedule(new SolarRadiationAbsorbedPartialCoverageJob()
				{
					SolarRadiationAbsorbed = solarRadiationIn[layerIndex],
					SolarRadiationIncoming = solarRadiation,
					SolarRadiationReflected = solarReflected[layerIndex],
					Coverage = dependent.WaterCoverage[layerIndex - _waterLayer0],
					Albedo = waterSlopeAlbedo,
				}, solarInJobHandle);
			}

			// flora
			solarInJobHandles[_floraLayer] = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedPartialCoverageConstantAlbedoJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_floraLayer],
				SolarRadiationIncoming = solarRadiation,
				SolarRadiationReflected = solarReflected[_floraLayer],
				Albedo = WorldData.AlbedoFlora,
				Coverage = dependent.FloraCoverage
			}, solarInJobHandle);


			solarInJobHandles[_terrainLayer] = solarInJobHandle =SimJob.Schedule(new SolarRadiationAbsorbedTerrainJob()
			{
				SolarRadiationAbsorbed = solarRadiationIn[_terrainLayer],
				SolarRadiationIncoming = solarRadiation,
				SolarRadiationReflected = solarReflected[_terrainLayer],
				worldData = worldData,
				SoilFertility = lastState.SoilFertility,
			}, solarInJobHandle);
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
					thermalInUpJobHandles[j] =SimJob.Schedule(new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[downIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[downIndex],
						AbsorptivityThermal = absorptivityThermal[airLayer],
						LayerElevation = dependent.LayerElevation[airLayer],
						LayerHeight = dependent.LayerHeight[airLayer],
						CloudElevation = dependent.CloudElevation,
						FromTop = false,
					}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex], absorptivityAirJobHandles[airLayer]));
				}
				else if (j == _iceLayer)
				{
					int downIndex = _waterLayer0 + _surfaceWaterLayer;
					thermalInUpJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[downIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[downIndex],
						Coverage = dependent.IceCoverage,

					}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex]));
				}
				else if (j == _floraLayer)
				{
					int downIndex = _terrainLayer;
					thermalInUpJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[downIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[downIndex],
						Coverage = dependent.FloraCoverage,

					}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex], thermalOutJobHandles[downIndex]));
				}
				else if (j > _waterLayer0 && j < _waterLayer0 + _waterLayers - 1)
				{
					int waterLayerIndex = j - _waterLayer0;
					int downIndex = (waterLayerIndex == 1) ? _floraLayer : (j - 1);
					thermalInUpJobHandles[j] =SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedUp[j],
						WindowRadiationTransmitted = windowRadiationTransmittedUp[j],

						WindowRadiationIncoming = windowRadiationTransmittedUp[downIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedUp[downIndex],
						Coverage = dependent.WaterCoverage[waterLayerIndex],
					}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex]));
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
			for (int j = _layerCount - 1; j >= 0; j--)
			{

				if (j == _terrainLayer)
				{
					// TERRAIN
					int upIndex = _floraLayer;
					var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
					thermalInDownJobHandles[j] =SimJob.Schedule(new ThermalEnergyAbsorbedTerrainJob()
					{
						ThermalRadiationAbsorbed = thermalRadiationDelta[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[upIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[upIndex],
					}, thermalInDependenciesHandle);
				}
				else if (j == _floraLayer)
				{
					// FLORA
					int upIndex = _waterLayer0 + _surfaceWaterLayer;
					var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
					thermalInDownJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[upIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[upIndex],
						Coverage = dependent.FloraCoverage,
					}, thermalInDependenciesHandle);
				}
				else if (j == _iceLayer)
				{
					// ICE
					int upIndex = _airLayer0 + 1;
					var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
					thermalInDownJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[upIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[upIndex],
						Coverage = dependent.IceCoverage,
					}, thermalInDependenciesHandle);
				}
				else if (j == _waterLayer0 + _surfaceWaterLayer)
				{
					// WATER
					int waterLayerIndex = j - _waterLayer0;
					int upIndex = (waterLayerIndex == _surfaceWaterLayer) ? _iceLayer : (j + 1);
					var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
					thermalInDownJobHandles[j] =SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
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
					thermalInDownJobHandles[j] =SimJob.Schedule(new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = thermalRadiationDelta[j],
						ThermalRadiationTransmitted = thermalRadiationTransmittedDown[j],
						WindowRadiationTransmitted = windowRadiationTransmittedDown[j],

						WindowRadiationIncoming = windowRadiationTransmittedDown[upIndex],
						ThermalRadiationIncoming = thermalRadiationTransmittedDown[upIndex],
						AbsorptivityThermal = absorptivityThermal[airLayer],
						LayerElevation = dependent.LayerElevation[airLayer],
						LayerHeight = dependent.LayerHeight[airLayer],
						CloudElevation = dependent.CloudElevation,
						FromTop = true,
					}, thermalInDependenciesHandle);
				}
			}
			#endregion

			// Conduction is calculated for each Surface that might touch another surface
			// Air to Cloud, Air to Ice, Air to Water, Air to Terrain, Ice to Water, Ice to Terrain, Water to Terrain
			#region Conduction
			// air to ice
			var conductionAirIceJobHandle =SimJob.Schedule(new ConductionBJob()
			{
				EnergyDelta = conductionAirIce,
				TemperatureA = dependent.SurfaceAirTemperatureAbsolute,
				TemperatureB = lastState.IceTemperature,
				EnergyB = dependent.IceEnergy,
				ConductionCoefficient = WorldData.ConductivityAirIce,
				SurfaceArea = dependent.SurfaceAreaAirIce,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// air to water
			var conductionAirWaterJobHandle =SimJob.Schedule(new ConductionBJob()
			{
				EnergyDelta = conductionAirWater,
				TemperatureA = dependent.SurfaceAirTemperatureAbsolute,
				TemperatureB = lastState.WaterTemperature[_surfaceWaterLayer],
				EnergyB = dependent.WaterPotentialEnergy[_surfaceWaterLayer],
				ConductionCoefficient = WorldData.ConductivityAirWater,
				SurfaceArea = dependent.SurfaceAreaAirWater,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// air to flora
			var conductionAirFloraJobHandle = SimJob.Schedule(new ConductionBJob()
			{
				EnergyDelta = conductionAirFlora,
				TemperatureA = dependent.SurfaceAirTemperatureAbsolute,
				TemperatureB = lastState.FloraTemperature,
				ConductionCoefficient = WorldData.ConductivityAirFlora,
				SurfaceArea = dependent.SurfaceAreaAirFlora,
				EnergyB = dependent.FloraEnergy,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// air to terrain
			var conductionAirTerrainJobHandle = SimJob.Schedule(new ConductionJob()
			{
				EnergyDelta = conductionAirTerrain,
				TemperatureA = dependent.SurfaceAirTemperatureAbsolute,
				TemperatureB = lastState.TerrainTemperature,
				ConductionCoefficient = WorldData.ConductivityAirTerrain,
				SurfaceArea = dependent.SurfaceAreaAirTerrain,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// ice to water
			var conductionIceWaterJobHandle = SimJob.Schedule(new ConductionABJob()
			{
				EnergyDelta = conductionIceWater,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.WaterTemperature[_surfaceWaterLayer],
				EnergyA = dependent.IceEnergy,
				EnergyB = dependent.WaterPotentialEnergy[_surfaceWaterLayer],
				ConductionCoefficient = WorldData.ConductivityIceWater,
				SurfaceArea = dependent.SurfaceAreaIceWater,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// ice to flora
			var conductionIceFloraJobHandle = SimJob.Schedule(new ConductionABJob()
			{
				EnergyDelta = conductionIceFlora,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.FloraTemperature,
				EnergyA = dependent.IceEnergy,
				EnergyB = dependent.FloraEnergy,
				ConductionCoefficient = WorldData.ConductivityIceFlora,
				SurfaceArea = dependent.SurfaceAreaIceFlora,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// ice to terrain
			var conductionIceTerrainJobHandle =SimJob.Schedule(new ConductionAJob()
			{
				EnergyDelta = conductionIceTerrain,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.TerrainTemperature,
				EnergyA = dependent.IceEnergy,
				ConductionCoefficient = WorldData.ConductivityIceTerrain,
				SurfaceArea = dependent.SurfaceAreaIceTerrain,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);

			// flora to terrain
			var conductionFloraTerrainJobHandle = SimJob.Schedule(new ConductionBJob()
			{
				EnergyDelta = conductionFloraTerrain,
				TemperatureA = lastState.FloraTemperature,
				TemperatureB = lastState.TerrainTemperature,
				EnergyB = dependent.FloraEnergy,
				ConductionCoefficient = WorldData.ConductivityFloraTerrain,
				SurfaceArea = dependent.SurfaceAreaFloraTerrain,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle);


			// water to terrain
			JobHandle conductionWaterTerrainJobHandle = default(JobHandle);
			JobHandle conductionWaterFloraJobHandle = default(JobHandle);
			for (int i = 1; i < _waterLayers-1; i++) {
				conductionWaterFloraJobHandle = JobHandle.CombineDependencies(conductionWaterTerrainJobHandle, SimJob.Schedule(new ConductionWaterBottomABJob()
				{
					EnergyDelta = conductionWaterFlora[i],
					EnergyDeltaTotal = conductionWaterFloraTotal,
					TemperatureA = lastState.WaterTemperature[i],
					TemperatureB = lastState.FloraTemperature,
					EnergyA = dependent.WaterPotentialEnergy[i],
					EnergyB = dependent.FloraEnergy,
					ConductionCoefficient = WorldData.ConductivityWaterFlora,
					SurfaceArea = dependent.SurfaceAreaWaterFlora,
					Coverage = dependent.WaterCoverage[i],
					CoverageBelow = dependent.WaterCoverage[i - 1],
					SecondsPerTick = worldData.SecondsPerTick
				}, conductionWaterFloraJobHandle));
				conductionWaterTerrainJobHandle = JobHandle.CombineDependencies(conductionWaterTerrainJobHandle, SimJob.Run(new ConductionWaterBottomAJob()
				{
					EnergyDelta = conductionWaterTerrain[i],
					EnergyDeltaTotal = conductionWaterTerrainTotal,
					TemperatureA = lastState.WaterTemperature[i],
					TemperatureB = lastState.TerrainTemperature,
					EnergyA = dependent.WaterPotentialEnergy[i],
					ConductionCoefficient = WorldData.ConductivityWaterTerrain,
					SurfaceArea = dependent.SurfaceAreaWaterTerrain,
					Coverage = dependent.WaterCoverage[i],
					CoverageBelow = dependent.WaterCoverage[i - 1],
					SecondsPerTick = worldData.SecondsPerTick
				}, conductionWaterTerrainJobHandle));
			}

			#endregion


			#region Change temperature due to energy flux

			var energyJobHandles = new NativeArray<JobHandle>(_layerCount, Allocator.TempJob);
			var terrainEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_terrainLayer],
				thermalOutJobHandles[_terrainLayer],
				thermalInDownJobHandles[_terrainLayer],
				thermalInUpJobHandles[_terrainLayer],
				conductionAirTerrainJobHandle,
				conductionIceTerrainJobHandle,
				conductionWaterTerrainJobHandle,
				conductionFloraTerrainJobHandle,
			};
			jobHandleDependencies.Add(terrainEnergyJobHandleDependencies);
			energyJobHandles[_terrainLayer] =SimJob.Schedule(new EnergyTerrainJob()
			{
				TerrainTemperature = nextState.TerrainTemperature,
				LastTemperature = lastState.TerrainTemperature,
				SoilFertility = lastState.SoilFertility,
				SolarRadiationIn = solarRadiationIn[_terrainLayer],
				ThermalRadiationDelta = thermalRadiationDelta[_terrainLayer],
				ConductionEnergyAir = conductionAirTerrain,
				ConductionEnergyIce = conductionIceTerrain,
				ConductionEnergyFlora = conductionFloraTerrain,
				ConductionEnergyWater = conductionWaterTerrainTotal,
				GeothermalEnergy = geothermalRadiation,
				HeatingDepth = worldData.SoilHeatDepth,
			}, JobHandle.CombineDependencies(terrainEnergyJobHandleDependencies));

			var energyIceJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_iceLayer],
				thermalOutJobHandles[_iceLayer],
				thermalInDownJobHandles[_iceLayer],
				thermalInUpJobHandles[_iceLayer],
				conductionAirIceJobHandle,
				conductionIceWaterJobHandle,
				conductionIceFloraJobHandle,
				conductionIceTerrainJobHandle,
			};
			jobHandleDependencies.Add(energyIceJobHandleDependencies);
			energyJobHandles[_iceLayer] = SimJob.Schedule(new EnergyIceJob()
			{
				Temperature = nextState.IceTemperature,
				LastTemperature = lastState.IceTemperature,
				LastMass = lastState.IceMass,
				SolarRadiationIn = solarRadiationIn[_iceLayer],
				ThermalRadiationDelta = thermalRadiationDelta[_iceLayer],
				ConductionEnergyAir = conductionAirIce,
				ConductionEnergyTerrain = conductionIceTerrain,
				ConductionEnergyWater = conductionIceWater,
				ConductionEnergyFlora = conductionIceFlora,
			}, JobHandle.CombineDependencies(energyIceJobHandleDependencies));

			var energyFloraJobHandleDependencies = new NativeList<JobHandle>(Allocator.TempJob)
			{
				solarInJobHandles[_floraLayer],
				thermalOutJobHandles[_floraLayer],
				thermalInDownJobHandles[_floraLayer],
				thermalInUpJobHandles[_floraLayer],
				conductionAirFloraJobHandle,
				conductionWaterFloraJobHandle,
				conductionIceFloraJobHandle,
				conductionFloraTerrainJobHandle,
			};
			jobHandleDependencies.Add(energyFloraJobHandleDependencies);
			energyJobHandles[_floraLayer] = SimJob.Schedule(new EnergyFloraJob()
			{
				FloraTemperature = nextState.FloraTemperature,
				LastTemperature = lastState.FloraTemperature,
				FloraMass = lastState.FloraMass,
				FloraWater = lastState.FloraWater,
				SolarRadiationIn = solarRadiationIn[_floraLayer],
				ThermalRadiationDelta = thermalRadiationDelta[_floraLayer],
				ConductionEnergyAir = conductionAirFlora,
				ConductionEnergyTerrain = conductionFloraTerrain,
				ConductionEnergyWater = conductionWaterFloraTotal,
				ConductionEnergyIce = conductionIceFlora,
			}, JobHandle.CombineDependencies(energyFloraJobHandleDependencies));

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
					airDependencies.Add(conductionAirFloraJobHandle);
					airDependencies.Add(conductionAirTerrainJobHandle);
					energyJobHandles[layerIndex] = SimJob.Schedule(new EnergyAirSurfaceJob()
					{
						AirTemperaturePotential = nextState.AirTemperaturePotential[j],
						LastTemperaturePotential = lastState.AirTemperaturePotential[j],
						LastVapor = lastState.AirVapor[j],
						AirMass = dependent.AirMass[j],
						ConductionEnergyWater = conductionAirWater,
						ConductionEnergyIce = conductionAirIce,
						ConductionEnergyFlora = conductionAirFlora,
						ConductionEnergyTerrain = conductionAirTerrain,
						SolarRadiationIn = solarRadiationIn[layerIndex],
						ThermalRadiationDelta = thermalRadiationDelta[layerIndex],
					}, JobHandle.CombineDependencies(airDependencies));
				}
				else
				{
					energyJobHandles[layerIndex] = SimJob.Schedule(new EnergyAirJob()
					{
						AirTemperaturePotential = nextState.AirTemperaturePotential[j],
						LastTemperaturePotential = lastState.AirTemperaturePotential[j],
						LastVapor = lastState.AirVapor[j],
						AirMass = dependent.AirMass[j],
						SolarRadiationIn = solarRadiationIn[layerIndex],
						ThermalRadiationDelta = thermalRadiationDelta[layerIndex],
					}, JobHandle.CombineDependencies(airDependencies));

				}
				jobHandleDependencies.Add(airDependencies);

			}

			for (int j = 1; j < _waterLayers - 1; j++)
			{
				int layerIndex = _waterLayer0 + j;
				var waterDependencies = new NativeList<JobHandle>(Allocator.TempJob)
				{
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
					conductionAirWaterJobHandle,
					conductionIceWaterJobHandle,
					conductionWaterFloraJobHandle,
					conductionWaterTerrainJobHandle,
				};
				jobHandleDependencies.Add(waterDependencies);

				energyJobHandles[layerIndex] = SimJob.Schedule(new EnergyWaterJob()
				{
					Temperature = nextState.WaterTemperature[j],
					LastMass = lastState.WaterMass[j],
					LastSaltMass = lastState.SaltMass[j],
					LastTemperature = lastState.WaterTemperature[j],
					SolarRadiationIn = solarRadiationIn[layerIndex],
					ThermalRadiationDelta = thermalRadiationDelta[layerIndex],
					CoverageUp = dependent.WaterCoverage[j + 1],
					CoverageDown = dependent.WaterCoverage[j - 1],
					ConductionEnergyAir = conductionAirWater,
					ConductionEnergyIce = conductionIceWater,
					ConductionEnergyFlora = conductionWaterFlora[j],
					ConductionEnergyTerrain = conductionWaterTerrain[j],
				}, JobHandle.CombineDependencies(waterDependencies));
			}
			JobHandle energyJobHandle = JobHandle.CombineDependencies(energyJobHandles);
			energyJobHandle.Complete();

			#endregion



			#region State Change (Flux)

			// surface water
			var fluxJobHandles = new NativeArray<JobHandle>(_layerCount, Allocator.TempJob);

			int surfaceWaterLayerIndex = _waterLayer0 + _surfaceWaterLayer;

			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layerIndex = j + _airLayer0;
				fluxJobHandles[layerIndex] =SimJob.Schedule(new FluxAirJob()
				{
					LatentHeat = latentHeat[layerIndex],
					CondensationCloudMass = condensationCloudMass[j],
					CondensationGroundMass = condensationGroundMass[j],

					TemperaturePotential = nextState.AirTemperaturePotential[j],
					LastVapor = lastState.AirVapor[j],
					AirMass = dependent.AirMass[j],
					AirPressure = dependent.AirPressure[j],
					CloudElevation = dependent.CloudElevation,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					LayerMiddle = dependent.LayerMiddle[j],

				});
			}


			fluxJobHandles[_waterLayer0+_surfaceWaterLayer] =SimJob.Schedule(new FluxWaterJob()
			{
				EvaporatedWaterMass = evaporationMassWater,
				EvaporatedWaterTemperaturePotential = evaporationTemperaturePotentialWater,
				FrozenMass = frozenMass,
				FrozenTemperature = frozenTemperature,
				LatentHeatWater = latentHeat[_waterLayer0 + _surfaceWaterLayer],
				LatentHeatAir = latentHeat[_airLayer0 + 1],
				SaltPlume = saltPlume,

				Temperature = nextState.WaterTemperature[_surfaceWaterLayer],
				AirTemperaturePotential = nextState.AirTemperaturePotential[1],
				LastMass = lastState.WaterMass[_surfaceWaterLayer],
				LastSaltMass = lastState.SaltMass[_surfaceWaterLayer],
				IceCoverage = dependent.IceCoverage,
				WaterCoverage = dependent.WaterCoverage[_surfaceWaterLayer],
				SurfaceWind = lastState.AirVelocity[1],
				AirMass = dependent.AirMass[1],
				AirPressure = dependent.AirPressure[1],
				AirVapor = lastState.AirVapor[1],
				LayerElevation = dependent.LayerElevation[1],
				LayerMiddle = dependent.LayerMiddle[1],
				WaterHeatingDepth = worldData.WaterHeatingDepth,
				FreezePointReductionPerSalinity = worldData.FreezePointReductionPerSalinity,
			}, fluxJobHandles[_airLayer0 + 1]);


			// CLOUD
			var fluxCloudJobHandle =SimJob.Schedule(new FluxCloudJob()
			{
				EvaporationMass = cloudEvaporationMass,
				PrecipitationMass = precipitationMass,
				PrecipitationTemperature = precipitationTemperature,
				DropletDelta = dropletDelta,

				SurfaceAirTemperaturePotential = nextState.AirTemperaturePotential[1],
				SurfaceLayerElevation = dependent.LayerElevation[1],
				SurfaceLayerMiddle = dependent.LayerMiddle[1],
				SurfaceSaltMass = lastState.SaltMass[_surfaceWaterLayer],
				LastCloudMass = lastState.CloudMass,
				LastVelocity = lastState.CloudVelocity,
				LastDropletMass = lastState.CloudDropletMass,
				CloudElevation = dependent.CloudElevation,
				DewPoint = dependent.DewPoint,
				AirDensityCloud = dependent.AirDensityCloud,
				Position = staticState.SphericalPosition,
				Gravity = lastState.PlanetState.Gravity,
				RainDropDragCoefficient = worldData.rainDropDragCoefficient,
				RainDropMaxSize = worldData.rainDropMaxSize,
				RainDropMinSize = worldData.rainDropMinSize,
				SecondsPerTick = worldData.SecondsPerTick,
				CloudDissapationRateDryAir = worldData.CloudDissapationRateDryAir,
				CloudDissapationRateWind = worldData.CloudDissapationRateWind,
			});

			fluxJobHandles[_floraLayer] = SimJob.Schedule(new FluxFloraJob()
			{
				LatentHeatAir = latentHeat[_airLayer0 + 1],
				EvaporatedWaterMass = evaporationMassFlora,
				EvaporatedWaterTemperaturePotential = evaporationTemperaturePotentialFlora,
				GroundWaterConsumed = groundWaterConsumed,
				FloraMassDelta = floraMassDelta,

				FloraTemperature = nextState.FloraTemperature,
				FloraMass = lastState.FloraMass,
				FloraWater = lastState.FloraWater,
				LayerElevation = dependent.LayerElevation[1],
				SurfaceWind = lastState.AirVelocity[1],
				AirMass = dependent.AirMass[1],
				AirPressure = dependent.AirPressure[1],
				AirVapor = lastState.AirVapor[1],
				SoilFertility = lastState.SoilFertility,
				GroundWater = lastState.GroundWater,
				GroundWaterMax = worldData.GroundWaterMax,
				FloraEvaporationRate = worldData.FloraEvaporationRate,
				FloraWaterConsumptionRate = worldData.FloraWaterConsumptionRate,
				FloraGrowthRate = worldData.FloraGrowthRate,
				FloraMax = worldData.FloraMax,
				FloraGrowthTemperatureRangeInverse = worldData.FloraGrowthTemperatureRangeInverse
			}, fluxJobHandles[_waterLayer0 + _surfaceWaterLayer]);

			fluxJobHandles[_iceLayer] = SimJob.Schedule(new FluxIceJob()
			{
				LatentHeatAir = latentHeat[_airLayer0 + 1],
				LatentHeatWater = latentHeat[_waterLayer0 + _surfaceWaterLayer],
				MeltedMass = iceMeltedMass,

				Temperature = nextState.IceTemperature,
				LastMass = lastState.IceMass,
				IceHeatingDepth = worldData.IceHeatingDepth,
				AirTemperaturePotential = nextState.AirTemperaturePotential[1],
				WaterTemperature = nextState.WaterTemperature[_surfaceWaterLayer],
				LayerElevation = dependent.LayerElevation[1],

			}, fluxJobHandles[_floraLayer]);

			JobHandle fluxJobHandle = JobHandle.CombineDependencies(fluxCloudJobHandle, JobHandle.CombineDependencies(fluxJobHandles));
			fluxJobHandle.Complete();

			#endregion

			#region Update Mass - Evaporation, Condensation, Melting, Rainfall

			JobHandle updateMassJobHandle = default(JobHandle);

			var updateMassWaterJobHandles = new JobHandle[_waterLayers];
			for (int j=1;j<_waterLayers-1;j++)
			{
				updateMassJobHandle = JobHandle.CombineDependencies(updateMassJobHandle,SimJob.Schedule(new UpdateMassWaterJob()
				{
					WaterMass = nextState.WaterMass[j],
					SaltMass = nextState.SaltMass[j],
					WaterTemperature = nextState.WaterTemperature[j],
					SaltPlume = saltPlume,
					SaltPlumeTemperature = frozenTemperature,
					LastSaltMass = lastState.SaltMass[j],
					LastWaterMass = lastState.WaterMass[j],
					DownLastWaterMass = lastState.WaterMass[j-1]
				}));
				updateMassWaterJobHandles[j] = updateMassJobHandle;
			}

			var surfaceWaterMassHandle = updateMassWaterJobHandles[_surfaceWaterLayer];
			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layerIndex = _airLayer0 + j;
				surfaceWaterMassHandle = JobHandle.CombineDependencies(surfaceWaterMassHandle,SimJob.Schedule(new UpdateMassCondensationGroundJob()
				{
					SurfaceWaterMass = nextState.WaterMass[_surfaceWaterLayer],
					SurfaceWaterTemperature = nextState.WaterTemperature[_surfaceWaterLayer],

					AirTemperaturePotential = nextState.AirTemperaturePotential[j],
					GroundCondensation = condensationGroundMass[j],
					SurfaceSaltMass = lastState.SaltMass[j],
					LayerMiddle = dependent.LayerMiddle[j],
				}, surfaceWaterMassHandle));
			}
			updateMassJobHandle = JobHandle.CombineDependencies(updateMassJobHandle, surfaceWaterMassHandle);

			surfaceWaterMassHandle =SimJob.Schedule(new UpdateMassWaterSurfaceJob()
			{
				WaterTemperature = nextState.WaterTemperature[_surfaceWaterLayer],
				WaterMass = nextState.WaterMass[_surfaceWaterLayer],
				SaltMass = nextState.SaltMass[_surfaceWaterLayer],
				SaltPlume = saltPlume,
				Evaporation = evaporationMassWater,
				IceMelted = iceMeltedMass,
				Precipitation = precipitationMass,
				PrecipitationTemperature = precipitationTemperature,
				WaterFrozen = frozenMass,
			}, surfaceWaterMassHandle);

			var updateCloudMassJobHandle =SimJob.Schedule(new UpdateMassCloudJob()
			{
				CloudMass = nextState.CloudMass,
				CloudDropletMass = nextState.CloudDropletMass,
				LastCloudMass = lastState.CloudMass,
				LastDropletMass = lastState.CloudDropletMass,
				CloudEvaporation = cloudEvaporationMass,
				PrecipitationMass = precipitationMass
			});
			updateMassJobHandle = JobHandle.CombineDependencies(updateMassJobHandle, updateCloudMassJobHandle);

			var updateMassAirJobHandles = new JobHandle[_airLayers];
			var updateAirMassJobHandle = default(JobHandle);
			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layerIndex = _airLayer0 + j;
				updateAirMassJobHandle =SimJob.Schedule(new UpdateMassAirJob()
				{
					VaporMass = nextState.AirVapor[j],
					CloudMass = nextState.CloudMass,
					CloudDropletMass = nextState.CloudDropletMass,

					CloudEvaporation = cloudEvaporationMass,
					CloudElevation = dependent.CloudElevation,
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					CloudCondensation = condensationCloudMass[j],
					GroundCondensation = condensationGroundMass[j],
					LastVaporMass = lastState.AirVapor[j],
					IsTop = j == _airLayers - 2,
					IsBottom = j == 1,
				}, JobHandle.CombineDependencies(updateCloudMassJobHandle, updateAirMassJobHandle));
				updateMassJobHandle = JobHandle.CombineDependencies(updateMassJobHandle, updateAirMassJobHandle);
				updateMassAirJobHandles[j] = updateMassJobHandle;
			}
			var updateMassEvaporationHandle =SimJob.Schedule(new UpdateMassEvaporationJob()
			{
				AirTemperaturePotential = nextState.AirTemperaturePotential[1],
				VaporMass = nextState.AirVapor[1],

				AirMass = dependent.AirMass[1],
				EvaporationWater = evaporationMassWater,
				EvaporationTemperaturePotentialWater = evaporationTemperaturePotentialWater,
				EvaporationFlora = evaporationMassFlora,
				EvaporationTemperaturePotentialFlora = evaporationTemperaturePotentialFlora,
			}, JobHandle.CombineDependencies(updateMassAirJobHandles[1], surfaceWaterMassHandle));
			updateMassJobHandle = JobHandle.CombineDependencies(updateMassJobHandle, updateMassEvaporationHandle);


			var updateIceMassJobHandle =SimJob.Schedule(new UpdateMassIceJob()
			{
				IceMass = nextState.IceMass,
				IceTemperature = nextState.IceTemperature,

				LastIceMass = lastState.IceMass,
				IceMelted = iceMeltedMass,
				WaterFrozen = frozenMass,
				WaterTemperature = frozenTemperature,
				Precipitation = precipitationMass,
				PrecipitationTemperature = precipitationTemperature,
			});
			updateMassJobHandle = JobHandle.CombineDependencies(updateMassJobHandle, updateIceMassJobHandle);

			updateMassJobHandle = JobHandle.CombineDependencies(updateMassJobHandle, SimJob.Schedule(new UpdateTerrainJob()
			{
				SoilFertility = nextState.SoilFertility,
				Roughness = nextState.Roughness,
				GroundWater = nextState.GroundWater,
				Elevation = nextState.Elevation,

				LastElevation = lastState.Elevation,
				LastRoughness = lastState.Roughness,
				LastSoilFertility = lastState.SoilFertility,
				LastGroundWater = lastState.GroundWater,
				GroundWaterConsumed = groundWaterConsumed
			}));

			updateMassJobHandle = JobHandle.CombineDependencies(updateMassJobHandle, SimJob.Schedule(new UpdateFloraJob()
			{
				FloraMass = nextState.FloraMass,
				FloraWater = nextState.FloraWater,
				FloraMassDelta = floraMassDelta,

				EvaporationMass = evaporationMassFlora,
				LastMass = lastState.FloraMass,
				LastWater = lastState.FloraWater,
				GroundWaterConsumed = groundWaterConsumed
			}));


			updateMassJobHandle.Complete();


			#endregion

			#region Apply Latent Heat

			var latentHeatJobHandle =SimJob.Schedule(new ApplyLatentHeatIceJob()
			{
				IceTemperature = nextState.IceTemperature,
				IceMass = nextState.IceMass,
				LatentHeat = latentHeat[_iceLayer]
			});

			//for (int i = 1; i < _airLayers - 1; i++)
			//{
			//	latentHeatJobHandle = JobHandle.CombineDependencies(latentHeatJobHandle, ApplyLatentHeatJob.Run(new ApplyLatentHeatAirJob()
			//	{
			//		AirTemperaturePotential = nextState.AirTemperaturePotential[i],
			//		AirMass = dependent.AirMass[i],
			//		VaporMass = nextState.AirVapor[i],
			//		LatentHeat = latentHeat[_airLayer0 + i]
			//	}));
			//}
			//for (int i = 1; i < _waterLayers - 1; i++)
			//{
			//	latentHeatJobHandle = JobHandle.CombineDependencies(latentHeatJobHandle, ApplyLatentHeatJob.Run(new ApplyLatentHeatWaterJob()
			//	{
			//		WaterTemperature = nextState.WaterTemperature[i],
			//		WaterMass = nextState.WaterMass[i],
			//		SaltMass = nextState.SaltMass[i],
			//		LatentHeat = latentHeat[_waterLayer0 + i]
			//	}));
			//}

			latentHeatJobHandle.Complete();
			#endregion

			#region Ground Water

			var groundWaterJob = default(JobHandle);
			groundWaterJob = SimJob.Schedule(new GroundWaterFlowJob()
			{
				GroundWater = nextState.GroundWater,
				GroundWaterTemperature = nextState.GroundWaterTemperature,

				LastGroundWater = lastState.GroundWater,
				LastGroundWaterTemperature = lastState.GroundWaterTemperature,
				Elevation = nextState.Elevation,
				Neighbors = staticState.Neighbors,
				NeighborDistInverse = staticState.NeighborDistInverse,
				FlowSpeed = worldData.GroundWaterFlowSpeed,
				GroundWaterMaxInverse = 1.0f / worldData.GroundWaterMax,
			}, groundWaterJob);

			groundWaterJob = SimJob.Schedule(new GroundWaterDiffusionJob()
			{
				GroundWater = groundWaterFlowMass,
				GroundWaterTemperature = groundWaterFlowTemperature,

				LastGroundWater = nextState.GroundWater,
				LastGroundWaterTemperature = nextState.GroundWaterTemperature,
				NeighborDist = staticState.NeighborDist,
				NeighborDistInverse = staticState.NeighborDistInverse,
				Neighbors = staticState.Neighbors,
				DiffusionCoefficient = worldData.GroundWaterDiffusionCoefficient
			}, groundWaterJob);


			for (int i = 1; i < _waterLayers - 1; i++) {
				groundWaterJob = SimJob.Schedule(new GroundWaterAbsorptionJob()
				{
					GroundWater = nextState.GroundWater,
					GroundWaterTemperature = nextState.GroundWaterTemperature,

					LastGroundWater = groundWaterFlowMass,
					LastGroundWaterTemperature = groundWaterFlowTemperature,
					WaterMass = nextState.WaterMass[i],
					WaterTemperature = nextState.WaterTemperature[i],
					WaterBelow = nextState.WaterMass[i - 1],
					GroundWaterAbsorptionRate = worldData.GroundWaterAbsorptionRate * worldData.SecondsPerTick,
					GroundWaterMaxInverse = 1.0f / worldData.GroundWaterMax,
				}, groundWaterJob);
			}

			groundWaterJob = SimJob.Schedule(new GroundWaterConductionJob()
			{
				GroundWaterTemperature = groundWaterFlowTemperature,
				TerrainTemperature = nextState.TerrainTemperature,

				GroundWater = nextState.GroundWater,
				LastGroundWaterTemperature = nextState.GroundWaterTemperature,
				SoilFertility = nextState.SoilFertility,
				GroundWaterConductionCoefficient = WorldData.ConductivityWaterTerrain,
				HeatingDepth = worldData.SoilHeatDepth,
				SecondsPerTick = worldData.SecondsPerTick
			}, groundWaterJob);


			groundWaterJob.Complete();

			#endregion

			// Buoyancy, Updrafts, and mixing occur across air layers and water layers
			// TODO: add an empty air layer on top and bottom so we can calculate up/down diffusion in a single step 
			// Temperature and trace elements diffuse into neighboring horizontal cells based on a diffusion constant
			// Air, Water, Cloud
			#region Update Velocity


			var airTerrainFrictionJobHandle =SimJob.Schedule(new AirTerrainFrictionJob()
			{
				Force = windFriction,
				IceCoverage = dependent.IceCoverage,
				WaterCoverage = dependent.WaterCoverage[_surfaceWaterLayer],
				FloraCoverage = dependent.FloraCoverage,
				Roughness = lastState.Roughness,
				IceFriction = worldData.WindIceFriction,
				TerrainFrictionMin = worldData.WindTerrainFrictionMin,
				TerrainFrictionMax = worldData.WindTerrainFrictionMax,
				FloraFriction = worldData.WindFloraFriction,
				WaterFriction = worldData.WindWaterFriction,
				MaxTerrainRoughness = worldData.MaxTerrainRoughnessForWindFriction
			});

			var waterFrictionJobHandle =SimJob.Schedule(new WaterSurfaceFrictionJob()
			{
				Force = waterFriction,

				Position = staticState.SphericalPosition,
				Current = lastState.WaterVelocity[_surfaceWaterLayer],
				AirVelocityUp = lastState.AirVelocity[1],
				AirVelocityDown = lastState.WaterVelocity[_surfaceWaterLayer - 1],
				LayerHeight = dependent.WaterLayerHeight[_surfaceWaterLayer],
				CoriolisMultiplier = staticState.CoriolisMultiplier,
				FrictionCoefficientUp = worldData.WindToWaterCurrentFrictionCoefficient,
				FrictionCoefficientDown = 0, // TODO: do we want to add a frictional force between layers of water?
				CoriolisTerm = coriolisTerm,
				WaterSurfaceFrictionDepth = worldData.WaterSurfaceFrictionDepth,
				SecondsPerTick = worldData.SecondsPerTick
			});


			var velocityJobHandle = default(JobHandle);
			JobHandle[] airAccelerationJobHandles = new JobHandle[_airLayers];
			for (int j = 1; j < _airLayers - 1; j++)
			{
				airAccelerationJobHandles[j] =SimJob.Schedule(new AccelerationAirJob()
				{
					Velocity = nextState.AirVelocity[j],
					Force = airAcceleration[j],

					Friction = windFriction,
					Pressure = dependent.AirPressure[j],
					AirMass = dependent.AirMass[j],
					TemperaturePotential = lastState.AirTemperaturePotential[j],
					VaporMass = lastState.AirVapor[j],
					LayerMiddle = dependent.LayerMiddle[j],
					Neighbors = staticState.Neighbors,
					NeighborDiffInverse = staticState.NeighborDiffInverse,
					Positions = staticState.SphericalPosition,
					PlanetRadius = staticState.PlanetRadius,
					Gravity = lastState.PlanetState.Gravity,
					GravityInverse = 1.0f / lastState.PlanetState.Gravity,
					UpTemperaturePotential = lastState.AirTemperaturePotential[j + 1],
					UpHumidity = lastState.AirVapor[j + 1],
					UpAirMass = dependent.AirMass[j + 1],
					UpLayerMiddle = dependent.LayerMiddle[j + 1],
					DownTemperaturePotential = lastState.AirTemperaturePotential[j - 1],
					DownHumidity = lastState.AirVapor[j - 1],
					DownAirMass = dependent.AirMass[j - 1],
					DownLayerMiddle = dependent.LayerMiddle[j - 1],
					IsTop = j == _airLayers - 2,
					IsBottom = j == 1,
					FrictionCoefficient = j == 1 ? 1 : 0,
					SecondsPerTick = worldData.SecondsPerTick

				}, JobHandle.CombineDependencies(waterFrictionJobHandle, airTerrainFrictionJobHandle));
				velocityJobHandle = JobHandle.CombineDependencies(velocityJobHandle, airAccelerationJobHandles[j]);
			}

			JobHandle updateCloudVelocityJobHandle = default(JobHandle);
			for (int j = 1; j < _airLayers - 1; j++)
			{

				updateCloudVelocityJobHandle = JobHandle.CombineDependencies(updateCloudVelocityJobHandle,SimJob.Schedule(new UpdateVelocityCloudJob()
				{
					Velocity = nextState.CloudVelocity,

					CloudElevation = dependent.CloudElevation,
					LayerVelocity = airAcceleration[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					DownVelocity = airAcceleration[j - 1],
					DownLayerElevation = dependent.LayerElevation[j - 1],
					UpVelocity = airAcceleration[j + 1],
					UpLayerElevation = dependent.LayerElevation[j + 1],
					UpLayerHeight = dependent.LayerHeight[j + 1],
					IsTop = j == _airLayers - 2,
					IsBottom = j == 1,
				}, JobHandle.CombineDependencies(updateCloudVelocityJobHandle, airTerrainFrictionJobHandle, JobHandle.CombineDependencies(airAccelerationJobHandles[j], airAccelerationJobHandles[j - 1], airAccelerationJobHandles[j + 1]))));
			}
			velocityJobHandle = JobHandle.CombineDependencies(velocityJobHandle, updateCloudVelocityJobHandle);


			JobHandle[] waterAccelerationJobHandles = new JobHandle[_airLayers];
			for (int j = 1; j < _waterLayers - 1; j++)
			{
				velocityJobHandle = JobHandle.CombineDependencies(velocityJobHandle,SimJob.Schedule(new AccelerationWaterJob()
				{
					Velocity = nextState.WaterVelocity[j],

					Positions = staticState.SphericalPosition,
					Neighbors = staticState.Neighbors,
					NeighborDiffInverse = staticState.NeighborDiffInverse,
					WaterDensity = dependent.WaterDensity[j],
					WaterPressure = dependent.WaterPressure[j],
					LayerDepth = dependent.WaterLayerDepth[j],
					LayerHeight = dependent.WaterLayerHeight[j],
					UpWaterDensity = dependent.WaterDensity[j + 1],
					UpWaterPressure = dependent.WaterPressure[j + 1],
					UpLayerDepth = dependent.WaterLayerDepth[j + 1],
					UpLayerHeight = dependent.WaterLayerHeight[j + 1],
					DownWaterDensity = dependent.WaterDensity[j - 1],
					DownWaterPressure = dependent.WaterPressure[j - 1],
					DownLayerDepth = dependent.WaterLayerDepth[j - 1],
					DownLayerHeight = dependent.WaterLayerHeight[j - 1],
					SurfaceElevation = dependent.LayerElevation[1],
					Friction = waterFriction,
					FrictionCoefficient = j == _surfaceWaterLayer ? 1 : 0,
					Gravity = lastState.PlanetState.Gravity,
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick

				}, waterFrictionJobHandle));
				waterAccelerationJobHandles[j] = velocityJobHandle;
			}

			velocityJobHandle.Complete();


			#endregion




			// Diffuse from last time step
			// Air, Water, Cloud
			#region Diffusion

			JobHandle[] diffusionJobHandles = new JobHandle[_layerCount];
			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layer = _airLayer0 + j;
				// TODO: is it a problem that we are using the dependent variables from last frame while referencing our newly calculated next frame values for temperature and such?
				diffusionJobHandles[layer] =SimJob.Schedule(new DiffusionAirJob()
				{
					Delta = diffusionAir[j],

					LastTemperature = nextState.AirTemperaturePotential[j],
					LastVapor = nextState.AirVapor[j],
					LastVelocity = nextState.AirVelocity[j],
					Neighbors = staticState.Neighbors,
					LayerHeight = dependent.LayerHeight[j],
					NeighborDistInverse = staticState.NeighborDistInverse,
					AirMass = dependent.AirMass[j],
					UpTemperature = nextState.AirTemperaturePotential[j + 1],
					UpVapor = nextState.AirVapor[j + 1],
					UpAirVelocity = nextState.AirVelocity[j + 1],
					UpAirMass = dependent.AirMass[j + 1],
					UpLayerElevation = dependent.LayerElevation[j + 1],
					UpLayerHeight = dependent.LayerHeight[j + 1],
					DownTemperature = nextState.AirTemperaturePotential[j - 1],
					DownVapor = nextState.AirVapor[j - 1],
					DownAirVelocity = nextState.AirVelocity[j - 1],
					DownAirMass = dependent.AirMass[j - 1],
					DownLayerElevation = dependent.LayerElevation[j - 1],
					DownLayerHeight = dependent.LayerHeight[j - 1],
					IsTop = j == _airLayers - 2,
					IsBottom = j == 1,
					DiffusionCoefficientHoriztonal = worldData.AirDiffusionCoefficientHorizontal,
					DiffusionCoefficientVertical = worldData.AirDiffusionCoefficientVertical,
				});
			}
			for (int j = 1; j < _waterLayers - 1; j++)
			{
				int layer = _waterLayer0 + j;
				diffusionJobHandles[layer] =SimJob.Schedule(new DiffusionWaterJob()
				{
					Delta = diffusionWater[j],

					LastTemperature = nextState.WaterTemperature[j],
					LastSalt = nextState.SaltMass[j],
					LastCurrent = nextState.WaterVelocity[j],
					LastMass = nextState.WaterMass[j],
					LayerHeight = dependent.LayerHeight[j],
					UpLayerHeight = dependent.LayerHeight[j + 1],
					DownLayerHeight = dependent.LayerHeight[j - 1],
					NeighborDistInverse = staticState.NeighborDistInverse,
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
				});
			}

			var diffusionCloudHandle =SimJob.Schedule(new DiffusionCloudJob()
			{
				Delta = diffusionCloud,

				LastMass = nextState.CloudMass,
				LastTemperature = nextState.CloudTemperature,
				LastDropletMass = nextState.CloudDropletMass,
				LastVelocity = nextState.CloudVelocity,
				Neighbors = staticState.Neighbors,
				DiffusionCoefficient = worldData.CloudDiffusionCoefficient,
			});

			#endregion

			#region Apply Diffusion

			JobHandle diffusionJobHandle = default(JobHandle);
			diffusionJobHandle = JobHandle.CombineDependencies(diffusionJobHandle,SimJob.Schedule(new ApplyAdvectionCloudJob()
			{
				Advection = diffusionCloud,
				CloudMass = nextState.CloudMass,
				Temperature = nextState.CloudTemperature,
				DropletMass = nextState.CloudDropletMass,
				Velocity = nextState.CloudVelocity
			}, diffusionCloudHandle));

			for (int i = 1; i < _waterLayers - 1; i++)
			{
				diffusionJobHandle = JobHandle.CombineDependencies(diffusionJobHandle,SimJob.Schedule(new ApplyAdvectionWaterJob()
				{
					Advection = diffusionWater[i],
					SaltMass = nextState.SaltMass[i],
					Temperature = nextState.WaterTemperature[i],
					Velocity = nextState.WaterVelocity[i],
					Mass = nextState.WaterMass[i]
				}, JobHandle.CombineDependencies(diffusionJobHandles[i + _waterLayer0], diffusionJobHandles[i + _waterLayer0 - 1], diffusionJobHandles[i + _waterLayer0 + 1])));
			}

			for (int i = 1; i < _airLayers - 1; i++)
			{
				diffusionJobHandle = JobHandle.CombineDependencies(diffusionJobHandle,SimJob.Schedule(new ApplyAdvectionAirJob()
				{
					Advection = diffusionAir[i],
					Vapor = nextState.AirVapor[i],
					Temperature = nextState.AirTemperaturePotential[i],
					AirVelocity = nextState.AirVelocity[i],
				}, JobHandle.CombineDependencies(diffusionJobHandles[i + _airLayer0], diffusionJobHandles[i + _airLayer0 - 1], diffusionJobHandles[i + _airLayer0 + 1])));
			}

			diffusionJobHandle.Complete();

			#endregion



			// Wind and currents move temperature and trace elements horizontally
			// Air, Water, Cloud
			// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
			#region Advection


			JobHandle[] advectionJobHandles = new JobHandle[_layerCount];
			for (int j = 1; j < _airLayers - 1; j++)
			{
				int layer = _airLayer0 + j;
				var airDestJob =SimJob.Schedule(new GetVectorDestCoordsJob()
				{
					Destination = destinationAir[j],
					DeflectedVelocity = dependent.DeflectedAirVelocity[j],
					Neighbors = staticState.Neighbors,
					Position = staticState.SphericalPosition,
					Velocity = nextState.AirVelocity[j],
					CoriolisMultiplier =staticState.CoriolisMultiplier,
					CoriolisTerm = coriolisTerm,
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick
				}, JobHandle.CombineDependencies(diffusionJobHandle, airAccelerationJobHandles[j]));
				advectionJobHandles[layer] =SimJob.Schedule(new AdvectionAirJob()
				{
					Delta = advectionAir[j],
					Temperature = nextState.AirTemperaturePotential[j],
					Vapor = nextState.AirVapor[j],
					Velocity = nextState.AirVelocity[j],
					DeflectedVelocity = dependent.DeflectedAirVelocity[j],
					Neighbors = staticState.Neighbors,
					Destination = destinationAir[j],
					LayerElevation = dependent.LayerElevation[j],
					LayerHeight = dependent.LayerHeight[j],
					UpTemperature = nextState.AirTemperaturePotential[j + 1],
					UpHumidity = nextState.AirVapor[j + 1],
					UpAirMass = dependent.AirMass[j + 1],
					UpLayerElevation = dependent.LayerElevation[j + 1],
					UpLayerHeight = dependent.LayerHeight[j + 1],
					DownTemperature = nextState.AirTemperaturePotential[j - 1],
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
				var waterDestJob =SimJob.Schedule(new GetVectorDestCoordsJob()
				{
					Destination = destinationWater[j],
					DeflectedVelocity = dependent.DeflectedWaterVelocity[j],
					Neighbors = staticState.Neighbors,
					Position = staticState.SphericalPosition,
					Velocity = nextState.WaterVelocity[j],
					CoriolisMultiplier = staticState.CoriolisMultiplier,
					CoriolisTerm = coriolisTerm,
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick
				}, JobHandle.CombineDependencies(diffusionJobHandle, waterAccelerationJobHandles[j]));
				advectionJobHandles[layer] =SimJob.Schedule(new AdvectionWaterJob()
				{
					Delta = advectionWater[j],
					Destination = destinationWater[j],
					Velocity = nextState.WaterVelocity[j],
					DeflectedVelocity = dependent.DeflectedWaterVelocity[j],
					Temperature = nextState.WaterTemperature[j],
					Mass = nextState.WaterMass[j],
					Salt = nextState.SaltMass[j],
					Neighbors = staticState.Neighbors,
				}, waterDestJob);
			}

			var cloudDestJob =SimJob.Schedule(new GetVectorDestCoordsJob()
			{
				Destination = destinationCloud,
				DeflectedVelocity = dependent.DeflectedCloudVelocity,
				Neighbors = staticState.Neighbors,
				Position = staticState.SphericalPosition,
				Velocity = nextState.CloudVelocity,
				CoriolisMultiplier = staticState.CoriolisMultiplier,
				CoriolisTerm = coriolisTerm,
				PlanetRadius = staticState.PlanetRadius,
				SecondsPerTick = worldData.SecondsPerTick
			}, JobHandle.CombineDependencies(diffusionJobHandle, updateCloudVelocityJobHandle));
			var advectionJobHandleCloud =SimJob.Schedule(new AdvectionCloudJob()
			{
				Delta = advectionCloud,
				Destination = destinationCloud,
				Mass = nextState.CloudMass,
				Temperature = nextState.CloudTemperature,
				DropletMass = nextState.CloudDropletMass,
				Velocity = nextState.CloudVelocity,
				DeflectedVelocity = dependent.DeflectedCloudVelocity,
				Neighbors = staticState.Neighbors,
			}, cloudDestJob);


			#endregion

			#region Apply Advection

			NativeList<JobHandle> applyAdvectionJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
			applyAdvectionJobHandles.Add(SimJob.Schedule(new ApplyAdvectionCloudJob()
			{
				Advection = advectionCloud,
				CloudMass = nextState.CloudMass,
				Temperature = nextState.CloudTemperature,
				DropletMass = nextState.CloudDropletMass,
				Velocity = nextState.CloudVelocity
			}, advectionJobHandleCloud));

			for (int i = 1; i < _waterLayers - 1; i++)
			{
				applyAdvectionJobHandles.Add(SimJob.Schedule(new ApplyAdvectionWaterJob()
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
				applyAdvectionJobHandles.Add(SimJob.Schedule(new ApplyAdvectionAirJob()
				{
					Advection = advectionAir[i],
					Vapor = nextState.AirVapor[i],
					Temperature = nextState.AirTemperaturePotential[i],
					AirVelocity = nextState.AirVelocity[i],
				}, JobHandle.CombineDependencies( advectionJobHandles[i+_airLayer0], advectionJobHandles[i + _airLayer0 - 1], advectionJobHandles[i + _airLayer0 + 1])));
			}
			JobHandle.CompleteAll(applyAdvectionJobHandles);

			#endregion

			#region Update dependent variables

			SimJobs.UpdateDependentVariables(SimJob, ref nextState, ref dependent, ref worldData, default(JobHandle), tempArrays).Complete();

			#endregion

			#region Debug
			if (checkForDegeneracy)
			{
				bool degen = false;
				SortedSet<int> degenIndices = new SortedSet<int>();
				List<string> degenVarNames = new List<string>();
				degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "TerrainTemperature", nextState.TerrainTemperature, 0, 400, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "FloraTemperature", nextState.FloraTemperature, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "FloraMass", nextState.FloraMass, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "FloraWater", nextState.FloraWater, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "GroundWater", nextState.GroundWater, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "CloudMass", nextState.CloudMass, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "CloudDropletMass", nextState.CloudDropletMass, degenVarNames);
				degen |= CheckDegen(_cellCount, degenIndices, "CloudVelocity", nextState.CloudVelocity, degenVarNames);
				degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "CloudElevation", dependent.CloudElevation, -100000, 100000, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "IceMass", nextState.IceMass, degenVarNames);
				degen |= CheckDegenPosValues(_cellCount, degenIndices, "IceTemperature", nextState.IceTemperature, degenVarNames);
				for (int i = 1; i < _airLayers - 1; i++) {
					degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "AirTemperature" + i, nextState.AirTemperaturePotential[i], 0, 400, degenVarNames);
					degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "AirVapor" + i, nextState.AirVapor[i], 0, 10000, degenVarNames);
					degen |= CheckDegen(_cellCount, degenIndices, "AirVelocity" + i, nextState.AirVelocity[i], degenVarNames);
				}
				for (int i=1;i<_waterLayers - 1;i++)
				{
					degen |= CheckDegenPosValues(_cellCount, degenIndices, "WaterMass" + i, nextState.WaterMass[i], degenVarNames);
					degen |= CheckDegenPosValues(_cellCount, degenIndices, "SaltMass" + i, nextState.SaltMass[i], degenVarNames);
					degen |= CheckDegenMinMaxValues(_cellCount, degenIndices, "WaterTemperature" + i, nextState.WaterTemperature[i], 0, 400, degenVarNames);
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

				var lastDisplay = display;
				display = new DisplayState();
				display.Init(_cellCount, _airLayers, _waterLayers, _layerCount);
				JobHandle initDisplayAirHandle = default(JobHandle);
				JobHandle initDisplayWaterHandle = default(JobHandle);
				for (int i = 1; i < _airLayers - 1; i++)
				{
					absorptivitySolar[i].CopyTo(display.AbsorptionSolar[i]);
					absorptivityThermal[i].CopyTo(display.AbsorptionThermal[i]);

					initDisplayAirHandle = JobHandle.CombineDependencies(initDisplayAirHandle, (new InitDisplayAirLayerJob()
					{
						DisplayPressure = display.Pressure[i],
						DisplayPressureGradientForce = display.PressureGradientForce[i],
						DisplayCondensationGround = display.CondensationGround,
						DisplayCondensationCloud = display.CondensationCloud,
						Enthalpy = display.EnthalpyAir[i],

						Gravity = curState.PlanetState.Gravity,
						AirTemperaturePotential = curState.AirTemperaturePotential[i],
						AirPressure = dependent.AirPressure[i],
						LayerMiddle = dependent.LayerMiddle[i],
						PressureGradientForce = airAcceleration[i],
						CondensationCloud = condensationCloudMass[i],
						CondensationGround = condensationGroundMass[i],
						AirMass = dependent.AirMass[i],
						VaporMass = nextState.AirVapor[i]
					}).Schedule(_cellCount, _batchCount, initDisplayAirHandle));
				}

				for (int i = 1; i < _waterLayers - 1; i++)
				{
					initDisplayWaterHandle = JobHandle.CombineDependencies(initDisplayWaterHandle, (new InitDisplayWaterLayerJob()
					{
						Enthalpy = display.EnthalpyWater[i],
						Salinity = display.Salinity[i],

						WaterTemperature = curState.WaterTemperature[i],
						SaltMass = curState.SaltMass[i],
						WaterMass = curState.WaterMass[i],
					}).Schedule(_cellCount, _batchCount, initDisplayWaterHandle));
				}
				for (int i = 0; i < _layerCount; i++)
				{
					solarRadiationIn[i].CopyTo(display.SolarDelta[i]);
					thermalRadiationDelta[i].CopyTo(display.ThermalDelta[i]);
				}

				var updateDisplayJob = new UpdateDisplayJob()
				{
					SolarRadiationAbsorbedSurface = display.SolarRadiationAbsorbedSurface,
					DisplayEvaporation = display.Evaporation,
					DisplayPrecipitation = display.Rainfall,
					Enthalpy = display.EnthalpyTerrain,
					
					SolarRadiationInTerrain = solarRadiationIn[_terrainLayer],
					SolarRadiationInIce = solarRadiationIn[_iceLayer],
					SolarRadiationInWaterSurface = solarRadiationIn[_waterLayer0 + _waterLayers - 2],
					Evaporation = evaporationMassWater,
					Precipitation = precipitationMass,
					SoilFertility = nextState.SoilFertility,
					TerrainTemperature = nextState.TerrainTemperature,
					Flora = nextState.FloraMass,
					FloraWater = nextState.FloraWater,
					HeatingDepth = worldData.SoilHeatDepth,
					CloudMass = nextState.CloudMass,
					IceMass = nextState.IceMass,
					IceTemperature = nextState.IceTemperature,
				};
				var updateDisplayJobHandle = updateDisplayJob.Schedule(_cellCount, _batchCount);
				var displayHandles = JobHandle.CombineDependencies(initDisplayAirHandle, initDisplayWaterHandle, updateDisplayJobHandle);
				displayHandles.Complete();

				if (displayGlobals)
				{
					float globalWaterMass = 0;
					float globalWaterSurfaceMass = 0;
					double globalAirMass = 0;
					double globalAirTemperature = 0;
					for (int i = 0; i < _cellCount; i++)
					{
						float waterMassSurface = curState.WaterMass[_surfaceWaterLayer][i];
						globalWaterSurfaceMass += waterMassSurface;
						display.GlobalOceanSurfaceTemperature += curState.WaterTemperature[_surfaceWaterLayer][i] * waterMassSurface;
						display.SolarRadiation += displaySolarRadiation[i];
						display.GeothermalRadiation += geothermalRadiation[i];
						display.GlobalCloudMass += curState.CloudMass[i];
						display.GlobalIceMass += curState.IceMass[i];
						display.GlobalOceanCoverage += dependent.WaterCoverage[_waterLayers - 2][i];
						display.GlobalSurfaceTemperature += dependent.SurfaceAirTemperatureAbsolute[i];
						display.GlobalOceanVolume += dependent.WaterLayerDepth[1][i];
						display.GlobalSeaLevel += dependent.LayerElevation[1][i];
						display.GlobalEvaporation += display.Evaporation[i];
						display.GlobalRainfall += display.Rainfall[i];
						display.GlobalCondensationCloud += display.CondensationCloud[i];
						display.GlobalCondensationGround += display.CondensationGround[i];
						display.GlobalEnthalpyTerrain += display.EnthalpyTerrain[i];
						display.GlobalTerrainTemperature += curState.TerrainTemperature[i];
						for (int j = 1; j < _airLayers - 1; j++)
						{
							globalAirTemperature += curState.AirTemperaturePotential[j][i] * (dependent.AirMass[j][i] + curState.AirVapor[j][i]);
							globalAirMass += dependent.AirMass[j][i];
							display.GlobalWaterVapor += curState.AirVapor[j][i];
							display.EnergySolarReflectedAtmosphere += solarReflected[j + _airLayer0][i];
							display.EnergySolarAbsorbedAtmosphere += solarRadiationIn[j + _airLayer0][i];
							display.GlobalEnthalpyAir += display.EnthalpyAir[j][i];
							display.GlobalCloudCoverage += absorptivitySolar[j][i].AbsorptivityCloud;
						}
						display.EnergySolarAbsorbedSurface += solarRadiationIn[_terrainLayer][i] + solarRadiationIn[_iceLayer][i];
						display.EnergySolarReflectedSurface += solarReflected[_terrainLayer][i] + solarReflected[_iceLayer][i];
						for (int j = 1; j < _waterLayers - 1; j++)
						{
							float waterMass = curState.WaterMass[j][i];
							globalWaterMass += waterMass;
							display.GlobalOceanTemperature += curState.WaterTemperature[j][i] * waterMass;

							float absorbed = solarRadiationIn[j + _waterLayer0][i];
							display.EnergySolarAbsorbedOcean += absorbed;
							display.EnergySolarAbsorbedSurface += absorbed;
							display.EnergySolarReflectedSurface += solarReflected[j + _waterLayer0][i];
							display.GlobalEnthalpyWater += display.EnthalpyWater[j][i];
						}
						display.EnergySurfaceConduction += conductionAirIce[i] + conductionAirTerrain[i] + conductionAirWater[i];
						display.EnergyOceanConduction += conductionAirWater[i];
						display.EnergyEvapotranspiration += (evaporationMassWater[i] + evaporationMassFlora[i]) * WorldData.LatentHeatWaterVapor;
						display.EnergyThermalBackRadiation += windowRadiationTransmittedDown[_airLayer0 + 1][i] + thermalRadiationTransmittedDown[_airLayer0 + 1][i];
						display.EnergyThermalOceanRadiation += (windowRadiationTransmittedUp[_waterLayer0 + _waterLayers - 2][i] + thermalRadiationTransmittedUp[_waterLayer0 + _surfaceWaterLayer][i]) * dependent.WaterCoverage[_waterLayers - 2][i];

						float surfaceRadiation = windowRadiationTransmittedUp[_iceLayer][i] + thermalRadiationTransmittedUp[_iceLayer][i];
						float surfaceRadiationOutWindow = windowRadiationTransmittedUp[_iceLayer][i];
						float radiationToSpace = thermalRadiationTransmittedUp[_airLayer0 + _airLayers - 2][i] + windowRadiationTransmittedUp[_airLayer0 + _airLayers - 2][i] - surfaceRadiationOutWindow;
						display.EnergyThermalAbsorbedAtmosphere += surfaceRadiation - surfaceRadiationOutWindow;
						display.EnergyThermalOutAtmosphere += radiationToSpace;
						display.EnergyThermalSurfaceOutAtmosphericWindow += surfaceRadiationOutWindow;
						display.EnergyThermalSurfaceRadiation += surfaceRadiation;

					}
					display.GlobalAirTemperaturePotential = globalAirTemperature / (globalAirMass + display.GlobalWaterVapor);
					display.GlobalOceanSurfaceTemperature /= globalWaterSurfaceMass;
					display.GlobalOceanTemperature /= globalWaterMass;
					display.GlobalTerrainTemperature /= _cellCount;
					display.GlobalEnthalpy = display.GlobalEnthalpyTerrain + display.GlobalEnthalpyAir + display.GlobalEnthalpyWater;
					display.GlobalEnthalpyDelta = display.GlobalEnthalpy - lastDisplay.GlobalEnthalpy;
					display.GlobalEnthalpyDeltaTerrain = display.GlobalEnthalpyTerrain - lastDisplay.GlobalEnthalpyTerrain;
					display.GlobalEnthalpyDeltaAir = display.GlobalEnthalpyAir - lastDisplay.GlobalEnthalpyAir;
					display.GlobalEnthalpyDeltaWater = display.GlobalEnthalpyWater - lastDisplay.GlobalEnthalpyWater;
				}
				lastDisplay.Dispose();

			}

			#endregion

			#region Dispose Temporary Arrays
			atmosphericWindowUp.Dispose();
			atmosphericWindowDown.Dispose();
			dropletDelta.Dispose();
			precipitationMass.Dispose();
			precipitationTemperature.Dispose();
			cloudEvaporationMass.Dispose();
			iceMeltedMass.Dispose();
			applyAdvectionJobHandles.Dispose();
			fluxJobHandles.Dispose();
			energyJobHandles.Dispose();
			foreach (var d in jobHandleDependencies)
			{
				d.Dispose();
			}
			foreach (var a in tempArrays)
			{
				a.Dispose();
			}
			conductionWaterTerrainTotal.Dispose();
			conductionWaterFloraTotal.Dispose();
			for (int i = 0; i < _layerCount; i++)
			{
				thermalRadiationDelta[i].Dispose();
				thermalRadiationTransmittedUp[i].Dispose();
				thermalRadiationTransmittedDown[i].Dispose();
				windowRadiationTransmittedUp[i].Dispose();
				windowRadiationTransmittedDown[i].Dispose();
				solarReflected[i].Dispose();
				latentHeat[i].Dispose();
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
		s.AppendFormat("Elevation: {0}\n", state.Elevation[i]);
		s.AppendFormat("Roughness: {0}\n", state.Roughness[i]);
		s.AppendFormat("SoilFertility: {0}\n", state.SoilFertility[i]);
		s.AppendFormat("Flora: {0}\n", state.FloraMass[i]);
		s.AppendFormat("Ground Water: {0} kg\n", state.GroundWater[i]);
		s.AppendFormat("TerrainTemperature: {0}\n", state.TerrainTemperature[i]);
		s.AppendFormat("CloudMass: {0}\n", state.CloudMass[i]);
		s.AppendFormat("CloudDropletMass: {0}\n", state.CloudDropletMass[i]);
		s.AppendFormat("CloudVelocity: {0}\n", state.CloudVelocity[i]);
		s.AppendFormat("IceMass: {0}\n", state.IceMass[i]);
		s.AppendFormat("IceTemperature: {0}\n", state.IceTemperature[i]);
		s.AppendFormat("FloraMass: {0}\n", state.FloraMass[i]);
		s.AppendFormat("FloraWater: {0}\n", state.FloraWater[i]);
		s.AppendFormat("FloraTemperature: {0}\n", state.FloraTemperature[i]);
		for (int j = 1; j < _waterLayers - 1; j++)
		{
			s.AppendFormat("WaterMass{0}: {1}\n", j, state.WaterMass[j][i]);
			s.AppendFormat("SaltMass{0}: {1}\n", j, state.SaltMass[j][i]);
			s.AppendFormat("WaterTemperature{0}: {1}\n", j, state.WaterTemperature[j][i]);
			s.AppendFormat("WaterVelocity{0}: {1}\n", j, state.WaterVelocity[j][i]);
		}
		for (int j = 1; j < _airLayers - 1; j++)
		{
			s.AppendFormat("AirTemperature{0}: {1}\n", j, state.AirTemperaturePotential[j][i]);
			s.AppendFormat("AirVapor{0}: {1}\n", j, state.AirVapor[j][i]);
			s.AppendFormat("AirVelocity{0}: {1}\n", j, state.AirVelocity[j][i]);
		}
		Debug.Log(s);
	}
	public void PrintDependentState(string title, int i, DependentState dependent)
	{
		StringBuilder s = new StringBuilder();
		s.AppendFormat("{0} Index: {1}\n", title, i);
		s.AppendFormat("Surface Elevation: {0}\n", dependent.LayerElevation[1][i]);
		s.AppendFormat("Water Depth: {0}\n", dependent.WaterLayerDepth[1][i]);
		s.AppendFormat("Ice Coverage: {0}\n", dependent.IceCoverage[i]);
		s.AppendFormat("Flora Coverage: {0}\n", dependent.FloraCoverage[i]);
		s.AppendFormat("Ice Terrain SA: {0}\n", dependent.SurfaceAreaIceTerrain[i]);
		s.AppendFormat("Ice Flora SA: {0}\n", dependent.SurfaceAreaIceFlora[i]);
		s.AppendFormat("Ice Water SA: {0}\n", dependent.SurfaceAreaIceWater[i]);
		s.AppendFormat("Air Ice SA: {0}\n", dependent.SurfaceAreaAirIce[i]);
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
