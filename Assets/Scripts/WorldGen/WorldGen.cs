using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;
using UnityEngine;
using Unity.Jobs;
using Unity.Collections;

public static class WorldGen {
	static System.Random _random;
	static FastNoise _noise;

	private static float GetPerlinMinMax(float x, float y, float frequency, float hash, float min, float max)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency) + 1.0f) * (max - min) / 2 + min;
	}
	private static float GetPerlinMinMax(float x, float y, float z, float frequency, float hash, float min, float max)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency, z * frequency) + 1.0f) * (max - min) / 2 + min;
	}
	private static float GetPerlinNormalized(float x, float y, float frequency, float hash)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency) + 1.0f) / 2;
	}
	private static float GetPerlinNormalized(float x, float y, float z, float frequency, float hash)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency, z * frequency) + 1.0f) / 2;
	}
	private static float GetPerlin(float x, float y, float frequency, float hash)
	{
		return _noise.GetPerlin(x * frequency + hash, y * frequency);
	}
	private static float GetPerlin(float x, float y, float z, float frequency, float hash)
	{
		return _noise.GetPerlin(x * frequency + hash, y * frequency, z * frequency);
	}
	public static void Generate(int seed, WorldGenData worldGenData, Icosphere icosphere, ref WorldData worldData, ref StaticState staticState, ref SimState state, ref DependentState dependentState)
	{

		float inversePI = 1.0f / math.PI;
		_noise = new FastNoise(seed);
		_noise.SetFrequency(10);
		_random = new System.Random(seed);

		staticState.Init(worldGenData.Radius, icosphere, ref worldData);

		state.PlanetState.Gravity = worldGenData.Gravity;
		state.PlanetState.DistanceToSun = worldGenData.DistanceToSun;
		state.PlanetState.Rotation = math.radians(math.float3(worldGenData.TiltAngle, 0, 0));
		state.PlanetState.Position = math.float3(1, 0, 0) * worldGenData.DistanceToSun;
		state.PlanetState.SpinSpeed = math.PI * 2 / (worldGenData.SpinTime * 60 * 60);
		state.PlanetState.OrbitSpeed = math.PI * 2 / worldGenData.OrbitTime;
		state.PlanetState.AngularSpeed = math.PI * 2 / (worldGenData.SpinTime * 60 * 60);
		state.PlanetState.GeothermalHeat = worldGenData.GeothermalHeat;
		state.PlanetState.SolarRadiation = worldGenData.SolarRadiation;
		state.PlanetState.StratosphereMass = worldGenData.StratosphereMass;
		state.PlanetState.CarbonDioxide = worldGenData.CarbonDioxide;

		float inverseDewPointTemperatureRange = 1.0f / worldData.DewPointTemperatureRange;

		var waterSaltMass = new NativeArray<WaterSaltMass>(staticState.Count, Allocator.TempJob);

		for (int i = 0; i < staticState.Count; i++)
		{
			var coord = staticState.Coordinate[i] * 2 * inversePI;
			var pos = staticState.SphericalPosition[i];
			var coordNormalized = (coord + 1) / 2;
			float roughness = math.max(1, math.pow(GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 660), 3) * worldGenData.MaxRoughness);
			float elevation =
				0.6f * GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 0, worldGenData.MinElevation, worldGenData.MaxElevation) +
				0.3f * GetPerlinMinMax(pos.x, pos.y, pos.z, 0.5f, 0, worldGenData.MinElevation, worldGenData.MaxElevation) +
				0.1f * GetPerlinMinMax(pos.x, pos.y, pos.z, 2.0f, 0, worldGenData.MinElevation, worldGenData.MaxElevation);
			float waterDepth = math.max(0, -elevation) + math.max(0, GetPerlinMinMax(pos.x, pos.y, pos.z, 1.0f, 25430, -100, 100));
			//cell.WaterDepth = GetPerlinNormalized(coordNormalized.x, coordNormalized.y, 2.5f, 110) * 2000;

			float surfaceElevation = elevation + waterDepth;

			float soilFertility =
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.2f, 6630) +
				0.3f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 435) +
				0.2f * GetPerlinNormalized(pos.x, pos.y, pos.z, 1f, 8740);

			float vegetation = elevation < 0 ? 0 : math.pow(math.sqrt(soilFertility) * math.sqrt(GetPerlinNormalized(pos.x, pos.y, pos.z, 1f, 6630423)), 3);

			float cloudCoverage =
				math.pow(
					0.75f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 30) +
					0.25f * GetPerlinNormalized(pos.x, pos.y, pos.z, 1.0f, 30),
					1);

			float relativeHumidity =
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 40) +
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 40);


			float troposphereColumnHeight = worldData.TropopauseElevation - surfaceElevation;

			float regionalTemperatureVariation =
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 15460, -5, 5) +
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.2f, 431, -10, 10) +
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.5f, 6952, -10, 10);
			float airTemperaturePotential =
				regionalTemperatureVariation + GetPerlinMinMax(pos.x, pos.y, pos.z, 0.15f, 80, -5, 5) +
				(1.0f - coord.y * coord.y) * (worldGenData.MaxTemperature - worldGenData.MinTemperature) + worldGenData.MinTemperature;
			float airTemperatureSurface = airTemperaturePotential + WorldData.TemperatureLapseRate * surfaceElevation;

			float airPressure = WorldData.StaticPressure - regionalTemperatureVariation * 1000;
			float cloudMass = Mathf.Pow(GetPerlinMinMax(pos.x, pos.y, pos.z, 1.0f, 2000, 0, 1), 0.85f) * Mathf.Pow(relativeHumidity, 0.5f);
			float cloudDropletMass = GetPerlinMinMax(pos.x, pos.y, pos.z, 1.0f, 2001, 0, 1) * relativeHumidity * Mathf.Pow(cloudMass, 2) * worldData.rainDropMaxSize;
			float cloudElevation = 5000;
			float cloudTemperature = airTemperaturePotential * WorldData.TemperatureLapseRate * cloudElevation;
			float totalAirMass = worldGenData.TroposphereMass * (1.0f - surfaceElevation / worldData.TropopauseElevation);

			float groundWaterDepth = GetPerlinMinMax(pos.x, pos.y, pos.z, 1f, 60423, worldGenData.GroundWaterDepthMin, worldGenData.GroundWaterDepthMax);
			float maxGroundWater = groundWaterDepth * soilFertility * WorldData.MassWater * worldData.MaxSoilPorousness;
			float groundWater;
			if (waterDepth > 0)
			{
				groundWater = maxGroundWater;
			}
			else
			{
				groundWater = maxGroundWater * 0.2f;
			}

			//float shallowDepth = Mathf.Min(data.DeepOceanDepth, depth);
			//float deepDepth = Mathf.Max(0, depth - data.DeepOceanDepth);

			float minOceanTemperature = WorldData.FreezingTemperature + 0.1f;
			//float oceanDepthMinTemperature = 2000;
			float waterTemperatureSurface = Mathf.Max(WorldData.FreezingTemperature, airTemperatureSurface + 2);
			float waterTemperatureBottom = waterTemperatureSurface;


			for (int j = 0; j < worldData.AirLayers; j++)
			{
				float2 wind = float2.zero;
				float windVertical = 0;
				float layerElevation = surfaceElevation;
				float layerHeight = j == 0 ? worldData.BoundaryZoneElevation : (worldData.TropopauseElevation - surfaceElevation - worldData.BoundaryZoneElevation) / (worldData.AirLayers - 1);
				float airTemperature = airTemperaturePotential + WorldData.TemperatureLapseRate * layerElevation;
				float airMass = 0;
				float vaporMass = Atmosphere.GetAbsoluteHumidity(ref worldData, airTemperature, relativeHumidity, totalAirMass, inverseDewPointTemperatureRange);
				state.AirVapor[j][i] = vaporMass;
				state.AirEnergy[j][i] = airTemperature * (WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWaterVapor * vaporMass);
				state.AirVelocity[j][i] = wind;
			}
			for (int j = 0; j < worldData.WaterLayers; j++)
			{
				float layerDepth = j == 0 ? waterDepth : 0;
				float waterTemperature = waterTemperatureSurface;
				float salinity = waterDepth == 0 ? 0 : (1.0f - Math.Abs(coord.y)) * (worldGenData.MaxSalinity - worldGenData.MinSalinity) + worldGenData.MinSalinity;
				//float deepSalinity = deepDepth == 0 ? 0 : Math.Abs(coord.y) * (worldGenData.MaxSalinity - worldGenData.MinSalinity) + worldGenData.MinSalinity;
				float waterAndSaltMass = GetWaterMass(worldData, layerDepth, waterTemperature, salinity);
				float waterMass = waterAndSaltMass * (1.0f - salinity);
				float saltMass = waterAndSaltMass * salinity;
				float2 current = float2.zero;
				float currentVertical = 0;

				state.WaterMass[j][i] = waterMass;
				state.WaterSaltMass[j][i] = saltMass;
				state.WaterEnergy[j][i] = waterTemperature;
				state.WaterVelocity[j][i] = current;
				waterSaltMass[i] = new WaterSaltMass()
				{
					WaterMass = waterSaltMass[i].WaterMass + waterMass,
					SaltMass = waterSaltMass[i].SaltMass + saltMass
				};
			}

			float2 cloudVelocity = float2.zero;

			float surfaceWaterMass = state.WaterEnergy[0][i];
			float iceMass = math.min(surfaceWaterMass, worldData.FullIceCoverage * WorldData.MassIce * Mathf.Clamp01(-(airTemperatureSurface - WorldData.FreezingTemperature) / 10));
			state.WaterEnergy[0][i] -= iceMass;
			float waterAndIceDepth = waterDepth + iceMass / WorldData.MassIce;

			////	float deepOceanTemperature = (state.ShallowWaterTemperature - minOceanTemperature) * (1.0f - Mathf.Pow(Mathf.Clamp01(deepDepth / oceanDepthMinTemperature), 2f)) + minOceanTemperature;
			//float deepOceanTemperature = minOceanTemperature;
			//float deepOceanMass = GetWaterMass(world, deepDepth, deepOceanTemperature, deepSalinity);
			//cell.DeepWaterMass = deepOceanMass * (1.0f - deepSalinity);
			//cell.DeepSaltMass = deepOceanMass * deepSalinity;
			//cell.DeepWaterEnergy = Atmosphere.GetWaterEnergy(world, deepOceanTemperature, state.DeepWaterMass, state.DeepSaltMass);
			//cell.DeepWaterDensity = Atmosphere.GetWaterDensity(world, state.DeepWaterEnergy, state.DeepSaltMass, state.DeepWaterMass);


			float waterCoverage = Mathf.Clamp01(waterDepth / worldData.FullWaterCoverage);
			float iceCoverage = Mathf.Clamp01(iceMass / (WorldData.MassIce * worldData.FullIceCoverage));

			if (waterDepth == 0)
			{
				vegetation = soilFertility * (groundWater) * (1.0f - waterCoverage) * (1.0f - iceCoverage) * Mathf.Clamp01((airTemperatureSurface - worldData.MinTemperatureCanopy) / (worldData.MaxTemperatureCanopy - worldData.MinTemperatureCanopy));
			}

			// TODO: ground water energy should probably be tracked independently
			float groundRadiationRate = Atmosphere.GetRadiationRate(airTemperatureSurface, WorldData.EmissivityDirt);
			float groundWaterEnergy = groundWater * WorldData.SpecificHeatWater * worldData.maxGroundWaterTemperature;
			float landMass = (WorldData.MassSand - WorldData.MassSoil) * soilFertility + WorldData.MassSoil;
			float heatingDepth = soilFertility * worldData.SoilHeatDepth;
			float terrainTemperature;
			if (waterCoverage >= 1)
			{
				terrainTemperature = waterTemperatureBottom;
			}
			else
			{
				terrainTemperature = airTemperatureSurface;
			}
			float soilEnergy = terrainTemperature * WorldData.SpecificHeatSoil * heatingDepth * landMass;
			float groundEnergy = groundWaterEnergy + soilEnergy;



			state.Terrain[i] = new CellTerrain()
			{
				Elevation = elevation,
				GroundWater = groundWater,
				GroundWaterDepth = groundWaterDepth,
				Roughness = roughness,
				SoilFertility = soilFertility,
				Vegetation = vegetation,
			};
			state.CloudDropletMass[i] = cloudDropletMass;
			state.CloudElevation[i] = cloudElevation;
			state.CloudMass[i] = cloudMass;
			state.CloudEnergy[i] = cloudTemperature;
			state.CloudVelocity[i] = cloudVelocity;
			state.IceMass[i] = iceMass;
			state.IceEnergy[i] = WorldData.FreezingTemperature;
			state.TerrainEnergy[i] = terrainTemperature;
		}

		dependentState.Init(staticState.Count, worldData.AirLayers, worldData.WaterLayers);

		var updateDependentStateJob = new UpdateDependentStateJob()
		{
			CloudCoverage = dependentState.CloudCoverage,
			IceCoverage = dependentState.IceCoverage,
			SurfaceElevation = dependentState.SurfaceElevation,
			VegetationCoverage = dependentState.VegetationCoverage,
			WaterCoverage = dependentState.WaterCoverage,
			WaterDepth = dependentState.WaterDepth,

			WaterSaltMass = waterSaltMass,
			CloudMass = state.CloudMass,
			IceMass = state.IceMass,
			Terrain = state.Terrain,
			worldData = worldData
		};
		var updateDependentStateJobHandle = updateDependentStateJob.Schedule(staticState.Count, 100);
		updateDependentStateJobHandle.Complete();

		waterSaltMass.Dispose();
	}

	static public float GetWaterMass(WorldData worldData, float depth, float temperature, float salinityPSU)
	{
		float density = worldData.waterDensity + worldData.OceanDensityPerDegree * (temperature - WorldData.FreezingTemperature) + worldData.OceanDensityPerSalinity * salinityPSU;
		return depth * density;
	}

}
