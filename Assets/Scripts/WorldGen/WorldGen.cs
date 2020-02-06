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
	public static void Generate(int seed, WorldGenData worldGenData, Icosphere icosphere, ref WorldData worldData, ref StaticState staticState, ref SimState state, ref DependentState dependent)
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
		dependent.Init(staticState.Count, worldData.AirLayers, worldData.WaterLayers);


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
			float waterDepth = math.max(0, -elevation);
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
			float cloudMass = Mathf.Pow(GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 2000, 0, 1), 1.0f) * Mathf.Pow(relativeHumidity, 1.0f);
			float cloudDropletMass = GetPerlinMinMax(pos.x, pos.y, pos.z, 1.0f, 2001, 0, 1) * relativeHumidity * Mathf.Pow(cloudMass, 2) * worldData.rainDropMaxSize;
			float cloudElevation = 5000;
			float cloudTemperature = airTemperaturePotential + WorldData.TemperatureLapseRate * cloudElevation;
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


			float minOceanTemperature = WorldData.FreezingTemperature + 0.1f;
			//float oceanDepthMinTemperature = 2000;
			float waterTemperatureSurface = Mathf.Max(WorldData.FreezingTemperature, airTemperaturePotential + 2);
			float waterTemperatureBottom = math.pow(waterTemperatureSurface - WorldData.FreezingTemperature, 1.0f / (1.0f + waterDepth / 500.0f)) + WorldData.FreezingTemperature;

			float layerElevation = surfaceElevation;
			float upperAirLayerHeight = (worldData.TropopauseElevation - surfaceElevation - worldData.BoundaryZoneElevation) / (worldData.AirLayers - 1);
			for (int j = 0; j < worldData.AirLayers; j++)
			{
				float2 wind = float2.zero;
				float windVertical = 0;
				float layerHeight = j == 0 ? worldData.BoundaryZoneElevation : upperAirLayerHeight;
				float layerMidHeight = layerElevation + layerHeight / 2;
				float airTemperature = airTemperaturePotential + WorldData.TemperatureLapseRate * layerMidHeight;
				float airMass = Atmosphere.GetAirMass(ref worldData, airPressure, layerMidHeight, airTemperature, WorldData.MolarMassAir, worldGenData.Gravity);
				float vaporMass = Atmosphere.GetAbsoluteHumidity(ref worldData, airTemperature, relativeHumidity, totalAirMass, inverseDewPointTemperatureRange);
				state.AirVapor[j][i] = vaporMass;
				state.AirEnergy[j][i] = airTemperature * (WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWaterVapor * vaporMass);
				state.AirVelocity[j][i] = wind;
				dependent.AirMass[j][i] = airMass;
				dependent.LayerElevation[j][i] = layerElevation;
				dependent.LayerHeight[j][i] = layerHeight;
				layerElevation += layerHeight;
			}

			float[] waterLayerDepth = new float[worldData.WaterLayers];
			int surfaceWaterLayer = worldData.WaterLayers - 1;
			float curDepth = math.min(waterDepth, worldData.ThermoclineDepth);
			waterLayerDepth[surfaceWaterLayer] = curDepth;
			dependent.WaterTemperature[surfaceWaterLayer][i] = waterTemperatureSurface;
			float deepWaterLayerDepths = (waterDepth - curDepth) / (worldData.WaterLayers - 1);
			for (int j = worldData.WaterLayers - 2; j >= 0; j--)
			{
				float layerDepth = math.min(deepWaterLayerDepths, waterDepth - curDepth);
				waterLayerDepth[j] = layerDepth;
				dependent.WaterTemperature[j][i] = math.pow(waterTemperatureSurface - WorldData.FreezingTemperature, 1.0f / ( 1.0f + curDepth / 500.0f)) + WorldData.FreezingTemperature;

				curDepth += layerDepth;
			}

			for (int j = 0; j < worldData.WaterLayers; j++)
			{
				float waterTemperature = dependent.WaterTemperature[j][i];
				float salinity = waterDepth == 0 ? 0 : (1.0f - Math.Abs(coord.y)) * (worldGenData.MaxSalinity - worldGenData.MinSalinity) + worldGenData.MinSalinity;
				//float deepSalinity = deepDepth == 0 ? 0 : Math.Abs(coord.y) * (worldGenData.MaxSalinity - worldGenData.MinSalinity) + worldGenData.MinSalinity;
				float waterAndSaltMass = GetWaterMass(worldData, waterLayerDepth[j], waterTemperature, salinity);
				float waterMass = waterAndSaltMass * (1.0f - salinity);
				float saltMass = waterAndSaltMass * salinity;
				float2 current = float2.zero;
				float currentVertical = 0;

				state.WaterMass[j][i] = waterMass;
				state.WaterSaltMass[j][i] = saltMass;
				state.WaterEnergy[j][i] = Atmosphere.GetWaterEnergy(waterTemperature, waterMass, saltMass);
				state.WaterVelocity[j][i] = current;
				waterSaltMass[i] = new WaterSaltMass()
				{
					WaterMass = waterSaltMass[i].WaterMass + waterMass,
					SaltMass = waterSaltMass[i].SaltMass + saltMass
				};
			}

			float2 cloudVelocity = float2.zero;

			float surfaceWaterMass = state.WaterMass[worldData.WaterLayers-1][i];
			float iceMass = math.min(surfaceWaterMass, worldData.FullIceCoverage * WorldData.MassIce * Mathf.Clamp01(-(airTemperaturePotential - WorldData.FreezingTemperature) / 10));
			state.IceEnergy[i] = iceMass * WorldData.FreezingTemperature * WorldData.SpecificHeatIce;
			state.WaterEnergy[0][i] -= iceMass;
			float waterAndIceDepth = waterDepth + iceMass / WorldData.MassIce;

			float waterCoverage = Mathf.Clamp01(waterDepth / worldData.FullWaterCoverage);
			float iceCoverage = Mathf.Clamp01(iceMass / (WorldData.MassIce * worldData.FullIceCoverage));

			if (waterDepth == 0)
			{
				vegetation = soilFertility * (groundWater) * (1.0f - waterCoverage) * (1.0f - iceCoverage) * Mathf.Clamp01((airTemperatureSurface - worldData.MinTemperatureCanopy) / (worldData.MaxTemperatureCanopy - worldData.MinTemperatureCanopy));
			}

			// TODO: ground water energy should probably be tracked independently
			float terrainTemperature;
			if (waterCoverage >= 1)
			{
				terrainTemperature = waterTemperatureBottom;
			}
			else
			{
				terrainTemperature = airTemperatureSurface;
			}
			float groundEnergy = Atmosphere.GetTerrainEnergy(ref worldData, terrainTemperature, groundWater, soilFertility, vegetation * worldData.inverseFullCanopyCoverage);



			state.Terrain[i] = new CellTerrain()
			{
				Elevation = elevation,
				GroundWaterDepth = groundWaterDepth,
				Roughness = roughness,
				SoilFertility = soilFertility,
				Vegetation = vegetation,
			};
			state.CloudDropletMass[i] = cloudDropletMass;
			state.CloudElevation[i] = cloudElevation;
			state.CloudMass[i] = cloudMass;
			state.CloudEnergy[i] = cloudTemperature * WorldData.SpecificHeatWater * cloudMass;
			state.CloudVelocity[i] = cloudVelocity;
			state.IceMass[i] = iceMass;
			state.TerrainEnergy[i] = groundEnergy;
			state.GroundWater[i] = groundWater;
		}

		int batchCount = 100;
		NativeList<JobHandle> updateDependenciesJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
		for (int j = 0; j < worldData.WaterLayers; j++)
		{
			var updateDependentWaterLayerJob = new UpdateDependentWaterLayerJob()
			{
				Temperature = dependent.WaterTemperature[j],
				Salinity = dependent.WaterSalinity[j],
				WaterCoverage = dependent.WaterCoverage[j],
				Energy = state.WaterEnergy[j],
				SaltMass = state.WaterSaltMass[j],
				WaterMass = state.WaterMass[j],
				Terrain = state.Terrain,
				worldData = worldData
			};
			updateDependenciesJobHandles.Add(updateDependentWaterLayerJob.Schedule(staticState.Count, batchCount));
		}
		JobHandle lowerAirHandle = default(JobHandle);
		for (int j = 0; j < worldData.AirLayers; j++)
		{
			var updateDependentAirLayerJob = new UpdateDependentAirLayerJob()
			{
				Temperature = dependent.AirTemperature[j],
				Pressure = dependent.AirPressure[j],
				RelativeHumidity = dependent.AirHumidityRelative[j],
				AbsoluteHumidity = dependent.AirHumidityAbsolute[j],
				LayerElevation = dependent.LayerElevation[j],
				LayerHeight = dependent.LayerHeight[j],
				AirMass = dependent.AirMass[j],
				AirEnergy = state.AirEnergy[j],
				AirVapor = state.AirVapor[j],
				VaporMass = state.AirVapor[j],
				IceMass = state.IceMass,
				Gravity = state.PlanetState.Gravity,
				WorldData = worldData,
			};
			var h = updateDependentAirLayerJob.Schedule(staticState.Count, batchCount);
			updateDependenciesJobHandles.Add(h);
			if (j == 0)
			{
				lowerAirHandle = h;
			}
		}

		var updateDependentStateJob = new UpdateDependentStateJob()
		{
			CloudCoverage = dependent.CloudCoverage,
			IceCoverage = dependent.IceCoverage,
			SurfaceElevation = dependent.SurfaceElevation,
			VegetationCoverage = dependent.VegetationCoverage,
			WaterDepth = dependent.WaterDepth,
			CloudTemperature = dependent.CloudTemperature,
			IceTemperature = dependent.IceTemperature,
			TerrainTemperature = dependent.TerrainTemperature,
			SurfaceAirTemperature = dependent.SurfaceAirTemperature,

			WaterSaltMass = waterSaltMass,
			CloudMass = state.CloudMass,
			CloudEnergy = state.CloudEnergy,
			IceMass = state.IceMass,
			IceEnergy = state.IceEnergy,
			Terrain = state.Terrain,
			TerrainEnergy = state.TerrainEnergy,
			GroundWater = state.GroundWater,
			worldData = worldData,
			LowerAirTemperature = dependent.AirTemperature[0],
			lowerAirHeight = dependent.LayerHeight[0]
		};
		updateDependenciesJobHandles.Add(updateDependentStateJob.Schedule(staticState.Count, batchCount, lowerAirHandle));
		JobHandle.CompleteAll(updateDependenciesJobHandles);
		updateDependenciesJobHandles.Dispose();
		waterSaltMass.Dispose();
	}

	static public float GetWaterMass(WorldData worldData, float depth, float temperature, float salinityPSU)
	{
		float density = worldData.waterDensity + worldData.OceanDensityPerDegree * (temperature - WorldData.FreezingTemperature) + worldData.OceanDensityPerSalinity * salinityPSU;
		return depth * density;
	}

}
