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
		var pressureGradient = new NativeArray<float2>(staticState.Count, Allocator.TempJob);
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
			float cloudMass = Mathf.Pow(GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 2000, 0, 1), 1.0f) * Mathf.Pow(relativeHumidity, 1.0f);
			float cloudDropletMass = WorldData.MassWater * math.pow(worldData.rainDropMinSize, 3) * 4 / 3 *math.PI;
			float cloudElevation = 5000;
			float cloudTemperature = airTemperaturePotential + WorldData.TemperatureLapseRate * cloudElevation;
//			float totalAirMass = worldGenData.TroposphereMass * (1.0f - surfaceElevation / worldData.TropopauseElevation);


			float waterTemperatureSurface = Mathf.Max(WorldData.FreezingTemperature, airTemperaturePotential + 2);
			float waterTemperatureBottom = waterDepth > 0 ? math.pow(waterTemperatureSurface - WorldData.FreezingTemperature, 1.0f / (1.0f + waterDepth / worldData.WaterTemperatureDepthFalloff)) + WorldData.FreezingTemperature : waterTemperatureSurface;

			float layerElevation = surfaceElevation;
			float upperAirLayerHeight = (worldData.TropopauseElevation - surfaceElevation - worldData.BoundaryZoneElevation) / (worldData.AirLayers - 1);
			for (int j = 0; j < worldData.AirLayers; j++)
			{
				//				float layerHeight = j == 0 ? worldData.BoundaryZoneElevation : upperAirLayerHeight;
				float layerHeight;
				if (j == 0)
				{
					layerHeight = worldData.BoundaryZoneElevation;
				} else if (j == 1)
				{
					layerHeight = worldData.TropopauseElevation - surfaceElevation - worldData.BoundaryZoneElevation;
				} else
				{
					layerHeight = worldData.stratosphereElevation - worldData.TropopauseElevation;
				}
				float layerMidHeight = layerElevation + layerHeight / 2;
				float airTemperature = airTemperaturePotential + WorldData.TemperatureLapseRate * layerMidHeight;
				float airPressure = GetStandardPressureAtElevation(layerMidHeight, airTemperature, state.PlanetState.Gravity);
				float airMass = Atmosphere.GetStandardAirMass(layerElevation, layerHeight, state.PlanetState.Gravity);
				float vaporMass = GetWaterVaporMass(ref worldData, airTemperature, relativeHumidity, airMass, inverseDewPointTemperatureRange);
				state.AirVapor[j][i] = vaporMass;
				state.Wind[j][i] = 0;
				state.AirTemperature[j][i] = airTemperature;
				dependent.AirMass[j][i] = airMass;
				dependent.LayerElevation[j][i] = layerElevation;
				dependent.LayerHeight[j][i] = layerHeight;
				layerElevation += layerHeight;
			}

			float[] waterLayerDepth = new float[worldData.WaterLayers];
			int surfaceWaterLayer = worldData.WaterLayers - 1;
			float curDepth = math.min(waterDepth, worldData.ThermoclineDepth);
			waterLayerDepth[surfaceWaterLayer] = curDepth;
			state.WaterTemperature[surfaceWaterLayer][i] = waterTemperatureSurface;
			float deepWaterLayerDepths = (waterDepth - curDepth) / (worldData.WaterLayers - 1);
			for (int j = worldData.WaterLayers - 2; j >= 0; j--)
			{
				float layerDepth = math.min(deepWaterLayerDepths, waterDepth - curDepth);
				waterLayerDepth[j] = layerDepth;
				state.WaterTemperature[j][i] = math.pow(waterTemperatureSurface - WorldData.FreezingTemperature, 1.0f / ( 1.0f + curDepth / 500.0f)) + WorldData.FreezingTemperature;

				curDepth += layerDepth;
			}

			for (int j = 0; j < worldData.WaterLayers; j++)
			{
				float waterTemperature = state.WaterTemperature[j][i];
				float salinity = waterDepth == 0 ? 0 : (1.0f - Math.Abs(coord.y)) * (worldGenData.MaxSalinity - worldGenData.MinSalinity) + worldGenData.MinSalinity;
				//float deepSalinity = deepDepth == 0 ? 0 : Math.Abs(coord.y) * (worldGenData.MaxSalinity - worldGenData.MinSalinity) + worldGenData.MinSalinity;
				float waterAndSaltMass = GetWaterMass(worldData, waterLayerDepth[j], waterTemperature, salinity);
				float waterMass = waterAndSaltMass * (1.0f - salinity);
				float saltMass = waterAndSaltMass * salinity;

				state.WaterMass[j][i] = waterMass;
				state.WaterSaltMass[j][i] = saltMass;
				state.WaterVelocity[j][i] = 0;
				if (waterMass == 0)
				{
					state.WaterTemperature[j][i] = 0;
				}
				waterSaltMass[i] = new WaterSaltMass()
				{
					WaterMass = waterSaltMass[i].WaterMass + waterMass,
					SaltMass = waterSaltMass[i].SaltMass + saltMass
				};
			}

			float2 cloudVelocity = float2.zero;

			float surfaceWaterMass = state.WaterMass[worldData.WaterLayers-1][i];
			float iceMass = math.min(surfaceWaterMass, worldData.FullIceCoverage * WorldData.MassIce * Mathf.Clamp01(-(airTemperaturePotential - WorldData.FreezingTemperature) / 10));
			state.IceTemperature[i] = iceMass > 0 ? airTemperaturePotential : 0;
			float waterAndIceDepth = waterDepth + iceMass / WorldData.MassIce;

			float waterCoverage = Mathf.Clamp01(waterDepth / worldData.FullWaterCoverage);
			float iceCoverage = Mathf.Clamp01(iceMass / (WorldData.MassIce * worldData.FullIceCoverage));

			if (waterDepth == 0)
			{
				vegetation = 
					worldData.FullVegetationCoverage 
					* soilFertility 
					* (1.0f - waterCoverage) 
					* (1.0f - iceCoverage) 
					* GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 410)
					* math.sin(math.PI* math.saturate((airTemperatureSurface - worldData.MinTemperatureCanopy) / (worldData.MaxTemperatureCanopy - worldData.MinTemperatureCanopy)));
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


			state.Terrain[i] = new CellTerrain()
			{
				Elevation = elevation,
				Roughness = roughness,
				SoilFertility = soilFertility,
				Vegetation = vegetation,
			};
			state.CloudDropletMass[i] = cloudDropletMass;
			state.CloudElevation[i] = cloudElevation;
			state.CloudMass[i] = cloudMass;
			state.CloudTemperature[i] = cloudTemperature;
			state.CloudVelocity[i] = cloudVelocity;
			state.IceMass[i] = iceMass;
			state.TerrainTemperature[i] = terrainTemperature;
		}

		int batchCount = 100;
		NativeList<JobHandle> updateDependenciesJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
		for (int j = 0; j < worldData.WaterLayers; j++)
		{
			var updateDependentWaterLayerJob = new UpdateDependentWaterLayerJob()
			{
				Salinity = dependent.WaterSalinity[j],
				WaterCoverage = dependent.WaterCoverage[j],
				PotentialEnergy = dependent.WaterPotentialEnergy[j],

				SaltMass = state.WaterSaltMass[j],
				WaterMass = state.WaterMass[j],
				Temperature = state.WaterTemperature[j],
				Terrain = state.Terrain,
				worldData = worldData
			};
			updateDependenciesJobHandles.Add(updateDependentWaterLayerJob.Schedule(staticState.Count, batchCount));
		}
		JobHandle lowerAirHandle = default(JobHandle);
		JobHandle updateDependentAirLayerJobHandle = default(JobHandle);
		for (int j = 0; j < worldData.AirLayers; j++)
		{
			var updateDependentAirLayerJob = new UpdateDependentAirLayerJob()
			{
				Pressure = dependent.AirPressure[j],
				RelativeHumidity = dependent.AirHumidityRelative[j],
				AbsoluteHumidity = dependent.AirHumidityAbsolute[j],
				LayerElevation = dependent.LayerElevation[j],
				LayerHeight = dependent.LayerHeight[j],
				AirMass = dependent.AirMass[j],
				PotentialEnergy = dependent.AirPotentialEnergy[j],
				AirHumidityRelativeCloud = dependent.AirHumidityRelativeCloud,
				AirLayerCloud = dependent.AirLayerCloud,
				AirMassCloud = dependent.AirMassCloud,
				AirPressureCloud = dependent.AirPressureCloud,
				AirTemperatureCloud = dependent.AirTemperatureCloud,
				AirVaporCloud = dependent.AirVaporCloud,
				PressureGradientAtCloudElevation = dependent.PressureGradientAtCloudElevation,

				PressureGradient = pressureGradient,
				CloudDropletMass = state.CloudDropletMass,
				CloudElevation = state.CloudElevation,
				CloudMass = state.CloudMass,
				CloudTemperature = state.CloudTemperature,
				AirTemperature = state.AirTemperature[j],
				VaporMass = state.AirVapor[j],
				IceMass = state.IceMass,
				Gravity = state.PlanetState.Gravity,
				DewPointZero = worldData.DewPointZero,
				InverseDewPointTemperatureRange = worldData.inverseDewPointTemperatureRange,
				WaterVaporMassToAirMassAtDewPoint = worldData.WaterVaporMassToAirMassAtDewPoint,
				LayerIndex = j,
			};
			updateDependentAirLayerJobHandle = updateDependentAirLayerJob.Schedule(staticState.Count, batchCount, updateDependentAirLayerJobHandle);
			updateDependenciesJobHandles.Add(updateDependentAirLayerJobHandle);
			if (j == 0)
			{
				lowerAirHandle = updateDependentAirLayerJobHandle;
			}
		}

		var updateDependentStateJob = new UpdateDependentStateJob()
		{
			CloudCoverage = dependent.CloudCoverage,
			IceCoverage = dependent.IceCoverage,
			SurfaceElevation = dependent.SurfaceElevation,
			VegetationCoverage = dependent.VegetationCoverage,
			WaterDepth = dependent.WaterDepth,
			SurfaceAirTemperature = dependent.SurfaceAirTemperature,
			CloudEnergy = dependent.CloudEnergy,
			IceEnergy = dependent.IceEnergy,

			WaterSaltMass = waterSaltMass,
			CloudMass = state.CloudMass,
			CloudTemperature = state.CloudTemperature,
			IceMass = state.IceMass,
			IceTemperature = state.IceTemperature,
			Terrain = state.Terrain,
			worldData = worldData,
			LowerAirTemperature = state.AirTemperature[0],
			lowerAirHeight = dependent.LayerHeight[0]
		};
		updateDependenciesJobHandles.Add(updateDependentStateJob.Schedule(staticState.Count, batchCount, lowerAirHandle));
		JobHandle.CompleteAll(updateDependenciesJobHandles);
		updateDependenciesJobHandles.Dispose();
		waterSaltMass.Dispose();
		pressureGradient.Dispose();
	}

	static public float GetWaterMass(WorldData worldData, float depth, float temperature, float salinityPSU)
	{
		float density = WorldData.DensityWater + worldData.WaterDensityPerDegree * (temperature - WorldData.FreezingTemperature) + worldData.WaterDensityPerSalinity * salinityPSU;
		return depth * density;
	}

	static public float GetWaterVaporMass(ref WorldData worldData, float temperature, float relativeHumidity, float airMass, float inverseDewPointTemperatureRange)
	{
		float maxWaterVaporPerKilogramAtmosphere = worldData.WaterVaporMassToAirMassAtDewPoint * Utils.Sqr(math.max(0, (temperature - worldData.DewPointZero) * inverseDewPointTemperatureRange));
		float maxHumidity = maxWaterVaporPerKilogramAtmosphere * airMass;
		if (maxHumidity <= 0)
		{
			return 0;
		}
		float humidity = relativeHumidity * maxHumidity;
		return humidity;
	}

	static public float GetStandardPressureAtElevation(float elevation, float temperature, float gravity)
	{
		float pressure = WorldData.StaticPressure * (float)Math.Pow(1.0f + WorldData.TemperatureLapseRate / WorldData.StdTemp * elevation, gravity * WorldData.PressureExponent * WorldData.MolarMassAir);
		return pressure;
	}



}
