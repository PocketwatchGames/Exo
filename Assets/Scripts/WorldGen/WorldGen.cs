﻿//#define ASYNC_WORLDGEN

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;
using UnityEngine;
using Unity.Jobs;
using Unity.Collections;
using Unity.Burst;

public static class WorldGen {
	static System.Random _random;
	static FastNoise _noise;

	private static void InitSync(WorldGenData worldGenData, Icosphere icosphere, ref WorldData worldData, ref StaticState staticState, ref SimState state, ref DependentState dependent)
	{
		staticState.Init(worldGenData.Radius, icosphere, worldGenData.StratosphereMass, ref worldData);

		state.PlanetState.Gravity = worldGenData.Gravity;
		state.PlanetState.DistanceToSun = worldGenData.DistanceToSun;
		state.PlanetState.Rotation = math.radians(math.float3(worldGenData.TiltAngle, 0, 0));
		state.PlanetState.Position = math.float3(1, 0, 0) * worldGenData.DistanceToSun;
		state.PlanetState.SpinSpeed = math.PI * 2 / (worldGenData.SpinTime * 60 * 60);
		state.PlanetState.OrbitSpeed = math.PI * 2 / worldGenData.OrbitTime;
		state.PlanetState.AngularSpeed = math.PI * 2 / (worldGenData.SpinTime * 60 * 60);
		state.PlanetState.GeothermalHeat = worldGenData.GeothermalHeat;
		state.PlanetState.SolarRadiation = worldGenData.SolarRadiation;
		state.PlanetState.CarbonDioxide = worldGenData.CarbonDioxide;

		float inverseDewPointTemperatureRange = 1.0f / worldData.DewPointTemperatureRange;

		dependent.Init(staticState.Count, worldData.AirLayers, worldData.WaterLayers);
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private struct WorldGenInitJob {

		public NativeArray<float> relativeHumidity;
		public NativeArray<float> potentialTemperature;
		public NativeArray<float> CloudMass;
		public NativeArray<CellTerrain> Terrain;
		public NativeArray<float> LayerElevationBase;

		[ReadOnly] public NativeArray<float2> Coordinate;
		[ReadOnly] public NativeArray<float3> SphericalPosition;
		[ReadOnly] public FastNoise noise;
		[ReadOnly] public float MaxRoughness;
		[ReadOnly] public float MinElevation;
		[ReadOnly] public float MaxElevation;
		[ReadOnly] public float MaxTemperature;
		[ReadOnly] public float MinTemperature;
		[ReadOnly] public float FullCoverageVegetation;
		[ReadOnly] public float MinTemperatureCanopy;
		[ReadOnly] public float MaxTemperatureCanopy;

		public void Run(int count)
		{
			for (int i=0; i<count;i++)
			{
				Execute(i);
			}
		}

		public void Execute(int i)
		{
			float inversePI = 1.0f / math.PI;

			var coord = Coordinate[i] * 2 * inversePI;
			var pos = SphericalPosition[i];
			float roughness = math.max(1, math.pow(GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 660), 3) * MaxRoughness);
			float elevation =
				0.6f * GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 0, MinElevation, MaxElevation) +
				0.3f * GetPerlinMinMax(pos.x, pos.y, pos.z, 0.5f, 0, MinElevation, MaxElevation) +
				0.1f * GetPerlinMinMax(pos.x, pos.y, pos.z, 2.0f, 0, MinElevation, MaxElevation);
			float soilFertility =
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.2f, 6630) +
				0.3f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 435) +
				0.2f * GetPerlinNormalized(pos.x, pos.y, pos.z, 1f, 8740);

			relativeHumidity[i] =
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 40) +
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 40);

			float surfaceElevation = math.max(0, elevation);

			float regionalTemperatureVariation =
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 15460, -5, 5) +
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.2f, 431, -10, 10) +
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.5f, 6952, -10, 10);
			potentialTemperature[i] =
				regionalTemperatureVariation + GetPerlinMinMax(pos.x, pos.y, pos.z, 0.15f, 80, -5, 5) +
				(1.0f - coord.y * coord.y) * (MaxTemperature - MinTemperature) + MinTemperature;

			float airTemperatureSurface = potentialTemperature[i] + WorldData.TemperatureLapseRate * surfaceElevation;
			float vegetation = 0;
			if (elevation > 0)
			{
				vegetation =
					FullCoverageVegetation
					* soilFertility
					* GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 410)
					* math.sin(math.PI * math.saturate((airTemperatureSurface - MinTemperatureCanopy) / (MaxTemperatureCanopy - MinTemperatureCanopy)));
			}
			float cloudMass = Mathf.Pow(GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 2000, 0, 1), 1.0f) * Mathf.Pow(relativeHumidity[i], 1.0f);

			LayerElevationBase[i] = surfaceElevation;
			Terrain[i] = new CellTerrain()
			{
				Elevation = elevation,
				Roughness = roughness,
				SoilFertility = soilFertility,
				Vegetation = vegetation,
			};

			CloudMass[i] = cloudMass;
		}
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private struct WorldGenJob : IJobParallelFor {

		public NativeArray<float> IceTemperature;
		public NativeArray<float> CloudDropletMass;
		public NativeArray<float3> CloudVelocity;
		public NativeArray<float> IceMass;
		public NativeArray<float> TerrainTemperature;
		public NativeArray<float> WaterTemperatureSurface;
		public NativeArray<float> WaterTemperatureBottom;
		
		[ReadOnly] public NativeArray<float> TemperaturePotential;
		[ReadOnly] public NativeArray<CellTerrain> terrain;
		[ReadOnly] public float inversePI;
		[ReadOnly] public float TropopauseElevation;
		[ReadOnly] public float rainDropMinSize;
		[ReadOnly] public float WaterTemperatureDepthFalloff;
		[ReadOnly] public float BoundaryZoneElevation;
		[ReadOnly] public float FullWaterCoverage;

		public void Execute(int i)
		{
			float elevation = terrain[i].Elevation;

			float waterDepth = math.max(0, -elevation);
			float surfaceElevation = elevation + waterDepth;

			float troposphereColumnHeight = TropopauseElevation - surfaceElevation;
			float airTemperatureSurface = TemperaturePotential[i] + WorldData.TemperatureLapseRate * surfaceElevation;

			WaterTemperatureSurface[i] = Mathf.Max(WorldData.FreezingTemperature, TemperaturePotential[i] + 2);
			WaterTemperatureBottom[i] = waterDepth > 0 ? GetWaterTemperatureBottom(WaterTemperatureSurface[i], waterDepth, WaterTemperatureDepthFalloff) : WaterTemperatureSurface[i];

			float waterCoverage = Mathf.Clamp01(waterDepth / FullWaterCoverage);

			// TODO: ground water energy should probably be tracked independently
			float terrainTemperature;
			if (waterCoverage >= 1)
			{
				terrainTemperature = WaterTemperatureBottom[i];
			}
			else
			{
				terrainTemperature = airTemperatureSurface;
			}


			CloudVelocity[i] = float3.zero;
			TerrainTemperature[i] = terrainTemperature;

		}
	}


#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private struct WorldGenAirLayerJob : IJobParallelFor {

		public NativeArray<float> AirVapor;
		public NativeArray<float> AirTemperaturePotential;
		public NativeArray<float> LayerElevation;
		public NativeArray<float> LayerHeight;

		[ReadOnly] public NativeArray<float> LayerElevationBelow;
		[ReadOnly] public NativeArray<float> LayerHeightBelow;
		[ReadOnly] public NativeArray<float> TemperaturePotential;
		[ReadOnly] public NativeArray<float> RelativeHumidity;
		[ReadOnly] public float InverseDewPointTemperatureRange;
		[ReadOnly] public float WaterVaporMassToAirMassAtDewPoint;
		[ReadOnly] public float DewPointZero;
		[ReadOnly] public float Gravity;
		[ReadOnly] public float SetLayerHeight;
		[ReadOnly] public float LayerCeiling;
		[ReadOnly] public float LayerCeilingCount;
		public void Execute(int i)
		{
			float layerElevation = LayerElevationBelow[i] + LayerHeightBelow[i];
			float layerHeight;
			if (SetLayerHeight > 0)
			{
				layerHeight = SetLayerHeight;
			} else
			{
				layerHeight = (LayerCeiling - layerElevation) / LayerCeilingCount;
			}
			float layerMidHeight = layerElevation + layerHeight / 2;
			float airTemperatureAbsolute = TemperaturePotential[i] + WorldData.TemperatureLapseRate * layerMidHeight;
			float airMass = Atmosphere.GetAirMass(layerElevation, layerHeight, TemperaturePotential[i], Gravity);
			float vaporMass = GetWaterVaporMass(airTemperatureAbsolute, RelativeHumidity[i], airMass, InverseDewPointTemperatureRange, WaterVaporMassToAirMassAtDewPoint, DewPointZero);
			AirVapor[i] = vaporMass;
			AirTemperaturePotential[i] = TemperaturePotential[i];
			LayerElevation[i] = layerElevation;
			LayerHeight[i] = layerHeight;
		}
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private struct WorldGenWaterLayerJob : IJobParallelFor {

		public NativeArray<float> WaterMass;
		public NativeArray<float> SaltMass;
		public NativeArray<float> WaterTemperature;
		public NativeArray<float> ElevationTop;

		[ReadOnly] public NativeArray<float> WaterTemperatureSurface;
		[ReadOnly] public NativeArray<float> WaterTemperatureBottom;
		[ReadOnly] public NativeArray<float2> coord;
		[ReadOnly] public NativeArray<CellTerrain> terrain;
		[ReadOnly] public float MinSalinity;
		[ReadOnly] public float MaxSalinity;
		[ReadOnly] public float WaterDensityPerDegree;
		[ReadOnly] public float WaterDensityPerSalinity;
		[ReadOnly] public float LayerDepthMax;
		[ReadOnly] public float LayerCount;


		public void Execute(int i)
		{
			float terrainElevation = terrain[i].Elevation;
			float layerDepth = math.clamp((ElevationTop[i] - terrainElevation) / LayerCount, 0, LayerDepthMax);

			if (layerDepth == 0)
			{
			} else
			{
				float salinity = (1.0f - Math.Abs(coord[i].y)) * (MaxSalinity - MinSalinity) + MinSalinity;
				//float deepSalinity = deepDepth == 0 ? 0 : Math.Abs(coord.y) * (worldGenData.MaxSalinity - worldGenData.MinSalinity) + worldGenData.MinSalinity;
				float waterAndSaltMass = GetWaterMass(layerDepth, WaterTemperature[i], salinity, WaterDensityPerDegree, WaterDensityPerSalinity);
				float waterMass = waterAndSaltMass * (1.0f - salinity);
				float saltMass = waterAndSaltMass * salinity;
				WaterTemperature[i] = math.pow(WaterTemperatureSurface[i] - WorldData.FreezingTemperature, 1.0f / (1.0f - ElevationTop[i] / 500.0f)) + WorldData.FreezingTemperature;
				WaterMass[i] = waterMass;
				SaltMass[i] = saltMass;
				ElevationTop[i] -= layerDepth;
			}
		}
	}


#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	public static void Generate(int seed, WorldGenData worldGenData, Icosphere icosphere, ref WorldData worldData, ref StaticState staticState, ref SimState state, ref DependentState dependent)
	{

		float inversePI = 1.0f / math.PI;
		_noise = new FastNoise(seed);
		_noise.SetFrequency(10);
		_random = new System.Random(seed);

		InitSync(worldGenData, icosphere, ref worldData, ref staticState, ref state, ref dependent);

		JobHelper worldGenJobHelper = new JobHelper(staticState.Count);
#if ASYNC_WORLDGEN
		worldGenJobHelper.Async = true;
#else
		worldGenJobHelper.Async = false;
#endif

		var airMassTotal = new NativeArray<float>(staticState.StratosphereMass, Allocator.TempJob);
		var waterMassTotal = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var waterDepthTotal = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var temperaturePotential = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var WaterTemperatureBottom = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var WaterTemperatureTop = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var RelativeHumidity = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var WaterLayerElevation = new NativeArray<float>(staticState.Count, Allocator.TempJob);

		var worldGenInitJob = new WorldGenInitJob()
		{
			Terrain = state.Terrain,
			CloudMass = state.CloudMass,
			potentialTemperature = temperaturePotential,
			relativeHumidity = RelativeHumidity,
			LayerElevationBase = dependent.LayerElevation[0],

			noise = _noise,
			SphericalPosition = staticState.SphericalPosition,
			Coordinate = staticState.Coordinate,
			FullCoverageVegetation =worldData.FullCoverageVegetation,
			MinTemperatureCanopy = worldData.MinTemperatureVegetation,
			MaxTemperatureCanopy = worldData.MaxTemperatureVegetation,
			MinElevation = worldGenData.MinElevation,
			MaxElevation = worldGenData.MaxElevation,
			MaxRoughness =worldGenData.MaxRoughness,
			MinTemperature = worldGenData.MinTemperature,
			MaxTemperature = worldGenData.MaxTemperature,
		};
		worldGenInitJob.Run(staticState.Count);

		var worldGenJobHandle = worldGenJobHelper.Run(new WorldGenJob()
		{
			terrain = state.Terrain,
			TerrainTemperature = state.TerrainTemperature,
			CloudDropletMass = state.CloudDropletMass,
			CloudVelocity = state.CloudVelocity,
			IceMass = state.IceMass,
			IceTemperature = state.IceTemperature,
			TemperaturePotential = temperaturePotential,
			WaterTemperatureBottom = WaterTemperatureBottom,
			WaterTemperatureSurface = WaterTemperatureTop,

			FullWaterCoverage = worldData.FullCoverageWater,
			inversePI = inversePI,
			rainDropMinSize = worldData.rainDropMinSize,
			TropopauseElevation = worldGenData.TropopauseElevation,
			WaterTemperatureDepthFalloff = worldGenData.WaterTemperatureDepthFalloff,
			BoundaryZoneElevation = worldGenData.BoundaryZoneElevation,
		});

		var worldGenWaterLayerJobHandle = worldGenJobHandle;
		for (int i = worldData.WaterLayers - 2; i >= 1; i--)
		{
			float layerDepthMax;
			float layerCount;
			if (i== worldData.WaterLayers - 2)
			{
				layerDepthMax = worldGenData.SurfaceWaterDepth;
				layerCount = 1;
			}
			else if (i == worldData.WaterLayers - 3)
			{
				layerDepthMax = worldGenData.ThermoclineDepth;
				layerCount = 1;
			}
			else
			{
				layerDepthMax = float.MaxValue;
				layerCount = i;
			}
		worldGenWaterLayerJobHandle = worldGenJobHelper.Run(new WorldGenWaterLayerJob()
			{
				WaterTemperature = state.WaterTemperature[i],
				SaltMass = state.SaltMass[i],
				WaterMass = state.WaterMass[i],

				ElevationTop= WaterLayerElevation,
				LayerCount = layerCount,
				LayerDepthMax = layerDepthMax,
				coord = staticState.Coordinate,
				MaxSalinity = worldGenData.MaxSalinity,
				MinSalinity = worldGenData.MinSalinity,
				terrain = state.Terrain,
				WaterDensityPerDegree = worldData.WaterDensityPerDegree,
				WaterDensityPerSalinity = worldData.WaterDensityPerSalinity,
				WaterTemperatureBottom = WaterTemperatureBottom,
				WaterTemperatureSurface = WaterTemperatureTop,
			}, worldGenWaterLayerJobHandle);
		}

		var worldGenAirLayerJobHandle = worldGenJobHandle;
		for (int i = 1; i < worldData.AirLayers-1; i++)
		{
			float setLayerHeight = 0;
			float layerCeiling = 0;
			int layerCeilingCount = 0;
			if (i == 1)
			{
				setLayerHeight = worldGenData.BoundaryZoneElevation;
			} else
			{
				layerCeiling = worldGenData.TropopauseElevation;
				layerCeilingCount = worldData.AirLayers - 1 - i;
			}

			worldGenAirLayerJobHandle = worldGenJobHelper.Run(new WorldGenAirLayerJob()
			{
				AirTemperaturePotential = state.AirTemperaturePotential[i],
				AirVapor = state.AirVapor[i],
				LayerElevation = dependent.LayerElevation[i],
				LayerHeight = dependent.LayerHeight[i],

				LayerElevationBelow = dependent.LayerElevation[i - 1],
				LayerHeightBelow = dependent.LayerHeight[i - 1],
				SetLayerHeight = setLayerHeight,
				LayerCeiling = layerCeiling,
				LayerCeilingCount = layerCeilingCount,
				TemperaturePotential = temperaturePotential,
				RelativeHumidity = RelativeHumidity,
				DewPointZero = worldData.DewPointZero,
				Gravity = state.PlanetState.Gravity,
				InverseDewPointTemperatureRange = worldData.inverseDewPointTemperatureRange,
				WaterVaporMassToAirMassAtDewPoint = worldData.WaterVaporMassToAirMassAtDewPoint,
			}, worldGenAirLayerJobHandle);
		}


		///////////////////////////////////
		// Update dependent variables

		NativeList<JobHandle> updateDependenciesJobHandles = new NativeList<JobHandle>(Allocator.TempJob);
		JobHandle updateDependentWaterLayerJobHandle = default(JobHandle);
		for (int j = worldData.WaterLayers - 2; j >= 1; j--)
		{
			updateDependentWaterLayerJobHandle = worldGenJobHelper.Run( new UpdateDependentWaterLayerJob()
			{
				Salinity = dependent.WaterSalinity[j],
				WaterCoverage = dependent.WaterCoverage[j],
				PotentialEnergy = dependent.WaterPotentialEnergy[j],
				Density = dependent.WaterDensity[j],
				Pressure = dependent.WaterPressure[j],
				LayerDepth = dependent.WaterLayerDepth[j],
				LayerHeight = dependent.WaterLayerHeight[j],
				WaterDepthTotal = waterDepthTotal,
				WaterMassTotal = waterMassTotal,
								
				SaltMass = state.SaltMass[j],
				WaterMass = state.WaterMass[j],
				Temperature = state.WaterTemperature[j],
				Terrain = state.Terrain,
				UpLayerDepth = dependent.WaterLayerDepth[j+1],
				UpLayerHeight = dependent.WaterLayerHeight[j+1],
				Gravity = state.PlanetState.Gravity,
				WaterDensityPerDegree = worldData.WaterDensityPerDegree,
				WaterDensityPerSalinity = worldData.WaterDensityPerSalinity
			}, worldGenWaterLayerJobHandle);
			updateDependenciesJobHandles.Add(updateDependentWaterLayerJobHandle);
		}

		var updateDependentStateJobHandle = worldGenJobHelper.Run(new UpdateDependentStateJob()
		{
			CloudCoverage = dependent.CloudCoverage,
			IceCoverage = dependent.IceCoverage,
			SurfaceElevation = dependent.SurfaceElevation,
			VegetationCoverage = dependent.VegetationCoverage,
			SurfaceAirTemperature = dependent.SurfaceAirTemperatureAbsolute,
			IceEnergy = dependent.IceEnergy,
			WaterDepth = dependent.WaterDepth,


			WaterDepthTotal = waterDepthTotal,
			CloudMass = state.CloudMass,
			IceMass = state.IceMass,
			IceTemperature = state.IceTemperature,
			Terrain = state.Terrain,
			worldData = worldData,
		}, updateDependentWaterLayerJobHandle);
		updateDependenciesJobHandles.Add(updateDependentStateJobHandle);


		JobHandle lowerAirHandle = default(JobHandle);
		JobHandle updateDependentAirLayerJobHandle = worldGenAirLayerJobHandle;
		// Top to bottom so we can add up air mass for pressure calculation
		for (int j = worldData.AirLayers-2; j >= 1; j--)
		{
			updateDependentAirLayerJobHandle = worldGenJobHelper.Run(new UpdateDependentAirLayerJob()
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
				AirVaporCloud = dependent.AirVaporCloud,
				CloudElevation = dependent.CloudElevation,
				DewPoint = dependent.DewPoint,
				AirMassTotal = airMassTotal,

				CloudDropletMass = state.CloudDropletMass,
				CloudMass = state.CloudMass,
				SurfaceElevation = dependent.SurfaceElevation,
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				VaporMass = state.AirVapor[j],
				IceMass = state.IceMass,
				Gravity = state.PlanetState.Gravity,
				DewPointZero = worldData.DewPointZero,
				InverseDewPointTemperatureRange = worldData.inverseDewPointTemperatureRange,
				WaterVaporMassToAirMassAtDewPoint = worldData.WaterVaporMassToAirMassAtDewPoint,
				LayerIndex = j,
			}, updateDependentAirLayerJobHandle);
			updateDependenciesJobHandles.Add(updateDependentAirLayerJobHandle);
			if (j == 1)
			{
				lowerAirHandle = updateDependentAirLayerJobHandle;
			}
		}

		var updateSurfaceDependentStateJobHandle = worldGenJobHelper.Run(new UpdateSurfaceDependentStateJob()
		{
			SurfaceAirTemperatureAbsolute = dependent.SurfaceAirTemperatureAbsolute,

			AirTemperaturePotential = state.AirTemperaturePotential[1],
			SurfaceElevation = dependent.LayerElevation[1]
		}, lowerAirHandle);
		updateDependenciesJobHandles.Add(updateSurfaceDependentStateJobHandle);

		JobHandle.CompleteAll(updateDependenciesJobHandles);
		updateDependenciesJobHandles.Dispose();
		airMassTotal.Dispose();
		waterMassTotal.Dispose();
		waterDepthTotal.Dispose();
		temperaturePotential.Dispose();
		WaterTemperatureBottom.Dispose();
		WaterTemperatureTop.Dispose();
		WaterLayerElevation.Dispose();
		RelativeHumidity.Dispose();
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	static public float GetWaterMass(float depth, float temperature, float salinityPSU, float WaterDensityPerDegree, float WaterDensityPerSalinity)
	{
		float density = WorldData.DensityWater + WaterDensityPerDegree * (temperature - WorldData.FreezingTemperature) + WaterDensityPerSalinity * salinityPSU;
		return depth * density;
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	static public float GetWaterVaporMass(float temperature, float relativeHumidity, float airMass, float inverseDewPointTemperatureRange, float WaterVaporMassToAirMassAtDewPoint, float DewPointZero)
	{
		float maxWaterVaporPerKilogramAtmosphere = WaterVaporMassToAirMassAtDewPoint * Utils.Sqr(math.max(0, (temperature - DewPointZero) * inverseDewPointTemperatureRange));
		float maxHumidity = maxWaterVaporPerKilogramAtmosphere * airMass;
		if (maxHumidity <= 0)
		{
			return 0;
		}
		float humidity = relativeHumidity * maxHumidity;
		return humidity;
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private static float GetWaterTemperatureBottom(float waterTemperatureSurface, float waterDepth, float WaterTemperatureDepthFalloff)
	{
		return math.pow(waterTemperatureSurface - WorldData.FreezingTemperature, 1.0f / (1.0f + waterDepth / WaterTemperatureDepthFalloff)) + WorldData.FreezingTemperature;
	}



	#region PerlinUtils
#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private static float GetPerlinMinMax(float x, float y, float frequency, float hash, float min, float max)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency) + 1.0f) * (max - min) / 2 + min;
	}
#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private static float GetPerlinMinMax(float x, float y, float z, float frequency, float hash, float min, float max)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency, z * frequency) + 1.0f) * (max - min) / 2 + min;
	}
#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private static float GetPerlinNormalized(float x, float y, float frequency, float hash)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency) + 1.0f) / 2;
	}
#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private static float GetPerlinNormalized(float x, float y, float z, float frequency, float hash)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency, z * frequency) + 1.0f) / 2;
	}
#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private static float GetPerlin(float x, float y, float frequency, float hash)
	{
		return _noise.GetPerlin(x * frequency + hash, y * frequency);
	}
#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private static float GetPerlin(float x, float y, float z, float frequency, float hash)
	{
		return _noise.GetPerlin(x * frequency + hash, y * frequency, z * frequency);
	}
#endregion


}
