#define ASYNC_WORLDGEN

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

	private static void InitSync(WorldGenData worldGenData, Icosphere icosphere, ref WorldData worldData, ref StaticState staticState, ref SimState state, ref TempState tempState)
	{
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
		state.PlanetState.Oxygen = worldGenData.Oxygen;
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private struct WorldGenInitJob {

		public NativeArray<float> relativeHumidity;
		public NativeArray<float> potentialTemperature;
		public NativeArray<float> CloudMass;
		public NativeArray<float> GroundWater;
		public NativeArray<float> Roughness;
		public NativeArray<float> SoilCarbon;
		public NativeArray<float> SurfaceElevation;
		public NativeArray<float> Elevation;
		public NativeArray<float> Flora;
		public NativeArray<float> FloraWater;
		public NativeArray<float> FloraTemperature;
		public NativeArray<float> FloraGlucose;
		public NativeArray<float> MagmaMass;
		public NativeArray<float> CrustDepth;
		public NativeArray<float> LavaMass;
		public NativeArray<float> LavaTemperature;

		[ReadOnly] public NativeArray<float2> Coordinate;
		[ReadOnly] public NativeArray<float3> SphericalPosition;
		[ReadOnly] public FastNoise noise;
		[ReadOnly] public float MaxRoughness;
		[ReadOnly] public float MinElevation;
		[ReadOnly] public float MaxElevation;
		[ReadOnly] public float MaxTemperature;
		[ReadOnly] public float MinTemperature;
		[ReadOnly] public float FullCoverageFlora;
		[ReadOnly] public float MinTemperatureFlora;
		[ReadOnly] public float MaxTemperatureFlora;
		[ReadOnly] public float MaxGroundWater;
		[ReadOnly] public float SoilCarbonMax;
		[ReadOnly] public float MagmaMin;
		[ReadOnly] public float MagmaMax;
		[ReadOnly] public float CrustMin;
		[ReadOnly] public float CrustMax;
		[ReadOnly] public float LavaEruptionTemperature;

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
			relativeHumidity[i] = 0.5f +
				0.25f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 40) +
				0.25f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 40);
			if (elevation > 0)
			{
				relativeHumidity[i] *= math.saturate(1.0f - elevation / MaxElevation);
			}
			float surfaceElevation = math.max(0, elevation);

			float regionalTemperatureVariation =
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 15460, -5, 5) +
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.2f, 431, -10, 10) +
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.5f, 6952, -10, 10);
			potentialTemperature[i] =
				regionalTemperatureVariation + GetPerlinMinMax(pos.x, pos.y, pos.z, 0.15f, 80, -5, 5) +
				(1.0f - coord.y * coord.y) * (MaxTemperature - MinTemperature) + MinTemperature;

			float soilFertilityPercent =
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.2f, 6630) +
				0.3f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 435) +
				0.2f * GetPerlinNormalized(pos.x, pos.y, pos.z, 1f, 8740);
			float soilFertility = soilFertilityPercent * SoilCarbonMax;

			float airTemperatureSurface = potentialTemperature[i] + WorldData.TemperatureLapseRate * surfaceElevation;
			float flora = 0;
			float floraWater = 0;
			float floraGlucose = 0;
			if (elevation > 0)
			{
				flora =
					FullCoverageFlora * math.pow(
					soilFertilityPercent
					* GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 410)
					* GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 61)
					* math.sin(math.PI * math.saturate((airTemperatureSurface - MinTemperatureFlora) / (MaxTemperatureFlora - MinTemperatureFlora))),
					2.5f);
				floraWater =
					flora * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 41630);
				floraGlucose = flora;
			}

			float cloudMass = Mathf.Pow(GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 2000, 0, 1), 1.0f) * Mathf.Pow(relativeHumidity[i], 2.0f);

			float groundWater = MaxGroundWater;
			if (elevation > 0)
			{
				groundWater *= 
					0.25f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 1560) +
					0.25f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.2f, 4381) +
					0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 692);
			}
			GroundWater[i] = groundWater;

			CrustDepth[i] = CrustMin + (CrustMax - CrustMin) * ((elevation - MinElevation) / (MaxElevation - MinElevation)) *
				(0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 4570) +
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 1430));

			MagmaMass[i] = math.max(0, (elevation - CrustDepth[i] + 2000)) * WorldData.MassLava;

			//LavaMass[i] = 10000 * math.max(0, (GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 1630) - 0.75f));
			//LavaMass[i] = 10000;
			LavaTemperature[i] = LavaEruptionTemperature;

			SurfaceElevation[i] = surfaceElevation;
			Elevation[i] = elevation;
			Flora[i] = flora;
			FloraTemperature[i] = airTemperatureSurface;
			FloraWater[i] = floraWater;
			FloraGlucose[i] = floraGlucose;
			Roughness[i] = roughness;
			SoilCarbon[i] = soilFertility;

			CloudMass[i] = cloudMass;
		}
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private struct WorldGenJob : IJobParallelFor {

		public NativeArray<float> IceTemperature;
		public NativeArray<float> CloudDropletMass;
		public NativeArray<float> CloudTemperature;
		public NativeArray<float> IceMass;
		public NativeArray<float> TerrainTemperature;
		public NativeArray<float> WaterTemperatureSurface;
		public NativeArray<float> WaterTemperatureBottom;
		public NativeArray<float> GroundWaterTemperature;

		[ReadOnly] public NativeArray<float> TemperaturePotential;
		[ReadOnly] public NativeArray<float> Elevation;
		[ReadOnly] public float inversePI;
		[ReadOnly] public float rainDropMinSize;
		[ReadOnly] public float WaterTemperatureDepthFalloff;
		[ReadOnly] public float FullWaterCoverage;

		public void Execute(int i)
		{
			float elevation = Elevation[i];

			float waterDepth = math.max(0, -elevation);
			float surfaceElevation = elevation + waterDepth;

			float airTemperatureSurface = TemperaturePotential[i] + WorldData.TemperatureLapseRate * surfaceElevation;

			WaterTemperatureSurface[i] = Mathf.Max(WorldData.FreezingTemperature, TemperaturePotential[i] + 0);
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


			TerrainTemperature[i] = terrainTemperature;
			GroundWaterTemperature[i] = terrainTemperature;
		}
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private struct WorldGenAirLayerJob : IJobParallelFor {

		public NativeSlice<float> AirTemperaturePotential;
		public NativeSlice<float> Dust;

		[ReadOnly] public NativeArray<float> TemperaturePotential;

		public void Execute(int i)
		{
			AirTemperaturePotential[i] = TemperaturePotential[i];
			Dust[i] = 0;
		}
	}


#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private struct WorldGenWaterVaporJob : IJobParallelFor {

		public NativeSlice<float> AirVapor;
		public NativeSlice<float> CarbonDioxide;

		[ReadOnly] public NativeSlice<float> AirMass;
		[ReadOnly] public NativeSlice<float> Pressure;
		[ReadOnly] public NativeSlice<float> LayerMiddle;
		[ReadOnly] public NativeSlice<float> TemperaturePotential;
		[ReadOnly] public NativeSlice<float> RelativeHumidity;
		[ReadOnly] public float CarbonDioxidePPM;
		public void Execute(int i)
		{
			float airTemperatureAbsolute = TemperaturePotential[i] + WorldData.TemperatureLapseRate * LayerMiddle[i];
			float airVapor = Atmosphere.GetMaxVaporAtTemperature(AirMass[i], airTemperatureAbsolute, Pressure[i]) * RelativeHumidity[i];
			AirVapor[i] = airVapor;
			CarbonDioxide[i] = CarbonDioxidePPM * AirMass[i];
		}
	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	private struct WorldGenWaterLayerJob : IJobParallelFor {

		public NativeArray<float> WaterMass;
		public NativeArray<float> CarbonMass;
		public NativeArray<float> SaltMass;
		public NativeArray<float> PlanktonMass;
		public NativeArray<float> PlanktonGlucose;
		public NativeArray<float> WaterTemperature;
		public NativeArray<float> ElevationTop;

		[ReadOnly] public NativeArray<float> WaterTemperatureSurface;
		[ReadOnly] public NativeArray<float> WaterTemperatureBottom;
		[ReadOnly] public NativeArray<float2> coord;
		[ReadOnly] public NativeArray<float> Elevation;
		[ReadOnly] public float CarbonPercent;
		[ReadOnly] public float MinSalinity;
		[ReadOnly] public float MaxSalinity;
		[ReadOnly] public float WaterDensityPerDegree;
		[ReadOnly] public float WaterDensityPerSalinity;
		[ReadOnly] public float LayerDepthMax;
		[ReadOnly] public float LayerCount;
		[ReadOnly] public float PlanktonInLayer;


		public void Execute(int i)
		{
			float layerDepth = math.clamp((ElevationTop[i] - Elevation[i]) / LayerCount, 0, LayerDepthMax);

			if (layerDepth == 0)
			{
			} else
			{
				WaterTemperature[i] = math.pow(WaterTemperatureSurface[i] - WorldData.FreezingTemperature, 1.0f / (1.0f - ElevationTop[i] / 500.0f)) + WorldData.FreezingTemperature;

				float salinity = Math.Abs(coord[i].y) * (MaxSalinity - MinSalinity) + MinSalinity;
				float waterDensity = Atmosphere.GetWaterDensity(salinity, WaterTemperature[i]);
				float waterAndSaltMass = layerDepth * waterDensity;
				float waterMass = waterAndSaltMass * (1.0f - salinity);
				float saltMass = waterAndSaltMass * salinity;
				float waterCarbonMass = CarbonPercent * waterAndSaltMass / (1.0f - CarbonPercent);

				float density = Atmosphere.GetWaterDensity(salinity, WaterTemperature[i]);
				float actualDepth = (waterMass + saltMass) / density;
				
				WaterMass[i] = waterMass;
				CarbonMass[i] = waterCarbonMass;
				SaltMass[i] = saltMass;
				ElevationTop[i] -= layerDepth;
				PlanktonMass[i] = PlanktonInLayer * layerDepth / LayerDepthMax;
				PlanktonGlucose[i] = PlanktonMass[i] / 2;

			}
		}
	}


#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	public static void Generate(int seed, WorldGenData worldGenData, Icosphere icosphere, ref WorldData worldData, ref StaticState staticState, ref SimState state, ref TempState tempState)
	{

		float inversePI = 1.0f / math.PI;
		_noise = new FastNoise(seed);
		_noise.SetFrequency(10);
		_random = new System.Random(seed);

		InitSync(worldGenData, icosphere, ref worldData, ref staticState, ref state, ref tempState);

		JobHelper worldGenJobHelper = new JobHelper(staticState.Count);

		var waterMassTotal = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var waterDepthTotal = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var temperaturePotential = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var WaterTemperatureBottom = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var WaterTemperatureTop = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var RelativeHumidity = new NativeArray<float>(staticState.Count, Allocator.TempJob);
		var WaterLayerElevation = new NativeArray<float>(staticState.Count, Allocator.TempJob);

		var worldGenInitJob = new WorldGenInitJob()
		{
			Roughness = state.Roughness,
			SoilCarbon = state.GroundCarbon,
			CloudMass = state.CloudMass,
			potentialTemperature = temperaturePotential,
			relativeHumidity = RelativeHumidity,
			SurfaceElevation = tempState.SurfaceElevation,
			GroundWater = state.GroundWater,
			Elevation = state.Elevation,
			Flora = state.FloraMass,
			FloraWater = state.FloraWater,
			FloraTemperature = state.FloraTemperature,
			FloraGlucose = state.FloraGlucose,
			CrustDepth = state.CrustDepth,
			MagmaMass = state.MagmaMass,
			LavaMass = state.LavaMass,
			LavaTemperature = state.LavaTemperature,

			noise = _noise,
			SphericalPosition = staticState.SphericalPosition,
			Coordinate = staticState.Coordinate,
			FullCoverageFlora =worldData.FullCoverageFlora,
			MinTemperatureFlora = WorldData.FreezingTemperature - 10,
			MaxTemperatureFlora = WorldData.FreezingTemperature + 40,
			MinElevation = worldGenData.MinElevation,
			MaxElevation = worldGenData.MaxElevation,
			MaxRoughness =worldGenData.MaxRoughness,
			MinTemperature = worldGenData.MinTemperature,
			MaxTemperature = worldGenData.MaxTemperature,
			MaxGroundWater = worldData.GroundWaterMax,
			SoilCarbonMax = worldGenData.SoilCarbonMass,
			MagmaMin = worldGenData.MagmaMin,
			MagmaMax = worldGenData.MagmaMax,
			CrustMin = worldGenData.CrustMin,
			CrustMax = worldGenData.CrustMax,
			LavaEruptionTemperature = worldData.MagmaTemperature
		};
		worldGenInitJob.Run(staticState.Count);

		var worldGenJobHandle = worldGenJobHelper.Schedule(new WorldGenJob()
		{
			Elevation = state.Elevation,
			TerrainTemperature = state.GroundTemperature,
			CloudDropletMass = state.CloudDropletMass,
			CloudTemperature = state.CloudTemperature,
			IceMass = state.IceMass,
			IceTemperature = state.IceTemperature,
			WaterTemperatureBottom = WaterTemperatureBottom,
			WaterTemperatureSurface = WaterTemperatureTop,
			GroundWaterTemperature = state.GroundWaterTemperature,

			TemperaturePotential = temperaturePotential,
			FullWaterCoverage = worldData.FullCoverageWater,
			inversePI = inversePI,
			rainDropMinSize = worldData.rainDropMinSize,
			WaterTemperatureDepthFalloff = worldGenData.WaterTemperatureDepthFalloff,
		});

		for (int i = worldData.WaterLayers - 2; i >= 1; i--)
		{
			float layerDepthMax;
			float layerCount;
			float plankton;
			if (i== worldData.WaterLayers - 2)
			{
				layerDepthMax = worldData.SurfaceWaterDepth;
				layerCount = 1;
				plankton = worldGenData.PlanktonMass;
			}
			else if (i == worldData.WaterLayers - 3)
			{
				layerDepthMax = worldData.ThermoclineDepth;
				layerCount = 1;
				plankton = 0;
			}
			else
			{
				layerDepthMax = float.MaxValue;
				layerCount = i;
				plankton = 0;
			}
			worldGenJobHandle = worldGenJobHelper.Schedule(new WorldGenWaterLayerJob()
			{
				WaterTemperature = state.WaterTemperature[i],
				SaltMass = state.SaltMass[i],
				WaterMass = state.WaterMass[i],
				CarbonMass = state.WaterCarbon[i],
				PlanktonMass = state.PlanktonMass[i],
				PlanktonGlucose = state.PlanktonGlucose[i],
				ElevationTop = WaterLayerElevation,

				LayerCount = layerCount,
				LayerDepthMax = layerDepthMax,
				coord = staticState.Coordinate,
				MaxSalinity = worldGenData.MaxSalinity,
				MinSalinity = worldGenData.MinSalinity,
				Elevation = state.Elevation,
				CarbonPercent = worldGenData.WaterCarbonPercent,
				WaterDensityPerDegree = worldData.WaterDensityPerDegree,
				WaterDensityPerSalinity = worldData.WaterDensityPerSalinity,
				WaterTemperatureBottom = WaterTemperatureBottom,
				WaterTemperatureSurface = WaterTemperatureTop,
				PlanktonInLayer = plankton
			}, worldGenJobHandle);
		}

		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			worldGenJobHandle = worldGenJobHelper.Schedule(new WorldGenAirLayerJob()
			{
				AirTemperaturePotential = staticState.GetSliceAirLayer(state.AirTemperaturePotential, i),
				Dust = staticState.GetSliceAirLayer(state.Dust, i),

				TemperaturePotential = temperaturePotential,
			}, worldGenJobHandle);
		}

		worldGenJobHandle.Complete();

		var tempArrays = new List<NativeArray<float>>();
		tempState.Update(ref state, ref worldData, ref staticState, worldGenJobHandle).Complete();

		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			worldGenJobHandle = worldGenJobHelper.Schedule(new WorldGenWaterVaporJob()
			{
				AirVapor = staticState.GetSliceAirLayer(state.AirVapor,i),
				CarbonDioxide = staticState.GetSliceAirLayer(state.AirCarbon,i),

				AirMass = staticState.GetSliceAirLayer(tempState.AirMass,i),
				Pressure = staticState.GetSliceAirLayer(tempState.AirPressure,i),
				LayerMiddle = staticState.GetSliceAirLayer(tempState.AirLayerMiddle,i),
				TemperaturePotential = staticState.GetSliceAirLayer(state.AirTemperaturePotential,i),
				RelativeHumidity = RelativeHumidity,
				CarbonDioxidePPM = worldGenData.AirCarbonPercent,

			}, worldGenJobHandle);
		}


		///////////////////////////////////
		// Update dependent variables

		tempState.Update(ref state, ref worldData, ref staticState, worldGenJobHandle).Complete();

		waterMassTotal.Dispose();
		waterDepthTotal.Dispose();
		temperaturePotential.Dispose();
		WaterTemperatureBottom.Dispose();
		WaterTemperatureTop.Dispose();
		WaterLayerElevation.Dispose();
		RelativeHumidity.Dispose();
		foreach (var a in tempArrays)
		{
			a.Dispose();
		}

		tempState.Clear(staticState.Count, ref worldData, default(JobHandle)).Complete();

	}

#if ASYNC_WORLDGEN
	[BurstCompile]
#endif
	static public float GetWaterVaporMass(float temperature, float relativeHumidity, float airMass, float pressure)
	{
		return math.max(0, relativeHumidity * Atmosphere.GetMaxVaporAtTemperature(airMass, temperature, pressure));
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
