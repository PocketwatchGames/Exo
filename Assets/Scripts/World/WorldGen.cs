﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;
using UnityEngine;

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
	public static void Generate(int seed, WorldGenData worldGenData, Icosphere icosphere, ref WorldData worldData, ref StaticState staticState, ref SimState state)
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

		for (int i = 0; i < state.CellDependents.Length; i++)
		{
			CellDependent dependent;
			CellState cell;
			CellTerrain terrain;
			var coord = staticState.Coordinate[i] * 2 * inversePI;
			var pos = staticState.SphericalPosition[i];
			var coordNormalized = (coord + 1) / 2;
			terrain.Roughness = math.max(1, math.pow(GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 660), 3) * worldGenData.MaxRoughness);
			float elevation =
				0.6f * GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 0, worldGenData.MinElevation, worldGenData.MaxElevation) +
				0.3f * GetPerlinMinMax(pos.x, pos.y, pos.z, 0.5f, 0, worldGenData.MinElevation, worldGenData.MaxElevation) +
				0.1f * GetPerlinMinMax(pos.x, pos.y, pos.z, 2.0f, 0, worldGenData.MinElevation, worldGenData.MaxElevation);
			terrain.Elevation = elevation;
			dependent.WaterDepth = math.max(0, -elevation) + math.max(0, GetPerlinMinMax(pos.x, pos.y, pos.z, 1.0f, 25430, -100, 100));
			//cell.WaterDepth = GetPerlinNormalized(coordNormalized.x, coordNormalized.y, 2.5f, 110) * 2000;

			float surfaceElevation = elevation + dependent.WaterDepth;

			float soil =
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.2f, 6630) +
				0.3f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 435) +
				0.2f * GetPerlinNormalized(pos.x, pos.y, pos.z, 1f, 8740);
			terrain.SoilFertility = soil;

			terrain.Vegetation = elevation < 0 ? 0 : math.pow(math.sqrt(soil) * math.sqrt(GetPerlinNormalized(pos.x, pos.y, pos.z, 1f, 6630423)), 3);

			dependent.CloudCoverage =
				math.pow(
					0.75f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 30) +
					0.25f * GetPerlinNormalized(pos.x, pos.y, pos.z, 1.0f, 30),
					1);

			dependent.RelativeHumidity =
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 40) +
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 40);

			dependent.WindSurface = float2.zero;
			dependent.WindTropopause = float2.zero;
			dependent.WindVertical = 0;
			dependent.CurrentDeep = float2.zero;
			dependent.CurrentSurface = float2.zero;
			dependent.CurrentVertical = 0;


			float troposphereColumnHeight = worldData.TropopauseElevation - surfaceElevation;

			float regionalTemperatureVariation =
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 15460, -5, 5) +
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.2f, 431, -10, 10) +
				GetPerlinMinMax(pos.x, pos.y, pos.z, 0.5f, 6952, -10, 10);
			cell.AirTemperature =
				regionalTemperatureVariation + GetPerlinMinMax(pos.x, pos.y, pos.z, 0.15f, 80, -5, 5) +
				(1.0f - coord.y * coord.y) * (worldGenData.MaxTemperature - worldGenData.MinTemperature) + worldGenData.MinTemperature + WorldData.TemperatureLapseRate * surfaceElevation;

			dependent.AirPressure = WorldData.StaticPressure - regionalTemperatureVariation * 1000;
			cell.CloudMass = Mathf.Pow(GetPerlinMinMax(pos.x, pos.y, pos.z, 1.0f, 2000, 0, 1), 0.85f) * Mathf.Pow(dependent.RelativeHumidity, 0.5f);
			cell.CloudDropletMass = GetPerlinMinMax(pos.x, pos.y, pos.z, 1.0f, 2001, 0, 1) * dependent.RelativeHumidity * Mathf.Pow(cell.CloudMass, 2) * worldData.rainDropMaxSize;
			float totalAirMass = Atmosphere.GetAirMass(ref worldData, dependent.AirPressure, surfaceElevation, cell.AirTemperature, WorldData.MolarMassAir, state.PlanetState.Gravity) - state.PlanetState.StratosphereMass;
			//float totalAirMassLower = Atmosphere.GetAirMass(world, cell.LowerAirPressure, surfaceElevation, cell.LowerAirTemperature, world.Data.MolarMassAir * (1.0f - relativeHumidity) + world.Data.MolarMassWater * relativeHumidity) - state.StratosphereMass - totalAirMassUpper;
			cell.AirWaterMass = Atmosphere.GetAbsoluteHumidity(ref worldData, cell.AirTemperature, dependent.RelativeHumidity, totalAirMass, inverseDewPointTemperatureRange);

			//cell.UpperAirMass = totalAirMassUpper - state.CloudMass;
			//cell.UpperAirPressure = Atmosphere.GetAirPressure(world, state.UpperAirMass + state.StratosphereMass + state.CloudMass, surfaceElevation + world.Data.BoundaryZoneElevation, state.UpperAirTemperature, Atmosphere.GetMolarMassAir(world, state.UpperAirMass + state.StratosphereMass, state.CloudMass));
			dependent.AirMass = totalAirMass - cell.AirWaterMass - cell.CloudMass;
			dependent.AirPressure = Atmosphere.GetAirPressure(
				ref worldData,
				dependent.AirMass + state.PlanetState.StratosphereMass + cell.AirWaterMass + cell.CloudMass,
				surfaceElevation,
				cell.AirTemperature,
				Atmosphere.GetMolarMassAir(dependent.AirMass + state.PlanetState.StratosphereMass, cell.AirWaterMass + cell.CloudMass),
				state.PlanetState.Gravity);

			cell.GroundWaterDepth = GetPerlinMinMax(pos.x, pos.y, pos.z, 1f, 60423, worldGenData.GroundWaterDepthMin, worldGenData.GroundWaterDepthMax);
			float maxGroundWater = cell.GroundWaterDepth * terrain.SoilFertility * WorldData.MassWater * worldData.MaxSoilPorousness;
			if (dependent.WaterDepth > 0)
			{
				cell.GroundWater = maxGroundWater;
			}
			else
			{
				cell.GroundWater = maxGroundWater * 0.2f;
			}

			//			cell.UpperAirEnergy = Atmosphere.GetAirEnergy(world, state.UpperAirTemperature, state.UpperAirMass, state.CloudMass, worldData.SpecificHeatWater);
			dependent.AirEnergy = Atmosphere.GetAirEnergy(cell.AirTemperature, dependent.AirMass, cell.CloudMass, cell.AirWaterMass);

			//float shallowDepth = Mathf.Min(data.DeepOceanDepth, depth);
			//float deepDepth = Mathf.Max(0, depth - data.DeepOceanDepth);
			float salinity = dependent.WaterDepth == 0 ? 0 : (1.0f - Math.Abs(coord.y)) * (worldGenData.MaxSalinity - worldGenData.MinSalinity) + worldGenData.MinSalinity;
			//float deepSalinity = deepDepth == 0 ? 0 : Math.Abs(coord.y) * (worldGenData.MaxSalinity - worldGenData.MinSalinity) + worldGenData.MinSalinity;

			float minOceanTemperature = WorldData.FreezingTemperature + 0.1f;
			//float oceanDepthMinTemperature = 2000;
			dependent.WaterTemperature = Mathf.Max(WorldData.FreezingTemperature, cell.AirTemperature + 2);
			float waterAndSaltMass = GetWaterMass(worldData, dependent.WaterDepth, dependent.WaterTemperature, salinity);
			cell.WaterMass = waterAndSaltMass * (1.0f - salinity);
			cell.SaltMass = waterAndSaltMass * salinity;
			cell.IceMass = math.min(cell.WaterMass, worldData.FullIceCoverage * WorldData.MassIce * Mathf.Clamp01(-(cell.AirTemperature - WorldData.FreezingTemperature) / 10));
			cell.WaterMass -= cell.IceMass;
			cell.WaterEnergy = Atmosphere.GetWaterEnergy(dependent.WaterTemperature, cell.WaterMass, cell.SaltMass);
			dependent.WaterDensity = Atmosphere.GetWaterDensity(ref worldData, cell.WaterEnergy, cell.SaltMass, cell.WaterMass);
			dependent.WaterAndIceDepth = dependent.WaterDepth + cell.IceMass / WorldData.MassIce;
			dependent.CloudElevation = Atmosphere.GetCloudElevation(ref worldData, cell.AirTemperature, Atmosphere.GetDewPoint(ref worldData, cell.AirTemperature, dependent.RelativeHumidity), surfaceElevation);

			////	float deepOceanTemperature = (state.ShallowWaterTemperature - minOceanTemperature) * (1.0f - Mathf.Pow(Mathf.Clamp01(deepDepth / oceanDepthMinTemperature), 2f)) + minOceanTemperature;
			//float deepOceanTemperature = minOceanTemperature;
			//float deepOceanMass = GetWaterMass(world, deepDepth, deepOceanTemperature, deepSalinity);
			//cell.DeepWaterMass = deepOceanMass * (1.0f - deepSalinity);
			//cell.DeepSaltMass = deepOceanMass * deepSalinity;
			//cell.DeepWaterEnergy = Atmosphere.GetWaterEnergy(world, deepOceanTemperature, state.DeepWaterMass, state.DeepSaltMass);
			//cell.DeepWaterDensity = Atmosphere.GetWaterDensity(world, state.DeepWaterEnergy, state.DeepSaltMass, state.DeepWaterMass);


			float waterCoverage = Mathf.Clamp01(dependent.WaterDepth / worldData.FullWaterCoverage);
			float iceCoverage = Mathf.Clamp01(cell.IceMass / (WorldData.MassIce * worldData.FullIceCoverage));

			if (dependent.WaterDepth == 0)
			{
				terrain.Vegetation = terrain.SoilFertility * (cell.GroundWater) * (1.0f - waterCoverage) * (1.0f - iceCoverage) * Mathf.Clamp01((cell.AirTemperature - worldData.MinTemperatureCanopy) / (worldData.MaxTemperatureCanopy - worldData.MinTemperatureCanopy));
			}

			// TODO: ground water energy should probably be tracked independently
			float groundRadiationRate = Atmosphere.GetRadiationRate(cell.AirTemperature, WorldData.EmissivityDirt);
			float groundWaterEnergy = cell.GroundWater * WorldData.SpecificHeatWater * worldData.maxGroundWaterTemperature;
			float landMass = (WorldData.MassSand - WorldData.MassSoil) * terrain.SoilFertility + WorldData.MassSoil;
			float heatingDepth = terrain.SoilFertility * worldData.SoilHeatDepth;
			float groundTemperature;
			if (waterCoverage >= 1)
			{
				groundTemperature = dependent.WaterTemperature;
			}
			else
			{
				groundTemperature = cell.AirTemperature;
			}
			float soilEnergy = groundTemperature * WorldData.SpecificHeatSoil * heatingDepth * landMass;
			cell.GroundEnergy = groundWaterEnergy + soilEnergy;


			state.CellDependents[i] = dependent;
			state.CellStates[i] = cell;
			state.CellTerrains[i] = terrain;
		}
	}

	static public float GetWaterMass(WorldData worldData, float depth, float temperature, float salinityPSU)
	{
		float density = worldData.waterDensity + worldData.OceanDensityPerDegree * (temperature - WorldData.FreezingTemperature) + worldData.OceanDensityPerSalinity * salinityPSU;
		return depth * density;
	}

}
