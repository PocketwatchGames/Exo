using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;
using Unity.Burst;
using UnityEngine;

public static class Atmosphere {

	[BurstCompile]
	static public float GetPotentialTemperature(float temperature, float elevation)
	{
		return temperature - WorldData.TemperatureLapseRate * elevation;
	}


	[BurstCompile]
	static public float GetPressureAtElevation(float elevation, float gravity, float referencePressure, float referenceTemperature, float referenceElevation)
	{
		float pressure = referencePressure * math.pow(referenceTemperature / (referenceTemperature + (elevation - referenceElevation) * WorldData.TemperatureLapseRate), gravity * WorldData.PressureExponent * WorldData.MolarMassAir);
		return pressure;
	}


	[BurstCompile]
	static public float GetElevationAtPressure(float pressure, float referenceTemperature, float referencePressure, float referenceElevation, float gravity)
	{
		float elevation = referenceElevation + referenceTemperature / WorldData.TemperatureLapseRate * (math.pow(pressure / referencePressure, -1.0f / (gravity * WorldData.PressureExponent * WorldData.MolarMassAir)) - 1);
		return elevation;
	}

	[BurstCompile]
	static public float GetStandardPressureAtElevation(float elevation, float temperature, float gravity)
	{
		float pressure = Atmosphere.GetPressureAtElevation(elevation, gravity, WorldData.StaticPressure, Atmosphere.GetPotentialTemperature(temperature, elevation), 0);
		return pressure;
	}

	[BurstCompile]
	static public float GetAirMass(float layerElevation, float layerHeight, float temperature, float gravity)
	{
		float layerMiddle = layerElevation + layerHeight / 2;
		float standardPressure = GetStandardPressureAtElevation(layerMiddle, temperature, gravity);
		return standardPressure* layerHeight * WorldData.MolarMassAir / (WorldData.UniversalGasConstant * temperature);
	}

	//[BurstCompile]
	//static public float GetStandardAirMass(float elevation, float layerHeight, float gravity)
	//{
	//	float temperatureLapseA = -WorldData.TemperatureLapseRate * elevation;
	//	float massA = WorldData.StaticPressure / (gravity * math.pow(1.0f - (temperatureLapseA) / (WorldData.StdTemp + temperatureLapseA), gravity * WorldData.PressureExponent * WorldData.MolarMassAir));
	//	float temperatureLapseB = -WorldData.TemperatureLapseRate * (elevation + layerHeight);
	//	float massB = WorldData.StaticPressure / (gravity * math.pow(1.0f - (temperatureLapseB) / (WorldData.StdTemp + temperatureLapseB), gravity * WorldData.PressureExponent * WorldData.MolarMassAir));
	//	//	Debug.Log("A: " + massA + " B: " + massB + " E: " + elevation + " H: " + layerHeight);
	//	return massA - massB;
	//}

	[BurstCompile]
	static public float GetMolarMassMoistAir(float airMass, float vaporMass)
	{
		return (WorldData.GasConstantAir * airMass + WorldData.GasConstantWaterVapor * vaporMass) / (airMass + vaporMass);
	}

	[BurstCompile]
	static public float GetAirDensity(float absolutePressure, float temperature, float airMass, float vaporMass)
	{
		return absolutePressure / (WorldData.UniversalGasConstant * temperature * GetMolarMassMoistAir(airMass, vaporMass));
	}

	[BurstCompile]
	static public float GetWaterDensityAtElevation(float temperature, float elevation)
	{
		// TODO: make this vary by temperature and elevation
		return WorldData.DensityWater;
	}



	[BurstCompile]
	static public float GetWaterDensity(float waterMass, float saltMass, float temperature, float waterDensityPerSalinity, float waterDensityPerDegree)
	{
		if (waterMass <= 0)
		{
			return 0;
		}
		return WorldData.DensityWater + (waterDensityPerSalinity * saltMass / (waterMass + saltMass) + waterDensityPerDegree * (temperature - WorldData.FreezingTemperature));
	}


	[BurstCompile]
	static public float GetWaterVolume(float waterMass, float saltMass, float temperature, float waterDensityPerSalinity, float waterDensityPerDegree)
	{
		if (waterMass <= 0)
		{
			return 0;
		}
		return (waterMass + saltMass) / GetWaterDensity(waterMass, saltMass, temperature, waterDensityPerSalinity, waterDensityPerDegree);
	}




	//	[BurstCompile]
	//static public float GetEvaporationRate(ref WorldData worldData, float iceCoverage, float temperature, float relativeHumidity, float inverseEvapTemperatureRange)
	//{
	//	float evapTemperature = math.saturate((temperature - worldData.EvapMinTemperature) * inverseEvapTemperatureRange);

	//	return math.saturate((1.0f - iceCoverage) * (1.0f - relativeHumidity) * Utils.Sqr(evapTemperature)) * worldData.EvaporationRate * WorldData.MassWater;
	//}

	//https://en.wikipedia.org/wiki/Dew_point
	[BurstCompile]
	static public float GetDewPoint(float relativeHumidity, float temperature)
	{
		return math.log(relativeHumidity) + (18.678f * (temperature - WorldData.FreezingTemperature)) / (temperature - 16.01f) + WorldData.FreezingTemperature;
	}

	[BurstCompile]
	static public float GetElevationAtDewPoint(float relativeHumidity, float temperature, float referenceElevation)
	{
		float dewPoint = GetDewPoint(relativeHumidity, temperature);
		return referenceElevation + (dewPoint - temperature) / WorldData.TemperatureLapseRate;
	}

	[BurstCompile]
	static public float GetRelativeHumidity(float airMass, float waterVaporMass, float temperature, float dewPointZero, float waterVaporMassToAirMassAtDewPoint, float inverseDewPointTemperatureRange)
	{
		float maxWaterVaporPerKilogramAir = waterVaporMassToAirMassAtDewPoint * Utils.Sqr(math.max(0, (temperature - dewPointZero) * inverseDewPointTemperatureRange));
		float maxWaterVapor = maxWaterVaporPerKilogramAir * airMass;
		if (maxWaterVapor <= 0)
		{
			return waterVaporMass > 0 ? 10000 : 0;
		}
		float relativeHumidity = waterVaporMass / maxWaterVapor;
		return relativeHumidity;
	}

	//	[BurstCompile]
	//static public float GetAirTemperature(float energy, float mass, float waterMass, float waterVaporMass)
	//{
	//	return energy / (mass * WorldData.SpecificHeatAtmosphere + waterMass * WorldData.SpecificHeatWater + waterVaporMass * WorldData.SpecificHeatWaterVapor);
	//}

	[BurstCompile]
	static public float GetWaterTemperature(float energy, float waterMass, float saltMass)
	{
		if (waterMass == 0)
		{
			return 0;
		}
		return math.max(0, energy / (waterMass * WorldData.SpecificHeatWater + saltMass * WorldData.SpecificHeatSalt));
	}

	[BurstCompile]
	static public float GetAlbedo(float surfaceAlbedo, float slope)
	{
		return surfaceAlbedo + (1.0f - surfaceAlbedo) * slope;
	}
	[BurstCompile]
	static public float GetDewPoint(float dewPointTemperaturePerRelativeHumidity, float airTemperature, float relativeHumidity)
	{
		return airTemperature - (1.0f - relativeHumidity) * dewPointTemperaturePerRelativeHumidity;
	}
	[BurstCompile]
	static public float GetCloudElevation(float dewPointElevationPerDegree, float airTemperature, float dewPoint, float elevationOrSeaLevel)
	{
		return elevationOrSeaLevel + math.max(0, (airTemperature - dewPoint) * dewPointElevationPerDegree);
	}

	[BurstCompile]
	static public float GetSpecificHeatOfWater(float waterMass, float saltMass)
	{
		return (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatSalt * saltMass) / (waterMass + saltMass);
	}

	[BurstCompile]
	static public float GetSpecificHeatOfAir(float airMass, float vaporMass)
	{
		return (WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWater * vaporMass) / (airMass + vaporMass);
	}

	[BurstCompile]
	static public float GetSpecificHeatTerrain(float heatingDepth, float soilFertility, float canopyCoverage)
	{
		float landMass = (WorldData.MassSand - WorldData.MassSoil) * soilFertility + WorldData.MassSoil;
		return (WorldData.SpecificHeatSoil * heatingDepth * soilFertility * landMass);
	}


	[BurstCompile]
	static public float GetRadiationRate(float temperature, float emissivity)
	{
		return temperature * temperature * temperature * temperature * emissivity * 0.001f * WorldData.StefanBoltzmannConstant;
	}

	[BurstCompile]
	static public float GetTerrainTemperature(float soilHeatingDepth, float terrainEnergy, float soilFertility, float canopyCoverage)
	{
		float landMass = (WorldData.MassSand - WorldData.MassSoil) * soilFertility + WorldData.MassSoil;
		float heatingDepth = soilFertility * soilHeatingDepth;
		return math.max(0, terrainEnergy / (WorldData.SpecificHeatSoil * heatingDepth * landMass));
	}

}
