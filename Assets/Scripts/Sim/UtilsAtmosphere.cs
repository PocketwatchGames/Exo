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
	static public float GetAbsoluteTemperature(float potentialTemperature, float elevation)
	{
		return potentialTemperature + WorldData.TemperatureLapseRate * elevation;
	}


	[BurstCompile]
	static public float GetPressureAtElevation(float elevation, float gravity, float referencePressure, float potentialTemperature, float referenceElevation)
	{

		float pressure = referencePressure * math.pow(GetAbsoluteTemperature(potentialTemperature, referenceElevation) / GetAbsoluteTemperature(potentialTemperature, elevation), gravity * WorldData.PressureExponent * WorldData.MolarMassAir);
		return pressure;
	}

	[BurstCompile]
	static public float GetAbsolutePressureAtElevation(float elevation, float gravity, float seaLevelPressure, float potentialTemperature)
	{
		float pressure = seaLevelPressure * math.pow(1 + WorldData.TemperatureLapseRate / potentialTemperature * elevation, -gravity * WorldData.PressureExponent * WorldData.MolarMassAir);
		return pressure;
	}


	[BurstCompile]
	static public float GetElevationAtPressure(float pressure, float temperaturePotential, float referencePressure, float referenceElevation, float gravity)
	{
		float elevation = referenceElevation + (temperaturePotential / WorldData.TemperatureLapseRate + referenceElevation) * (math.pow(pressure / referencePressure, -1.0f / (gravity * WorldData.PressureExponent * WorldData.MolarMassAir)) - 1);
		return elevation;
	}

	[BurstCompile]
	static public float GetStandardPressureAtElevation(float elevation, float potentialTemperature, float gravity)
	{
		float pressure = Atmosphere.GetPressureAtElevation(elevation, gravity, WorldData.StandardPressure, potentialTemperature, 0);
		return pressure;
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
		return (airMass * WorldData.MolarMassAir + vaporMass * WorldData.MolarMassWater) / (airMass + vaporMass);
	}

	[BurstCompile]
	static public float GetAirDensity(float absolutePressure, float temperature, float airMass, float vaporMass)
	{
		return absolutePressure * GetMolarMassMoistAir(airMass, vaporMass) / (WorldData.UniversalGasConstant * temperature);
	}

	[BurstCompile]
	static public float GetInverseAirDensity(float absolutePressure, float temperature, float airMass, float vaporMass)
	{
		return (WorldData.UniversalGasConstant * temperature) / (absolutePressure * GetMolarMassMoistAir(airMass, vaporMass));
	}

	[BurstCompile]
	static public float GetWaterDensityAtElevation(float temperature, float elevation)
	{
		return WorldData.DensityWater;
	}

	[BurstCompile]
	static public float GetDepthAtPressure(float pressure, float referencePressure, float referenceDepth, float referenceDensity, float gravity)
	{
		return referenceDepth + (pressure - referencePressure) / (gravity * referenceDensity);
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
	static public float GetWaterSalinity(float waterMass, float saltMass)
	{
		if (waterMass <= 0)
		{
			return 0;
		}
		return saltMass / (waterMass + saltMass);
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

	[BurstCompile]
	static public float GetFreezingPoint(float salinity, float freezePointReductionPerSalinity)
	{
		return WorldData.FreezingTemperature - salinity * freezePointReductionPerSalinity;
	}




	//	[BurstCompile]
	//static public float GetEvaporationRate(ref WorldData worldData, float iceCoverage, float temperature, float relativeHumidity, float inverseEvapTemperatureRange)
	//{
	//	float evapTemperature = math.saturate((temperature - worldData.EvapMinTemperature) * inverseEvapTemperatureRange);

	//	return math.saturate((1.0f - iceCoverage) * (1.0f - relativeHumidity) * Utils.Sqr(evapTemperature)) * worldData.EvaporationRate * WorldData.MassWater;
	//}

	//https://en.wikipedia.org/wiki/Dew_point
	[BurstCompile]
	static public float GetDewPoint(float relativeHumidity, float referenceTemperature)
	{
		if (relativeHumidity == 0)
		{
			return referenceTemperature;
		}
		return math.log(math.min(1, relativeHumidity)) + (18.678f * (referenceTemperature - WorldData.FreezingTemperature)) / (referenceTemperature - 16.01f) + WorldData.FreezingTemperature;
	}

	[BurstCompile]
	static public float GetElevationAtDewPoint(float dewPoint, float potentialTemperature)
	{
		return (dewPoint - potentialTemperature) / WorldData.TemperatureLapseRate;
	}

	[BurstCompile]
	static public float GetRelativeHumidity(float airMass, float waterVaporMass, float temperature, float pressure)
	{
		float maxWaterVapor = GetMaxVaporAtTemperature(airMass, temperature, pressure);
		if (maxWaterVapor <= 0)
		{
			return waterVaporMass > 0 ? 10000 : 0;
		}
		float relativeHumidity = waterVaporMass / maxWaterVapor;
		return relativeHumidity;
	}

	[BurstCompile]
	static public float GetMaxVaporAtTemperature(float airMass, float temperatureAbsolute, float pressure)
	{
		// https://www.engineeringtoolbox.com/water-vapor-saturation-pressure-air-d_689.html
		float saturationPressureOfWaterVapor = math.exp(77.345f + 0.0057f * temperatureAbsolute - 7235 / temperatureAbsolute) / math.pow(temperatureAbsolute, 8.2f);

		//https://www.engineeringtoolbox.com/humidity-ratio-air-d_686.html
		return airMass * 0.62198f * saturationPressureOfWaterVapor / (pressure - saturationPressureOfWaterVapor);
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
	static public float GetSpecificHeatTerrain(float heatingDepth, float soilFertility, float vegetationMass)
	{
		float landMass = (WorldData.MassSand * (1.0f - soilFertility) + WorldData.MassSoil * soilFertility) * heatingDepth;
		return (WorldData.SpecificHeatSoil * landMass + WorldData.SpecificHeatVegetation * vegetationMass);
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

	[BurstCompile]
	static public float GetDropletRadius(float mass, float waterDensity)
	{
		const float ThreeQuarterInversePi = 3 / (4 * math.PI);
		return math.pow(mass / waterDensity * ThreeQuarterInversePi, 0.333f);
	}

}
