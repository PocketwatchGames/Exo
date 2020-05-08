using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;
using Unity.Burst;
using UnityEngine;

public static class Atmosphere {

	static public float GetPotentialTemperature(float temperature, float elevation)
	{
		return temperature - WorldData.TemperatureLapseRate * elevation;
	}

	static public float GetAbsoluteTemperature(float potentialTemperature, float elevation)
	{
		return potentialTemperature + WorldData.TemperatureLapseRate * elevation;
	}


	static public float GetPressureAtElevation(float elevation, float gravity, float referencePressure, float potentialTemperature, float referenceElevation)
	{

		float pressure = referencePressure * math.pow(GetAbsoluteTemperature(potentialTemperature, referenceElevation) / GetAbsoluteTemperature(potentialTemperature, elevation), gravity * WorldData.PressureExponent * WorldData.MolarMassAir);
		return pressure;
	}

	static public float GetAbsolutePressureAtElevation(float elevation, float gravity, float seaLevelPressure, float potentialTemperature)
	{
		float pressure = seaLevelPressure * math.pow(1 + WorldData.TemperatureLapseRate / potentialTemperature * elevation, -gravity * WorldData.PressureExponent * WorldData.MolarMassAir);
		return pressure;
	}


	static public float GetElevationAtPressure(float pressure, float temperaturePotential, float referencePressure, float referenceElevation, float inverseGravity)
	{
		float elevation = referenceElevation + (temperaturePotential * WorldData.TemperatureLapseRateInverse + referenceElevation) * (math.pow(pressure / referencePressure, -inverseGravity * WorldData.UniversalGasConstant * WorldData.TemperatureLapseRate * WorldData.MolarMassAirInverse) - 1);
		return elevation;
	}

	static public float GetStandardPressureAtElevation(float elevation, float potentialTemperature, float gravity)
	{
		float pressure = Atmosphere.GetPressureAtElevation(elevation, gravity, WorldData.StandardPressure, potentialTemperature, 0);
		return pressure;
	}

	//static public float GetStandardAirMass(float elevation, float layerHeight, float gravity)
	//{
	//	float temperatureLapseA = -WorldData.TemperatureLapseRate * elevation;
	//	float massA = WorldData.StaticPressure / (gravity * math.pow(1.0f - (temperatureLapseA) / (WorldData.StdTemp + temperatureLapseA), gravity * WorldData.PressureExponent * WorldData.MolarMassAir));
	//	float temperatureLapseB = -WorldData.TemperatureLapseRate * (elevation + layerHeight);
	//	float massB = WorldData.StaticPressure / (gravity * math.pow(1.0f - (temperatureLapseB) / (WorldData.StdTemp + temperatureLapseB), gravity * WorldData.PressureExponent * WorldData.MolarMassAir));
	//	//	Debug.Log("A: " + massA + " B: " + massB + " E: " + elevation + " H: " + layerHeight);
	//	return massA - massB;
	//}

	static public float GetMolarMassMoistAir(float airMass, float vaporMass)
	{
		return (airMass * WorldData.MolarMassAir + vaporMass * WorldData.MolarMassWater) / (airMass + vaporMass);
	}

	static public float GetAirDensity(float absolutePressure, float absoluteTemperature, float airMass, float vaporMass)
	{
		return absolutePressure * GetMolarMassMoistAir(airMass, vaporMass) / (WorldData.UniversalGasConstant * absoluteTemperature);
	}

	static public float GetInverseAirDensity(float absolutePressure, float temperature, float airMass, float vaporMass)
	{
		return (WorldData.UniversalGasConstant * temperature) / (absolutePressure * GetMolarMassMoistAir(airMass, vaporMass));
	}

	static public float GetWaterDensityAtElevation(float temperature, float elevation)
	{
		return WorldData.DensityWater;
	}

	static public float GetDepthAtPressure(float pressure, float referencePressure, float referenceDepth, float referenceDensity, float gravity)
	{
		return referenceDepth + (pressure - referencePressure) / (gravity * referenceDensity);
	}



	static public float GetWaterDensity(float salinity, float temperature)
	{
		// https://link.springer.com/content/pdf/bbm%3A978-3-319-18908-6%2F1.pdf
		double tempCelsius = temperature - WorldData.FreezingTemperature;
		double standardMeanOceanWaterDensity =
			+999.842594
			+ 0.06793953 * tempCelsius
			- 0.009095290 * tempCelsius * tempCelsius
			+ 0.0001001685 * tempCelsius * tempCelsius * tempCelsius
			- 0.000001120083 * tempCelsius * tempCelsius * tempCelsius * tempCelsius
			+ 0.000000006536332 * tempCelsius * tempCelsius * tempCelsius * tempCelsius * tempCelsius;

		double b1 =
			+0.82449
			- 0.0040899 * tempCelsius
			+ 0.000076438 * tempCelsius * tempCelsius
			- 0.00000082467 * tempCelsius * tempCelsius * tempCelsius
			+ 0.0000000053875 * tempCelsius * tempCelsius * tempCelsius * tempCelsius;

		double c1 =
			-0.0057246
			+ 0.00010227 * tempCelsius
			- 0.0000016546 * tempCelsius * tempCelsius;

		double d0 = 0.00048314;

		double salinityPSU = salinity * 10000;
		double density = standardMeanOceanWaterDensity + b1 * salinityPSU + c1 * math.pow(salinityPSU, 1.5) + d0 * salinityPSU * salinityPSU;
		return (float)density;


		//	return WorldData.DensityWater + (waterDensityPerSalinity * saltMass / (waterMass + saltMass) + waterDensityPerDegree * (temperature - WorldData.FreezingTemperature));
	}

	static public float GetWaterSalinity(float waterMass, float saltMass)
	{
		if (waterMass <= 0)
		{
			return 0;
		}
		return saltMass / (waterMass + saltMass);
	}


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
	static public float GetDewPoint(float relativeHumidity, float referenceTemperature)
	{
		if (relativeHumidity == 0)
		{
			return referenceTemperature;
		}
		return math.log(math.min(1, relativeHumidity)) + (18.678f * (referenceTemperature - WorldData.FreezingTemperature)) / (referenceTemperature - 16.01f) + WorldData.FreezingTemperature;
	}

	static public float GetElevationAtDewPoint(float dewPoint, float potentialTemperature)
	{
		return (dewPoint - potentialTemperature) * WorldData.TemperatureLapseRateInverse;
	}

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

	static public float GetMaxVaporAtTemperature(float airMass, float temperatureAbsolute, float pressure)
	{

		// https://en.wikipedia.org/wiki/Vapour_pressure_of_water
		float saturationPressureOfWaterVaporTentens = 610.78f * math.exp(17.27f * (temperatureAbsolute - 273.15f) / (temperatureAbsolute - 35.85f));

		//https://www.engineeringtoolbox.com/humidity-ratio-air-d_686.html
		return airMass * 0.62198f * saturationPressureOfWaterVaporTentens / math.max(1, pressure - saturationPressureOfWaterVaporTentens);
	}

	static public float GetEvaporationMass(float airMass, float airPressure, float airVapor, float3 wind, float waterTemperature, float maxMass)
	{
		// evap formula from here:
		// https://www.engineeringtoolbox.com/evaporation-water-surface-d_690.html
		// NOTE: I've made adjustments to this because my finite differencing sometimes means that the water surface and air temperature are a bit out of sync
		// so i'm using the air temperature instead of the water temperature, which means the the formula just reduces to (1-RH)*WindCoefficient
		float evaporationCoefficient = (25 + 19 * math.length(wind));
		return math.clamp(evaporationCoefficient * (Atmosphere.GetMaxVaporAtTemperature(airMass, waterTemperature, airPressure) - airVapor) / airMass, 0, maxMass);
	}

	//static public float GetAirTemperature(float energy, float mass, float waterMass, float waterVaporMass)
	//{
	//	return energy / (mass * WorldData.SpecificHeatAtmosphere + waterMass * WorldData.SpecificHeatWater + waterVaporMass * WorldData.SpecificHeatWaterVapor);
	//}

	static public float GetWaterTemperature(float energy, float waterMass, float saltMass)
	{
		if (waterMass == 0)
		{
			return 0;
		}
		return math.max(0, energy / (waterMass * WorldData.SpecificHeatWater + saltMass * WorldData.SpecificHeatSalt));
	}

	static public float GetAlbedo(float surfaceAlbedo, float slope)
	{
		return surfaceAlbedo + (1.0f - surfaceAlbedo) * slope;
	}

	static public float GetDewPoint(float dewPointTemperaturePerRelativeHumidity, float airTemperature, float relativeHumidity)
	{
		return airTemperature - (1.0f - relativeHumidity) * dewPointTemperaturePerRelativeHumidity;
	}

	static public float GetCloudElevation(float dewPointElevationPerDegree, float airTemperature, float dewPoint, float elevationOrSeaLevel)
	{
		return elevationOrSeaLevel + math.max(0, (airTemperature - dewPoint) * dewPointElevationPerDegree);
	}

	static public float GetSpecificHeatWater(float waterMass, float saltMass)
	{
		return WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatSalt * saltMass;
	}

	static public float GetSpecificHeatAir(float airMass, float vaporMass, float cloudMass)
	{
		return (airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor + cloudMass * WorldData.SpecificHeatWater);
	}

	static public float GetCloudMassInLayer(float cloudMass, float cloudElevation, float layerElevation, float layerHeight)
	{
		if (cloudElevation >= layerElevation && cloudElevation < layerElevation + layerHeight)
		{
			return cloudMass;
		}
		return 0;
	}

	static public float GetSpecificHeatTerrain(float heatingDepth, float soilFertility, float floraMass, float floraWater)
	{
		float sandPercent = soilFertility;
		float massSand = heatingDepth * WorldData.MassSand * sandPercent;
		float massSoil = heatingDepth * WorldData.MassSand * (1.0f - sandPercent);
		return WorldData.SpecificHeatSoil * massSoil + WorldData.SpecificHeatSand * massSand + WorldData.SpecificHeatFlora * floraMass + WorldData.SpecificHeatWater * floraWater;
	}


	static public float GetRadiationRate(float temperature, float emissivity)
	{
		return temperature * temperature * temperature * temperature * emissivity * 0.001f * WorldData.StefanBoltzmannConstant;
	}

	static public float GetDropletRadius(float mass, float waterDensity)
	{
		const float ThreeQuarterInversePi = 3 / (4 * math.PI);
		return math.pow(mass / waterDensity * ThreeQuarterInversePi, 0.333f);
	}

	public static float GetDiffusionAmount(float massA, float massB, float diffusionCoefficient, float surfaceAreaA, float surfaceAreaB, float distInverse)
	{
		return massB * diffusionCoefficient * (surfaceAreaA + surfaceAreaB) * distInverse / (2 * (massA + massB));
	}

	public static float GetDiffusionAmount(float massA, float massB, float diffusionCoefficient, float dist)
	{
		return massB * diffusionCoefficient / (dist * (massA + massB));
	}

}
