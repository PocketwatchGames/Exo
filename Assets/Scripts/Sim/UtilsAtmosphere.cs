using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;

public static class Atmosphere {
	static public float GetAirPressure(ref WorldData worldData, float mass, float elevation, float temperature, float molarMass, float gravity)
	{
		float temperatureLapse = -WorldData.TemperatureLapseRate * elevation;
		float pressure = mass * gravity * math.pow(1.0f - (temperatureLapse) / (temperature + temperatureLapse), gravity * worldData.PressureExponent * molarMass);
		return pressure;
	}

	static public float GetAirMass(ref WorldData worldData, float pressure, float elevation, float temperature, float molarMass, float gravity)
	{
		float temperatureLapse = -WorldData.TemperatureLapseRate * elevation;
		float mass = pressure / (gravity * math.pow(1.0f - (temperatureLapse) / (temperature + temperatureLapse), gravity * worldData.PressureExponent * molarMass));
		return mass;
	}

	static public float GetMolarMassAir(float airMass, float waterMass)
	{
		return (airMass * WorldData.MolarMassAir + waterMass * WorldData.MolarMassWater) / (airMass + waterMass);
	}

	static public float GetAirDensity(float absolutePressure, float temperature, float molarMassAir)
	{
		return absolutePressure * molarMassAir / (WorldData.UniversalGasConstant * temperature);
	}

	static public float GetWaterDensity(ref WorldData worldData, float oceanEnergy, float saltMass, float mass)
	{
		if (mass <= 0)
		{
			return 0;
		}
		return worldData.waterDensity + (worldData.OceanDensityPerSalinity * saltMass / (mass + saltMass) + worldData.OceanDensityPerDegree * (GetWaterTemperature(oceanEnergy, mass, saltMass) - WorldData.FreezingTemperature));
	}

	static public float GetWaterDensityByTemperature(ref WorldData worldData, float temperature, float saltMass, float mass)
	{
		if (mass <= 0)
		{
			return 0;
		}
		return worldData.waterDensity + (worldData.OceanDensityPerSalinity * saltMass / (mass + saltMass) + worldData.OceanDensityPerDegree * (temperature - WorldData.FreezingTemperature));
	}

	static public float GetWaterVolume(ref WorldData worldData, float mass, float salt, float energy)
	{
		if (mass <= 0)
		{
			return 0;
		}
		return (mass + salt) / GetWaterDensity(ref worldData, energy, salt, mass);
	}

	static public float GetWaterVolumeByTemperature(ref WorldData worldData, float mass, float salt, float temperature)
	{
		if (mass <= 0)
		{
			return 0;
		}
		return (mass + salt) / GetWaterDensityByTemperature(ref worldData, temperature, salt, mass);
	}



	static public float GetPressureAtElevation(ref WorldData worldData, int index, float elevation, float molarMass, float gravity)
	{
		// Units: Pascals
		// Barometric Formula
		// Pressure = StaticPressure * (StdTemp / (StdTemp + StdTempLapseRate * (Elevation - ElevationAtBottomOfAtmLayer)) ^ (GravitationalAcceleration * MolarMassOfEarthAir / (UniversalGasConstant * StdTempLapseRate))
		// https://en.wikipedia.org/wiki/Barometric_formula
		// For the bottom layer of atmosphere ( < 11000 meters), ElevationAtBottomOfAtmLayer == 0)

		//	float standardPressure = Data.StaticPressure * (float)Math.Pow(Data.StdTemp / (Data.StdTemp + Data.StdTempLapseRate * elevation), Data.PressureExponent);
		float pressure = WorldData.StaticPressure * (float)Math.Pow(WorldData.StdTemp / (WorldData.StdTemp + WorldData.TemperatureLapseRate * elevation), gravity * worldData.PressureExponent * molarMass);
		return pressure;
	}

	static public float GetEvaporationRate(ref WorldData worldData, float iceCoverage, float temperature, float relativeHumidity, float inverseEvapTemperatureRange)
	{
		float evapTemperature = math.saturate((temperature - worldData.EvapMinTemperature) * inverseEvapTemperatureRange);

		return math.saturate((1.0f - iceCoverage) * (1.0f - relativeHumidity) * Utils.Sqr(evapTemperature)) * worldData.EvaporationRate * WorldData.MassWater;
	}

	static public float GetTemperatureAtElevation(ref WorldData worldData, float elevation, float lowerTemperature, float upperTemperature, float elevationOrSeaLevel)
	{
		float temperatureLapseRate = (upperTemperature - lowerTemperature) / worldData.BoundaryZoneElevation;
		return (elevation - elevationOrSeaLevel) * temperatureLapseRate + lowerTemperature;
	}

	static public float GetRelativeHumidity(ref WorldData worldData, float temperature, float humidity, float airMass, float inverseDewPointTemperatureRange)
	{
		float maxWaterVaporPerKilogramAir = worldData.WaterVaporMassToAirMassAtDewPoint * Utils.Sqr(math.max(0, (temperature - worldData.DewPointZero) * inverseDewPointTemperatureRange));
		float maxHumidity = maxWaterVaporPerKilogramAir * airMass;
		if (maxHumidity <= 0)
		{
			return humidity > 0 ? 10000 : 0;
		}
		float relativeHumidity = humidity / maxHumidity;
		return relativeHumidity;
	}

	static public float GetAbsoluteHumidity(ref WorldData worldData, float temperature, float relativeHumidity, float totalAtmosphericMass, float inverseDewPointTemperatureRange)
	{
		float maxWaterVaporPerKilogramAtmosphere = worldData.WaterVaporMassToAirMassAtDewPoint * Utils.Sqr(math.max(0, (temperature - worldData.DewPointZero) * inverseDewPointTemperatureRange)) / (1.0f + worldData.WaterVaporMassToAirMassAtDewPoint);
		float maxHumidity = maxWaterVaporPerKilogramAtmosphere * totalAtmosphericMass;
		if (maxHumidity <= 0)
		{
			return 0;
		}
		float humidity = relativeHumidity * maxHumidity;
		return humidity;
	}

	static public float GetAirTemperature(float energy, float mass, float waterMass, float waterVaporMass)
	{
		return energy / (mass * WorldData.SpecificHeatAtmosphere + waterMass * WorldData.SpecificHeatWater + waterVaporMass * WorldData.SpecificHeatWaterVapor);
	}
	static public float GetAirEnergy(float temperature, float mass, float waterMass, float waterVaporMass)
	{
		return temperature * (mass * WorldData.SpecificHeatAtmosphere + waterMass * WorldData.SpecificHeatWater + waterVaporMass * WorldData.SpecificHeatWaterVapor);
	}

	static public float GetWaterTemperature(float energy, float waterMass, float saltMass)
	{
		if (waterMass == 0)
		{
			return 0;
		}
		return math.max(0, energy / (waterMass * WorldData.SpecificHeatWater + saltMass * WorldData.SpecificHeatSalt));
	}
	static public float GetWaterEnergy(float temperature, float waterMass, float saltMass)
	{
		return temperature * (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatSalt * saltMass);
	}
	static public float GetAlbedo(float surfaceAlbedo, float slope)
	{
		return surfaceAlbedo + (1.0f - surfaceAlbedo) * slope;
	}
	static public float GetDewPoint(ref WorldData worldData, float lowerAirTemperature, float relativeHumidity)
	{
		return lowerAirTemperature - (1.0f - relativeHumidity) * worldData.DewPointTemperaturePerRelativeHumidity;
	}
	static public float GetCloudElevation(ref WorldData worldData, float airTemperature, float dewPoint, float elevationOrSeaLevel)
	{
		return elevationOrSeaLevel + math.max(0, (airTemperature - dewPoint) * worldData.DewPointElevationPerDegree);
	}

	static public float GetSpecificHeatOfWater(float waterMass, float saltMass)
	{
		return (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatSalt * saltMass) / (waterMass + saltMass);
	}

	static public float GetAtmosphericEmissivity(ref WorldData worldData, float airMass, float greenhouseGasMass, float humidity, float cloudMass)
	{
		return airMass * worldData.AbsorptivityAir +
			greenhouseGasMass * worldData.AbsorptivityCarbonDioxide +
			humidity * worldData.AbsorptivityWaterVapor +
			cloudMass * worldData.AbsorptivityWaterLiquid;
	}

	static public float GetLandRadiationRate(ref WorldData worldData, float landEnergy, float groundWater, float soilFertility, float canopyCoverage)
	{
		float soilTemperature = GetLandTemperature(ref worldData, landEnergy, groundWater, soilFertility, canopyCoverage);
		return GetRadiationRate(soilTemperature, WorldData.EmissivityDirt) * (1.0f - canopyCoverage * 0.5f);
	}

	static public float GetRadiationRate(float temperature, float emissivity)
	{
		return temperature * temperature * temperature * temperature * emissivity * 0.001f * WorldData.StefanBoltzmannConstant;
	}

	static public float GetLandTemperature(ref WorldData worldData, float landEnergy, float groundWater, float soilFertility, float canopyCoverage)
	{
		float soilEnergy = landEnergy - groundWater * worldData.maxGroundWaterTemperature * WorldData.SpecificHeatWater;
		float landMass = (WorldData.MassSand - WorldData.MassSoil) * soilFertility + WorldData.MassSoil;
		float heatingDepth = soilFertility * worldData.SoilHeatDepth;
		return math.max(0, soilEnergy / (WorldData.SpecificHeatSoil * heatingDepth * landMass));
	}

	static public float RepeatExclusive(float x, float y)
	{
		while (x < 0)
		{
			x += y;
		}
		while (x >= y)
		{
			x -= y;
		}
		return x;
	}
}
