using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;

public static class Atmosphere {
	static public float GetAirPressure(WorldData worldData, float mass, float elevation, float temperature, float molarMass)
	{
		float temperatureLapse = -worldData.TemperatureLapseRate * elevation;
		float pressure = mass * worldData.GravitationalAcceleration * math.pow(1.0f - (temperatureLapse) / (temperature + temperatureLapse), worldData.PressureExponent * molarMass);
		return pressure;
	}

	static public float GetAirMass(WorldData worldData, float pressure, float elevation, float temperature, float molarMass)
	{
		float temperatureLapse = -worldData.TemperatureLapseRate * elevation;
		float mass = pressure / (worldData.GravitationalAcceleration * math.pow(1.0f - (temperatureLapse) / (temperature + temperatureLapse), worldData.PressureExponent * molarMass));
		return mass;
	}

	static public float GetMolarMassAir(WorldData worldData, float airMass, float waterMass)
	{
		return (airMass * worldData.MolarMassAir + waterMass * worldData.MolarMassWater) / (airMass + waterMass);
	}

	static public float GetAirDensity(WorldData worldData, float absolutePressure, float temperature, float molarMassAir)
	{
		return absolutePressure * molarMassAir / (worldData.UniversalGasConstant * temperature);
	}

	static public float GetWaterDensity(WorldData worldData, float oceanEnergy, float saltMass, float mass)
	{
		if (mass <= 0)
		{
			return 0;
		}
		return worldData.waterDensity + (worldData.OceanDensityPerSalinity * saltMass / (mass + saltMass) + worldData.OceanDensityPerDegree * (GetWaterTemperature(worldData, oceanEnergy, mass, saltMass) - worldData.FreezingTemperature));
	}

	static public float GetWaterDensityByTemperature(WorldData worldData, float temperature, float saltMass, float mass)
	{
		if (mass <= 0)
		{
			return 0;
		}
		return worldData.waterDensity + (worldData.OceanDensityPerSalinity * saltMass / (mass + saltMass) + worldData.OceanDensityPerDegree * (temperature - worldData.FreezingTemperature));
	}

	static public float GetWaterVolume(WorldData worldData, float mass, float salt, float energy)
	{
		if (mass <= 0)
		{
			return 0;
		}
		return (mass + salt) / GetWaterDensity(worldData, energy, salt, mass);
	}

	static public float GetWaterVolumeByTemperature(WorldData worldData, float mass, float salt, float temperature)
	{
		if (mass <= 0)
		{
			return 0;
		}
		return (mass + salt) / GetWaterDensityByTemperature(worldData, temperature, salt, mass);
	}



	static public float GetPressureAtElevation(WorldData worldData, int index, float elevation, float molarMass)
	{
		// Units: Pascals
		// Barometric Formula
		// Pressure = StaticPressure * (StdTemp / (StdTemp + StdTempLapseRate * (Elevation - ElevationAtBottomOfAtmLayer)) ^ (GravitationalAcceleration * MolarMassOfEarthAir / (UniversalGasConstant * StdTempLapseRate))
		// https://en.wikipedia.org/wiki/Barometric_formula
		// For the bottom layer of atmosphere ( < 11000 meters), ElevationAtBottomOfAtmLayer == 0)

		//	float standardPressure = Data.StaticPressure * (float)Math.Pow(Data.StdTemp / (Data.StdTemp + Data.StdTempLapseRate * elevation), Data.PressureExponent);
		float pressure = worldData.StaticPressure * (float)Math.Pow(worldData.StdTemp / (worldData.StdTemp + worldData.TemperatureLapseRate * elevation), worldData.PressureExponent * molarMass);
		return pressure;
	}

	static private float GetEvaporationRate(WorldData worldData, float iceCoverage, float temperature, float relativeHumidity, float inverseEvapTemperatureRange)
	{
		float evapTemperature = math.clamp((temperature - worldData.EvapMinTemperature) * inverseEvapTemperatureRange, 0, 1);

		return math.clamp((1.0f - iceCoverage) * (1.0f - relativeHumidity) * Utils.Sqr(evapTemperature), 0, 1) * worldData.EvaporationRate * worldData.MassWater;
	}

	static public float GetTemperatureAtElevation(WorldData worldData, float elevation, float lowerTemperature, float upperTemperature, float elevationOrSeaLevel)
	{
		float temperatureLapseRate = (upperTemperature - lowerTemperature) / worldData.BoundaryZoneElevation;
		return (elevation - elevationOrSeaLevel) * temperatureLapseRate + lowerTemperature;
	}

	static public float GetRelativeHumidity(WorldData worldData, float temperature, float humidity, float airMass, float inverseDewPointTemperatureRange)
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

	static public float GetAbsoluteHumidity(WorldData worldData, float temperature, float relativeHumidity, float totalAtmosphericMass, float inverseDewPointTemperatureRange)
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

	static public float GetAirTemperature(WorldData worldData, float energy, float mass, float waterMass, float waterSpecificHeat)
	{
		return energy / (mass * worldData.SpecificHeatAtmosphere + waterMass * waterSpecificHeat);
	}
	static public float GetAirEnergy(WorldData worldData, float temperature, float mass, float waterMass, float waterSpecificHeat)
	{
		return temperature * (mass * worldData.SpecificHeatAtmosphere + waterMass * waterSpecificHeat);
	}

	static public float GetWaterTemperature(WorldData worldData, float energy, float waterMass, float saltMass)
	{
		if (waterMass == 0)
		{
			return 0;
		}
		return math.max(0, energy / (waterMass * worldData.SpecificHeatWater + saltMass * worldData.SpecificHeatSalt));
	}
	static public float GetWaterEnergy(WorldData worldData, float temperature, float waterMass, float saltMass)
	{
		return temperature * (worldData.SpecificHeatWater * waterMass + worldData.SpecificHeatSalt * saltMass);
	}
	static public float GetAlbedo(float surfaceAlbedo, float slope)
	{
		return surfaceAlbedo + (1.0f - surfaceAlbedo) * slope;
	}
	static public float GetDewPoint(WorldData worldData, float lowerAirTemperature, float relativeHumidity)
	{
		return lowerAirTemperature - (1.0f - relativeHumidity) * worldData.DewPointTemperaturePerRelativeHumidity;
	}
	static public float GetCloudElevation(WorldData worldData, float upperAirTemperature, float dewPoint, float elevationOrSeaLevel)
	{
		return elevationOrSeaLevel + math.max(0, worldData.BoundaryZoneElevation + (upperAirTemperature - dewPoint) * worldData.DewPointElevationPerDegree);
	}

	static public float GetSpecificHeatOfWater(WorldData worldData, float waterMass, float saltMass)
	{
		return (worldData.SpecificHeatWater * waterMass + worldData.SpecificHeatSalt * saltMass) / (waterMass + saltMass);
	}

	static public float GetAtmosphericEmissivity(WorldData worldData, float airMass, float greenhouseGasMass, float humidity, float cloudMass)
	{
		return airMass * worldData.AbsorptivityAir +
			greenhouseGasMass * worldData.AbsorptivityCarbonDioxide +
			humidity * worldData.AbsorptivityWaterVapor +
			cloudMass * worldData.AbsorptivityWaterLiquid;
	}

	static public float GetLandRadiationRate(WorldData worldData, float landEnergy, float groundWater, float soilFertility, float canopyCoverage)
	{
		float soilTemperature = GetLandTemperature(worldData, landEnergy, groundWater, soilFertility, canopyCoverage);
		return GetRadiationRate(worldData, soilTemperature, worldData.EmissivityDirt) * (1.0f - canopyCoverage * 0.5f);
	}

	static public float GetRadiationRate(WorldData worldData, float temperature, float emissivity)
	{
		return temperature * temperature * temperature * temperature * emissivity * 0.001f * worldData.StefanBoltzmannConstant;
	}

	static public float GetLandTemperature(WorldData worldData, float landEnergy, float groundWater, float soilFertility, float canopyCoverage)
	{
		float soilEnergy = landEnergy - groundWater * worldData.maxGroundWaterTemperature * worldData.SpecificHeatWater;
		float landMass = (worldData.MassSand - worldData.MassSoil) * soilFertility + worldData.MassSoil;
		float heatingDepth = soilFertility * worldData.SoilHeatDepth;
		return math.max(0, soilEnergy / (worldData.SpecificHeatSoil * heatingDepth * landMass));
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
