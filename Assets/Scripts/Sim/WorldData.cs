﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

[Serializable]
public struct WorldData {
	public float SecondsPerTick;
	public int AirLayers;
	public int WaterLayers;
	public float FullIceCoverage;
	public float FullWaterCoverage;
	public float FullVegetationCoverage;

	[Header("Fresh Water")]
	public float MaxSoilPorousness;
	public float GroundWaterReplenishmentSpeed;
	public float GroundWaterFlowSpeed;

	[Header("Pressure and Wind")]
	public float AirDiffusionCoefficientHorizontal;
	public float AirDiffusionCoefficientVertical;
	public float WaterDiffusionCoefficientHorizontal;
	public float WaterDiffusionCoefficientVertical;
	public float WindWaterFriction;
	public float WindIceFriction;
	public float WindTerrainFrictionMin;
	public float WindTerrainFrictionMax;
	public float WindVegetationFriction;
	public float CloudDiffusionCoefficient;
	public float MaxTerrainRoughnessForWindFriction;
	public float WaterSurfaceFrictionDepth;
	public float WindToWaterCurrentFrictionCoefficient;

	[Header("Atmospheric Energy Cycle")]
	// atmospheric heat balance https://energyeducation.ca/encyclopedia/Earth%27s_heat_balance
	// https://en.wikipedia.org/wiki/Earth%27s_energy_budget
	// https://en.wikipedia.org/wiki/Electromagnetic_absorption_by_water
	// Water vapor is responsible for 70% of solar absorption and about 60% of absorption of thermal radiation.
	public float SolarAbsorptivityAir; // total absorbed by atmosphere AFTER reflection about 30%
	public float SolarAbsorptivityWaterVapor; // total absorbed by atmosphere AFTER reflection about 30%
	public float SolarAbsorptivityCloud; // 6% absorbed by clouds
	public float SolarReflectivityAir; // 7% is reflected due to atmospheric scattering 
	public float SolarReflectivityWater; // 7% is reflected due to atmospheric scattering 
	public float SolarReflectivityCloud;

	//public float EvaporativeHeatLoss = 0.6f; // global average = 78 watts
	// Net Back Radiation: The ocean transmits electromagnetic radiation into the atmosphere in proportion to the fourth power of the sea surface temperature(black-body radiation)
	// https://eesc.columbia.edu/courses/ees/climate/lectures/o_atm.html
	public float AirWaterConductionPositive; // global avg = 16 watts per degree delta between air and ocean (global avg = 24 watts per m^2 of ocean)

	public float maxGroundWaterTemperature;
	public float SoilHeatDepth;
	public float AlbedoReductionSoilQuality;

	// TODO: tune these to match the science
	public float ThermalReflectivityCloud;
	public float CloudMassFullAbsorption; // how much heat gain/loss is caused by cloud cover (cumulus cloud is 0.3g/cubic meter, and about 3 kilometers high)
	public float EnergyLostThroughAtmosphereWindow; // AKA Atmospheric window global average = 40 watts = 6.7% of all surface and atmospheric radiation
	public float minCloudFreezingTemperature;
	public float maxCloudFreezingTemperature;
	public float maxCloudSlopeAlbedo;
	public float rainDropSizeAlbedoMin;
	public float rainDropSizeAlbedoMax;

	// https://en.wikipedia.org/wiki/Electromagnetic_absorption_by_water
	// Water vapor is responsible for 70% of solar absorption and about 60% of absorption of thermal radiation.
	// carbon dioxide accounts for just 26% of the greenhouse effect.
	public float AbsorptivityCarbonDioxide;
	public float AbsorptivityAir;
	public float AbsorptivityWaterLiquid;
	public float AbsorptivityWaterVapor;

	[Header("Evap, Humidity and Clouds")]
	public float DewPointTemperatureRange;
	public float DewPointZero;
	public float WaterVaporMassToAirMassAtDewPoint;
	public float EvapMinTemperature; // -30 celsius
	public float EvapMaxTemperature; // 70 celsius
	public float EvaporationRate; // TODO: evaporation on earth maxes out around 2.5M per year 
	public float rainDropDragCoefficient;
	public float rainDropMaxSize;
	public float rainDropMinSize;
	public float CloudDissapationRateWind;
	public float CloudDissapationRateDryAir;
	public float DewPointElevationPerDegree;
	public float DewPointTemperaturePerRelativeHumidity;
	public float IceHeatingDepth;
	public float WaterHeatingDepth;

	[Header("Water")]
	public float WaterDensityPerSalinity;
	public float WaterDensityPerDegree;
	public float WaterDensityCurrentSpeed;

	[Header("Ecology")]
	public float MinTemperatureCanopy;
	public float MaxTemperatureCanopy;

	[Header("Physical Constants")]
	public const float TemperatureLapseRate = -0.0065f;
	public const float AdiabaticLapseRate = 0.0098f;
	public const float StaticPressure = 101325;
	public const float StdTemp = 288.15f;
	public const float MolarMassAir = 0.0289647f;
	public const float MolarMassWater = 0.01802f;
	public const float UniversalGasConstant = 8.3144598f;
	public const float FreezingTemperature = 273.15f;
	public const float StefanBoltzmannConstant = 0.00000005670373f;
	// specific heat is joules to raise one degree (kJ/kgK)
	public const float AlbedoWater = 0.06f; // How much heat is reflected back by the water
	public const float AlbedoIce = 0.5f; // How much heat is reflected back by the water
	public const float AlbedoLand = 0.4f;
	public const float AlbedoFoliage = 0.1f;
	//public const float AlbedoCloud = 0.05f; // 24% incoming  reflected back to space by clouds (avg, globally)
	public const float SpecificHeatIce = 2.108f;
	public const float SpecificHeatWater = 4.187f;
	public const float SpecificHeatWaterVapor = 1.996f;
	public const float SpecificHeatSalt = 0.85f;
	public const float SpecificHeatAtmosphere = 1.158f;
	public const float SpecificHeatSoil = 0.84f;
	public const float LatentHeatWaterLiquid = 334.0f;
	public const float LatentHeatWaterVapor = 2264.705f;
	// emissivity values obtained here: https://www.thermoworks.com/emissivity-table
	// and here https://www.aspen-electronics.com/uploads/3/7/1/2/37123419/emissivity-table.pdf
	public const float EmissivityWater = 0.96f;
	public const float EmissivitySalt = 0.34f;
	public const float EmissivityIce = 0.97f;
	public const float EmissivityDirt = 0.92f;
	public const float EmissivitySand = 0.76f;
	public const float EmissivityAir = 0.8f;
	public const float EmissivityWaterVapor = 0.4f;
	public const float EmissivityVegetation = 0.95f;
	public const float MassEarthAir = 1.29f;
	public const float MassWater = 1000f;
	public const float MassSalt = 2170f;
	public const float MassIce = 919f;
	public const float MassSoil = 1200f;
	public const float MassSand = 1600f;
	public const float DensityWater = 997f;
	public const float DensityAir = 1.21f;
	public const float ConductivityAir = 0.0262f;
	public const float ConductivityWater = 0.606f;
	public const float ConductivityIce = 2.18f;
	public const float ConductivityTerrain = 0.2f;
	public const float ThermalContactResistance = 0.00005f;
	public const float ConductivityAirWater = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityWater + ThermalContactResistance);
	public const float ConductivityAirIce = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityIce + ThermalContactResistance);
	public const float ConductivityAirTerrain = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityTerrain + ThermalContactResistance);
	public const float ConductivityIceWater = 1.0f / (1.0f / ConductivityWater + 1.0f / ConductivityIce + ThermalContactResistance);
	public const float ConductivityIceTerrain = 1.0f / (1.0f / ConductivityTerrain + 1.0f / ConductivityIce + ThermalContactResistance);
	public const float ConductivityWaterTerrain = 1.0f / (1.0f / ConductivityTerrain + 1.0f / ConductivityWater + ThermalContactResistance);
	public const float GasConstantAir = UniversalGasConstant / MolarMassAir * 1000;
	public const float GasConstantWaterVapor = UniversalGasConstant / MolarMassWater * 1000;
	public const float PressureExponent = 1.0f / (UniversalGasConstant * TemperatureLapseRate);
	public const float DryAirAdiabaticLapseRate = AdiabaticLapseRate / SpecificHeatAtmosphere;
	public const float inverseSpecificHeatIce = 1.0f / SpecificHeatIce;
	public const float InverseDensityAir = 1.0f / DensityAir;

	[NonSerialized]	public float EvapTemperatureRange;
	[NonSerialized] public float TicksPerSecond;
	[NonSerialized] public float TicksPerYear;
	[NonSerialized] public float inverseFullCanopyCoverage;
	[NonSerialized] public float inverseFullWaterCoverage;
	[NonSerialized] public float inverseFullIceCoverage;
	[NonSerialized] public float inverseCloudMassFullAbsorption;
	[NonSerialized] public float wattsToKJPerTick;
	[NonSerialized] public float declinationOfSun;
	[NonSerialized] public float sunHitsAtmosphereBelowHorizonAmount;
	[NonSerialized] public float inverseSunAtmosphereAmount;
	[NonSerialized] public float inverseDewPointTemperatureRange;
	[NonSerialized] public float inverseEvapTemperatureRange;

	public void Init()
	{
		EvapTemperatureRange = EvapMaxTemperature - EvapMinTemperature;

		TicksPerSecond = 1.0f / SecondsPerTick;
		TicksPerYear = 60 * 60 * 24 * 365 / SecondsPerTick;

		inverseFullCanopyCoverage = 1.0f / FullVegetationCoverage;
		inverseFullWaterCoverage = 1.0f / FullWaterCoverage;
		inverseFullIceCoverage = 1.0f / (MassIce * FullIceCoverage);
		inverseCloudMassFullAbsorption = 1.0f / CloudMassFullAbsorption;
		wattsToKJPerTick = SecondsPerTick * 1000;
		sunHitsAtmosphereBelowHorizonAmount = 0.055f;
		inverseSunAtmosphereAmount = 1.0f / (1.0f + sunHitsAtmosphereBelowHorizonAmount);
		inverseDewPointTemperatureRange = 1.0f / DewPointTemperatureRange;
		inverseEvapTemperatureRange = 1.0f / EvapTemperatureRange;

	}

}
