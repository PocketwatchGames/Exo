using System;
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

	[Header("Atmosphere")]
	public float TropopauseElevation;
	public float BoundaryZoneElevation;
	public float stratosphereElevation;
	public float MaxTropopauseElevation;
	public float MinTropopauseElevation;
	public float TropopauseElevationSeason;


	[Header("Pressure and Wind")]
	//public float tradeWindSpeed = 12.0f; // average wind speeds around trade winds around 12 m/s
	//	public float pressureDifferentialWindSpeed = 70.0f; // hurricane wind speeds 70 m/s
	public float PressureToVerticalWindSpeed;
	public float DestinationPressureDifferentialToVerticalWindSpeed;
	public float MountainUpdraftWindSpeed;
	public float MaxTerrainNormalForFriction;
	public float AirMassDiffusionSpeedHorizontal;
	public float AirMassDiffusionSpeedVertical;
	public float WindOceanFriction;
	public float WindIceFriction;
	public float WindLandFriction;
	public float WindAirMovementHorizontal;
	public float WindAirMovementVertical;
	public float HumidityToCloudPercent;
	// http://tornado.sfsu.edu/geosciences/classes/e260/Coriolis_rdg/GeostrophicApproximation.html
	// states that geostrophic wind is only realistic at middle altitudes (500M), less so at 10K and at the surface, so we reduce overall coriolis effect by 25% to account
	public float GlobalCoriolisInfluenceWindUpper;
	public float GlobalCoriolisInfluenceWindLower;
	public float GlobalCoriolisInfluenceOcean;

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
	public float RainfallRate;
	public float EvapMinTemperature; // -30 celsius
	public float EvapMaxTemperature; // 70 celsius
	public float EvaporationRate; // TODO: evaporation on earth maxes out around 2.5M per year 
	public float RainDropFormationSpeedTemperature;
	public float RainDropCoalescenceWind;
	public float rainDropDragCoefficient;
	public float rainDropMaxSize;
	public float rainDropMinSize;
	public float airDensity;
	public float waterDensity;
	public float CloudDissapationRateWind;
	public float CloudDissapationRateDryAir;
	public float DewPointElevationPerDegree;
	public float DewPointTemperaturePerRelativeHumidity;

	[Header("Ocean")]
	public float ThermoclineDepth;
	public float WindToOceanCurrentFactor;
	public float WaterDiffuseSpeed;
	public float OceanCurrentSpeed;
	public float OceanHorizontalMixingSpeed;
	public float OceanUpwellingSpeed;
	public float OceanTemperatureVerticalMixingSpeed;
	public float SalinityVerticalMixingSpeed;
	public float OceanDensityPerSalinity;
	public float OceanDensityPerDegree;
	public float OceanDensityCurrentSpeed;

	[Header("Ecology")]
	public float MinTemperatureCanopy;
	public float MaxTemperatureCanopy;

	[Header("Physical Constants")]
	public const float TemperatureLapseRate = -0.0065f;
	public const float AdiabaticLapseRate = 0.0098f;
	public const float StaticPressure = 101325;
	public const float StdTemp = 288.15f;
	public const float MolarMassEarthAir = 0.0289644f;
	public const float MolarMassAir = 0.02857f;
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
	public const float EmissivityWater = 0.96f;
	public const float EmissivitySalt = 0.34f;
	public const float EmissivityIce = 0.96f;
	public const float EmissivityDirt = 0.92f;
	public const float EmissivitySand = 0.76f;
	public const float EmissivityAir = 0.8f;
	public const float EmissivityVegetation = 0.25f;
	public const float MassEarthAir = 1.29f;
	public const float MassWater = 1000f;
	public const float MassSalt = 2170f;
	public const float MassIce = 919f;
	public const float MassSoil = 1200f;
	public const float MassSand = 1600f;
	public const float ConductivityAir = 0.0262f;
	public const float ConductivityWater = 0.606f;
	public const float ConductivityIce = 2.18f;
	public const float ConductivityTerrain = 0.2f;
	public const float ConductivityAirWater = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityWater);
	public const float ConductivityAirIce = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityIce);
	public const float ConductivityAirTerrain = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityTerrain);
	public const float ConductivityIceWater = 1.0f / (1.0f / ConductivityWater + 1.0f / ConductivityIce);
	public const float ConductivityIceTerrain = 1.0f / (1.0f / ConductivityTerrain + 1.0f / ConductivityIce);
	public const float ConductivityWaterTerrain = 1.0f / (1.0f / ConductivityTerrain + 1.0f / ConductivityWater);



	[NonSerialized]	public float SpecificGasConstantDryAir;
	[NonSerialized]	public float DryAirAdiabaticLapseRate;
	[NonSerialized]	public float EvapTemperatureRange;
	[NonSerialized] public float TicksPerSecond;
	[NonSerialized] public float TicksPerYear;
	[NonSerialized] public float inverseFullCanopyCoverage;
	[NonSerialized] public float inverseFullWaterCoverage;
	[NonSerialized] public float inverseFullIceCoverage;
	[NonSerialized] public float inverseSpecificHeatIce;
	[NonSerialized] public float inverseCloudMassFullAbsorption;
	[NonSerialized] public float inverseBoundaryZoneElevation;
	[NonSerialized] public float wattsToKJPerTick;
	[NonSerialized] public float declinationOfSun;
	[NonSerialized] public float sunHitsAtmosphereBelowHorizonAmount;
	[NonSerialized] public float inverseSunAtmosphereAmount;
	[NonSerialized] public float inverseDewPointTemperatureRange;
	[NonSerialized] public float inverseEvapTemperatureRange;
	[NonSerialized] public float PressureExponent;

	public void Init()
	{
		EvapTemperatureRange = EvapMaxTemperature - EvapMinTemperature;
		SpecificGasConstantDryAir = UniversalGasConstant / MolarMassEarthAir;
		PressureExponent = 1.0f / (UniversalGasConstant * TemperatureLapseRate);

		DryAirAdiabaticLapseRate = AdiabaticLapseRate / SpecificHeatAtmosphere;
		TicksPerSecond = 1.0f / SecondsPerTick;
		TicksPerYear = 60 * 60 * 24 * 365 / SecondsPerTick;

		inverseFullCanopyCoverage = 1.0f / FullVegetationCoverage;
		inverseFullWaterCoverage = 1.0f / FullWaterCoverage;
		inverseFullIceCoverage = 1.0f / (MassIce * FullIceCoverage);
		inverseSpecificHeatIce = 1.0f / SpecificHeatIce;
		inverseCloudMassFullAbsorption = 1.0f / CloudMassFullAbsorption;
		inverseBoundaryZoneElevation = 1.0f / BoundaryZoneElevation;
		wattsToKJPerTick = SecondsPerTick * 1000;
		sunHitsAtmosphereBelowHorizonAmount = 0.055f;
		inverseSunAtmosphereAmount = 1.0f / (1.0f + sunHitsAtmosphereBelowHorizonAmount);
		inverseDewPointTemperatureRange = 1.0f / DewPointTemperatureRange;
		inverseEvapTemperatureRange = 1.0f / EvapTemperatureRange;

	}

}
