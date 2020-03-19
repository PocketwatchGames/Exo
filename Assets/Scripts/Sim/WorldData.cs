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
	public float TropopauseElevation;
	public float BoundaryZoneElevation;

	[Header("Solar Energy")]
	// atmospheric heat balance https://energyeducation.ca/encyclopedia/Earth%27s_heat_balance
	// https://en.wikipedia.org/wiki/Earth%27s_energy_budget
	// https://en.wikipedia.org/wiki/Electromagnetic_absorption_by_water
	// Water vapor is responsible for 70% of solar absorption and about 60% of absorption of thermal radiation.
	public float SolarAbsorptivityAir; // total absorbed by atmosphere AFTER reflection about 30%
	public float SolarAbsorptivityWaterVapor; // total absorbed by atmosphere AFTER reflection about 30%
	public float SolarAbsorptivityDust; // total absorbed by atmosphere AFTER reflection about 30%
	public float SolarAbsorptivityCloud; // 6% absorbed by clouds
	public float AlbedoAir; // 7% is reflected due to atmospheric scattering 
	public float AlbedoWaterVapor;
	public float AlbedoDust;
	public float AlbedoCloud;
	public float AlbedoReductionSoilQuality;
	public float minCloudFreezingTemperature;
	public float maxCloudFreezingTemperature;
	public float maxCloudSlopeAlbedo;
	public float rainDropSizeAlbedoMin;
	public float rainDropSizeAlbedoMax;

	[Header("Thermal Energy")]
	//public float EvaporativeHeatLoss = 0.6f; // global average = 78 watts
	// Net Back Radiation: The ocean transmits electromagnetic radiation into the atmosphere in proportion to the fourth power of the sea surface temperature(black-body radiation)
	// https://eesc.columbia.edu/courses/ees/climate/lectures/o_atm.html


	// TODO: tune these to match the science
	public float EnergyLostThroughAtmosphereWindow; // AKA Atmospheric window global average = 40 watts = 6.7% of all surface and atmospheric radiation

	// https://en.wikipedia.org/wiki/Electromagnetic_absorption_by_water
	// Water vapor is responsible for 70% of solar absorption and about 60% of absorption of thermal radiation.
	// carbon dioxide accounts for just 26% of the greenhouse effect.
	public float ThermalAbsorptivityAir;
	public float ThermalAbsorptivityWaterVapor;
	public float ThermalAbsorptivityDust;
	public float ThermalAbsorptivityCloud;

	// emissivity values obtained here: https://www.thermoworks.com/emissivity-table
	// and here https://www.aspen-electronics.com/uploads/3/7/1/2/37123419/emissivity-table.pdf
	public float ThermalEmissivityWater;
	public float ThermalEmissivitySalt;
	public float ThermalEmissivityIce;
	public float ThermalEmissivityAir;
	public float ThermalEmissivityCarbonDioxide;
	public float ThermalEmissivityWaterVapor;
	public float ThermalEmissivityDust;
	public float ThermalEmissivityDirt;
	public float ThermalEmissivitySand;
	public float ThermalEmissivityFlora;

	// TODO: should we parameterize the micro-conduction that allows for water to heat the air faster than it can cool it?
	//public float AirWaterConductionPositive; // global avg = 16 watts per degree delta between air and ocean (global avg = 24 watts per m^2 of ocean)

	[Header("Evaporation")] // evaporation on earth maxes out around 2.5M per year
	public float WaterHeatingDepth;

	[Header("Ice")]
	public float IceHeatingDepth;
	public float FreezePointReductionPerSalinity;

	[Header("Rain and Clouds")]
	public float rainDropDragCoefficient;
	public float rainDropMaxSize;
	public float rainDropMinSize;
	public float CloudDissapationRateWind;
	public float CloudDissapationRateDryAir;

	[Header("Diffusion")]
	public float CloudDiffusionCoefficient;
	public float AirDiffusionCoefficientHorizontal;
	public float AirDiffusionCoefficientVertical;
	public float WaterDiffusionCoefficientHorizontal;
	public float WaterDiffusionCoefficientVertical;

	[Header("Wind")]
	public float WindWaterFriction;
	public float WindIceFriction;
	public float WindTerrainFrictionMin;
	public float WindTerrainFrictionMax;
	public float WindFloraFriction;
	public float MaxTerrainRoughnessForWindFriction;
	public float WaterSurfaceFrictionDepth;

	[Header("Water Current")]
	public float WindToWaterCurrentFrictionCoefficient;
	public float WaterDensityPerSalinity;
	public float WaterDensityPerDegree;
	public float WaterDensityCurrentSpeed;

	[Header("Terrain and Ground Water")]
	public float SoilHeatDepth;
	public float GroundWaterFlowSpeed;
	public float GroundWaterMax;
	public float GroundWaterMaxDepth;
	public float GroundWaterAbsorptionRate;
	public float GroundWaterDiffusionCoefficient;
	public float FullCoverageIce;
	public float FullCoverageWater;
	public float FullCoverageFlora;

	[Header("Flora")]
	public float FloraMax;
	public float FloraGrowthRate;
	public float FloraDeathRateWater;
	public float FloraDeathRateCrowding;
	public float FloraDeathRateTemperature;
	public float FloraDeathRateAge;
	public float FloraEvaporationRate;
	public float FloraWaterConsumptionRate;
	public float FloraGrowthTemperatureRangeInverse;
	public float FloraAirSurfaceArea;

	[Header("Lava")]
	public float LavaSolidificationTemperature;
	public float CrustDepthForEruption;


	#region Constants
	public const float TemperatureLapseRate = -0.0065f;
	public const float TemperatureLapseRateInverse = 1.0f / TemperatureLapseRate;
	public const float AdiabaticLapseRate = 0.0098f;
	public const float StandardPressure = 101325;
	public const float StandardTemperature = 288.15f;
	public const float MolarMassAir = 0.0289647f;
	public const float MolarMassWater = 0.01802f;
	public const float MolarMassAirInverse = 1.0f / MolarMassAir;
	public const float UniversalGasConstant = 8.3144598f;
	public const float FreezingTemperature = 273.15f;
	public const float StefanBoltzmannConstant = 0.00000005670373f;
	// specific heat is joules to raise one degree (kJ/kgK)
	public const float AlbedoWater = 0.06f; // How much heat is reflected back by the water
	public const float AlbedoIce = 0.5f; // How much heat is reflected back by the water
	public const float AlbedoLand = 0.4f;
	public const float AlbedoFlora = 0.1f;
	//public const float AlbedoCloud = 0.05f; // 24% incoming  reflected back to space by clouds (avg, globally)
	public const float SpecificHeatIce = 2.108f;
	public const float SpecificHeatFlora = 1.76f;
	public const float SpecificHeatWater = 4.187f;
	public const float SpecificHeatWaterVapor = 1.996f;
	public const float SpecificHeatSalt = 0.85f;
	public const float SpecificHeatAtmosphere = 1.158f;
	public const float SpecificHeatSoil = 0.84f;
	public const float LatentHeatWaterLiquid = 334.0f;
	public const float LatentHeatWaterVapor = 2264.705f;
	public const float MassEarthAir = 1.29f;
	public const float MassWater = 1000f;
	public const float MassSalt = 2170f;
	public const float MassIce = 919f;
	public const float MassSoil = 1200f;
	public const float MassSand = 1600f;
	public const float MassLava = 3100f;
	public const float DensityWater = 997f;
	public const float DensityAir = 1.21f;
	public const float ConductivityAir = 0.0262f;
	public const float ConductivityWater = 0.606f;
	public const float ConductivityIce = 2.18f;
	public const float ConductivityTerrain = 0.2f;
	public const float ConductivityFlora = 0.1f;
	public const float ThermalContactResistance = 0.00005f;
	public const float ConductivityAirWater = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityWater + ThermalContactResistance);
	public const float ConductivityAirIce = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityIce + ThermalContactResistance);
	public const float ConductivityAirFlora = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityFlora + ThermalContactResistance);
	public const float ConductivityAirTerrain = 1.0f / (1.0f / ConductivityAir + 1.0f / ConductivityTerrain + ThermalContactResistance);
	public const float ConductivityIceWater = 1.0f / (1.0f / ConductivityWater + 1.0f / ConductivityIce + ThermalContactResistance);
	public const float ConductivityIceFlora = 1.0f / (1.0f / ConductivityFlora + 1.0f / ConductivityIce + ThermalContactResistance);
	public const float ConductivityIceTerrain = 1.0f / (1.0f / ConductivityTerrain + 1.0f / ConductivityIce + ThermalContactResistance);
	public const float ConductivityWaterFlora = 1.0f / (1.0f / ConductivityFlora + 1.0f / ConductivityWater + ThermalContactResistance);
	public const float ConductivityWaterTerrain = 1.0f / (1.0f / ConductivityTerrain + 1.0f / ConductivityWater + ThermalContactResistance);
	public const float ConductivityFloraTerrain = 1.0f / (1.0f / ConductivityTerrain + 1.0f / ConductivityFlora + ThermalContactResistance);
	public const float GasConstantAir = UniversalGasConstant / MolarMassAir * 1000;
	public const float GasConstantWaterVapor = UniversalGasConstant / MolarMassWater * 1000;
	public const float PressureExponent = 1.0f / (UniversalGasConstant * TemperatureLapseRate);
	public const float DryAirAdiabaticLapseRate = AdiabaticLapseRate / SpecificHeatAtmosphere;
	public const float inverseSpecificHeatIce = 1.0f / SpecificHeatIce;
	public const float InverseDensityAir = 1.0f / DensityAir;
	#endregion

	#region Nonserialized
	[NonSerialized] public float TicksPerSecond;
	[NonSerialized] public float TicksPerYear;
	[NonSerialized] public float inverseFullCoverageFlora;
	[NonSerialized] public float inverseFullCoverageWater;
	[NonSerialized] public float inverseFullCoverageIce;
	[NonSerialized] public float wattsToKJPerTick;
	[NonSerialized] public float declinationOfSun;
	[NonSerialized] public float sunHitsAtmosphereBelowHorizonAmount;
	[NonSerialized] public float inverseSunAtmosphereAmount;
	#endregion

	public void Init()
	{

		TicksPerSecond = 1.0f / SecondsPerTick;
		TicksPerYear = 60 * 60 * 24 * 365 / SecondsPerTick;

		inverseFullCoverageFlora = 1.0f / FullCoverageFlora;
		inverseFullCoverageWater = 1.0f / FullCoverageWater;
		inverseFullCoverageIce = 1.0f / FullCoverageIce;
		wattsToKJPerTick = SecondsPerTick * 1000;
		sunHitsAtmosphereBelowHorizonAmount = 0.055f;
		inverseSunAtmosphereAmount = 1.0f / (1.0f + sunHitsAtmosphereBelowHorizonAmount);

	}

}
