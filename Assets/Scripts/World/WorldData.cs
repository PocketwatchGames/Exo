using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

[Serializable]
public class WorldData {
	public float WaterDiffuseSpeed = 0.5f;
	public float FullIceCoverage = 3.0f;
	public float FullWaterCoverage = 50.0f;

	[Header("Fresh Water")]
	public float MaxSoilPorousness = 0.1f;
	public float GroundWaterReplenishmentSpeed = 10.0f;
	public float GroundWaterFlowSpeed = 0.5f;

	[Header("Atmosphere")]
	public float TropopauseElevation = 10000;
	public float BoundaryZoneElevation = 1000;
	public float stratosphereElevation = 50000;
	public float MaxTropopauseElevation = 17000f;
	public float MinTropopauseElevation = 9000f;
	public float TropopauseElevationSeason = 1000f;


	[Header("Pressure and Wind")]
	//public float tradeWindSpeed = 12.0f; // average wind speeds around trade winds around 12 m/s
	//	public float pressureDifferentialWindSpeed = 70.0f; // hurricane wind speeds 70 m/s
	public float PressureToVerticalWindSpeed = 0.0001f;
	public float DestinationPressureDifferentialToVerticalWindSpeed = 0.0001f;
	public float MountainUpdraftWindSpeed = 0.1f;
	public float MaxTerrainNormalForFriction = 0.25f;
	public float AirMassDiffusionSpeedHorizontal = 0.1f;
	public float AirMassDiffusionSpeedVertical = 0.1f;
	public float AirEnergyDiffusionSpeedHorizontal = 0.1f;
	public float AirEnergyDiffusionSpeedVertical = 0.1f;
	public float WindOceanFriction = 0.2f;
	public float WindIceFriction = 0.1f;
	public float WindLandFriction = 0.5f;
	public float WindAirMovementHorizontal = 0.001f;
	public float WindAirMovementVertical = 0.01f;
	public float WindHumidityMovement = 0.1f;
	public float HumidityToCloudPercent = 0.1f;
	public float PressureGradientWindMultiplier = 4000;
	// http://tornado.sfsu.edu/geosciences/classes/e260/Coriolis_rdg/GeostrophicApproximation.html
	// states that geostrophic wind is only realistic at middle altitudes (500M), less so at 10K and at the surface, so we reduce overall coriolis effect by 25% to account
	public float GlobalCoriolisInfluenceWindUpper = 0.75f;
	public float GlobalCoriolisInfluenceWindLower = 0.25f;
	public float GlobalCoriolisInfluenceOcean = 0.5f;

	[Header("Atmospheric Energy Cycle")]
	// atmospheric heat balance https://energyeducation.ca/encyclopedia/Earth%27s_heat_balance
	// https://en.wikipedia.org/wiki/Earth%27s_energy_budget
	public float AtmosphericHeatAbsorption = 0.297f; // total absorbed by atmosphere AFTER reflection about 30%
	public float AtmosphericHeatReflection = 0.07f; // 7% is reflected due to atmospheric scattering 
	public float CloudSolarAbsorptivity = 0.06f; // 6% absorbed by clouds
	public float CloudOutgoingReflectionRate = 0.90f;
	//public float EvaporativeHeatLoss = 0.6f; // global average = 78 watts
	// Net Back Radiation: The ocean transmits electromagnetic radiation into the atmosphere in proportion to the fourth power of the sea surface temperature(black-body radiation)
	// https://eesc.columbia.edu/courses/ees/climate/lectures/o_atm.html
	public float OceanAirConductionWarming = 0.016f; // global avg = 16 watts per degree delta between air and ocean (global avg = 24 watts per m^2 of ocean)
	public float OceanAirConductionCooling = 0.008f; // 
	public float WaterAirConductionDepth = 10.0f;
	public float OceanIceConduction = 0.00001f; // small
	public float IceAirConductionCooling = 0.008f;
	public float maxGroundWaterTemperature = 283;
	public float SoilHeatDepth = 2;

	// TODO: tune these to match the science
	public float CloudMassFullAbsorption = 0.5f; // how much heat gain/loss is caused by cloud cover (cumulus cloud is 0.3g/cubic meter, and about 3 kilometers high)
	public float EnergyEmittedByAtmosphere = 0.000000199f; // how fast a cell loses heat an min elevation, no cloud cover, global average = 199 watts
	public float EnergyLostThroughAtmosphereWindow = 0.067f; // AKA Atmospheric window global average = 40 watts = 6.7% of all surface and atmospheric radiation

	[Header("Evap, Humidity and Clouds")]
	public float DewPointTemperatureRange = 100.0f;
	public float DewPointZero = 213.0f;
	public float WaterVaporMassToAirMassAtDewPoint = 0.2f;
	public float RainfallRate = 0.01f;
	public float EvapMinTemperature = 243; // -30 celsius
	public float EvapMaxTemperature = 343; // 70 celsius
	public float EvaporationRate = 0.002f; // TODO: evaporation on earth maxes out around 2.5M per year 
	public float RainDropFormationSpeedTemperature = 10f;
	public float RainDropCoalescenceWind = 0.00001f;
	public float rainDropDragCoefficient = 0.5f;
	public float rainDropMaxSize = 0.004f;
	public float rainDropMinSize = 0.00001f;
	public float airDensity = 1.21f;
	public float waterDensity = 997;
	public float CloudDissapationRateWind = 0.1f;
	public float CloudDissapationRateDryAir = 0.1f;

	[Header("Ocean")]
	public float DeepOceanDepth = 500;
	public float WindToOceanCurrentFactor = 0.1f;
	public float OceanCurrentSpeed = 0.00001f;
	public float OceanHorizontalMixingSpeed = 0.01f;
	public float OceanUpwellingSpeed = 0.0001f;
	public float OceanTemperatureVerticalMixingSpeed = 0.0001f;
	public float SalinityVerticalMixingSpeed = 0.001f;
	public float OceanDensityPerSalinity = 800f;
	public float OceanDensityPerDegree = 0.2f;
	public float OceanDensityCurrentSpeed = 0.1f;

	[Header("Ecology")]
	public float MinTemperatureCanopy = 223.15f;
	public float MaxTemperatureCanopy = 323.15f;

	[Header("Physical Constants")]
	public float TemperatureLapseRate = -0.0065f;
	public float AdiabaticLapseRate = 0.0098f;
	public float GravitationalAcceleration = 9.80665f;
	public float StaticPressure = 101325;
	public float StdTemp = 288.15f;
	public float MolarMassEarthAir = 0.0289644f;
	public float MolarMassAir = 0.02857f;
	public float MolarMassWater = 0.01802f;
	public float UniversalGasConstant = 8.3144598f;
	public float FreezingTemperature = 273.15f;
	public float DewPointElevationPerDegree = 67.73f;
	public float DewPointTemperaturePerRelativeHumidity = 20;
	public float StefanBoltzmannConstant = 0.00000005670373f;

	[Header("Substance Properties")]
	// specific heat is joules to raise one degree (kJ/kgK)
	public float AlbedoWater = 0.06f; // How much heat is reflected back by the water
	public float AlbedoIce = 0.5f; // How much heat is reflected back by the water
	public float AlbedoLand = 0.4f;
	public float AlbedoReductionSoilQuality = 0.15f;
	public float AlbedoFoliage = 0.1f;
	//public float AlbedoCloud = 0.05f; // 24% incoming  reflected back to space by clouds (avg, globally)
	public float SpecificHeatIce = 2.108f;
	public float SpecificHeatWater = 4.187f;
	public float SpecificHeatWaterVapor = 1.996f;
	public float SpecificHeatSalt = 0.85f;
	public float SpecificHeatAtmosphere = 1.158f;
	public float SpecificHeatSoil = 0.84f;
	public float LatentHeatWaterLiquid = 334.0f;
	public float LatentHeatWaterVapor = 2264.705f;
	public float EmissivityWater = 0.96f;
	public float EmissivityIce = 0.96f;
	public float EmissivityDirt = 0.92f;
	public float EmissivitySand = 0.76f;
	public float EmissivityAir = 0.8f;
	public float MassEarthAir = 1.29f;
	public float MassWater = 1000f;
	public float MassSalt = 2170f;
	public float MassIce = 919f;
	public float MassSoil = 1200f;
	public float MassSand = 1600f;
	public float AbsorptivityWaterLiquid = 0.05f;
	public float AbsorptivityCarbonDioxide = 0.25f;
	public float AbsorptivityWaterVapor = 0.00155f;
	public float AbsorptivityAir = 0.0005f;


	[NonSerialized]
	public float PressureExponent;
	[NonSerialized]
	public float SpecificGasConstantDryAir;
	[NonSerialized]
	public float DryAirAdiabaticLapseRate;
	[NonSerialized]
	public float EvapTemperatureRange;


	public void Init(int size)
	{

		EvapTemperatureRange = EvapMaxTemperature - EvapMinTemperature;
		SpecificGasConstantDryAir = UniversalGasConstant / MolarMassEarthAir;
		PressureExponent = GravitationalAcceleration / (UniversalGasConstant * TemperatureLapseRate);

		DryAirAdiabaticLapseRate = AdiabaticLapseRate / SpecificHeatAtmosphere;

		for (int y = 0; y < size; y++)
		{
			float latitude = ((float)y / size) * 2 - 1.0f;
			float yaw = (float)(latitude * Math.PI * 1.5f);
			float pitch = (float)(latitude * Math.PI * 3f);
			float absSinPitch = (float)(Math.Abs(Math.Sin(pitch)));
			float cosYaw = (float)Math.Cos(yaw);
			float cosPitch = (float)Math.Cos(pitch);

			float tropopauseElevation = (1.0f - Math.Abs(latitude)) * (MaxTropopauseElevation - MinTropopauseElevation) + MinTropopauseElevation + TropopauseElevationSeason * latitude;
		}

	}

}
