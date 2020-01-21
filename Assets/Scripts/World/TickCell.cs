using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

[BurstCompile]
public struct TickCellJob : IJobParallelFor {

	public NativeArray<SimCell> Cells;
	public NativeArray<DisplayCell> DisplayCells;

	[ReadOnly] public SimPlanetState PlanetState;
	[ReadOnly] public NativeArray<CellDiffusion> Diffusion;
	[ReadOnly] public NativeArray<CellDiffusion> DiffusionLimit;
	[ReadOnly] public NativeArray<SimCell> Last;
	[ReadOnly] public NativeArray<SimWind> LastWind;
	[ReadOnly] public WorldData worldData;
	[ReadOnly] public StaticState staticState;

	public void Execute(int i)
	{
		var last = Last[i];
		var lastWind = LastWind[i];
		var next = new SimCell();
		var display = new DisplayCell();

		next.CloudCoverage = last.CloudCoverage;
		next.RelativeHumidity = last.RelativeHumidity;
		next.WaterDepth = last.WaterDepth;
		next.CloudElevation = last.CloudElevation;
		next.Elevation = last.Elevation;
		next.Roughness = last.Roughness;
		next.IceMass = last.IceMass;
		next.Vegetation = last.Vegetation;
		next.SoilFertility = last.SoilFertility;
		next.WaterMass = last.WaterMass;
		next.WaterEnergy = last.WaterEnergy;
		next.SaltMass = last.SaltMass;
		next.GroundEnergy = last.GroundEnergy;
		next.GroundWater = last.GroundWater;
		next.GroundWaterDepth = last.GroundWaterDepth;
		next.AirMass = last.AirMass;
		next.AirEnergy = last.AirEnergy;
		next.CloudMass = last.CloudMass;
		next.CloudDropletMass = last.CloudDropletMass;
		next.AirWaterMass = last.AirWaterMass;
		next.AirTemperature = last.AirTemperature;
		next.AirPressure = last.AirPressure;
		next.WaterTemperature = last.WaterTemperature;

		float iceCoverage = math.min(1.0f, math.pow(last.IceMass * worldData.inverseFullIceCoverage, 0.6667f));
		float surfaceElevation = last.Elevation + last.WaterAndIceDepth;
		float evaporationRate = Atmosphere.GetEvaporationRate(ref worldData, iceCoverage, last.AirTemperature, last.RelativeHumidity, worldData.inverseEvapTemperatureRange);
		float dewPoint = Atmosphere.GetDewPoint(ref worldData, last.AirTemperature, last.RelativeHumidity);

		DoEnergyCycle(i, ref next, ref display, surfaceElevation, evaporationRate, dewPoint, iceCoverage);
		DoVerticalWaterMovement(i, ref last, ref lastWind, ref next, ref display, surfaceElevation, dewPoint, evaporationRate);
		DoDiffusion(i, ref next);



		{

			//float shallowOceanVolume = GetWaterVolume(world, shallowOceanMass, oceanSalinityShallow, oceanEnergyShallow);
			//float moveFromShallow = shallowOceanVolume - world.Data.DeepOceanDepth;
			//if (moveFromShallow < 0)
			//{
			//	if (deepOceanMass > 0)
			//	{
			//		float deepOceanVolume = GetWaterVolume(world, deepOceanMass, oceanSalinityDeep, oceanEnergyDeep);
			//		moveFromShallow = Mathf.Max(moveFromShallow, -deepOceanVolume);
			//		float movePercent = moveFromShallow / deepOceanVolume;
			//		nextState.ShallowWaterMass[index] -= deepOceanMass * movePercent;
			//		nextState.ShallowWaterEnergy[index] -= oceanEnergyDeep * movePercent;
			//		nextState.ShallowSaltMass[index] -= oceanSalinityDeep * movePercent;
			//		nextState.DeepWaterMass[index] += deepOceanMass * movePercent;
			//		nextState.DeepWaterEnergy[index] += oceanEnergyDeep * movePercent;
			//		nextState.DeepSaltMass[index] += oceanSalinityDeep * movePercent;
			//	}
			//}
			//else if (shallowOceanVolume > 0)
			//{
			//	float movePercent = moveFromShallow / shallowOceanVolume;
			//	nextState.ShallowWaterMass[index] -= shallowOceanMass * movePercent;
			//	nextState.ShallowWaterEnergy[index] -= oceanEnergyShallow * movePercent;
			//	nextState.ShallowSaltMass[index] -= oceanSalinityShallow * movePercent;
			//	nextState.DeepWaterMass[index] += shallowOceanMass * movePercent;
			//	nextState.DeepWaterEnergy[index] += oceanEnergyShallow * movePercent;
			//	nextState.DeepSaltMass[index] += oceanSalinityShallow * movePercent;
			//}


			float waterTemperature = Atmosphere.GetWaterTemperature(next.WaterEnergy, next.WaterMass, next.SaltMass);
			//float deepWaterTemperature = GetWaterTemperature(world, nextState.DeepWaterEnergy[index], nextState.DeepWaterMass[index], nextState.DeepSaltMass[index]);
			next.WaterTemperature = waterTemperature;
			float waterDepth = math.max(0, Atmosphere.GetWaterVolumeByTemperature(ref worldData, next.WaterMass, next.SaltMass, waterTemperature));
			next.WaterDepth = waterDepth;
			next.WaterAndIceDepth = waterDepth + next.IceMass / WorldData.MassIce;
			next.WaterDensity = Atmosphere.GetWaterDensity(ref worldData, next.WaterEnergy, next.SaltMass, next.WaterMass);

		}
		float newSurfaceElevation = next.Elevation + next.WaterAndIceDepth;
		next.AirTemperature = Atmosphere.GetAirTemperature(next.AirEnergy, next.AirMass, next.CloudMass, next.AirWaterMass);
		next.AirPressure = math.max(0, Atmosphere.GetAirPressure(
			ref worldData, 
			next.AirMass + PlanetState.StratosphereMass + next.AirWaterMass + next.CloudMass, 
			surfaceElevation, 
			next.AirTemperature, 
			Atmosphere.GetMolarMassAir(next.AirMass + PlanetState.StratosphereMass, next.CloudMass + next.AirWaterMass),
			PlanetState.Gravity));
		next.RelativeHumidity = Atmosphere.GetRelativeHumidity(ref worldData, next.AirTemperature, next.AirWaterMass, next.AirMass, worldData.inverseDewPointTemperatureRange);
		//next.UpperAirTemperature[index] = GetAirTemperature(world, nextState.UpperAirEnergy[index], nextState.UpperAirMass[index], nextState.CloudMass[index], world.Data.SpecificHeatWater);
		//next.UpperAirPressure[index] = Mathf.Max(0, GetAirPressure(world, nextState.UpperAirMass[index] + nextState.StratosphereMass + nextState.CloudMass[index], elevationOrSeaLevel + world.Data.BoundaryZoneElevation, nextState.UpperAirTemperature[index], GetMolarMassAir(world, nextState.LowerAirMass[index] + nextState.UpperAirMass[index] + nextState.StratosphereMass, nextState.Humidity[index] + nextState.CloudMass[index])));

		float newDewPoint = Atmosphere.GetDewPoint(ref worldData, next.AirTemperature, next.RelativeHumidity);
		next.CloudElevation = Atmosphere.GetCloudElevation(ref worldData, next.AirTemperature, newDewPoint, newSurfaceElevation);

		if (math.isnan(next.CloudElevation) || math.isinf(next.CloudElevation) ||next.CloudElevation < 0 || next.CloudElevation > 10000)
		{
			surfaceElevation = 0;
		}

		Cells[i] = next;
		DisplayCells[i] = display;

	}

	private void DoEnergyCycle(int i, ref SimCell next, ref DisplayCell displayCell, float surfaceElevation, float evaporationRate, float dewPoint, float iceCoverage)
	{
		var last = Last[i];

		float waterCoverage = math.min(1.0f, math.pow(last.WaterDepth * worldData.inverseFullWaterCoverage, 0.6667f));
		float canopyCoverage = math.min(1.0f, math.pow(last.Vegetation * worldData.inverseFullCanopyCoverage, 0.6667f));
		float cloudCoverage = math.min(1.0f, math.pow(last.CloudMass * worldData.inverseCloudMassFullAbsorption, 0.6667f)); // bottom surface of volume

		float sunDotSurface = math.dot(math.normalize(PlanetState.Position), math.rotate(UnityEngine.Quaternion.Euler(math.degrees(PlanetState.Rotation)), -staticState.SphericalPosition[i]));

		float waterSlopeAlbedo = math.pow(1.0f - math.max(0, sunDotSurface), 9);
		//float groundWaterSaturation = Animals.GetGroundWaterSaturation(state.GroundWater[index], state.WaterTableDepth[index], soilFertility * world.Data.MaxSoilPorousness);


		// SOLAR RADIATION

		float solarRadiationAbsorbed = 0;
		float solarRadiation = PlanetState.SolarRadiation * worldData.SecondsPerTick * math.max(0, (sunDotSurface + worldData.sunHitsAtmosphereBelowHorizonAmount) * worldData.inverseSunAtmosphereAmount);

		// get the actual atmospheric depth here based on radius of earth plus atmosphere
		//float inverseSunAngle = PIOver2 + sunAngle;
		//float angleFromSunToLatitudeAndAtmophereEdge = math.Asin(state.PlanetRadius * math.Sin(inverseSunAngle) / (state.PlanetRadius + world.Data.TropopauseElevation));
		//float angleFromPlanetCenterToLatitudeAndAtmosphereEdge = math.PI - inverseSunAngle - angleFromSunToLatitudeAndAtmophereEdge;
		//float atmosphericDepthInMeters = math.Sin(angleFromPlanetCenterToLatitudeAndAtmosphereEdge) * state.PlanetRadius / math.Sin(angleFromSunToLatitudeAndAtmophereEdge);
		//float atmosphericDepth = math.max(1.0f, atmosphericDepthInMeters / world.Data.TropopauseElevation);

		float atmosphericDepth = 1.0f + (1.0f - sunDotSurface);

		// These constants obtained here, dunno if I've interpreted them correctly
		// https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass

		////// MAJOR TODO:
		///// USE THIS LINK: https://www.ftexploring.com/solar-energy/sun-angle-and-insolation2.htm
		/// With the sun 90 degrees above the horizon (SEA° = 90°), the air mass lowers the intensity of the sunlight from the 1,367 W / m2 that it is in outerspace down to about 1040 W / m2.
		//				float consumedByAtmosphere = 1.0f - math.pow(0.7f, math.pow(atmosphericDepth, 0.678f));

		if (solarRadiation > 0)
		{
			displayCell.EnergyIncoming += solarRadiation;

			// TODO: reflect/absorb more in the atmosphere with a lower sun angle

			// reflect some rads off atmosphere and clouds
			// TODO: this process feels a little broken -- are we giving too much priority to reflecting/absorbing in certain layers?
			float energyReflectedAtmosphere = solarRadiation * math.min(1, worldData.AtmosphericHeatReflection * (last.AirMass + last.AirWaterMass));
			solarRadiation -= energyReflectedAtmosphere;
			displayCell.EnergySolarReflectedAtmosphere += energyReflectedAtmosphere;

			if (last.CloudMass > 0)
			{
				float cloudTemperatureAlbedo = WorldData.AlbedoIce + (WorldData.AlbedoWater - WorldData.AlbedoIce) * math.saturate((dewPoint - worldData.minCloudFreezingTemperature) / (worldData.maxCloudFreezingTemperature - worldData.minCloudFreezingTemperature));
				float rainDropSizeAlbedo = math.clamp(1.0f - last.CloudDropletMass / last.CloudMass, 0,1) * (worldData.rainDropSizeAlbedoMax - worldData.rainDropSizeAlbedoMin) + worldData.rainDropSizeAlbedoMin;
				float cloudReflectivity = math.min(1.0f, worldData.cloudAlbedo * cloudTemperatureAlbedo * last.CloudMass * rainDropSizeAlbedo / math.max(worldData.maxCloudSlopeAlbedo, 1.0f - waterSlopeAlbedo));
				float energyReflectedClouds = solarRadiation * cloudReflectivity;
				solarRadiation -= energyReflectedClouds;

				float absorbedByCloudsIncoming = solarRadiation * math.min(1.0f, worldData.AtmosphericHeatAbsorption * last.CloudMass);
				solarRadiation -= absorbedByCloudsIncoming;
				next.AirEnergy += absorbedByCloudsIncoming;
				displayCell.EnergySolarAbsorbedCloud += absorbedByCloudsIncoming;
				displayCell.EnergySolarAbsorbedAtmosphere += absorbedByCloudsIncoming;
				displayCell.EnergySolarReflectedCloud += energyReflectedClouds;
			}

			// Absorbed by atmosphere
			// stratosphere accounts for about a quarter of atmospheric mass
			//	float absorbedByStratosphere = incomingRadiation * world.Data.AtmosphericHeatAbsorption * (state.StratosphereMass / massOfAtmosphericColumn);

			float atmosphereAbsorptionRate = math.min(1, worldData.AtmosphericHeatAbsorption * (last.AirMass + last.AirWaterMass));
			float absorbedByAtmosphereIncoming = solarRadiation * atmosphereAbsorptionRate * atmosphericDepth;

			next.AirEnergy += absorbedByAtmosphereIncoming;
			solarRadiation -= absorbedByAtmosphereIncoming;
			displayCell.EnergySolarAbsorbedAtmosphere += absorbedByAtmosphereIncoming;

			// reflection off surface
			float energyReflected = 0;
			{
				if (iceCoverage > 0)
				{
					energyReflected += solarRadiation * iceCoverage * Atmosphere.GetAlbedo(WorldData.AlbedoIce, 0);
				}
				if (waterCoverage > 0)
				{
					energyReflected += waterCoverage * solarRadiation * Atmosphere.GetAlbedo(WorldData.AlbedoWater, waterSlopeAlbedo) * (1.0f - iceCoverage);
				}
				if (waterCoverage < 1 && iceCoverage < 1)
				{
					// reflect some incoming radiation
					float slopeAlbedo = 0;
					float soilReflectivity = Atmosphere.GetAlbedo(WorldData.AlbedoLand - worldData.AlbedoReductionSoilQuality * last.SoilFertility, slopeAlbedo);
					float heatReflectedLand = canopyCoverage * WorldData.AlbedoFoliage + math.max(0, 1.0f - canopyCoverage) * soilReflectivity;
					energyReflected += solarRadiation * math.saturate(heatReflectedLand) * (1.0f - iceCoverage) * (1.0f - waterCoverage);
				}
				solarRadiation -= energyReflected;

				// TODO: do we absorb some of this energy on the way back out of the atmosphere?
				displayCell.EnergySolarReflectedSurface += energyReflected;
			}

			solarRadiationAbsorbed += solarRadiation;

		}

		// THERMAL RADIATION

		float evaporation = 0;
		float backRadiation = 0;
		float reflected = 0;

		// radiate heat from land
		// TODO: deal with the fact that this also incorporates ground water
		float shallowWaterRadiation = Atmosphere.GetRadiationRate(last.WaterTemperature, WorldData.EmissivityWater) * worldData.SecondsPerTick * waterCoverage;
		float soilEnergy = last.GroundEnergy - last.GroundWater * worldData.maxGroundWaterTemperature * WorldData.SpecificHeatWater;
		float radiationRate = Atmosphere.GetLandRadiationRate(ref worldData, last.GroundEnergy, last.GroundWater, last.SoilFertility, canopyCoverage);
		float thermalEnergyRadiatedLand = math.min(soilEnergy, radiationRate * worldData.SecondsPerTick);
		next.GroundEnergy += PlanetState.GeothermalHeat * worldData.SecondsPerTick - thermalEnergyRadiatedLand;

		float thermalEnergyRadiatedToIce = 0;
		float thermalEnergyRadiatedToShallowWater = 0;
		float thermalEnergyRadiatedToAir = 0;

		next.WaterEnergy -= shallowWaterRadiation;
		next.GroundEnergy += shallowWaterRadiation;

		thermalEnergyRadiatedToShallowWater += thermalEnergyRadiatedLand * waterCoverage;
		thermalEnergyRadiatedToIce += iceCoverage * (thermalEnergyRadiatedLand - thermalEnergyRadiatedToShallowWater);
		thermalEnergyRadiatedToAir += thermalEnergyRadiatedLand - thermalEnergyRadiatedToShallowWater - thermalEnergyRadiatedToIce;
		next.WaterEnergy += thermalEnergyRadiatedToShallowWater;

		// radiate from ice to air above
		if (iceCoverage > 0)
		{
			float thermalEnergyRadiatedFromIce = iceCoverage * Atmosphere.GetRadiationRate(WorldData.FreezingTemperature, WorldData.EmissivityIce) * worldData.SecondsPerTick;
			thermalEnergyRadiatedToAir += thermalEnergyRadiatedFromIce;
			thermalEnergyRadiatedToIce -= thermalEnergyRadiatedFromIce;
		}

		// lose heat to air via conduction AND radiation
		if (iceCoverage < 1)
		{
			// radiate heat, will be absorbed by air
			// Net Back Radiation: The ocean transmits electromagnetic radiation into the atmosphere in proportion to the fourth power of the sea surface temperature(black-body radiation)
			// https://eesc.columbia.edu/courses/ees/climate/lectures/o_atm.html
			float oceanRadiation = math.min(last.WaterEnergy, shallowWaterRadiation);
			next.WaterEnergy -= oceanRadiation;
			thermalEnergyRadiatedToIce += oceanRadiation * iceCoverage;
			thermalEnergyRadiatedToAir += oceanRadiation * (1.0f - iceCoverage);
			displayCell.EnergyThermalOceanRadiation += oceanRadiation;
		}

		// TODO: track and emit heat from ice

		float atmosphereEmissivity = Atmosphere.GetAtmosphericEmissivity(last.AirMass, last.AirMass * PlanetState.CarbonDioxide, last.AirWaterMass, last.CloudMass);
		float emittedByAtmosphere = 0;
		float surfaceEnergyReflected = 0;

		// Thermal energy from surface to air, space, reflected off clouds
		{
			displayCell.EnergyThermalSurfaceRadiation += thermalEnergyRadiatedToAir;
			float energyThroughAtmosphericWindow = thermalEnergyRadiatedToAir * worldData.EnergyLostThroughAtmosphereWindow;
			thermalEnergyRadiatedToAir -= energyThroughAtmosphericWindow;
			displayCell.EnergyThermalOutAtmosphericWindow += energyThroughAtmosphericWindow;

			float absorbed = thermalEnergyRadiatedToAir * atmosphereEmissivity;
			thermalEnergyRadiatedToAir -= absorbed;
			next.AirEnergy += absorbed;
			displayCell.EnergyThermalAbsorbedAtmosphere += absorbed;

			surfaceEnergyReflected = thermalEnergyRadiatedToAir * cloudCoverage * worldData.CloudOutgoingReflectionRate;
			reflected += surfaceEnergyReflected;
			displayCell.EnergyThermalOutAtmosphere += thermalEnergyRadiatedToAir - surfaceEnergyReflected;
		}

		// atmosphere radiation
		{
			float energyEmitted = Atmosphere.GetRadiationRate(last.AirTemperature, atmosphereEmissivity) * worldData.SecondsPerTick;
			next.AirEnergy -= 2 * energyEmitted;
			emittedByAtmosphere += 2 * energyEmitted;

			backRadiation += energyEmitted;

			float energyThroughAtmosphericWindow = energyEmitted * worldData.EnergyLostThroughAtmosphereWindow;
			energyEmitted -= energyThroughAtmosphericWindow;
			displayCell.EnergyThermalOutAtmosphericWindow += energyThroughAtmosphericWindow;

			float energyReflected = energyEmitted * cloudCoverage * worldData.CloudOutgoingReflectionRate;
			reflected += energyReflected;
			displayCell.EnergyThermalOutAtmosphere += energyEmitted - energyReflected;
		}

		// reflected thermal radiation
		{
			float absorbed = reflected * atmosphereEmissivity;
			next.AirEnergy += absorbed;
			reflected -= absorbed;

			displayCell.EnergyThermalAbsorbedAtmosphere += absorbed;

			backRadiation += reflected;
		}

		displayCell.EnergyThermalBackRadiation += backRadiation;
		displayCell.EnergySolarAbsorbedSurface += solarRadiationAbsorbed;
		displayCell.Heat = solarRadiationAbsorbed;

		float radiationToSurface = solarRadiationAbsorbed + backRadiation;

		// ice
		float remainingIceMass = last.IceMass;
		{
			// melt ice at surface from air temp and incoming radiation
			if (remainingIceMass > 0)
			{
				float radiationAbsorbedByIce = radiationToSurface * iceCoverage;
				radiationToSurface -= radiationAbsorbedByIce;
				radiationAbsorbedByIce += thermalEnergyRadiatedToIce;

				// world.Data.SpecificHeatIce * world.Data.MassIce == KJ required to raise one cubic meter by 1 degree
				if (radiationAbsorbedByIce > 0)
				{
					// Remove the latent heat from the incoming energy
					float iceMelted = math.min(remainingIceMass, radiationAbsorbedByIce / WorldData.LatentHeatWaterLiquid);
					next.IceMass -= iceMelted;
					next.WaterMass += iceMelted;
					next.WaterEnergy += iceMelted * (WorldData.SpecificHeatWater * WorldData.FreezingTemperature);
					remainingIceMass -= iceMelted;
				}
				else
				{
					next.WaterEnergy += radiationAbsorbedByIce * waterCoverage;
					next.GroundEnergy += radiationAbsorbedByIce * (1.0f - waterCoverage);
				}
				if (last.AirTemperature > WorldData.FreezingTemperature)
				{
					// Remove the latent heat of water from the air
					float temperatureDiff = last.AirTemperature - WorldData.FreezingTemperature;
					float energyTransfer = math.min(last.AirEnergy, temperatureDiff * worldData.SecondsPerTick * iceCoverage * worldData.IceAirConductionCooling);
					float iceMeltedFromConduction = remainingIceMass * math.saturate(energyTransfer / WorldData.LatentHeatWaterLiquid);
					next.AirEnergy -= energyTransfer;
					next.IceMass -= iceMeltedFromConduction;
					next.WaterMass += iceMeltedFromConduction;
					next.WaterEnergy += iceMeltedFromConduction * (WorldData.SpecificHeatWater * WorldData.FreezingTemperature);
					remainingIceMass -= iceMeltedFromConduction;
					displayCell.EnergySurfaceConduction -= energyTransfer;
				}

			}

			// freeze the top meter based on surface temperature (plus incoming radiation)
			if (iceCoverage < 1)
			{
				if (last.WaterMass > 0)
				{
					// world.Data.SpecificHeatIce * world.Data.MassIce == KJ required to raise one cubic meter by 1 degree
					float specificHeatWater = Atmosphere.GetSpecificHeatOfWater(last.WaterMass, last.SaltMass);
					float seaWaterHeatingRate = WorldData.MassWater / specificHeatWater;
					//float surfaceTemp = lowerAirTemperature + incomingRadiation * seaWaterHeatingRate;
					float localHeating = 0; // TODO: add in some local heating
					float surfaceTemp = (last.AirTemperature + localHeating) * (1.0f - iceCoverage) + last.WaterTemperature * iceCoverage;
					if (surfaceTemp < WorldData.FreezingTemperature)
					{
						float iceMassFrozen = math.min(last.WaterMass, math.min(math.max(0, worldData.FullIceCoverage * WorldData.MassIce - last.IceMass), (WorldData.FreezingTemperature - surfaceTemp) * seaWaterHeatingRate));
						next.IceMass += iceMassFrozen;
						next.WaterMass -= iceMassFrozen;
						// TODO: shouldnt the latent heat be added to the air, not the water?
						next.WaterEnergy -= iceMassFrozen * (WorldData.SpecificHeatWater * surfaceTemp - WorldData.LatentHeatWaterLiquid);
					}

					// TODO this should be using absolute pressure not barometric
					float inverseLowerAirPressure = 1.0f / last.AirPressure;
					// evaporation
					if (evaporationRate > 0)
					{
						float evapotranspiration;
						// TODO: absorb incoming radiation as latent heat (rather than from the surrounding air)
						EvaporateWater(
							waterCoverage,
							evaporationRate,
							last.GroundWaterDepth,
							last.WaterMass,
							last.SaltMass,
							last.WaterEnergy,
							last.WaterTemperature,
							ref next.AirWaterMass,
							ref next.AirEnergy,
							ref next.WaterEnergy,
							ref next.WaterMass,
							out evaporation,
							out evapotranspiration);
						displayCell.EnergyEvapotranspiration += evapotranspiration;
						displayCell.Evaporation = evaporation;
					}
				}
			}
		}


		// absorbed by surface
		{
			// absorb the remainder and radiate heat
			float absorbedByLand = (1.0f - waterCoverage) * radiationToSurface;
			next.GroundEnergy += absorbedByLand;
			if (waterCoverage > 0)
			{
				// absorb remaining incoming radiation (we've already absorbed radiation in surface ice above)
				float absorbedByWater = waterCoverage * radiationToSurface;
				next.WaterEnergy += absorbedByWater;
				//
				// heat transfer (both ways) based on temperature differential
				// conduction to ice from below
				if (last.IceMass > 0 && last.WaterTemperature > WorldData.FreezingTemperature)
				{
					float oceanConductionRate = (last.WaterTemperature - WorldData.FreezingTemperature) * worldData.OceanIceConduction * worldData.SecondsPerTick * iceCoverage;

					float energyToIce = math.max(0, oceanConductionRate * iceCoverage);
					float iceMelted = math.min(remainingIceMass, energyToIce * worldData.inverseSpecificHeatIce);
					next.IceMass -= iceMelted;
					next.WaterMass += iceMelted;
					next.WaterEnergy += iceMelted * (WorldData.SpecificHeatWater * WorldData.FreezingTemperature - WorldData.LatentHeatWaterLiquid);
					remainingIceMass -= iceMelted;
				}
				// lose heat to air via conduction AND radiation
				if (iceCoverage < 1)
				{
					// when ocean is warmer than air, it creates a convection current, which makes conduction more efficient)
					float oceanConduction = (last.WaterTemperature - last.AirTemperature) * worldData.SecondsPerTick * (1.0f - iceCoverage) * math.min(1.0f, last.WaterDepth / worldData.WaterAirConductionDepth);
					if (oceanConduction > 0)
					{
						oceanConduction *= worldData.OceanAirConductionWarming;
					}
					else
					{
						oceanConduction *= worldData.OceanAirConductionCooling;
					}
					next.AirEnergy += oceanConduction;
					next.WaterEnergy -= oceanConduction;
					displayCell.EnergyOceanConduction += oceanConduction;
					displayCell.EnergySurfaceConduction += oceanConduction;
				}

				if (last.WaterTemperature < WorldData.FreezingTemperature)
				{
					float specificHeatSaltWater = (WorldData.SpecificHeatWater * last.WaterMass + WorldData.SpecificHeatWater * last.SaltMass);
					float massFrozen = math.min(last.WaterMass,
						specificHeatSaltWater * (last.WaterTemperature - WorldData.FreezingTemperature) /
						(WorldData.LatentHeatWaterLiquid - WorldData.FreezingTemperature * (WorldData.SpecificHeatWater + WorldData.SpecificHeatIce)));

					next.IceMass += massFrozen;
					next.WaterMass -= massFrozen;
					next.WaterEnergy -= massFrozen * (WorldData.SpecificHeatWater * WorldData.FreezingTemperature - WorldData.LatentHeatWaterLiquid);
				}
			}
		}

		if (Atmosphere.GetAirTemperature(next.AirEnergy, next.AirMass, next.CloudMass, next.AirWaterMass) > 350)
		{
			surfaceElevation = 0;
		}


	}

	static private void EvaporateWater(
		float waterCoverage,
		float evapRate,
		float waterTableDepth,
		float shallowWaterMass,
		float shallowSaltMass,
		float shallowWaterEnergy,
		float shallowWaterTemperature,
		ref float newHumidity,
		ref float newLowerAirEnergy,
		ref float newShallowWaterEnergy,
		ref float newShallowWaterMass,
		out float evaporation,
		out float evapotranspiration)
	{
		evaporation = 0;
		evapotranspiration = 0;


		if (waterCoverage > 0)
		{
			float evapMass = math.min(shallowWaterMass, waterCoverage * evapRate);
			newHumidity += evapMass;
			newShallowWaterMass -= evapMass;
			evaporation += evapMass;
			// this sucks energy out of the lower atmosphere since it uses up some energy to fill up the latent heat of water vapor
			newShallowWaterEnergy -= evapMass * WorldData.SpecificHeatWater * shallowWaterTemperature;
			newLowerAirEnergy += evapMass * (WorldData.SpecificHeatWaterVapor * shallowWaterTemperature - WorldData.LatentHeatWaterVapor);
			evapotranspiration = evapMass * (WorldData.LatentHeatWaterVapor + WorldData.SpecificHeatWaterVapor * shallowWaterTemperature);
		}
	}



	private void DoVerticalWaterMovement(int i, ref SimCell last, ref SimWind lastWind, ref SimCell next, ref DisplayCell display, float surfaceElevation, float dewPoint, float evaporationRate)
	{


		// condensation
		if (last.RelativeHumidity > 1)
		{
			float condensationMass = next.AirWaterMass * (last.RelativeHumidity - 1.0f) / last.RelativeHumidity;
			next.AirWaterMass -= condensationMass;
			next.AirEnergy -= condensationMass * (WorldData.SpecificHeatWaterVapor * last.AirTemperature - WorldData.LatentHeatWaterVapor);
			if (last.AirTemperature <= WorldData.FreezingTemperature)
			{
				next.IceMass += condensationMass;
			}
			else
			{
				next.WaterMass += condensationMass;
				next.WaterEnergy += condensationMass * (WorldData.SpecificHeatWater * last.AirTemperature);
			}
		}

		if (lastWind.WindVertical > 0)
		{
			float humidityToCloud = math.min(1.0f, lastWind.WindVertical * worldData.SecondsPerTick / (last.CloudElevation - surfaceElevation)) * next.AirWaterMass * worldData.HumidityToCloudPercent;
			next.CloudMass += humidityToCloud;
			next.AirWaterMass -= humidityToCloud;

			// TODO: figure out what to do about the 2 layers of atmosphere
			// We're moving the latent heat of water vapor here since we want it to heat up the upper air around the cloud
			next.AirEnergy -= humidityToCloud * WorldData.SpecificHeatWaterVapor * last.AirTemperature;
			next.AirEnergy += humidityToCloud * (WorldData.SpecificHeatWater * last.AirTemperature + WorldData.LatentHeatWaterVapor);
		}

		if (last.CloudMass > 0)
		{

			// TODO: airDesntiy and rainDensity should probably be cleaned up (derived from other data?)
			float rainDropVolume = math.max(worldData.rainDropMinSize, last.CloudDropletMass / (last.CloudMass * worldData.waterDensity));
			float rainDropRadius = math.min(math.pow(rainDropVolume, 0.333f), worldData.rainDropMaxSize);
			float rainDropVelocity = lastWind.WindVertical - math.sqrt(8 * rainDropRadius * worldData.waterDensity * PlanetState.Gravity / (3 * worldData.airDensity * worldData.rainDropDragCoefficient));

			next.CloudDropletMass = math.max(0, next.CloudDropletMass + last.CloudMass * (worldData.RainDropFormationSpeedTemperature / dewPoint * math.pow(math.max(0, -rainDropVelocity) * worldData.RainDropCoalescenceWind, 2)));

			if (rainDropVelocity < 0 && last.CloudDropletMass > 0)
			{
				float rainfall = last.CloudMass * math.saturate(-rainDropVelocity * worldData.RainfallRate);

				next.CloudMass -= rainfall;
				next.CloudDropletMass = math.max(0, next.CloudDropletMass - rainfall / last.CloudMass);

				{
					float rainDropFallTime = -last.CloudElevation / rainDropVelocity;
					// evap rate is based on full tile surface coverage, an occurs in the top millimeter
					float rainDropSurfaceArea = 4 * math.PI * rainDropRadius * rainDropRadius;
					float totalRainSurfaceArea = rainfall / (worldData.waterDensity * rainDropVolume) * rainDropSurfaceArea;
					float rainEvapRate = evaporationRate * totalRainSurfaceArea * 1000 * staticState.InverseCellDiameter * staticState.InverseCellDiameter;
					float rainDropMassToHumidity = math.min(rainfall, rainDropFallTime * rainEvapRate * worldData.TicksPerSecond);
					rainfall -= rainDropMassToHumidity;
					next.AirWaterMass += rainDropMassToHumidity;
					// This sucks heat out of the lower atmosphere in the form of latent heat of water vapor
					next.AirEnergy -= rainDropMassToHumidity * last.AirTemperature * WorldData.SpecificHeatWater;
					next.AirEnergy += rainDropMassToHumidity * (last.AirTemperature * WorldData.SpecificHeatWaterVapor - WorldData.LatentHeatWaterVapor);
				}
				if (rainfall > 0)
				{
					display.Rainfall = rainfall;
					next.WaterMass += rainfall;
					// No real state change here
					float energyTransfer = rainfall * last.AirTemperature * WorldData.SpecificHeatWater;
					next.WaterEnergy += energyTransfer;
					next.AirEnergy -= energyTransfer;
				}
			}

			// dissapation
			float dissapationSpeed = math.min(1.0f, worldData.CloudDissapationRateWind * math.max(0, -lastWind.WindVertical) + worldData.CloudDissapationRateDryAir) * (1.0f - last.RelativeHumidity);
			float dissapationMass = last.CloudMass * dissapationSpeed;
			next.CloudDropletMass = math.max(0, next.CloudDropletMass - dissapationSpeed);
			next.CloudMass -= dissapationMass;
			next.AirWaterMass += dissapationMass;
			next.AirEnergy -= dissapationMass * (last.AirTemperature * WorldData.SpecificHeatWater + WorldData.LatentHeatWaterVapor);
			next.AirEnergy += dissapationMass * WorldData.SpecificHeatWaterVapor * last.AirTemperature;
		}

	}


	//static private void SeepWaterIntoGround(World world, float groundWater, float groundEnergy, float shallowWaterMass, float soilFertility, float waterTableDepth, float shallowWaterTemperature, ref float newGroundWater, ref float newShallowWater, ref float newGroundEnergy, ref float newShallowEnergy)
	//{
	//	float maxGroundWater = soilFertility * waterTableDepth * world.Data.MaxSoilPorousness * world.Data.MassWater;
	//	if (groundWater >= maxGroundWater && groundWater > 0)
	//	{
	//		float massTransfer = groundWater - maxGroundWater;
	//		newShallowWater += massTransfer;
	//		newGroundWater -= massTransfer;
	//		float energyTransfer = massTransfer / groundWater * groundEnergy; // TODO: this isn't great, some of that ground energy is in the terrain, not just in the water
	//		newShallowEnergy += energyTransfer;
	//		newGroundEnergy -= energyTransfer;
	//	}
	//	else if (shallowWaterMass > 0)
	//	{
	//		float massTransfer = Mathf.Min(shallowWaterMass, Math.Min(soilFertility * world.Data.GroundWaterReplenishmentSpeed * world.Data.SecondsPerTick, maxGroundWater - groundWater));
	//		newGroundWater += massTransfer;
	//		newShallowWater -= massTransfer;
	//		float energyTransfer = massTransfer * shallowWaterTemperature * world.Data.SpecificHeatWater;
	//		newShallowEnergy -= energyTransfer;
	//		newGroundEnergy += energyTransfer;
	//	}
	//}



	private void DoDiffusion(int i, ref SimCell cell)
	{
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = staticState.Neighbors[neighborIndex];
			if (n >= 0)
			{
				cell.CloudCoverage -= Diffusion[neighborIndex].Cloud * DiffusionLimit[i].Cloud;
				cell.RelativeHumidity -= Diffusion[neighborIndex].Humidity * DiffusionLimit[i].Humidity;
				cell.WaterDepth -= Diffusion[neighborIndex].Water * DiffusionLimit[i].Water;
				for (int k = 0; k < 6; k++)
				{
					int nToMeIndex = n * 6 + k;
					if (staticState.Neighbors[nToMeIndex] == i)
					{
						cell.CloudCoverage += Diffusion[nToMeIndex].Cloud * DiffusionLimit[n].Cloud;
						cell.RelativeHumidity += Diffusion[nToMeIndex].Humidity * DiffusionLimit[n].Humidity;
						cell.WaterDepth += Diffusion[nToMeIndex].Water * DiffusionLimit[n].Water;
						break;
					}
				}
			}
		}
	}
}
