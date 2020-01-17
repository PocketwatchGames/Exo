using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

[BurstCompile]
public struct TickCellJob : IJobParallelFor {

	public NativeArray<SimStateCell> Cells;

	[ReadOnly] public SimPlanetState PlanetState;
	[ReadOnly] public NativeArray<CellDiffusion> Diffusion;
	[ReadOnly] public NativeArray<CellDiffusion> DiffusionLimit;
	[ReadOnly] public NativeArray<SimStateCell> Last;
	[ReadOnly] public WorldData worldData;
	[ReadOnly] public StaticState staticState;

	public void Execute(int i)
	{
		var curCell = Last[i];
		var cell = new SimStateCell();

		cell.CloudCoverage = curCell.CloudCoverage;
		cell.RelativeHumidity = curCell.RelativeHumidity;
		cell.WaterDepth = curCell.WaterDepth;
		cell.CloudElevation = curCell.CloudElevation;
		cell.Elevation = curCell.Elevation;
		cell.Roughness = curCell.Roughness;
		cell.IceMass = curCell.IceMass;
		cell.Vegetation = curCell.Vegetation;
		cell.SoilFertility = curCell.SoilFertility;
		cell.WaterMass = curCell.WaterMass;
		cell.WaterEnergy = curCell.WaterEnergy;
		cell.SaltMass = curCell.SaltMass;
		cell.GroundEnergy = curCell.GroundEnergy;
		cell.GroundWater = curCell.GroundWater;
		cell.GroundWaterDepth = curCell.GroundWaterDepth;
		cell.AirMass = curCell.AirMass;
		cell.AirEnergy = curCell.AirEnergy;
		cell.CloudMass = curCell.CloudMass;
		cell.CloudDropletMass = curCell.CloudDropletMass;
		cell.AirWaterMass = curCell.AirWaterMass;
		cell.AirTemperature = curCell.AirTemperature;
		cell.AirPressure = curCell.AirPressure;
		cell.WaterTemperature = curCell.WaterTemperature;
		cell.WindVertical = curCell.WindVertical;
		cell.WindSurface = curCell.WindSurface;
		cell.WindTropopause = curCell.WindTropopause;

		DoEnergyCycle(i, ref cell);
		DoVerticalWaterMovement(i, ref cell);
		DoDiffusion(i, ref cell);

		Cells[i] = cell;

	}

	private void DoEnergyCycle(int i, ref SimStateCell next)
	{
		var last = Last[i];

		// TODO: account for ice here
		float surfaceElevation = last.Elevation + last.WaterDepth;
		float iceCoverage = math.min(1.0f, math.pow(last.IceMass * worldData.inverseFullIceCoverage, 0.6667f));
		float waterCoverage = math.min(1.0f, math.pow(last.WaterDepth * worldData.inverseFullWaterCoverage, 0.6667f));
		float canopyCoverage = math.min(1.0f, math.pow(last.Vegetation * worldData.inverseFullCanopyCoverage, 0.6667f));
		float relativeHumidity = Atmosphere.GetRelativeHumidity(ref worldData, last.AirTemperature, last.AirWaterMass, last.AirMass, worldData.inverseDewPointTemperatureRange);
		float evapRate = Atmosphere.GetEvaporationRate(ref worldData, iceCoverage, last.AirTemperature, relativeHumidity, worldData.inverseEvapTemperatureRange);
		float dewPoint = Atmosphere.GetDewPoint(ref worldData, last.AirTemperature, relativeHumidity);
		float cloudElevation = Atmosphere.GetCloudElevation(ref worldData, last.AirTemperature, dewPoint, surfaceElevation);
		float cloudCoverage = math.min(1.0f, math.pow(last.CloudMass * worldData.inverseCloudMassFullAbsorption, 0.6667f)); // bottom surface of volume

		float3 sunVector = -(math.rotate(PlanetState.Rotation, staticState.SphericalPosition[i]) + PlanetState.Position);

		float waterSlopeAlbedo = math.pow(1.0f - math.max(0, sunVector.z), 9);
		//float groundWaterSaturation = Animals.GetGroundWaterSaturation(state.GroundWater[index], state.WaterTableDepth[index], soilFertility * world.Data.MaxSoilPorousness);


		// SOLAR RADIATION

		float solarRadiationAbsorbed = 0;
		float solarRadiation = PlanetState.SolarRadiation * worldData.SecondsPerTick * math.max(0, (sunVector.z + worldData.sunHitsAtmosphereBelowHorizonAmount) * worldData.inverseSunAtmosphereAmount);

		// get the actual atmospheric depth here based on radius of earth plus atmosphere
		//float inverseSunAngle = PIOver2 + sunAngle;
		//float angleFromSunToLatitudeAndAtmophereEdge = math.Asin(state.PlanetRadius * math.Sin(inverseSunAngle) / (state.PlanetRadius + world.Data.TropopauseElevation));
		//float angleFromPlanetCenterToLatitudeAndAtmosphereEdge = math.PI - inverseSunAngle - angleFromSunToLatitudeAndAtmophereEdge;
		//float atmosphericDepthInMeters = math.Sin(angleFromPlanetCenterToLatitudeAndAtmosphereEdge) * state.PlanetRadius / math.Sin(angleFromSunToLatitudeAndAtmophereEdge);
		//float atmosphericDepth = math.max(1.0f, atmosphericDepthInMeters / world.Data.TropopauseElevation);

		float atmosphericDepth = 1.0f + sunVector.y;

		// These constants obtained here, dunno if I've interpreted them correctly
		// https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass

		////// MAJOR TODO:
		///// USE THIS LINK: https://www.ftexploring.com/solar-energy/sun-angle-and-insolation2.htm
		/// With the sun 90 degrees above the horizon (SEA° = 90°), the air mass lowers the intensity of the sunlight from the 1,367 W / m2 that it is in outerspace down to about 1040 W / m2.
		//				float consumedByAtmosphere = 1.0f - math.pow(0.7f, math.pow(atmosphericDepth, 0.678f));


		if (solarRadiation > 0)
		{
//			globalEnergyIncoming += solarRadiation;

			// TODO: reflect/absorb more in the atmosphere with a lower sun angle

			// reflect some rads off atmosphere and clouds
			// TODO: this process feels a little broken -- are we giving too much priority to reflecting/absorbing in certain layers?
			float energyReflectedAtmosphere = solarRadiation * math.min(1, worldData.AtmosphericHeatReflection * (last.AirMass + last.AirWaterMass));
			solarRadiation -= energyReflectedAtmosphere;
//			globalEnergySolarReflectedAtmosphere += energyReflectedAtmosphere;

			if (last.CloudMass > 0)
			{
				float cloudTemperatureAlbedo = WorldData.AlbedoIce + (WorldData.AlbedoWater - WorldData.AlbedoIce) * math.clamp((dewPoint - worldData.minCloudFreezingTemperature) / (worldData.maxCloudFreezingTemperature - worldData.minCloudFreezingTemperature),0,1);
				float rainDropSizeAlbedo = math.clamp(1.0f - last.CloudDropletMass / last.CloudMass, 0,1) * (worldData.rainDropSizeAlbedoMax - worldData.rainDropSizeAlbedoMin) + worldData.rainDropSizeAlbedoMin;
				float cloudReflectivity = math.min(1.0f, worldData.cloudAlbedo * cloudTemperatureAlbedo * last.CloudMass * rainDropSizeAlbedo / math.max(worldData.maxCloudSlopeAlbedo, 1.0f - waterSlopeAlbedo));
				float energyReflectedClouds = solarRadiation * cloudReflectivity;
				solarRadiation -= energyReflectedClouds;

				float absorbedByCloudsIncoming = solarRadiation * math.min(1.0f, worldData.AtmosphericHeatAbsorption * last.CloudMass);
				solarRadiation -= absorbedByCloudsIncoming;
				next.AirEnergy += absorbedByCloudsIncoming;
				//globalEnergySolarAbsorbedClouds += absorbedByCloudsIncoming;
				//globalEnergySolarAbsorbedAtmosphere += absorbedByCloudsIncoming;
				//globalSolarEnergyReflectedClouds += energyReflectedClouds;
			}

			// Absorbed by atmosphere
			// stratosphere accounts for about a quarter of atmospheric mass
			//	float absorbedByStratosphere = incomingRadiation * world.Data.AtmosphericHeatAbsorption * (state.StratosphereMass / massOfAtmosphericColumn);

			float atmosphereAbsorptionRate = math.min(1, worldData.AtmosphericHeatAbsorption * (last.AirMass + last.AirWaterMass));
			float absorbedByLowerAtmosphereIncoming = solarRadiation * atmosphereAbsorptionRate * atmosphericDepth;

			next.AirEnergy += absorbedByLowerAtmosphereIncoming;
			solarRadiation -= absorbedByLowerAtmosphereIncoming;
			//globalEnergySolarAbsorbedAtmosphere += absorbedByLowerAtmosphereIncoming + absorbedByUpperAtmosphereIncoming;

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
					energyReflected += solarRadiation * math.clamp(heatReflectedLand,0,1) * (1.0f - iceCoverage) * (1.0f - waterCoverage);
				}
				solarRadiation -= energyReflected;

				// TODO: do we absorb some of this energy on the way back out of the atmosphere?
				//globalEnergySolarReflectedSurface += energyReflected;
			}

			//solarRadiationAbsorbed += solarRadiation;

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
			thermalEnergyRadiatedToAir += oceanRadiation - thermalEnergyRadiatedToIce;
			//globalEnergyThermalOceanRadiation += oceanRadiation;
		}

		// TODO: track and emit heat from ice

		float atmosphereEmissivity = Atmosphere.GetAtmosphericEmissivity(last.AirMass, last.AirMass * PlanetState.CarbonDioxide, last.AirWaterMass, last.CloudMass);
		float emittedByAtmosphere = 0;
		float surfaceEnergyReflected = 0;

		// Thermal energy from surface to air, space, reflected off clouds
		{
			//globalEnergyThermalSurfaceRadiation += thermalEnergyRadiatedToAir;
			float energyThroughAtmosphericWindow = thermalEnergyRadiatedToAir * worldData.EnergyLostThroughAtmosphereWindow;
			thermalEnergyRadiatedToAir -= energyThroughAtmosphericWindow;
			//globalEnergyThermalOutAtmosphericWindow += energyThroughAtmosphericWindow;

			float absorbed = thermalEnergyRadiatedToAir * atmosphereEmissivity;
			thermalEnergyRadiatedToAir -= absorbed;
			next.AirEnergy += absorbed;
			//globalEnergyThermalAbsorbedAtmosphere += absorbed;

			surfaceEnergyReflected = thermalEnergyRadiatedToAir * cloudCoverage * worldData.CloudOutgoingReflectionRate;
			reflected += surfaceEnergyReflected;
			//globalEnergyThermalOutAtmosphere += thermalEnergyRadiatedToAir - surfaceEnergyReflected;
		}

		// atmosphere radiation
		{
			float energyEmitted = Atmosphere.GetRadiationRate(last.AirTemperature, atmosphereEmissivity) * worldData.SecondsPerTick;
			next.AirEnergy -= 2 * energyEmitted;
			emittedByAtmosphere += 2 * energyEmitted;

			backRadiation += energyEmitted;

			float energyThroughAtmosphericWindow = energyEmitted * worldData.EnergyLostThroughAtmosphereWindow;
			energyEmitted -= energyThroughAtmosphericWindow;
			//globalEnergyThermalOutAtmosphericWindow += energyThroughAtmosphericWindow;

			float lowerEnergyReflected = energyEmitted * cloudCoverage * worldData.CloudOutgoingReflectionRate;
			reflected += lowerEnergyReflected;
			//globalEnergyThermalOutAtmosphere += lowerEnergyEmitted - lowerEnergyReflected;
		}

		// reflected thermal radiation
		{
			float absorbed = reflected * atmosphereEmissivity;
			next.AirEnergy += absorbed;
			reflected -= absorbed;

			//globalEnergyThermalAbsorbedAtmosphere += surfaceEnergyReflected * upperAtmosphereInfraredAbsorption + surfaceEnergyReflected * (1.0f - upperAtmosphereInfraredAbsorption) * lowerAtmosphereInfraredAbsorption;

			backRadiation += reflected;
		}

		//globalEnergyThermalBackRadiation += backRadiation;
		//globalEnergySolarAbsorbedSurface += solarRadiationAbsorbed;

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
					float iceMeltedFromConduction = remainingIceMass * math.clamp(energyTransfer / WorldData.LatentHeatWaterLiquid,0,1);
					next.AirEnergy -= energyTransfer;
					next.IceMass -= iceMeltedFromConduction;
					next.WaterMass += iceMeltedFromConduction;
					next.WaterEnergy += iceMeltedFromConduction * (WorldData.SpecificHeatWater * WorldData.FreezingTemperature);
					remainingIceMass -= iceMeltedFromConduction;
					//globalEnergySurfaceConduction -= energyTransfer;
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
					if (evapRate > 0)
					{
						float evapotranspiration;
						// TODO: absorb incoming radiation as latent heat (rather than from the surrounding air)
						EvaporateWater(
							waterCoverage,
							evapRate,
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
						//globalEnergyEvapotranspiration += evapotranspiration;
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
					//globalEnergyOceanConduction += oceanConduction;
					//globalEnergySurfaceConduction += oceanConduction;
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



	private void DoVerticalWaterMovement(int i, ref SimStateCell cell)
	{
		var curCell = Last[i];
		float evap = math.min(curCell.WaterDepth, 1f) * math.clamp(1.0f - curCell.RelativeHumidity / 100, 0, 1);
		cell.WaterDepth -= evap;
		cell.RelativeHumidity += evap;
		float hToC = curCell.RelativeHumidity / 100;
		cell.RelativeHumidity -= hToC;
		cell.CloudCoverage += hToC;
		if (cell.CloudCoverage > 1000)
		{
			cell.WaterDepth += curCell.CloudCoverage;
			cell.CloudCoverage = 0;
		}
	}


	private void DoDiffusion(int i, ref SimStateCell cell)
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
