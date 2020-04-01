﻿//#define DISABLE_FREEZE_TOP
//#define DISABLE_FREEZE_BOTTOM
//#define DISABLE_EVAPORATION
//#define DISABLE_CONDENSATION
//#define DISABLE_CLOUD_DISSAPATION
//#define DISABLE_RAINFALL
//#define DISABLE_MELTING_TOP
//#define DISABLE_MELTING_BOTTOM
#define DISABLE_ERUPTION
//#define DISABLE_EVAPORATION_LATENT_HEAT
//#define DISABLE_FLORA_WATER_ABSORPTION

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

[BurstCompile]
public struct FluxWaterJob : IJobParallelFor {
	public NativeArray<float> EvaporatedWaterMass;
	public NativeArray<float> FrozenMass;
	public NativeArray<float> FrozenTemperature;
	public NativeArray<float> SaltPlume;
	public NativeArray<float> LatentHeatWater;
	public NativeArray<float> LatentHeatAir;
	public NativeArray<float> PlanktonMassDelta;
	public NativeArray<float> PlanktonGlucoseDelta;
	public NativeArray<float> PlanktonDeath;
	public NativeArray<float> WaterCarbonDelta;
	[ReadOnly] public NativeArray<float> SolarRadiation;
	[ReadOnly] public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> WaterCarbon;
	[ReadOnly] public NativeArray<float> PlanktonMass;
	[ReadOnly] public NativeArray<float> PlanktonGlucoseMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirVapor;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> AirLayerElevation;
	[ReadOnly] public NativeArray<float3> SurfaceWind;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public float WaterHeatingDepth;
	[ReadOnly] public float FreezePointReductionPerSalinity;
	[ReadOnly] public float PlanktonEnergyForPhotosynthesis;
	[ReadOnly] public float PlanktonDensityMax;
	[ReadOnly] public float PlanktonPhotosynthesisSpeed;
	[ReadOnly] public float PlanktonCarbonDioxideExtractionEfficiency;
	[ReadOnly] public float PlanktonRespirationSpeed;
	[ReadOnly] public float PlanktonGrowthRate;
	[ReadOnly] public float PlanktonDeathRate;
	[ReadOnly] public float PlanktonRespirationPerDegree;
	public void Execute(int i)
	{
		float temperature = WaterTemperature[i];
		float waterMass = WaterMass[i];
		float saltMass = SaltMass[i];

		float latentHeatFromAir = 0;
		float evapMass = 0;
		float frozenMass = 0;
		float freezingTemperature = 0;
		float solarRadiation = SolarRadiation[i];
		float energyFlux = solarRadiation;
		float saltPlume = 0;
		float glucoseDelta = 0;
		float planktonMassDelta = 0;
		float planktonDeath = 0;
		float waterCarbonDelta = 0;

		if (waterMass > 0)
		{

			float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], AirLayerElevation[i]);

#if !DISABLE_EVAPORATION
			evapMass = Atmosphere.GetEvaporationMass(AirMass[i], AirPressure[i], AirVapor[i], SurfaceWind[i], temperature, waterMass);
			waterMass -= evapMass;

#if !DISABLE_EVAPORATION_LATENT_HEAT
			const bool latentHeatFromAirOrWater = false;
			if (latentHeatFromAirOrWater)
			{
				latentHeatFromAir = evapMass * WorldData.LatentHeatWaterVapor;
			} else
			{
				energyFlux -= evapMass * WorldData.LatentHeatWaterVapor;
			}
			//energyTop -= evapMass * WorldData.LatentHeatWaterVapor;
#endif

#endif

			float specificHeatSaltWater = (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatSalt * saltMass) / (waterMass + saltMass);
			float heatingMass = math.min(waterMass, WaterHeatingDepth * WorldData.MassWater);
			freezingTemperature = Atmosphere.GetFreezingPoint(Atmosphere.GetWaterSalinity(waterMass, saltMass), FreezePointReductionPerSalinity);

#if !DISABLE_FREEZE_TOP
			if (waterMass > 0)
			{
				// TODO: this is only applicable for water exposed to the air... what about water under ice?
				if (airTemperatureAbsolute < freezingTemperature)
				{
					// TODO: some fo the energy in the ehating mass should also be accounted for
					float energyToRelease = (freezingTemperature - airTemperatureAbsolute) * specificHeatSaltWater * heatingMass;
					frozenMass = math.min(waterMass, energyToRelease / WorldData.LatentHeatWaterLiquid);
					saltPlume += frozenMass / waterMass * saltMass;
					energyFlux += frozenMass * WorldData.LatentHeatWaterLiquid;
					waterMass -= frozenMass;
				}
			}
#endif

#if !DISABLE_MARINE_PHOTOSYNTHESIS

			float planktonMass = PlanktonMass[i];
			float waterCarbon = WaterCarbon[i];
			if (planktonMass > 0)
			{
				float glucose = PlanktonGlucoseMass[i];
				float inversePlanktonMass = 1.0f / planktonMass;
				float waterSupply = math.min(1, waterMass * inversePlanktonMass * PlanktonDensityMax);

				// Photosynthesis: Consume solarRadiation, carbon dioxide, to produce glucose (water and oxygen are ignored)
				// TODO: carbon dioxide extraction should curve to a limit rather than cap
				float photosynthesis =
					PlanktonPhotosynthesisSpeed * planktonMass
					* waterSupply// minimum water mass requirements
					* (1.0f - math.pow(math.min(1, glucose * inversePlanktonMass), 3)) // glucose saturation
					* math.min(1, solarRadiation / PlanktonEnergyForPhotosynthesis * inversePlanktonMass) // energy
					* math.min(1, waterCarbon * PlanktonCarbonDioxideExtractionEfficiency * inversePlanktonMass); // carbon dioxide 

				energyFlux -= photosynthesis * PlanktonEnergyForPhotosynthesis;
				waterCarbonDelta -= photosynthesis;
				glucoseDelta += photosynthesis;

				// Respiration: Consume glucose, oxygen, produce water and carbon dioxide
				float respirationExpenditure = math.min(waterCarbon, math.min(glucose,
					PlanktonRespirationSpeed 
					* math.max(0, 1.0f + (temperature - WorldData.FreezingTemperature) * PlanktonRespirationPerDegree) // temperature
					* waterSupply //water
					* glucose)); // energy (glucose)

				float respiration = respirationExpenditure / PlanktonRespirationSpeed;

				energyFlux += respiration * PlanktonEnergyForPhotosynthesis;
				glucoseDelta -= respirationExpenditure;
				waterCarbonDelta += respirationExpenditure;

				// Growth: Consume glucose, produce plant growth
				float growth =
					PlanktonGrowthRate
					* planktonMass
					* respiration * inversePlanktonMass
					* Utils.Sqr(glucose * inversePlanktonMass);

				glucoseDelta -= growth;
				planktonMassDelta += growth;

				// Death: Produce carbon dioxide, water
				float deathPercent = Utils.Sqr(math.max(0, planktonMass - respiration) * inversePlanktonMass) * PlanktonDeathRate;
				float death = deathPercent * planktonMass;
				// convert mass back through the respiration process, and dump any stored glucose and water
				glucoseDelta -= deathPercent * glucose;
				planktonMassDelta -= death;
				planktonDeath += death + deathPercent * glucose;
			}
#endif

		}

		EvaporatedWaterMass[i] = evapMass;
		FrozenMass[i] = frozenMass;
		FrozenTemperature[i] = freezingTemperature;
		SaltPlume[i] = saltPlume;
		LatentHeatWater[i] = energyFlux;
		LatentHeatAir[i] += -latentHeatFromAir;
		WaterCarbonDelta[i] = waterCarbonDelta;
		PlanktonGlucoseDelta[i] = glucoseDelta;
		PlanktonMassDelta[i] = planktonMassDelta;
		PlanktonDeath[i] = planktonDeath;
	}
}


[BurstCompile]
public struct FluxCloudJob : IJobParallelFor {
	public NativeArray<float> DropletDelta;
	public NativeArray<float> PrecipitationMass;
	public NativeArray<float> PrecipitationTemperature;
	public NativeArray<float> EvaporationMass;
	[ReadOnly] public NativeArray<float3> LastVelocity;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> DewPoint;
	[ReadOnly] public NativeArray<float> LastDropletMass;
	[ReadOnly] public NativeArray<float> LastCloudMass;
	[ReadOnly] public NativeArray<float> AirDensityCloud;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	[ReadOnly] public NativeArray<float> SurfaceAirTemperaturePotential;
	[ReadOnly] public NativeArray<float> SurfaceLayerElevation;
	[ReadOnly] public NativeArray<float> SurfaceLayerMiddle;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float RainDropMinSize;
	[ReadOnly] public float RainDropMaxSize;
	[ReadOnly] public float RainDropDragCoefficient;
	[ReadOnly] public float CloudDissapationRateDryAir;
	[ReadOnly] public float CloudDissapationRateWind;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float surfaceElevation = SurfaceLayerElevation[i];

		float cloudMass = LastCloudMass[i];
		float dropletMass = LastDropletMass[i];
		float cloudElevation = CloudElevation[i];
		float3 lastVelocity = LastVelocity[i];
		float rainfallWaterMass = 0;
		float cloudEvaporationMass = 0;
		float dewPoint = DewPoint[i];

	//	var velocityUp = math.dot(velocity, Position[i]) * Position[i];

		float precipitationMass = 0;
#if !DISABLE_RAINFALL

		if (cloudMass > 0)
		{
			// TODO: improve this somehow
			dropletMass += cloudMass * 0.00000001f;


			float waterDensityAtElevation = Atmosphere.GetWaterDensityAtElevation(dewPoint, cloudElevation);
			float rainDropRadius = math.clamp(Atmosphere.GetDropletRadius(dropletMass, waterDensityAtElevation), RainDropMinSize, RainDropMaxSize);
			float rainDropVolume = 4 / 3 * math.PI * math.pow(rainDropRadius, 3);

			// TODO: See wikipedia Terminal velocity:
			//https://en.wikipedia.org/wiki/Terminal_velocity
			// We shouldn't be using the Air's buoyancy force as a stand in for the water droplets buoyancy
			// We should instead just be adding the vertical velocity of the air parcel
			float windVertical = math.dot(lastVelocity, Position[i]);
			float terminalVelocity = windVertical - math.sqrt(8 * rainDropRadius * waterDensityAtElevation * Gravity / (3 * AirDensityCloud[i] * RainDropDragCoefficient));
			if (terminalVelocity < 0 && dropletMass > 0)
			{

				float rainDropFallTime = (cloudElevation - surfaceElevation) / -terminalVelocity;
				if (rainDropFallTime < SecondsPerTick)
				{
					// TODO: account for drop size variance so we don't just dump it all at once
					rainfallWaterMass = cloudMass * math.saturate(1.0f - rainDropFallTime / SecondsPerTick);
					dropletMass *= 1.0f - rainfallWaterMass / cloudMass;
					cloudMass -= rainfallWaterMass;
					precipitationMass += rainfallWaterMass;
				}
			}
		}
#endif

#if !DISABLE_CLOUD_DISSAPATION
		// dissapation
		//if (cloudMass > 0)
		//{
		//	float dissapationSpeed = math.saturate((1.0f - RelativeHumidityCloud[i]) * EvaporationRate * (DewPoint[i] - EvapTemperatureMin) / (EvapTemperatureMax - EvapTemperatureMin));
		////	float dissapationSpeed = math.min(1.0f, CloudDissapationRateDryAir) * math.max(0, 1.0f - RelativeHumidityCloud[i]);
		//	cloudEvaporationMass = cloudMass * dissapationSpeed;
		//	dropletMass *= 1.0f - dissapationSpeed;
		//	cloudMass -= cloudEvaporationMass;
		//}
#endif

		if (cloudMass <= 0)
		{
			dropletMass = 0;
		}

		if (precipitationMass > 0)
		{
			PrecipitationTemperature[i] = Atmosphere.GetAbsoluteTemperature(SurfaceAirTemperaturePotential[i], SurfaceLayerMiddle[i]);
		}
		PrecipitationMass[i] = precipitationMass;
		DropletDelta[i] = dropletMass - LastDropletMass[i];
		EvaporationMass[i] = cloudEvaporationMass;
	}
}


[BurstCompile]
public struct FluxAirJob : IJobParallelFor {
	public NativeArray<float> LatentHeat;
	public NativeArray<float> CondensationGroundMass;
	public NativeArray<float> CondensationCloudMass;
	public NativeArray<float> DustUp;
	public NativeArray<float> DustDown;
	[ReadOnly] public NativeArray<float> TemperaturePotential;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastDust;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerMiddle;
	[ReadOnly] public NativeArray<float3> AirVelocity;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public float DustVerticalVelocity;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float condensationGroundMass = 0;
		float condensationCloudMass = 0;
		float energyFlux = 0;

#if !DISABLE_CONDENSATION
		float temperatureAbsolute = Atmosphere.GetAbsoluteTemperature(TemperaturePotential[i], LayerMiddle[i]);
		var relativeHumidity = Atmosphere.GetRelativeHumidity(AirMass[i], LastVapor[i], temperatureAbsolute, AirPressure[i]);
		if (relativeHumidity > 1.0f)
		{
			float aboveCloud = math.saturate((LayerElevation[i] - CloudElevation[i]) / LayerHeight[i]);
			float vaporToCondense = (relativeHumidity - 1.0f) / relativeHumidity * LastVapor[i];
			condensationCloudMass = aboveCloud * vaporToCondense;
			condensationGroundMass = (1.0f - aboveCloud) * vaporToCondense;
			energyFlux += vaporToCondense * WorldData.LatentHeatWaterVapor;
		}
#endif

		float layerHeight = LayerHeight[i];
		var velVertical = (math.dot(AirVelocity[i], Positions[i]) + DustVerticalVelocity) * SecondsPerTick;
		if (velVertical > 0)
		{
			DustDown[i] = 0;
			DustUp[i] = math.min(1, velVertical / LayerHeight[i]) * LastDust[i];
		}
		else
		{
			DustUp[i] = 0;
			DustDown[i] = math.min(1, -velVertical / LayerHeight[i]) * LastDust[i];
		}


		LatentHeat[i] = energyFlux;
		CondensationGroundMass[i] = condensationGroundMass;
		CondensationCloudMass[i] = condensationCloudMass;
	}
}



[BurstCompile]
public struct FluxIceJob : IJobParallelFor {
	public NativeArray<float> LatentHeatAir;
	public NativeArray<float> LatentHeatWater;
	public NativeArray<float> LatentHeatTerrain;
	public NativeArray<float> LatentHeatIce;
	public NativeArray<float> MeltedMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> WaterIceSurfaceArea;
	[ReadOnly] public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public float IceHeatingDepth;
	public void Execute(int i)
	{
		float meltedMass = 0;
		float iceMass = LastMass[i];
		float iceTemperature = Temperature[i];

#if !DISABLE_MELTING_TOP
		if (iceMass > 0)
		{
			if (iceTemperature > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (iceTemperature - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * iceMass;
				float melted = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				LatentHeatIce[i] -= melted * WorldData.LatentHeatWaterLiquid;
				iceMass -= melted;
				meltedMass += melted;
			}
		}
		if (iceMass > 0)
		{
			float iceDepth = iceMass / WorldData.MassIce;
			float heatingDepth = math.min(IceHeatingDepth, iceDepth);
			float heatingMass = heatingDepth * WorldData.MassIce;
			float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], LayerElevation[i]);
			if (airTemperatureAbsolute > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (airTemperatureAbsolute - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * heatingMass;
				float melted = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				LatentHeatAir[i] -= melted * WorldData.LatentHeatWaterLiquid;
				iceMass -= melted;
				meltedMass += melted;
			}
		}
#endif

#if !DISABLE_MELTING_BOTTOM
		float waterCoverage = WaterIceSurfaceArea[i];
		if (iceMass > 0 && waterCoverage > 0)
		{
			float iceDepth = iceMass / WorldData.MassIce;
			float heatingDepth = math.min(IceHeatingDepth, iceDepth);
			float heatingMass = heatingDepth * WorldData.MassIce;
			float waterTemperature = WaterTemperature[i];
			if (waterTemperature > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (waterTemperature - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * heatingMass * waterCoverage;
				float melted = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				float latentEnergy = melted * WorldData.LatentHeatWaterLiquid;
				LatentHeatWater[i] -= latentEnergy;
				iceMass -= melted;
				meltedMass += melted;
			}
		}
		float terrainCoverage = (1.0f - waterCoverage); // This ignores foliage and lava and conducts directly with the terrain for simplicity
		if (iceMass > 0 && waterCoverage < 1)
		{
			float iceDepth = iceMass / WorldData.MassIce;
			float heatingDepth = math.min(IceHeatingDepth, iceDepth);
			float heatingMass = heatingDepth * WorldData.MassIce;
			float terrainTemperature = TerrainTemperature[i];
			if (terrainTemperature > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (terrainTemperature - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * heatingMass * (1.0f - waterCoverage);
				float melted = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				float latentEnergy = melted * WorldData.LatentHeatWaterLiquid;
				LatentHeatTerrain[i] -= latentEnergy;
				iceMass -= melted;
				meltedMass += melted;
			}
		}
#endif

		MeltedMass[i] = meltedMass;
	}
}

[BurstCompile]
public struct FluxFloraJob : IJobParallelFor {
	public NativeArray<float> EvaporatedWaterMass;
	public NativeArray<float> LatentHeatAir;
	public NativeArray<float> LatentHeatFlora;
	public NativeArray<float> GroundWaterConsumed;
	public NativeArray<float> CarbonDioxideDelta;
	public NativeArray<float> OxygenDelta;
	public NativeArray<float> FloraMassDelta;
	public NativeArray<float> FloraWaterDelta;
	public NativeArray<float> FloraGlucoseDelta;
	public NativeArray<float> SurfaceWaterDelta;
	public NativeArray<float> FloraDeath;
	[ReadOnly] public NativeArray<float> FloraTemperature;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> CarbonDioxide;
	[ReadOnly] public NativeArray<float> FloraMass;
	[ReadOnly] public NativeArray<float> FloraGlucose;
	[ReadOnly] public NativeArray<float> FloraWater;
	[ReadOnly] public NativeArray<float> FloraCoverage;
	[ReadOnly] public NativeArray<float> GroundWater;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirVapor;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float3> SurfaceWind;
	[ReadOnly] public float GroundWaterMax;
	[ReadOnly] public float FloraWaterConsumptionRate;
	[ReadOnly] public float FloraGrowthTemperatureRangeInverse;
	[ReadOnly] public float FloraGrowthRate;
	[ReadOnly] public float FloraDeathRate;
	[ReadOnly] public float FloraEnergyForPhotosynthesis;
	[ReadOnly] public float FloraCarbonDioxideExtractionEfficiency;
	[ReadOnly] public float FloraOxygenExtractionEfficiency;
	[ReadOnly] public float FloraPhotosynthesisSpeed;
	[ReadOnly] public float FloraRespirationSpeed;
	[ReadOnly] public float FloraRespirationPerDegree;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float OxygenPercent;
	public void Execute(int i)
	{
		float floraMass = FloraMass[i];

		float latentHeatFromAir = 0;
		float latentHeatFromFlora = 0;
		float evapMass = 0;
		float groundWaterConsumed = 0;
		float floraMassDelta = 0;
		float carbonDioxideDelta = 0;
		float oxygenDelta = 0;
		float glucoseDelta = 0;
		float floraWaterDelta = 0;
		float surfaceWaterDelta = 0;
		float floraDeath = 0;

		if (floraMass > 0)
		{
			float temperature = FloraTemperature[i];
			float floraWater = FloraWater[i];
			float energy = SolarRadiationIn[i];


#if !DISABLE_PHOTOSYNTHESIS
			float glucose = FloraGlucose[i];
			float carbonDioxide = CarbonDioxide[i];
			float floraCoverage = FloraCoverage[i];
			float soilFertility = SoilFertility[i];
			float inverseFloraMass = 1.0f / floraMass;

			float airDensity = Atmosphere.GetAirDensity(Atmosphere.GetPressureAtElevation(LayerElevation[i],Gravity,AirPressure[i],AirTemperaturePotential[i],LayerElevation[i]+LayerHeight[i]/2), Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], LayerElevation[i]), AirMass[i], AirVapor[i]);

			// Photosynthesis: Consume solarRadiation, water, carbon dioxide, to produce glucose and oxygen
			float carbonDioxideExtractedMax = math.min(carbonDioxide, airDensity * carbonDioxide / AirMass[i] * FloraCarbonDioxideExtractionEfficiency);
			float photosynthesis = 
				FloraPhotosynthesisSpeed * floraWater // water
				* (1.0f - math.pow(math.min(1, glucose * inverseFloraMass), 3)) // glucose saturation
				* math.min(1, energy / FloraEnergyForPhotosynthesis * inverseFloraMass) // energy
				* math.min(1, carbonDioxideExtractedMax * inverseFloraMass); // carbon dioxide 
			// TODO: carbon dioxide extraction should curve to a limit rather than cap

			energy -= photosynthesis * FloraEnergyForPhotosynthesis;
			floraWaterDelta -= photosynthesis;
			carbonDioxideDelta -= photosynthesis;
			glucoseDelta += photosynthesis;
			oxygenDelta += photosynthesis;

			// Respiration: Consume glucose, oxygen, produce water and carbon dioxide
			float oxygenExtractedMax = floraCoverage * OxygenPercent * AirMass[i] * FloraOxygenExtractionEfficiency;
			float respirationExpenditure = math.min(floraWater, math.min(glucose, math.min(oxygenExtractedMax, 
				FloraRespirationSpeed
				* math.max(0, 1.0f + FloraRespirationPerDegree * (FloraTemperature[i] - WorldData.FreezingTemperature))
				* floraWater //water
				* soilFertility
				* glucose * inverseFloraMass // energy (glucose)
				* math.min(1, oxygenExtractedMax * inverseFloraMass)))); // oxygen
			// TODO: oxygen extraction should curve to a limit rather than cap

			float respiration = respirationExpenditure / FloraRespirationSpeed;

			energy += respiration * FloraEnergyForPhotosynthesis;
			glucoseDelta -= respirationExpenditure;
			floraWaterDelta -= respirationExpenditure;
			oxygenDelta -= respirationExpenditure;
			surfaceWaterDelta += 2 * respirationExpenditure;
			carbonDioxideDelta += respirationExpenditure;

			// Growth: Consume glucose, produce plant growth
			float growth = 
				FloraGrowthRate
				* soilFertility
				* floraMass
				* respiration * inverseFloraMass
				* Utils.Sqr(glucose * inverseFloraMass);

			glucoseDelta -= growth;
			floraMassDelta += growth;

			// Death: Produce carbon dioxide, water
			float deathPercent = Utils.Sqr(math.max(0, floraMass - respiration) * inverseFloraMass) * FloraDeathRate;
			float death = deathPercent * floraMass;
			// dump any stored glucose and water, put glucose in the ground
			glucoseDelta -= deathPercent * glucose;
			surfaceWaterDelta += death + deathPercent * (glucose + floraWater);
			floraMassDelta -= death;
			floraWaterDelta -= deathPercent * floraWater;
			floraDeath += death + deathPercent * glucose;

			evapMass = Atmosphere.GetEvaporationMass(AirMass[i], AirPressure[i], AirVapor[i], SurfaceWind[i], FloraTemperature[i], surfaceWaterDelta);
			latentHeatFromFlora = evapMass * WorldData.LatentHeatWaterVapor;
			surfaceWaterDelta -= evapMass;

			latentHeatFromFlora = -evapMass * WorldData.LatentHeatWaterVapor;

#endif

#if !DISABLE_FLORA_WATER_ABSORPTION
			groundWaterConsumed = math.min(GroundWater[i], floraMass * (GroundWater[i] / GroundWaterMax) * math.max(0, 1.0f - floraWater * inverseFloraMass) * FloraWaterConsumptionRate);
			floraWaterDelta += groundWaterConsumed;
#endif

			// Convert leftover solar radiation into heat, since we didn't do it in the energy step
			latentHeatFromFlora += energy;

		}

		EvaporatedWaterMass[i] = evapMass;
		LatentHeatAir[i] += -latentHeatFromAir;
		LatentHeatFlora[i] += latentHeatFromFlora;
		GroundWaterConsumed[i] = groundWaterConsumed;
		CarbonDioxideDelta[i] = carbonDioxideDelta;
		OxygenDelta[i] = oxygenDelta;
		FloraMassDelta[i] = floraMassDelta;
		FloraWaterDelta[i] = floraWaterDelta;
		SurfaceWaterDelta[i] = surfaceWaterDelta;
		FloraGlucoseDelta[i] = glucoseDelta;
		FloraDeath[i] = floraDeath;
	}
}

[BurstCompile]
public struct FluxLavaJob : IJobParallelFor {
	public NativeArray<float> LatentHeat;
	public NativeArray<float> CrystalizedMass;
	public NativeArray<float> LavaEjected;
	public NativeArray<float> DustEjected;
	public NativeArray<float> CrustDelta;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> LavaTemperature;
	[ReadOnly] public NativeArray<float> LavaMass;
	[ReadOnly] public NativeArray<float> MagmaMass;
	[ReadOnly] public NativeArray<float> CrustDepth;
	[ReadOnly] public NativeArray<float> Elevation;
	[ReadOnly] public float LavaCrystalizationTemperature;
	[ReadOnly] public float CrustEruptionDepth;
	[ReadOnly] public float DustPerLavaEjected;
	[ReadOnly] public float LavaEruptionSpeed;
	[ReadOnly] public float MagmaPressureCrustReductionSpeed;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float crystalizedMass = 0;
		float latentHeat = 0;
		float lavaEjected = 0;
		float dustEjected = 0;
		float crustDelta = 0;
		float mass = LavaMass[i];
		float temperature = LavaTemperature[i];
		float crystalizationTempDelta = LavaCrystalizationTemperature - temperature;

		if (crystalizationTempDelta > 0)
		{
			//			float crystalized = math.min(mass, (crystalizationTempDelta * mass * WorldData.SpecificHeatLava) / WorldData.LatentHeatLava);
			crystalizedMass = mass;
		}

#if !DISABLE_ERUPTION
		float magmaPressure = MagmaMass[i] / (WorldData.MassLava * (Elevation[i] - CrustDepth[i] + 20000));
		if (magmaPressure > 1)
		{
			crustDelta = -CrustDepth[i] * math.min(1, SecondsPerTick * (magmaPressure - 1) * MagmaPressureCrustReductionSpeed);
			if (CrustDepth[i] < CrustEruptionDepth)
			{
				lavaEjected = math.min(MagmaMass[i], (magmaPressure - 1.0f) * (1.0f - CrustDepth[i] / CrustEruptionDepth) * LavaEruptionSpeed * SecondsPerTick);
				dustEjected = lavaEjected * DustPerLavaEjected * (1.0f - WaterCoverage[i]);
			}
		}
#endif

		CrystalizedMass[i] = crystalizedMass;
		LatentHeat[i] = latentHeat;
		LavaEjected[i] = lavaEjected;
		DustEjected[i] = dustEjected;
		CrustDelta[i] = crustDelta;
	}
}

[BurstCompile]
public struct FluxTerrainJob : IJobParallelFor {
	public NativeArray<float> SoilRespiration;
	[ReadOnly] public NativeArray<float> SoilCarbon;
	[ReadOnly] public float SoilRespirationSpeed;
	[ReadOnly] public float OxygenPercent;
	public void Execute(int i)
	{
		float soilRespiration = 0;
#if !DISABLE_SOIL_RESPIRATION

		soilRespiration = SoilCarbon[i] * OxygenPercent * SoilRespirationSpeed;

#endif

		SoilRespiration[i] = soilRespiration;
	}
}
