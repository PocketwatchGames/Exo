//#define EnergyAirJobDebug
//#define EnergyWaterSurfaceJobDebug
//#define DISABLE_CORIOLIS

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !EnergyAirJobDebug
[BurstCompile]
#endif
public struct EnergyAirJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Vapor;
	public NativeArray<float3> Wind;
	public NativeArray<float> CondensationGroundMass;
	public NativeArray<float> CondensationCloudMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float3> LastWind;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> WindFriction;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float WindFrictionMultiplier;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float DewPointZero;
	[ReadOnly] public float WaterVaporMassToAirMassAtDewPoint;
	[ReadOnly] public float InverseDewPointTemperatureRange;
	public void Execute(int i)
	{
		float airMass = AirMass[i];
		float conductionDelta =
			+ConductionEnergyIce[i]
			+ ConductionEnergyTerrain[i]
			+ ConductionEnergyWater[i];
		float specificHeat = WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWaterVapor * LastVapor[i];

		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta;
		float temperature = LastTemperature[i];
		float vapor = LastVapor[i];

		float3 lastWind = LastWind[i];

		float3 wind = lastWind;
		wind *= (1.0f - WindFriction[i] * WindFrictionMultiplier);
		wind += PressureGradientForce[i] * SecondsPerTick;

#if !DISABLE_CORIOLIS
//		var windUp = math.dot(Position[i], lastWind) * Position[i];
		var windRight = math.cross(Position[i], lastWind);
		wind += windRight * math.clamp(CoriolisMultiplier[i] * CoriolisTerm * SecondsPerTick, -1, 1);
#endif

		// TODO: deal with byouancy here, but rather than storing a vertical velocity, let's just move to neutral buoyancy every time step

		//float lastWindVertical = WindVertical[i];
		//WindVertical[i] = lastWindVertical;
		//// TODO: this can overshoot
		//WindVertical[i] += Buoyancy[i] /* * SecondsPerTick */;

		float condensationGroundMass = 0;
		float condensationCloudMass = 0;

#if !DISABLE_CONDENSATION
		var relativeHumidity = Atmosphere.GetRelativeHumidity(airMass, vapor, temperature, DewPointZero, WaterVaporMassToAirMassAtDewPoint, InverseDewPointTemperatureRange);
		if (relativeHumidity > 1.0f)
		{
			float aboveCloud = math.saturate((LayerElevation[i] - CloudElevation[i]) / LayerHeight[i]);
			float vaporToCondense = (relativeHumidity - 1.0f) / relativeHumidity * vapor;
			condensationCloudMass = aboveCloud * vaporToCondense;
			condensationGroundMass = (1.0f - aboveCloud) * vaporToCondense;
			vapor -= vaporToCondense;
			energy += vaporToCondense * WorldData.LatentHeatWaterVapor;
		}
#endif
		temperature += energy / specificHeat;

		CondensationGroundMass[i] = condensationGroundMass;
		CondensationCloudMass[i] = condensationCloudMass;

		Temperature[i] = temperature;
		Vapor[i] = vapor;
		Wind[i] = wind;
	}
}


#if !EnergyWaterJobDebug
[BurstCompile]
#endif
public struct EnergyWaterJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> SaltMass;
	public NativeArray<float3> Velocity;
	public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastWaterMass;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float3> LastVelocity;
	[ReadOnly] public NativeArray<float3> Force;
	[ReadOnly] public NativeArray<float> FrozenMass;
	[ReadOnly] public NativeArray<float> EvaporatedMass;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public NativeArray<float> StateFluxEnergy;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float specificHeat = WorldData.SpecificHeatWater * LastWaterMass[i] + WorldData.SpecificHeatSalt * LastSaltMass[i];
		float energySources = 0;
		if (specificHeat > 0)
		{
			energySources = (SolarRadiationIn[i] + ThermalRadiationDelta[i] + ConductionEnergyTerrain[i] + StateFluxEnergy[i]) / specificHeat;
		}
		Mass[i] = LastWaterMass[i] - FrozenMass[i] - EvaporatedMass[i];
		SaltMass[i] = LastSaltMass[i];
		if (Mass[i] > 0)
		{
			Temperature[i] = LastTemperature[i] + energySources;
			Velocity[i] = LastVelocity[i] + Force[i] * SecondsPerTick;
		} else
		{
			Temperature[i] = 0;
			Velocity[i] = 0;
		}
	}
}

#if !EnergyWaterSurfaceJobDebug
[BurstCompile]
#endif
public struct EnergyWaterJobSurface : IJobParallelFor {
	public NativeArray<float> EvaporatedWaterMass;
	public NativeArray<float> FrozenMass;
	public NativeArray<float> StateFluxEnergy;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaTop;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaBottom;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public NativeArray<float> RelativeHumidity;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public float EvaporationRate;
	[ReadOnly] public float EvapTemperatureMin;
	[ReadOnly] public float EvapTemperatureMax;
	[ReadOnly] public float WaterHeatingDepth;
	public void Execute(int i)
	{
		float temperature = LastTemperature[i];
		float waterMass = LastMass[i];
		float saltMass = LastSaltMass[i];
		float energyTop = -ConductionEnergyAir[i] - ConductionEnergyIce[i] + ThermalRadiationDeltaTop[i] + SolarRadiationIn[i];
		float energyBottom = ConductionEnergyTerrain[i] + ThermalRadiationDeltaBottom[i];

		float stateFluxEnergy = 0;
		float evapMass = 0;
		float frozenTopMass = 0;
		float frozenBottomMass = 0;
		if (waterMass > 0)
		{

#if !DISABLE_EVAPORATION
			if (energyTop > 0)
			{
				float evapRate = math.saturate((1.0f - RelativeHumidity[i]) * math.max(0, WaterCoverage[i] - IceCoverage[i]) * EvaporationRate * (LastTemperature[i] - EvapTemperatureMin) / (EvapTemperatureMax - EvapTemperatureMin));
				evapMass = math.clamp(evapRate * energyTop / WorldData.LatentHeatWaterVapor, 0, waterMass);
				if (evapRate > 0)
				{
					float evapEnergy = evapMass * WorldData.LatentHeatWaterVapor;
					energyTop -= evapEnergy;
					stateFluxEnergy -= evapEnergy;
				}
				waterMass -= evapMass;
			}
#endif

#if !DISABLE_FREEZE_TOP
			if (waterMass > 0)
			{
				float heatingMass = math.min(waterMass, WaterHeatingDepth * WorldData.MassWater);
				float specificHeatSaltWater = (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatWater * saltMass) / (waterMass + saltMass);
				float newTempTop = temperature + energyTop / (specificHeatSaltWater * heatingMass);
				if (newTempTop < WorldData.FreezingTemperature)
				{
					float energyToRelease = (WorldData.FreezingTemperature - newTempTop) * specificHeatSaltWater * heatingMass;
					frozenTopMass = math.min(waterMass, energyToRelease / WorldData.LatentHeatWaterLiquid);
					float latentEnergy = frozenTopMass * WorldData.LatentHeatWaterLiquid;
					energyTop += latentEnergy;
					stateFluxEnergy += latentEnergy;
					waterMass -= frozenTopMass;
				}
			}
#endif

#if !DISABLE_FREEZE_BOTTOM
			if (waterMass > 0)
			{
				float heatingMass = math.min(waterMass, WaterHeatingDepth * WorldData.MassWater);
				float specificHeatSaltWater = (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatWater * saltMass) / (waterMass + saltMass);
				float newTempBottom = temperature + energyBottom / (specificHeatSaltWater * heatingMass);
				if (newTempBottom < WorldData.FreezingTemperature)
				{
					float energyToRelease = (WorldData.FreezingTemperature - newTempBottom) * specificHeatSaltWater * heatingMass;
					frozenBottomMass = math.min(waterMass, energyToRelease / WorldData.LatentHeatWaterLiquid);
					float latentEnergy = frozenBottomMass * WorldData.LatentHeatWaterLiquid;
					energyBottom += latentEnergy;
					stateFluxEnergy += latentEnergy;
					waterMass -= frozenBottomMass;
				}
			}
#endif
		}


		EvaporatedWaterMass[i] = evapMass;
		FrozenMass[i] = frozenTopMass + frozenBottomMass;

		StateFluxEnergy[i] = stateFluxEnergy;
	}
}

#if !EnergyCloudJobDebug
[BurstCompile]
#endif
public struct EnergyCloudJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> SurfaceWaterMass;
	public NativeArray<float> SurfaceWaterTemperature;
	public NativeArray<float> IceMass;
	public NativeArray<float> IceTemperature;
	public NativeArray<float> DropletMass;
	public NativeArray<float3> Velocity;
	public NativeArray<float> RainfallWaterMass;
	public NativeArray<float> CloudEvaporationMass;
	[ReadOnly] public NativeArray<float3> LastVelocity;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> DewPoint;
	[ReadOnly] public NativeArray<float> LastDropletMass;
	[ReadOnly] public NativeArray<float> LastCloudMass;
	[ReadOnly] public NativeArray<float> AirMassCloud;
	[ReadOnly] public NativeArray<float> WaterVaporCloud;
	[ReadOnly] public NativeArray<float> AirPressureCloud;
	[ReadOnly] public NativeArray<float> RelativeHumidityCloud;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float> WindFriction;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	[ReadOnly] public NativeArray<float> SurfaceAirTemperature;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float WindFrictionMultiplier;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float RainDropMinSize;
	[ReadOnly] public float RainDropMaxSize;
	[ReadOnly] public float RainDropDragCoefficient;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float TicksPerSecond;
	[ReadOnly] public float CloudDissapationRateDryAir;
	[ReadOnly] public float CloudDissapationRateWind;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float EvapTemperatureMax;
	[ReadOnly] public float EvapTemperatureMin;
	[ReadOnly] public float EvaporationRate;
	public void Execute(int i)
	{
		float surfaceElevation = SurfaceElevation[i];

		float cloudMass = LastCloudMass[i];
		float dropletMass = LastDropletMass[i];
		float cloudElevation = CloudElevation[i];
		float3 lastVelocity = LastVelocity[i];
		float rainfallWaterMass = 0;
		float cloudEvaporationMass = 0;
		float dewPoint = DewPoint[i];

		float3 velocity = lastVelocity * (1.0f - WindFriction[i] * WindFrictionMultiplier) + PressureGradientForce[i] * SecondsPerTick;

		var velocityUp =  velocity * Position[i];
		var velocityRight = math.cross(Position[i], velocity);
		velocity += velocityRight * CoriolisMultiplier[i] * CoriolisTerm * SecondsPerTick;

		float precipitationMass = 0;
		if (cloudElevation <= surfaceElevation)
		{
			precipitationMass += cloudMass;
			cloudMass = 0;
		}

#if !DISABLE_RAINFALL
		if (cloudMass > 0)
		{
			// TODO: improve this somehow
			dropletMass += cloudMass * 0.00000001f;


			float airDensityAtElevation = Atmosphere.GetAirDensity(AirPressureCloud[i], dewPoint, AirMassCloud[i], WaterVaporCloud[i]);
			float waterDensityAtElevation = Atmosphere.GetWaterDensityAtElevation(dewPoint, cloudElevation);
			float rainDropRadius = math.clamp(Atmosphere.GetDropletRadius(dropletMass, waterDensityAtElevation), RainDropMinSize, RainDropMaxSize);
			float rainDropVolume = 4 / 3 * math.PI * math.pow(rainDropRadius, 3);

			// TODO: See wikipedia Terminal velocity:
			//https://en.wikipedia.org/wiki/Terminal_velocity
			// We shouldn't be using the Air's buoyancy force as a stand in for the water droplets buoyancy
			// We should instead just be adding the vertical velocity of the air parcel
			float windVertical = math.dot(lastVelocity, Position[i]);
			float terminalVelocity = windVertical - math.sqrt(8 * rainDropRadius * waterDensityAtElevation * Gravity / (3 * airDensityAtElevation * RainDropDragCoefficient));
			if (terminalVelocity < 0 && dropletMass > 0)
			{

				// TODO: use this to detemine temperature when it hits the ground (snow or rain)
				float rainDropFallTime = (cloudElevation - surfaceElevation) / -terminalVelocity;
				if (rainDropFallTime < SecondsPerTick)
				{
					// TODO: account for drop size variance so we don't just dump it all at once
					rainfallWaterMass = cloudMass * (1.0f - rainDropFallTime / SecondsPerTick);
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
			velocity = 0;
		}

		if (precipitationMass > 0)
		{
			if (SurfaceAirTemperature[i] < WorldData.FreezingTemperature)
			{
				IceTemperature[i] = (IceTemperature[i] * IceMass[i] + SurfaceAirTemperature[i] * precipitationMass) / (IceMass[i] + precipitationMass);
				IceMass[i] += precipitationMass;
			}
			else
			{
				SurfaceWaterTemperature[i] = (SurfaceWaterTemperature[i] * SurfaceWaterMass[i] + SurfaceAirTemperature[i] * precipitationMass) / (SurfaceWaterMass[i] + precipitationMass);
				SurfaceWaterMass[i] += precipitationMass;
			}
		}

		DropletMass[i] = dropletMass;
		CloudMass[i] = cloudMass;
		Velocity[i] = velocity;

		RainfallWaterMass[i] = rainfallWaterMass;
		CloudEvaporationMass[i] = cloudEvaporationMass;
	}
}

#if !EnergyIceJobDebug
[BurstCompile]
#endif
public struct EnergyIceJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Mass;
	public NativeArray<float> MeltedMass;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaBottom;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaTop;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public float IceHeatingDepth;
	public void Execute(int i)
	{
		float energyTop = -ConductionEnergyAir[i] + ThermalRadiationDeltaTop[i] + SolarRadiationIn[i];
		float energyBottom = ConductionEnergyWater[i] + ConductionEnergyTerrain[i] + ThermalRadiationDeltaBottom[i];

		float meltedTopMass = 0;
		float meltedBottomMass = 0;
		float iceMass = LastMass[i];
		float iceTemperature = LastTemperature[i];

#if !DISABLE_MELTING_TOP
		if (iceMass > 0)
		{
			float iceDepth = iceMass / WorldData.MassIce;
			float heatingDepth = math.min(IceHeatingDepth, iceDepth);
			float heatingMass = heatingDepth * WorldData.MassIce;
			float newTempTop = iceTemperature + energyTop / (WorldData.SpecificHeatIce * heatingMass);
			if (newTempTop > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (newTempTop - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * heatingMass;
				meltedTopMass = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				energyTop -= meltedTopMass * WorldData.LatentHeatWaterLiquid;
				iceMass -= meltedTopMass;
			}
		}
#endif

#if !DISABLE_MELTING_BOTTOM
		if (iceMass > 0)
		{
			float iceDepth = iceMass / WorldData.MassIce;
			float heatingDepth = math.min(IceHeatingDepth, iceDepth);
			float heatingMass = heatingDepth * WorldData.MassIce;
			float newTempBottom = iceTemperature + energyBottom / (WorldData.SpecificHeatIce * heatingMass);
			if (newTempBottom > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (newTempBottom - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * heatingMass;
				meltedBottomMass = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				energyBottom -= meltedBottomMass * WorldData.LatentHeatWaterLiquid;
				iceMass -= meltedBottomMass;
			}
		}
#endif

		Mass[i] = iceMass;
		MeltedMass[i] = meltedTopMass + meltedBottomMass;
		Temperature[i] = iceMass > 0 ? iceTemperature + (energyTop + energyBottom) / (iceMass * WorldData.SpecificHeatIce) : 0;
	}
}

#if !EnergyTerrainJobDebug
[BurstCompile]
#endif
public struct EnergyTerrainJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> VegetationCoverage;
	[ReadOnly] public float GeothermalEnergy;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyAir[i]
			- ConductionEnergyWater[i]
			- ConductionEnergyIce[i];
		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta + GeothermalEnergy;
		float specificHeat = Atmosphere.GetSpecificHeatTerrain(HeatingDepth, Terrain[i].SoilFertility, VegetationCoverage[i]);
		Temperature[i] = LastTemperature[i] + energy / specificHeat;
	}
}

