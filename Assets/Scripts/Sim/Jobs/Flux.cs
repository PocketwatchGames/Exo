//#define DISABLE_RAINFALL
//#define FluxCloudJobDebug
//#define FluxWaterJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

#if !FluxWaterJobDebug
[BurstCompile]
#endif
public struct FluxWaterJob : IJobParallelFor {
	public NativeArray<float> EvaporatedWaterMass;
	public NativeArray<float> FrozenMass;
	public NativeArray<float> FrozenTemperature;
	public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaTop;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaBottom;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public NativeArray<float> SurfaceAirMass;
	[ReadOnly] public NativeArray<float> SurfaceVaporMass;
	[ReadOnly] public NativeArray<float> SurfaceAirPressure;
	[ReadOnly] public NativeArray<float3> SurfaceWind;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> Salinity;
	[ReadOnly] public float WaterHeatingDepth;
	[ReadOnly] public float FreezePointReductionPerSalinity;
	public void Execute(int i)
	{
		float temperature = LastTemperature[i];
		float waterMass = LastMass[i];
		float saltMass = LastSaltMass[i];
		float energyTop = -ConductionEnergyAir[i] - ConductionEnergyIce[i] + ThermalRadiationDeltaTop[i] + SolarRadiationIn[i];
		float energyBottom = ConductionEnergyTerrain[i] + ThermalRadiationDeltaBottom[i];

		float evapMass = 0;
		float frozenTopMass = 0;
		float frozenBottomMass = 0;
		float freezingTemperature = Atmosphere.GetFreezingPoint(Salinity[i], FreezePointReductionPerSalinity);

		if (waterMass > 0)
		{
			float specificHeatSaltWater = (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatWater * saltMass) / (waterMass + saltMass);
			float heatingMass = math.min(waterMass, WaterHeatingDepth * WorldData.MassWater);
			float newTempTop = temperature + energyTop / (specificHeatSaltWater * heatingMass);

#if !DISABLE_EVAPORATION
			// evap formula from here:
			// https://www.engineeringtoolbox.com/evaporation-water-surface-d_690.html
			float evaporationCoefficient = 25 + 19 * math.length(SurfaceWind[i]);
			float maxVaporAtWaterTemperature = Atmosphere.GetMaxVaporAtTemperature(SurfaceAirMass[i], newTempTop, SurfaceAirPressure[i]);
			evapMass = math.max(0, evaporationCoefficient * (maxVaporAtWaterTemperature - SurfaceVaporMass[i]) / SurfaceAirMass[i]);

			if (evapMass > 0)
			{
				energyTop -= evapMass * WorldData.LatentHeatWaterVapor;
				waterMass -= evapMass;
			}
#endif

#if !DISABLE_FREEZE_TOP
			if (waterMass > 0)
			{
				if (newTempTop < freezingTemperature)
				{
					float energyToRelease = (freezingTemperature - newTempTop) * specificHeatSaltWater * heatingMass;
					frozenTopMass = math.min(waterMass, energyToRelease / WorldData.LatentHeatWaterLiquid);
					energyTop += frozenTopMass * WorldData.LatentHeatWaterLiquid;
					waterMass -= frozenTopMass;
				}
			}
#endif

#if !DISABLE_FREEZE_BOTTOM
			if (waterMass > 0)
			{
				float newTempBottom = temperature + energyBottom / (specificHeatSaltWater * heatingMass);
				if (newTempBottom < freezingTemperature)
				{
					float energyToRelease = (freezingTemperature - newTempBottom) * specificHeatSaltWater * heatingMass;
					frozenBottomMass = math.min(waterMass, energyToRelease / WorldData.LatentHeatWaterLiquid);
					energyBottom += frozenBottomMass * WorldData.LatentHeatWaterLiquid;
					waterMass -= frozenBottomMass;
				}
			}
#endif
		}


		EvaporatedWaterMass[i] = evapMass;
		FrozenMass[i] = frozenTopMass + frozenBottomMass;
		FrozenTemperature[i] = freezingTemperature;
		Energy[i] = energyTop + energyBottom;
	}
}

#if !FluxCloudJobDebug
[BurstCompile]
#endif
public struct FluxCloudJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> DropletMass;
	public NativeArray<float3> Velocity;
	public NativeArray<float> PrecipitationMass;
	public NativeArray<float> PrecipitationTemperature;
	public NativeArray<float> EvaporationMass;
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
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	[ReadOnly] public NativeArray<float> SurfaceAirTemperature;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
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
		float surfaceElevation = SurfaceElevation[i];

		float cloudMass = LastCloudMass[i];
		float dropletMass = LastDropletMass[i];
		float cloudElevation = CloudElevation[i];
		float3 lastVelocity = LastVelocity[i];
		float rainfallWaterMass = 0;
		float cloudEvaporationMass = 0;
		float dewPoint = DewPoint[i];

		// I already multiplied PressureGradientForce by SecondsPerTick during the calculation of the force
		float3 velocity = lastVelocity + PressureGradientForce[i];

	//	var velocityUp = math.dot(velocity, Position[i]) * Position[i];

		float precipitationMass = 0;
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
			PrecipitationTemperature[i] = SurfaceAirTemperature[i];
			PrecipitationMass[i] = precipitationMass;
		}

		DropletMass[i] = dropletMass;
		CloudMass[i] = cloudMass;
		Velocity[i] = velocity;

		EvaporationMass[i] = cloudEvaporationMass;
	}
}


#if !FluxAirJobDebug
[BurstCompile]
#endif
public struct FluxAirJob : IJobParallelFor {
	public NativeArray<float> FluxEnergy;
	public NativeArray<float> CondensationGroundMass;
	public NativeArray<float> CondensationCloudMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastTemperaturePotential;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		float condensationGroundMass = 0;
		float condensationCloudMass = 0;

		float energy =
			SolarRadiationIn[i] 
			+ ThermalRadiationDelta[i]
			+ ConductionEnergyIce[i]
			+ ConductionEnergyTerrain[i]
			+ ConductionEnergyWater[i];

#if !DISABLE_CONDENSATION
		float temperatureAbsolute = Atmosphere.GetAbsoluteTemperature(LastTemperaturePotential[i], LayerElevation[i] + LayerHeight[i] / 2);
		var relativeHumidity = Atmosphere.GetRelativeHumidity(AirMass[i], LastVapor[i], temperatureAbsolute, AirPressure[i]);
		if (relativeHumidity > 1.0f)
		{
			float aboveCloud = math.saturate((LayerElevation[i] - CloudElevation[i]) / LayerHeight[i]);
			float vaporToCondense = (relativeHumidity - 1.0f) / relativeHumidity * LastVapor[i];
			condensationCloudMass = aboveCloud * vaporToCondense;
			condensationGroundMass = (1.0f - aboveCloud) * vaporToCondense;
			energy += vaporToCondense * WorldData.LatentHeatWaterVapor;
		}
#endif

		FluxEnergy[i] = energy;
		CondensationGroundMass[i] = condensationGroundMass;
		CondensationCloudMass[i] = condensationCloudMass;
	}
}



#if !FluxIceJobDebug
[BurstCompile]
#endif
public struct FluxIceJob : IJobParallelFor {
	public NativeArray<float> Energy;
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
				float latentEnergy = meltedBottomMass * WorldData.LatentHeatWaterLiquid;
				energyBottom -= latentEnergy;
				iceMass -= meltedBottomMass;
			}
		}
#endif

		Energy[i] = energyTop + energyBottom;
		MeltedMass[i] = meltedTopMass + meltedBottomMass;
	}
}

