#define DISABLE_FREEZE_TOP
#define DISABLE_EVAPORATION
#define DISABLE_CONDENSATION
#define DISABLE_CLOUD_DISSAPATION
#define DISABLE_RAINFALL
#define DISABLE_MELTING_TOP
#define DISABLE_MELTING_BOTTOM

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
	public NativeArray<float> EvaporatedWaterTemperaturePotential;
	public NativeArray<float> FrozenMass;
	public NativeArray<float> FrozenTemperature;
	public NativeArray<float> LatentHeatWater;
	public NativeArray<float> LatentHeatAir;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirVapor;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float3> SurfaceWind;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public float WaterHeatingDepth;
	[ReadOnly] public float FreezePointReductionPerSalinity;
	public void Execute(int i)
	{
		float temperature = Temperature[i];
		float waterMass = LastMass[i];
		float saltMass = LastSaltMass[i];

		float latentHeatFromAir = 0;
		float evapMass = 0;
		float frozenMass = 0;
		float freezingTemperature = 0;
		float energyFlux = 0;
		float evapTemperaturePotential = 0;

		if (waterMass > 0)
		{

			float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], LayerElevation[i] + LayerHeight[i] / 2);

			#if !DISABLE_EVAPORATION
			// evap formula from here:
			// https://www.engineeringtoolbox.com/evaporation-water-surface-d_690.html
			// NOTE: I've made adjustments to this because my finite differencing sometimes means that the water surface and air temperature are a bit out of sync
			// so i'm using the air temperature instead of the water temperature, which means the the formula just reduces to (1-RH)*WindCoefficient
			float evaporationCoefficient = 25 + 19 * math.length(SurfaceWind[i]);
			evapMass = math.clamp(evaporationCoefficient * (Atmosphere.GetMaxVaporAtTemperature(AirMass[i], airTemperatureAbsolute, AirPressure[i]) - AirVapor[i]) / AirMass[i], 0, waterMass);
			waterMass -= evapMass;
			latentHeatFromAir = evapMass * WorldData.LatentHeatWaterVapor;
			evapTemperaturePotential = Atmosphere.GetPotentialTemperature(temperature, LayerElevation[i]);
			//energyTop -= evapMass * WorldData.LatentHeatWaterVapor;

#endif

			float specificHeatSaltWater = (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatWater * saltMass) / (waterMass + saltMass);
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
					energyFlux += frozenMass * WorldData.LatentHeatWaterLiquid;
					waterMass -= frozenMass;
				}
			}
#endif
		}


		EvaporatedWaterTemperaturePotential[i] = evapTemperaturePotential;
		EvaporatedWaterMass[i] = evapMass;
		FrozenMass[i] = frozenMass;
		FrozenTemperature[i] = freezingTemperature;
		LatentHeatWater[i] = energyFlux;
		LatentHeatAir[i] += -latentHeatFromAir;

	}
}

#if !FluxCloudJobDebug
[BurstCompile]
#endif
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
	[ReadOnly] public NativeArray<float> SurfaceLayerHeight;
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

				// TODO: use this to detemine temperature when it hits the ground (snow or rain)
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
			PrecipitationTemperature[i] = Atmosphere.GetAbsoluteTemperature(SurfaceAirTemperaturePotential[i], SurfaceLayerElevation[i] + SurfaceLayerHeight[i] / 2);
		}
		PrecipitationMass[i] = precipitationMass;
		DropletDelta[i] = dropletMass - LastDropletMass[i];
		EvaporationMass[i] = cloudEvaporationMass;
	}
}


#if !FluxAirJobDebug
[BurstCompile]
#endif
public struct FluxAirJob : IJobParallelFor {
	public NativeArray<float> LatentHeat;
	public NativeArray<float> CondensationGroundMass;
	public NativeArray<float> CondensationCloudMass;
	[ReadOnly] public NativeArray<float> TemperaturePotential;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	public void Execute(int i)
	{
		float condensationGroundMass = 0;
		float condensationCloudMass = 0;
		float energyFlux = 0;

#if !DISABLE_CONDENSATION
		float temperatureAbsolute = Atmosphere.GetAbsoluteTemperature(TemperaturePotential[i], LayerElevation[i] + LayerHeight[i] / 2);
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

		LatentHeat[i] = energyFlux;
		CondensationGroundMass[i] = condensationGroundMass;
		CondensationCloudMass[i] = condensationCloudMass;
	}
}



#if !FluxIceJobDebug
[BurstCompile]
#endif
public struct FluxIceJob : IJobParallelFor {
	public NativeArray<float> LatentHeatAir;
	public NativeArray<float> LatentHeatWater;
	public NativeArray<float> MeltedMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public float IceHeatingDepth;
	public void Execute(int i)
	{
		float meltedTopMass = 0;
		float meltedBottomMass = 0;
		float iceMass = LastMass[i];
		float iceTemperature = Temperature[i];

#if !DISABLE_MELTING_TOP
		if (iceMass > 0)
		{
			float iceDepth = iceMass / WorldData.MassIce;
			float heatingDepth = math.min(IceHeatingDepth, iceDepth);
			float heatingMass = heatingDepth * WorldData.MassIce;
			float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], LayerElevation[i]);
			if (airTemperatureAbsolute > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (airTemperatureAbsolute - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * heatingMass;
				meltedTopMass = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				LatentHeatAir[i] -= meltedTopMass * WorldData.LatentHeatWaterLiquid;
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
			float waterTemperature = WaterTemperature[i];
			if (waterTemperature > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (WaterTemperature[i] - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * heatingMass;
				meltedBottomMass = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				float latentEnergy = meltedBottomMass * WorldData.LatentHeatWaterLiquid;
				LatentHeatWater[i] -= latentEnergy;
				iceMass -= meltedBottomMass;
			}
		}
#endif

		MeltedMass[i] = meltedTopMass + meltedBottomMass;
	}
}

