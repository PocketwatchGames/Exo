using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

#if !StateChangeJobDebug
[BurstCompile]
#endif
public struct StateChangeJob : IJobParallelFor {
	public NativeArray<float> SurfaceAirTemperature;
	public NativeArray<float> SurfaceAirVapor;
	public NativeArray<float> SurfaceWaterMass;
	public NativeArray<float> SurfaceWaterTemperature;
	public NativeArray<float> IceMass;
	public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> IceMeltedMass;
	[ReadOnly] public NativeArray<float> WaterEvaporatedMass;
	[ReadOnly] public NativeArray<float> WaterFrozenMass;
	[ReadOnly] public NativeArray<float> RainfallWaterMass;
	[ReadOnly] public NativeArray<float> RainfallTemperature;
	[ReadOnly] public NativeArray<float> SurfaceAirMass;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	public void Execute(int i)
	{
		float vaporMass = SurfaceAirVapor[i];
		float airMass = SurfaceAirMass[i];
		float evaporatedMass = WaterEvaporatedMass[i];
		float waterMass = SurfaceWaterMass[i];
		float rainfallMass = RainfallWaterMass[i];
		float waterTemperature = SurfaceWaterTemperature[i];
		float rainfallTemperature = RainfallTemperature[i];
		float airTemperature = SurfaceAirTemperature[i];
		float iceMeltedMass = IceMeltedMass[i];
		float waterFrozenMass = WaterFrozenMass[i];
		float iceMass = IceMass[i];
		float saltMass = SurfaceSaltMass[i];

		float newWaterMass = waterMass + rainfallMass + iceMeltedMass;
		float newVaporMass = vaporMass + evaporatedMass;
		float newAirTemperature = airTemperature;

		float specificHeatAir = WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWater * vaporMass;
		newAirTemperature = (newAirTemperature * (airMass + vaporMass) + evaporatedMass * waterTemperature) / (airMass + vaporMass + evaporatedMass);

		float newWaterTemperature;
		if (newWaterMass > 0)
		{
			newWaterTemperature = (
				waterTemperature * (saltMass + waterMass)
				+ rainfallTemperature * rainfallMass
				+ WorldData.FreezingTemperature * iceMeltedMass)
				/ (waterMass + saltMass + rainfallMass + iceMeltedMass);
		}
		else
		{
			newWaterTemperature = 0;
		}

		float newIceMass = iceMass + waterFrozenMass;
		float newIceTemperature;
		if (newIceMass > 0)
		{
			newIceTemperature = (IceTemperature[i] * iceMass + WorldData.FreezingTemperature * waterFrozenMass) / newIceMass;
		}
		else
		{
			newIceTemperature = 0;
		}

		SurfaceWaterTemperature[i] = newWaterTemperature;
		SurfaceAirTemperature[i] = newAirTemperature;
		SurfaceWaterMass[i] = newWaterMass;
		SurfaceAirVapor[i] = newVaporMass;
		IceTemperature[i] = newIceTemperature;
		IceMass[i] = newIceMass;


	}
}

#if !StateChangeAirLayerJobDebug
[BurstCompile]
#endif
public struct StateChangeAirLayerJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudDropletMass;
	public NativeArray<float> SurfaceWaterMass;
	public NativeArray<float> SurfaceWaterTemperature;
	public NativeArray<float> AirTemperature;
	[ReadOnly] public NativeArray<float> GroundCondensationMass;
	[ReadOnly] public NativeArray<float> CloudCondensationMass;
	[ReadOnly] public NativeArray<float> CloudEvaporationMass;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	// TODO: implement
	//[ReadOnly] public NativeArray<float> AirTemperatureTop;
	//[ReadOnly] public NativeArray<float> AirTemperatureBottom;
	[ReadOnly] public int LayerIndex;
	public void Execute(int i)
	{

		//float layerElevation = LayerElevation[LayerIndex];
		//float layerHeight = LayerHeight[LayerIndex];
		//bool isCloudLayer = newCloudElevation >= layerElevation && newCloudElevation < layerElevation + layerHeight;

		float cloudMass = CloudMass[i];
		float cloudDropletMass = CloudDropletMass[i];
		float groundCondensationMass = GroundCondensationMass[i];
		float cloudCondensationMass = CloudCondensationMass[i];
		float cloudEvaporationMass = CloudEvaporationMass[i];
		float airTemperature = AirTemperature[i];


		float cloudMassDelta = math.max(-cloudMass, cloudCondensationMass - cloudEvaporationMass);
		float newCloudMass = cloudMass + cloudMassDelta;

		float newDropletSize = 0;
		if (newCloudMass > 0)
		{
			newDropletSize = cloudDropletMass * (cloudMass - cloudEvaporationMass) / (cloudMass + cloudCondensationMass);
		}
		else
		{
			newDropletSize = 0;
		}

		CloudMass[i] = newCloudMass;
		CloudDropletMass[i] = newDropletSize;

		float newWaterMass = SurfaceWaterMass[i] + groundCondensationMass;
		if (newWaterMass == 0)
		{
			SurfaceWaterMass[i] = 0;
			SurfaceWaterTemperature[i] = 0;
		}
		else
		{
			SurfaceWaterTemperature[i] = ((SurfaceWaterMass[i] + SurfaceSaltMass[i]) * SurfaceWaterTemperature[i] + groundCondensationMass * AirTemperature[i]) / (newWaterMass + SurfaceSaltMass[i]);
			SurfaceWaterMass[i] = newWaterMass;
		}
	}
}
