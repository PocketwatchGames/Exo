//#define StateChangeJobDebug
//#define StateChangeAirLayerJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

#if !StateChangeJobDebug
[BurstCompile]
#endif
public struct StateChangeJob : IJobParallelFor {
	public NativeArray<float> LowerAirTemperaturePotential;
	public NativeArray<float> LowerAirVapor;
	public NativeArray<float> SurfaceWaterMass;
	public NativeArray<float> SurfaceWaterTemperature;
	public NativeArray<float> IceMass;
	public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> LastIceMass;
	[ReadOnly] public NativeArray<float> LastIceTemperature;
	[ReadOnly] public NativeArray<float> IceMeltedMass;
	[ReadOnly] public NativeArray<float> WaterEvaporatedMass;
	[ReadOnly] public NativeArray<float> WaterFrozenMass;
	[ReadOnly] public NativeArray<float> WaterFrozenTemperature;
	[ReadOnly] public NativeArray<float> PrecipitationMass;
	[ReadOnly] public NativeArray<float> PrecipitationTemperature;
	[ReadOnly] public NativeArray<float> SurfaceAirMass;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	[ReadOnly] public NativeArray<float> LowerAirLayerElevation;
	[ReadOnly] public NativeArray<float> LowerAirLayerHeight;
	public void Execute(int i)
	{
		float vaporMass = LowerAirVapor[i];
		float airMass = SurfaceAirMass[i];
		float evaporatedMass = WaterEvaporatedMass[i];
		float waterMass = SurfaceWaterMass[i];
		float waterTemperature = SurfaceWaterTemperature[i];
		float iceMeltedMass = IceMeltedMass[i];
		float waterFrozenMass = WaterFrozenMass[i];
		float saltMass = SurfaceSaltMass[i];
		float precipitationTemperature = PrecipitationTemperature[i];

		float rainfallMass = 0;
		float snowMass = 0;
		if (precipitationTemperature >= WorldData.FreezingTemperature)
		{
			rainfallMass = PrecipitationMass[i];
		}
		else
		{
			snowMass = PrecipitationMass[i];
		}

		float newWaterMass = waterMass + iceMeltedMass + rainfallMass - evaporatedMass - waterFrozenMass;
		float newVaporMass = vaporMass + evaporatedMass;
		float newIceMass = LastIceMass[i] + waterFrozenMass + snowMass - iceMeltedMass;
		float elevation = LowerAirLayerElevation[i] + LowerAirLayerHeight[i] / 2;
		float newAirTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(LowerAirTemperaturePotential[i], elevation);

		float specificHeatAir = WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWater * vaporMass;
		newAirTemperatureAbsolute = (
			newAirTemperatureAbsolute * (airMass + vaporMass) 
			+ evaporatedMass * waterTemperature)
			/ (airMass + vaporMass + evaporatedMass);

		float newWaterTemperature;
		if (newWaterMass > 0)
		{
			newWaterTemperature = (
				waterTemperature * (saltMass + waterMass - evaporatedMass)
				+ precipitationTemperature * rainfallMass
				+ WorldData.FreezingTemperature * iceMeltedMass)
				/ (waterMass + saltMass - evaporatedMass + rainfallMass + iceMeltedMass);
		}
		else
		{
			newWaterTemperature = 0;
		}

		float newIceTemperature;
		if (newIceMass > 0)
		{
			newIceTemperature = (
				LastIceTemperature[i] * (LastIceMass[i] - iceMeltedMass)
				+ WaterFrozenTemperature[i] * waterFrozenMass
				+ precipitationTemperature * snowMass
				) / newIceMass;
		}
		else
		{
			newIceTemperature = 0;
		}

		SurfaceWaterTemperature[i] = newWaterTemperature;
		LowerAirTemperaturePotential[i] = Atmosphere.GetPotentialTemperature(newAirTemperatureAbsolute, elevation);
		SurfaceWaterMass[i] = newWaterMass;
		LowerAirVapor[i] = newVaporMass;
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
	public NativeArray<float> AirTemperaturePotential;
	public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> LastAirTemperaturePotential;
	[ReadOnly] public NativeArray<float> LastVaporMass;
	[ReadOnly] public NativeArray<float> GroundCondensationMass;
	[ReadOnly] public NativeArray<float> CloudCondensationMass;
	[ReadOnly] public NativeArray<float> CloudEvaporationMass;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float> AirMass;
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


		float newCloudMass = math.max(0, cloudMass + cloudCondensationMass - cloudEvaporationMass);

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
		VaporMass[i] = LastVaporMass[i] - cloudCondensationMass - groundCondensationMass + cloudEvaporationMass;
		AirTemperaturePotential[i] = LastAirTemperaturePotential[i];

		float newWaterMass = SurfaceWaterMass[i] + groundCondensationMass;
		if (newWaterMass == 0)
		{
			SurfaceWaterMass[i] = 0;
			SurfaceWaterTemperature[i] = 0;
		}
		else
		{
			float surfaceAirTemperature = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], SurfaceElevation[i]);
			SurfaceWaterTemperature[i] = ((SurfaceWaterMass[i] + SurfaceSaltMass[i]) * SurfaceWaterTemperature[i] + groundCondensationMass * surfaceAirTemperature) / (newWaterMass + SurfaceSaltMass[i]);
			SurfaceWaterMass[i] = newWaterMass;
		}
	}
}
