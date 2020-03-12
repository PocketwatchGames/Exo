//#define StateChangeAirLayerJobDebug
#define UpdateMassEvaporationJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;



#if !UpdateMassWaterJobDebug
[BurstCompile]
#endif
public struct UpdateMassWaterJob : IJobParallelFor {
	public NativeArray<float> WaterMass;
	public NativeArray<float> SaltMass;
	public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> SaltPlume;
	[ReadOnly] public NativeArray<float> SaltPlumeTemperature;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float> LastWaterMass;
	[ReadOnly] public NativeArray<float> DownLastWaterMass;
	public void Execute(int i)
	{
		WaterMass[i] = LastWaterMass[i];
		SaltMass[i] = LastSaltMass[i];
		if (DownLastWaterMass[i] == 0 && LastWaterMass[i] > 0)
		{
			WaterTemperature[i] = (WaterTemperature[i] * (LastWaterMass[i] * WorldData.SpecificHeatWater + LastSaltMass[i] * WorldData.SpecificHeatSalt) + SaltPlumeTemperature[i] * SaltPlume[i] * WorldData.SpecificHeatSalt) / (LastWaterMass[i] * WorldData.SpecificHeatWater + (LastSaltMass[i] + SaltPlume[i]) * WorldData.SpecificHeatSalt);
			SaltMass[i] += SaltPlume[i];
		}
	}
}
#if !UpdateMassWaterSurfaceJobDebug
[BurstCompile]
#endif
public struct UpdateMassWaterSurfaceJob : IJobParallelFor {
	public NativeArray<float> WaterTemperature;
	public NativeArray<float> WaterMass;
	public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Evaporation;
	[ReadOnly] public NativeArray<float> IceMelted;
	[ReadOnly] public NativeArray<float> Precipitation;
	[ReadOnly] public NativeArray<float> PrecipitationTemperature;
	[ReadOnly] public NativeArray<float> WaterFrozen;
	[ReadOnly] public NativeArray<float> SaltPlume;
	public void Execute(int i)
	{
		float precipitationTemperature = PrecipitationTemperature[i];
		float rainMass = precipitationTemperature > WorldData.FreezingTemperature ? Precipitation[i] : 0;

		float waterMass = WaterMass[i];
		float newMass = waterMass + IceMelted[i] + rainMass - Evaporation[i] - WaterFrozen[i];
		WaterMass[i] = newMass;
		SaltMass[i] -= SaltPlume[i];

		if (newMass <= 0)
		{
			WaterTemperature[i] = 0;
		}
		else
		{
			float remainingWater = waterMass - Evaporation[i] - WaterFrozen[i];
			float newTemperature = 
				WaterTemperature[i] * (SaltMass[i] * WorldData.SpecificHeatSalt + remainingWater * WorldData.SpecificHeatWater) 
				+ precipitationTemperature * rainMass * WorldData.SpecificHeatWater
				+ WorldData.FreezingTemperature * IceMelted[i] * WorldData.SpecificHeatWater;
			newTemperature /= (newMass * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt);
		}
	}
}


#if !UpdateMassAirJobDebug
[BurstCompile]
#endif
public struct UpdateMassCondensationGroundJob : IJobParallelFor {
	public NativeArray<float> SurfaceWaterMass;
	public NativeArray<float> SurfaceWaterTemperature;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> GroundCondensation;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	public void Execute(int i)
	{
		float waterMass = SurfaceWaterMass[i];
		float newWaterMass = waterMass + GroundCondensation[i];
		if (newWaterMass <= 0)
		{
			SurfaceWaterTemperature[i] = 0;
		}
		else
		{
			float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], LayerElevation[i] + LayerHeight[i] / 2);
			SurfaceWaterTemperature[i] = 
				((waterMass * WorldData.SpecificHeatWater + SurfaceSaltMass[i] * WorldData.SpecificHeatSalt) * SurfaceWaterTemperature[i] 
				+ GroundCondensation[i] * WorldData.SpecificHeatWater * airTemperatureAbsolute) 
				/ (newWaterMass * WorldData.SpecificHeatWater + SurfaceSaltMass[i] * WorldData.SpecificHeatSalt);
		}
		SurfaceWaterMass[i] = newWaterMass;
	}
}

#if !UpdateMassCloudJobDebug
[BurstCompile]
#endif
public struct UpdateMassCloudJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> LastCloudMass;
	[ReadOnly] public NativeArray<float> LastDropletMass;
	[ReadOnly] public NativeArray<float> CloudEvaporation;
	[ReadOnly] public NativeArray<float> PrecipitationMass;
	public void Execute(int i)
	{
		CloudMass[i] = LastCloudMass[i] - CloudEvaporation[i] - PrecipitationMass[i];
		CloudDropletMass[i] = LastDropletMass[i];
	}
}


#if !UpdateMassAirJobDebug
[BurstCompile]
#endif
public struct UpdateMassAirJob : IJobParallelFor {
	public NativeArray<float> VaporMass;
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> LastVaporMass;
	[ReadOnly] public NativeArray<float> CloudEvaporation;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CloudCondensation;
	[ReadOnly] public NativeArray<float> GroundCondensation;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	public void Execute(int i)
	{

		float cloudMass = CloudMass[i];
		float cloudEvaporationInLayer = 0;
		if ((CloudElevation[i] >= LayerElevation[i] || IsBottom) && (CloudElevation[i] < LayerElevation[i] + LayerHeight[i] || IsTop))
		{
			cloudEvaporationInLayer = CloudEvaporation[i];
		}
		float newVaporMass = LastVaporMass[i] - CloudCondensation[i] - GroundCondensation[i];
		float newCloudMass = cloudMass + CloudCondensation[i];

		float newDropletSize = 0;
		if (newCloudMass > 0)
		{
			newDropletSize = CloudDropletMass[i] * cloudMass / (cloudMass + CloudCondensation[i]);
		}
		else
		{
			newDropletSize = 0;
		}


		CloudMass[i] = newCloudMass;
		CloudDropletMass[i] = newDropletSize;
		VaporMass[i] = newVaporMass;
	}
}

#if !UpdateMassEvaporationJobDebug
[BurstCompile]
#endif
public struct UpdateMassEvaporationJob : IJobParallelFor {
	public NativeArray<float> AirTemperaturePotential;
	public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> Evaporation;
	[ReadOnly] public NativeArray<float> EvaporationTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
 	public void Execute(int i)
	{
		float vaporMass = VaporMass[i];
		float airMass = AirMass[i];
		float elevation = LayerElevation[i] + LayerHeight[i] / 2;

		AirTemperaturePotential[i] =
			((airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * AirTemperaturePotential[i] +
			Evaporation[i] * EvaporationTemperaturePotential[i] * WorldData.SpecificHeatWaterVapor) /
			(airMass * WorldData.SpecificHeatAtmosphere + (vaporMass + Evaporation[i]) * WorldData.SpecificHeatWaterVapor);

		VaporMass[i] = vaporMass + Evaporation[i];
	}
}

#if !UpdateMassIceJobDebug
[BurstCompile]
#endif
public struct UpdateMassIceJob : IJobParallelFor {
	public NativeArray<float> IceMass;
	public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> LastIceMass;
	[ReadOnly] public NativeArray<float> IceMelted;
	[ReadOnly] public NativeArray<float> Precipitation;
	[ReadOnly] public NativeArray<float> PrecipitationTemperature;
	[ReadOnly] public NativeArray<float> WaterFrozen;
	[ReadOnly] public NativeArray<float> WaterTemperature;
	public void Execute(int i)
	{
		float iceMass = LastIceMass[i];
		float precipitationTemperature = PrecipitationTemperature[i];
		float snowMass = precipitationTemperature <= WorldData.FreezingTemperature ? Precipitation[i] : 0;
		float newIceMass = iceMass + snowMass - IceMelted[i] + WaterFrozen[i];
		if (newIceMass > 0)
		{
			IceTemperature[i] = (
				IceTemperature[i] * (iceMass - IceMelted[i])
				+ precipitationTemperature * snowMass
				+ WaterTemperature[i] * WaterFrozen[i]
				) / newIceMass;
		} else
		{
			IceTemperature[i] = 0;
		}
		IceMass[i] = newIceMass;
	}
}


#if !ApplyLatentHeatIceJobDebug
[BurstCompile]
#endif
public struct ApplyLatentHeatIceJob : IJobParallelFor {
	public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> LatentHeat;
	public void Execute(int i)
	{
		if (IceMass[i] > 0)
		{
			IceTemperature[i] += LatentHeat[i] / (WorldData.SpecificHeatIce * IceMass[i]);
		}
	}
}

#if !ApplyLatentHeatAirJobDebug
[BurstCompile]
#endif
public struct ApplyLatentHeatAirJob : IJobParallelFor {
	public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> LatentHeat;
	public void Execute(int i)
	{
		AirTemperaturePotential[i] += LatentHeat[i] / (AirMass[i] * WorldData.SpecificHeatAtmosphere + VaporMass[i] * WorldData.SpecificHeatWaterVapor);
	}
}

#if !ApplyLatentHeatAirJobDebug
[BurstCompile]
#endif
public struct ApplyLatentHeatWaterJob : IJobParallelFor {
	public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> LatentHeat;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			WaterTemperature[i] += LatentHeat[i] / (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt);
		}
	}
}

