using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;


[BurstCompile]
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
[BurstCompile]
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
				+ precipitationTemperature * rainMass * WorldData.SpecificHeatWater
				+ WorldData.FreezingTemperature * IceMelted[i] * WorldData.SpecificHeatWater;
			float divisor =
				+ rainMass * WorldData.SpecificHeatWater
				+ IceMelted[i] * WorldData.SpecificHeatWater;
			if (remainingWater > 0)
			{
				newTemperature += WaterTemperature[i] * (SaltMass[i] * WorldData.SpecificHeatSalt + remainingWater * WorldData.SpecificHeatWater);
				divisor += (SaltMass[i] * WorldData.SpecificHeatSalt + remainingWater * WorldData.SpecificHeatWater);
			}
			newTemperature /= divisor;
			WaterTemperature[i] = newTemperature;
		}
	}
}


[BurstCompile]
public struct UpdateMassCondensationGroundJob : IJobParallelFor {
	public NativeArray<float> SurfaceWaterMass;
	public NativeArray<float> SurfaceWaterTemperature;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> GroundCondensation;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	[ReadOnly] public NativeArray<float> LayerMiddle;
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
			float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], LayerMiddle[i]);
			SurfaceWaterTemperature[i] = 
				((waterMass * WorldData.SpecificHeatWater + SurfaceSaltMass[i] * WorldData.SpecificHeatSalt) * SurfaceWaterTemperature[i] 
				+ GroundCondensation[i] * WorldData.SpecificHeatWater * airTemperatureAbsolute) 
				/ (newWaterMass * WorldData.SpecificHeatWater + SurfaceSaltMass[i] * WorldData.SpecificHeatSalt);
		}
		SurfaceWaterMass[i] = newWaterMass;
	}
}

[BurstCompile]
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


[BurstCompile]
public struct UpdateMassAirJob : IJobParallelFor {
	public NativeArray<float> VaporMass;
	public NativeArray<float> CarbonDioxideMass;
	public NativeArray<float> DustMass;
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> LastVaporMass;
	[ReadOnly] public NativeArray<float> LastCarbonDioxideMass;
	[ReadOnly] public NativeArray<float> LastDustMass;
	[ReadOnly] public NativeArray<float> CloudEvaporation;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CloudCondensation;
	[ReadOnly] public NativeArray<float> GroundCondensation;
	[ReadOnly] public NativeArray<float> DustFromAbove;
	[ReadOnly] public NativeArray<float> DustFromBelow;
	[ReadOnly] public NativeArray<float> DustUp;
	[ReadOnly] public NativeArray<float> DustDown;
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
		float newDustMass = LastDustMass[i] + DustFromAbove[i] + DustFromBelow[i] - DustDown[i];
		float newCarbonDioxide = LastCarbonDioxideMass[i];
		if (!IsTop)
		{
			newDustMass -= DustUp[i];
		}

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
		DustMass[i] = newDustMass;
		CarbonDioxideMass[i] = newCarbonDioxide;
	}
}

[BurstCompile]
public struct UpdateMassEvaporationJob : IJobParallelFor {
	public NativeArray<float> AirTemperaturePotential;
	public NativeArray<float> VaporMass;
	public NativeArray<float> DustMass;
	[ReadOnly] public NativeArray<float> EvaporationWater;
	[ReadOnly] public NativeArray<float> EvaporationTemperaturePotentialWater;
	[ReadOnly] public NativeArray<float> EvaporationFlora;
	[ReadOnly] public NativeArray<float> EvaporationTemperaturePotentialFlora;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> DustEjected;
	public void Execute(int i)
	{
		float vaporMass = VaporMass[i];
		float airMass = AirMass[i];

		AirTemperaturePotential[i] =
			((airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * AirTemperaturePotential[i] +
			EvaporationWater[i] * EvaporationTemperaturePotentialWater[i] * WorldData.SpecificHeatWaterVapor +
			EvaporationFlora[i] * EvaporationTemperaturePotentialFlora[i] * WorldData.SpecificHeatWaterVapor) /
			(airMass * WorldData.SpecificHeatAtmosphere + (vaporMass + EvaporationWater[i] + EvaporationFlora[i]) * WorldData.SpecificHeatWaterVapor);

		DustMass[i] = DustMass[i] + DustEjected[i];
		VaporMass[i] = vaporMass + EvaporationWater[i] + EvaporationFlora[i];
	}
}

[BurstCompile]
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
		float newIceMass = math.max(0, iceMass + snowMass - IceMelted[i] + WaterFrozen[i]);
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


[BurstCompile]
public struct UpdateTerrainJob : IJobParallelFor {
	public NativeArray<float> SoilFertility;
	public NativeArray<float> Roughness;
	public NativeArray<float> Elevation;
	public NativeArray<float> GroundWater;
	public NativeArray<float> LavaMass;
	public NativeArray<float> LavaTemperature;
	public NativeArray<float> MagmaMass;
	public NativeArray<float> CrustDepth;

	[ReadOnly] public NativeArray<float> LastSoilFertility;
	[ReadOnly] public NativeArray<float> LastRoughness;
	[ReadOnly] public NativeArray<float> LastElevation;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> LastLavaMass;
	[ReadOnly] public NativeArray<float> LastMagmaMass;
	[ReadOnly] public NativeArray<float> LastCrustDepth;
	[ReadOnly] public NativeArray<float> GroundWaterConsumed;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> DustSettled;
	[ReadOnly] public NativeArray<float> LavaCrystalized;
	[ReadOnly] public NativeArray<float> LavaEjected;
	[ReadOnly] public NativeArray<float> CrustDelta;
	[ReadOnly] public float MagmaTemperature;
	[ReadOnly] public float LavaDensityAdjustment;
	public void Execute(int i)
	{
		float lavaCrystalized = LavaCrystalized[i];
		float lavaEjected = LavaEjected[i];
		float lavaMass = LastLavaMass[i];
		float newLavaMass = lavaMass - lavaCrystalized + lavaEjected;

		float lavaCrystalizedDepth = lavaCrystalized / (WorldData.MassLava * LavaDensityAdjustment);
		Elevation[i] = LastElevation[i] + lavaCrystalizedDepth;
		// TODO: improve soil fertility when dust settles
		SoilFertility[i] = LastSoilFertility[i];
		Roughness[i] = LastRoughness[i];
		GroundWater[i] = LastGroundWater[i] - GroundWaterConsumed[i];
		CrustDepth[i] = LastCrustDepth[i] + lavaCrystalizedDepth + CrustDelta[i];
		MagmaMass[i] = LastMagmaMass[i] - lavaEjected;
		LavaMass[i] = newLavaMass;
		if (newLavaMass > 0)
		{
			LavaTemperature[i] = ((lavaMass - lavaCrystalized) * LavaTemperature[i] + lavaEjected * MagmaTemperature) / newLavaMass;
		} else
		{
			LavaTemperature[i] = 0;
		}
	}

}


[BurstCompile]
public struct UpdateFloraJob : IJobParallelFor {
	public NativeArray<float> FloraMass;
	public NativeArray<float> FloraWater;

	[ReadOnly] public NativeArray<float> FloraMassDelta;
	[ReadOnly] public NativeArray<float> EvaporationMass;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastWater;
	[ReadOnly] public NativeArray<float> GroundWaterConsumed;
	public void Execute(int i)
	{
		FloraMass[i] = math.max(0, LastMass[i] + FloraMassDelta[i]);
		FloraWater[i] = math.max(0, LastWater[i] - EvaporationMass[i] + GroundWaterConsumed[i]);
	}

}



