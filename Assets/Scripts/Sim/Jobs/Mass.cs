using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;


[BurstCompile]
public struct UpdateMassWaterJob : IJobParallelFor {
	public NativeArray<float> WaterMass;
	public NativeArray<float> SaltMass;
	public NativeArray<float> CarbonMass;
	public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> SoilRespiration;
	[ReadOnly] public NativeArray<float> SaltPlume;
	[ReadOnly] public NativeArray<float> SaltPlumeTemperature;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float> LastWaterMass;
	[ReadOnly] public NativeArray<float> LastCarbonMass;
	[ReadOnly] public NativeArray<float> DownLastWaterMass;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverageBelow;
	public void Execute(int i)
	{
		WaterMass[i] = LastWaterMass[i];
		SaltMass[i] = LastSaltMass[i];
		CarbonMass[i] = LastCarbonMass[i] + SoilRespiration[i] * math.max(0, WaterCoverage[i] - WaterCoverageBelow[i]);
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
	public NativeArray<float> PlanktonMass;
	public NativeArray<float> PlanktonGlucose;
	public NativeArray<float> CarbonMass;
	[ReadOnly] public NativeArray<float> LastPlanktonMass;
	[ReadOnly] public NativeArray<float> LastPlanktonGlucose;
	[ReadOnly] public NativeArray<float> Evaporation;
	[ReadOnly] public NativeArray<float> IceMelted;
	[ReadOnly] public NativeArray<float> Precipitation;
	[ReadOnly] public NativeArray<float> PrecipitationTemperature;
	[ReadOnly] public NativeArray<float> WaterFrozen;
	[ReadOnly] public NativeArray<float> SaltPlume;
	[ReadOnly] public NativeArray<float> FloraRespirationWater;
	[ReadOnly] public NativeArray<float> FloraTemperature;
	[ReadOnly] public NativeArray<float> PlanktonMassDelta;
	[ReadOnly] public NativeArray<float> PlanktonGlucoseDelta;
	[ReadOnly] public NativeArray<float> WaterCarbonDelta;
	public void Execute(int i)
	{
		float precipitationTemperature = PrecipitationTemperature[i];
		float rainMass = precipitationTemperature > WorldData.FreezingTemperature ? Precipitation[i] : 0;

		float waterMass = WaterMass[i];
		float newMass = waterMass + IceMelted[i] + rainMass - Evaporation[i] - WaterFrozen[i] + FloraRespirationWater[i];
		WaterMass[i] = newMass;
		SaltMass[i] -= SaltPlume[i];
		PlanktonMass[i] = LastPlanktonMass[i] + PlanktonMassDelta[i];
		PlanktonGlucose[i] = LastPlanktonGlucose[i] + PlanktonGlucoseDelta[i];
		CarbonMass[i] += WaterCarbonDelta[i];

		if (newMass <= 0)
		{
			WaterTemperature[i] = 0;
			PlanktonMass[i] = 0;
			PlanktonGlucose[i] = 0;
		}
		else
		{
			float remainingWater = waterMass - Evaporation[i] - WaterFrozen[i];
			float newTemperature = WorldData.SpecificHeatWater * (
				+ precipitationTemperature * rainMass
				+ WorldData.FreezingTemperature * IceMelted[i]
				+ FloraTemperature[i] * FloraRespirationWater[i]);
			float divisor = WorldData.SpecificHeatWater * (
				+ rainMass
				+ IceMelted[i]
				+ FloraRespirationWater[i]);
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
public struct UpdateMassAirSurfaceJob : IJobParallelFor {
	public NativeArray<float> AirTemperaturePotential;
	public NativeArray<float> VaporMass;
	public NativeArray<float> DustMass;
	public NativeArray<float> CarbonDioxide;
	[ReadOnly] public NativeArray<float> EvaporationWater;
	[ReadOnly] public NativeArray<float> EvaporationTemperatureWater;
	[ReadOnly] public NativeArray<float> EvaporationFlora;
	[ReadOnly] public NativeArray<float> EvaporationTemperatureFlora;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> DustEjected;
	[ReadOnly] public NativeArray<float> AirCarbonDelta;
	[ReadOnly] public NativeArray<float> SoilRespiration;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> Elevation;
	public void Execute(int i)
	{
		float vaporMass = VaporMass[i];
		float airMass = AirMass[i];
		float elevation = Elevation[i];

		AirTemperaturePotential[i] =
			((airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * AirTemperaturePotential[i] +
			EvaporationWater[i] * Atmosphere.GetPotentialTemperature(EvaporationTemperatureWater[i], elevation) * WorldData.SpecificHeatWaterVapor +
			EvaporationFlora[i] * Atmosphere.GetPotentialTemperature(EvaporationTemperatureFlora[i], elevation) * WorldData.SpecificHeatWaterVapor) /
			(airMass * WorldData.SpecificHeatAtmosphere + (vaporMass + EvaporationWater[i] + EvaporationFlora[i]) * WorldData.SpecificHeatWaterVapor);

		DustMass[i] = DustMass[i] + DustEjected[i];
		VaporMass[i] = vaporMass + EvaporationWater[i] + EvaporationFlora[i];
		CarbonDioxide[i] += AirCarbonDelta[i] + SoilRespiration[i] * WaterCoverage[i];
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
	public NativeArray<float> SoilCarbon;
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
	[ReadOnly] public NativeArray<float> SoilRespiration;
	[ReadOnly] public NativeArray<float> FloraDeath;
	[ReadOnly] public NativeArray<float> PlanktonDeath;
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
		SoilCarbon[i] = LastSoilFertility[i] - SoilRespiration[i] + FloraDeath[i] + PlanktonDeath[i];
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
	public NativeArray<float> FloraGlucose;

	[ReadOnly] public NativeArray<float> FloraGlucoseDelta;
	[ReadOnly] public NativeArray<float> FloraMassDelta;
	[ReadOnly] public NativeArray<float> FloraWaterDelta;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastGlucose;
	[ReadOnly] public NativeArray<float> LastWater;
	public void Execute(int i)
	{
		FloraMass[i] = math.max(0, LastMass[i] + FloraMassDelta[i]);
		FloraGlucose[i] = math.max(0, LastGlucose[i] + FloraGlucoseDelta[i]);
		FloraWater[i] = math.max(0, LastWater[i] + FloraWaterDelta[i]);
	}

}




[BurstCompile]
public struct UpdateWaterAirDiffusionJob : IJobParallelFor {
	public NativeArray<float> AirCarbon;
	public NativeArray<float> WaterCarbon;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public float WaterAirCarbonDiffusionCoefficient;
	public void Execute(int i)
	{
#if !DISABLE_WATER_AIR_CARBON_TRANSFER

		float airCarbon = AirCarbon[i];
		float waterCarbon = WaterCarbon[i];
		if ((waterCarbon > 0 || airCarbon > 0) && WaterMass[i] > 0)
		{
			float waterLayerMass = WaterMass[i] + SaltMass[i] + waterCarbon;
			float airLayerMass = AirMass[i] + airCarbon;
			float desiredCarbonDensity = (airCarbon + waterCarbon) / (waterLayerMass + airLayerMass);
			float diffusion = (desiredCarbonDensity - waterCarbon / waterLayerMass) * math.min(1, WaterAirCarbonDiffusionCoefficient / math.min(waterLayerMass, airLayerMass));

			AirCarbon[i] = airCarbon - diffusion;
			WaterCarbon[i] = waterCarbon + diffusion;
		}
#endif
	}

}



