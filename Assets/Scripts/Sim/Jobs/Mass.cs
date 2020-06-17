using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;


[BurstCompile]
public struct UpdateMassWaterJob : IJobParallelFor {
	public NativeSlice<float> WaterMass;
	public NativeSlice<float> SaltMass;
	public NativeSlice<float> CarbonMass;
	public NativeSlice<float> OxygenMass;
	public NativeSlice<float> NitrogenMass;
	public NativeSlice<float> GlucoseMass;
	public NativeSlice<float> MineralsMass;
	public NativeSlice<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> SaltPlume;
	[ReadOnly] public NativeArray<float> SaltPlumeTemperature;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float> LastWaterMass;
	[ReadOnly] public NativeArray<float> LastCarbonMass;
	[ReadOnly] public NativeArray<float> LastOxygenMass;
	[ReadOnly] public NativeArray<float> LastNitrogenMass;
	[ReadOnly] public NativeArray<float> LastGlucoseMass;
	[ReadOnly] public NativeArray<float> LastMineralsMass;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		int columnIndex = i % Count;
		int index = i + Count;
		int downIndex = i;
		float waterMass = LastWaterMass[index];
		float saltMass = LastSaltMass[index];
		WaterMass[i] = waterMass;
		SaltMass[i] = saltMass;
		CarbonMass[i] = LastCarbonMass[index];
		OxygenMass[i] = LastOxygenMass[index];
		NitrogenMass[i] = LastNitrogenMass[index];
		GlucoseMass[i] = LastGlucoseMass[index];
		MineralsMass[i] = LastMineralsMass[index];
		if (LastWaterMass[downIndex] == 0 && waterMass > 0)
		{
			float saltPlume = SaltPlume[columnIndex];
			float newTemp = 
				(WaterTemperature[i] * (waterMass * WorldData.SpecificHeatWater + saltMass * WorldData.SpecificHeatSalt) + 
				SaltPlumeTemperature[columnIndex] * saltPlume * WorldData.SpecificHeatSalt)
				/ (waterMass * WorldData.SpecificHeatWater + (saltMass + saltPlume) * WorldData.SpecificHeatSalt);
			WaterTemperature[i] = newTemp;
			SaltMass[i] += saltPlume;
		}
	}
}
[BurstCompile]
public struct UpdateMassWaterSurfaceJob : IJobParallelFor {
	public NativeSlice<float> WaterTemperature;
	public NativeSlice<float> WaterMass;
	public NativeSlice<float> SaltMass;
	public NativeSlice<float> CarbonMass;
	public NativeSlice<float> OxygenMass;
	public NativeSlice<float> NitrogenMass;
	public NativeSlice<float> GlucoseMass;
	public NativeSlice<float> MineralsMass;
	[ReadOnly] public NativeArray<float> Evaporation;
	[ReadOnly] public NativeArray<float> IceMelted;
	[ReadOnly] public NativeArray<float> Precipitation;
	[ReadOnly] public NativeArray<float> PrecipitationTemperature;
	[ReadOnly] public NativeArray<float> WaterFrozen;
	[ReadOnly] public NativeArray<float> SaltPlume;
	[ReadOnly] public NativeArray<float> FloraRespirationWater;
	[ReadOnly] public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeArray<float> WaterCarbonDelta;
	[ReadOnly] public NativeArray<float> WaterOxygenDelta;
	[ReadOnly] public NativeArray<float> WaterNitrogenDelta;
	[ReadOnly] public NativeArray<float> WaterGlucoseDelta;
	[ReadOnly] public NativeArray<float> WaterMineralsDelta;
	public void Execute(int i)
	{
		float precipitationTemperature = PrecipitationTemperature[i];
		float rainMass = precipitationTemperature > WorldData.FreezingTemperature ? Precipitation[i] : 0;

		float waterMass = WaterMass[i];
		float newMass = waterMass + IceMelted[i] + rainMass - Evaporation[i] - WaterFrozen[i] + FloraRespirationWater[i];
		WaterMass[i] = newMass;
		SaltMass[i] -= SaltPlume[i];
		CarbonMass[i] += WaterCarbonDelta[i];
		OxygenMass[i] += WaterOxygenDelta[i];
		NitrogenMass[i] += WaterNitrogenDelta[i];
		GlucoseMass[i] += WaterGlucoseDelta[i];
		MineralsMass[i] += WaterMineralsDelta[i];

		if (newMass <= 0)
		{
			WaterTemperature[i] = 0;
		}
		else
		{
			float remainingWater = waterMass - Evaporation[i] - WaterFrozen[i];
			float newTemperature = WorldData.SpecificHeatWater * (
				+ precipitationTemperature * rainMass
				+ WorldData.FreezingTemperature * IceMelted[i]
				+ TerrainTemperature[i] * FloraRespirationWater[i]);
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
	public NativeSlice<float> SurfaceWaterMass;
	public NativeSlice<float> SurfaceWaterTemperature;
	[ReadOnly] public NativeSlice<float> SurfaceSaltMass;
	[ReadOnly] public NativeSlice<float> AirTemperaturePotential;
	[ReadOnly] public NativeSlice<float> GroundCondensation;
	[ReadOnly] public NativeSlice<float> LayerMiddle;
	[ReadOnly] public int LayerCount;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		float waterMass = SurfaceWaterMass[i];
		float newWaterMass = waterMass;
		float newWaterTemperature = SurfaceWaterTemperature[i];
		for (int j = 0; j < LayerCount; j++)
		{
			int layerIndex = j * Count + i;
			float condensation = GroundCondensation[layerIndex];
			if (condensation > 0)
			{
				newWaterMass += condensation;
				float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[layerIndex], LayerMiddle[layerIndex]);
				float specificHeatSalt = SurfaceSaltMass[i] * WorldData.SpecificHeatSalt;
				newWaterTemperature =
					((waterMass * WorldData.SpecificHeatWater + specificHeatSalt) * newWaterTemperature
					+ condensation * WorldData.SpecificHeatWater * airTemperatureAbsolute)
					/ (newWaterMass * WorldData.SpecificHeatWater + specificHeatSalt);
				waterMass = newWaterMass;
			}
		}
		SurfaceWaterTemperature[i] = newWaterTemperature;
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
	[ReadOnly] public NativeArray<float> DropletDelta;
	[ReadOnly] public NativeSlice<float> CloudCondensation;
	[ReadOnly] public int LayerCount;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		float cloudMass = LastCloudMass[i] - CloudEvaporation[i] - PrecipitationMass[i];
		float dropletMass = LastDropletMass[i] + DropletDelta[i];

		for (int j = 0; j < LayerCount; j++)
		{
			int layerIndex = j * Count + i;
			float condensation = CloudCondensation[layerIndex];
			if (cloudMass > 0)
			{
				dropletMass = dropletMass * cloudMass / (cloudMass + condensation);
			}
			cloudMass += condensation;
		}

		CloudMass[i] = cloudMass;
		CloudDropletMass[i] = dropletMass;

	}
}


[BurstCompile]
public struct UpdateMassAirJob : IJobParallelFor {
	public NativeSlice<float> VaporMass;
	public NativeSlice<float> CarbonDioxideMass;
	public NativeSlice<float> DustMass;
	public NativeSlice<float> TemperaturePotential;
	[ReadOnly] public NativeSlice<float> LastVaporMass;
	[ReadOnly] public NativeSlice<float> LastCarbonDioxideMass;
	[ReadOnly] public NativeSlice<float> LastDustMass;
	[ReadOnly] public NativeSlice<float> CloudCondensation;
	[ReadOnly] public NativeSlice<float> GroundCondensation;
	[ReadOnly] public NativeSlice<float> LayerElevation;
	[ReadOnly] public NativeSlice<float> LayerHeight;
	[ReadOnly] public NativeSlice<float> CloudElevation;
	[ReadOnly] public NativeSlice<float> CloudMass;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeArray<float> DustUp;
	[ReadOnly] public NativeArray<float> DustDown;
	[ReadOnly] public int LayerCount;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		int layerIndex = i / Count;
		int columnIndex = i - layerIndex * Count;
		int dustLayerIndex = i + Count;
		int dustLayerIndexDown = i;
		int dustLayerIndexUp = dustLayerIndex + Count;
		float newVaporMass = LastVaporMass[i] - CloudCondensation[i] - GroundCondensation[i];
		float newDustMass = LastDustMass[i] + DustDown[dustLayerIndexUp] + DustUp[dustLayerIndexDown] - DustDown[dustLayerIndex];
		float newCarbonDioxide = LastCarbonDioxideMass[i];
		bool isTop = layerIndex == LayerCount - 1;
		if (!isTop)
		{
			newDustMass -= DustUp[dustLayerIndex];
		}

		float newTemperature = TemperaturePotential[i];
		float cloudCondensationInLayer = 0;
		float cloudMassInLayer = 0;
		if (CloudElevation[columnIndex] >= LayerElevation[i] && CloudElevation[columnIndex] < LayerElevation[i] + LayerHeight[i])
		{
			cloudCondensationInLayer = CloudCondensation[i];
			cloudMassInLayer = CloudMass[columnIndex];
		}


		float specificHeatAir = Atmosphere.GetSpecificHeatAir(AirMass[i], VaporMass[i], cloudMassInLayer);
		float specificHeatAirNew = Atmosphere.GetSpecificHeatAir(AirMass[i], VaporMass[i], cloudMassInLayer + cloudCondensationInLayer);
		newTemperature = (newTemperature * specificHeatAir + TemperaturePotential[i] * cloudCondensationInLayer * WorldData.SpecificHeatWater) / specificHeatAirNew;


		TemperaturePotential[i] = newTemperature;
		VaporMass[i] = newVaporMass;
		DustMass[i] = newDustMass;
		CarbonDioxideMass[i] = newCarbonDioxide;
	}
}


[BurstCompile]
public struct UpdateMassAirSurfaceJob : IJobParallelFor {
	public NativeSlice<float> AirTemperaturePotential;
	public NativeSlice<float> VaporMass;
	public NativeSlice<float> DustMass;
	public NativeSlice<float> CarbonDioxide;
	public NativeSlice<float> Oxygen;
	public NativeSlice<float> Nitrogen;
	public NativeSlice<float> Methane;
	public NativeSlice<float> Minerals;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeArray<float> EvaporationWater;
	[ReadOnly] public NativeSlice<float> EvaporationTemperatureWater;
	[ReadOnly] public NativeArray<float> EvaporationFlora;
	[ReadOnly] public NativeArray<float> EvaporationTemperatureFlora;
	[ReadOnly] public NativeArray<float> DustEjected;
	[ReadOnly] public NativeArray<float> CarbonDioxideDelta;
	[ReadOnly] public NativeArray<float> OxygenDelta;
	[ReadOnly] public NativeArray<float> NitrogenDelta;
	[ReadOnly] public NativeArray<float> MethaneDelta;
	[ReadOnly] public NativeArray<float> MineralsDelta;
	[ReadOnly] public NativeSlice<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> Elevation;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeSlice<float> LayerElevation;
	[ReadOnly] public NativeSlice<float> LayerHeight;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		float vaporMass = VaporMass[i];
		float airMass = AirMass[i];
		float elevation = Elevation[i];
		int columnIndex = i % Count;
		float cloudMassInLayer = Atmosphere.GetCloudMassInLayer(CloudMass[columnIndex], CloudElevation[columnIndex], LayerElevation[i], LayerHeight[i]);
		float specificHeatAir = Atmosphere.GetSpecificHeatAir(AirMass[i], VaporMass[i], cloudMassInLayer);
		float specificHeatAirNew = Atmosphere.GetSpecificHeatAir(AirMass[i], VaporMass[i] + EvaporationWater[i] + EvaporationFlora[i], cloudMassInLayer);
		AirTemperaturePotential[i] =
			(specificHeatAir * AirTemperaturePotential[i] +
			EvaporationWater[i] * Atmosphere.GetPotentialTemperature(EvaporationTemperatureWater[i], elevation) * WorldData.SpecificHeatWaterVapor +
			EvaporationFlora[i] * Atmosphere.GetPotentialTemperature(EvaporationTemperatureFlora[i], elevation) * WorldData.SpecificHeatWaterVapor) /
			specificHeatAirNew;

		DustMass[i] = DustMass[i] + DustEjected[i];
		VaporMass[i] = vaporMass + EvaporationWater[i] + EvaporationFlora[i];
		CarbonDioxide[i] += CarbonDioxideDelta[i];
		Oxygen[i] += OxygenDelta[i];
		Nitrogen[i] += NitrogenDelta[i];
		Methane[i] += MethaneDelta[i];
		Minerals[i] += MineralsDelta[i];
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
	public NativeArray<float> Roughness;
	public NativeArray<float> Elevation;
	public NativeArray<float> GroundWater;
	public NativeArray<float> LavaMass;
	public NativeArray<float> LavaTemperature;
	public NativeArray<float> MagmaMass;
	public NativeArray<float> CrustDepth;

	[ReadOnly] public NativeArray<float> LastRoughness;
	[ReadOnly] public NativeArray<float> LastElevation;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> LastLavaMass;
	[ReadOnly] public NativeArray<float> LastMagmaMass;
	[ReadOnly] public NativeArray<float> LastCrustDepth;
	[ReadOnly] public NativeArray<float> GroundWaterConsumed;
	[ReadOnly] public NativeSlice<float> WaterCoverage;
	[ReadOnly] public NativeSlice<float> DustSettled;
	[ReadOnly] public NativeArray<float> LavaCrystalized;
	[ReadOnly] public NativeArray<float> LavaEjected;
	[ReadOnly] public NativeArray<float> CrustDelta;
	[ReadOnly] public float MagmaTemperature;
	[ReadOnly] public float LavaToRockMassAdjustment;
	public void Execute(int i)
	{
		float lavaCrystalized = LavaCrystalized[i];
		float lavaEjected = LavaEjected[i];
		float lavaMass = LastLavaMass[i];
		float newLavaMass = lavaMass - lavaCrystalized + lavaEjected;

		float lavaCrystalizedDepth = lavaCrystalized * LavaToRockMassAdjustment / WorldData.MassLava;
		Elevation[i] = LastElevation[i] + lavaCrystalizedDepth;
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
	[ReadOnly] public NativeArray<float> FloraWaterConsumed;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastGlucose;
	[ReadOnly] public NativeArray<float> LastWater;
	public void Execute(int i)
	{
		FloraMass[i] = math.max(0, LastMass[i] + FloraMassDelta[i]);
		FloraGlucose[i] = math.max(0, LastGlucose[i] + FloraGlucoseDelta[i]);
		FloraWater[i] = math.max(0, LastWater[i] + FloraWaterDelta[i] + FloraWaterConsumed[i]);
	}

}




[BurstCompile]
public struct UpdateWaterAirDiffusionJob : IJobParallelFor {
	public NativeSlice<float> AirCarbon;
	public NativeSlice<float> AirOxygen;
	public NativeSlice<float> WaterCarbon;
	public NativeSlice<float> WaterOxygen;
	[ReadOnly] public NativeSlice<float> WaterMass;
	[ReadOnly] public NativeSlice<float> SaltMass;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeSlice<float> WaterDepth;
	[ReadOnly] public float WaterAirDiffusionCoefficient;
	[ReadOnly] public float WaterAirCarbonDepth;
	public void Execute(int i)
	{

		if (WaterMass[i] > 0)
		{
			float airCarbon = AirCarbon[i];
			float airOxygen = AirOxygen[i];
			float waterCarbon = WaterCarbon[i];
			float waterOxygen = WaterOxygen[i];
			float waterLayerMass = WaterMass[i] + SaltMass[i] + waterCarbon + waterOxygen;
			if (waterCarbon > 0 || airCarbon > 0)
			{
				float desiredCarbonDensity = (airCarbon + waterCarbon) / (waterLayerMass + AirMass[i]);
				float diffusionDepth = math.min(1, WaterAirCarbonDepth / WaterDepth[i]);
				float diffusionCarbon = (desiredCarbonDensity * waterLayerMass - waterCarbon) * WaterAirDiffusionCoefficient;
				diffusionCarbon *= diffusionDepth;

				AirCarbon[i] = airCarbon - diffusionCarbon;
				WaterCarbon[i] = waterCarbon + diffusionCarbon;
			}
			if (waterOxygen > 0 || airOxygen > 0)
			{
				float desiredOxygenDensity = (airOxygen + waterOxygen) / (waterLayerMass + AirMass[i]);
				float diffusionDepth = math.min(1, WaterAirCarbonDepth / WaterDepth[i]);
				float diffusionOxygen = (desiredOxygenDensity * waterLayerMass - waterOxygen) * WaterAirDiffusionCoefficient;
				diffusionOxygen *= diffusionDepth;

				AirOxygen[i] = airOxygen - diffusionOxygen;
				WaterOxygen[i] = waterOxygen + diffusionOxygen;
			}

		}
	}

}



