﻿using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


[BurstCompile]
public struct ApplyLatentHeatIceJob : IJobParallelFor {
	public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> LatentHeat;
	public void Execute(int i)
	{
		if (IceMass[i] > 0)
		{
			IceTemperature[i] = math.max(0, IceTemperature[i] + LatentHeat[i] / (WorldData.SpecificHeatIce * IceMass[i]));
		}
	}
}
[BurstCompile]
public struct ApplyLatentHeatTerrainJob : IJobParallelFor {
	public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeArray<float> LatentHeat;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		TerrainTemperature[i] = math.max(0, TerrainTemperature[i] + LatentHeat[i] / Atmosphere.GetSpecificHeatTerrain(HeatingDepth, SoilFertility[i]));
	}
}
[BurstCompile]
public struct ApplyLatentHeatFloraJob : IJobParallelFor {
	public NativeArray<float> FloraTemperature;
	[ReadOnly] public NativeArray<float> LatentHeat;
	[ReadOnly] public NativeArray<float> FloraWater;
	[ReadOnly] public NativeArray<float> FloraMass;
	public void Execute(int i)
	{
		float specificHeat = (FloraMass[i] * WorldData.SpecificHeatFlora + FloraWater[i] * WorldData.SpecificHeatWater);
		FloraTemperature[i] = math.max(0, FloraTemperature[i] + LatentHeat[i] / specificHeat);
	}
}

[BurstCompile]
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

[BurstCompile]
public struct ApplyLatentHeatWaterJob : IJobParallelFor {
	public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> LatentHeat;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			WaterTemperature[i] = math.max(0, WaterTemperature[i] + LatentHeat[i] / (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt));
		}
	}
}

