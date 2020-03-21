
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !UpdateTerrainJobDebug
[BurstCompile]
#endif
public struct UpdateTerrainJob : IJobParallelFor {
	public NativeArray<float> SoilFertility;
	public NativeArray<float> Roughness;
	public NativeArray<float> Elevation;
	public NativeArray<float> GroundWater;
	public NativeArray<float> LavaMass;
	public NativeArray<float> MagmaMass;
	public NativeArray<float> CrustDepth;

	[ReadOnly] public NativeArray<float> LastSoilFertility;
	[ReadOnly] public NativeArray<float> LastRoughness;
	[ReadOnly] public NativeArray<float> LastElevation;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> LastLavaMass;
	[ReadOnly] public NativeArray<float> LavaTemperature;
	[ReadOnly] public NativeArray<float> LastMagmaMass;
	[ReadOnly] public NativeArray<float> LastCrustDepth;
	[ReadOnly] public NativeArray<float> GroundWaterConsumed;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> DustSettled;
	[ReadOnly] public NativeArray<float> LavaCrystalized;
	public void Execute(int i)
	{
		Elevation[i] = LastElevation[i] + LavaCrystalized[i] / WorldData.MassLava;
		// TODO: improve soil fertility when dust settles
		SoilFertility[i] = LastSoilFertility[i];
		Roughness[i] = LastRoughness[i];
		GroundWater[i] = LastGroundWater[i] - GroundWaterConsumed[i];
		LavaMass[i] = LastLavaMass[i] - LavaCrystalized[i];
		MagmaMass[i] = LastMagmaMass[i];
		CrustDepth[i] = LastCrustDepth[i];
	}

}


#if !UpdateFloraJobDebug
[BurstCompile]
#endif
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



