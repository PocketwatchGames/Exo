#define TerrainGradientJobDebug

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

	[ReadOnly] public NativeArray<float> LastSoilFertility;
	[ReadOnly] public NativeArray<float> LastRoughness;
	[ReadOnly] public NativeArray<float> LastElevation;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> GroundWaterConsumed;
	public void Execute(int i)
	{
		Elevation[i] = LastElevation[i];
		SoilFertility[i] = LastSoilFertility[i];
		Roughness[i] = LastRoughness[i];
		GroundWater[i] = LastGroundWater[i] - GroundWaterConsumed[i];
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
		FloraMass[i] = LastMass[i] + FloraMassDelta[i];
		FloraWater[i] = LastWater[i] - EvaporationMass[i] + GroundWaterConsumed[i];
	}

}



