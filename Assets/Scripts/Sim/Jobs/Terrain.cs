#define TerrainGradientJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !UpdateTerrainJobDebug
[BurstCompile]
#endif
public struct UpdateTerrainJob : IJobParallelFor {
	public NativeArray<CellTerrain> Terrain;
	public NativeArray<float> Elevation;
	public NativeArray<float> GroundWater;

	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;
	[ReadOnly] public NativeArray<float> LastElevation;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> GroundWaterConsumed;
	public void Execute(int i)
	{
		Elevation[i] = LastElevation[i];
		Terrain[i] = LastTerrain[i];
		GroundWater[i] = LastGroundWater[i] - GroundWaterConsumed[i];
	}

}


#if !UpdateFloraJobDebug
[BurstCompile]
#endif
public struct UpdateFloraJob : IJobParallelFor {
	public NativeArray<float> FloraMass;
	public NativeArray<float> FloraWater;

	[ReadOnly] public NativeArray<float> EvaporationMass;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastWater;
	[ReadOnly] public NativeArray<float> GroundWaterConsumed;
	public void Execute(int i)
	{
		FloraMass[i] = LastMass[i];
		FloraWater[i] = LastWater[i] - EvaporationMass[i] + GroundWaterConsumed[i];
	}

}



