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
	public void Execute(int i)
	{
		Elevation[i] = LastElevation[i];
		Terrain[i] = LastTerrain[i];
		GroundWater[i] = LastGroundWater[i];
	}

}


#if !UpdateFloraJobDebug
[BurstCompile]
#endif
public struct UpdateFloraJob : IJobParallelFor {
	public NativeArray<float> FloraMass;
	public NativeArray<float> FloraWater;

	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastWater;
	public void Execute(int i)
	{
		FloraMass[i] = LastMass[i];
		FloraWater[i] = LastWater[i];
	}

}



