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



