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
	public NativeArray<float> GroundWater;

	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	public void Execute(int i)
	{
		Terrain[i] = LastTerrain[i];
		GroundWater[i] = LastGroundWater[i];
	}

}



#if !TerrainGradientJobDebug
[BurstCompile]
#endif
public struct TerrainGradientJob : IJobParallelFor {
	public NativeArray<float> Gradient;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	public void Execute(int i)
	{
		float elevation = Terrain[i / 6].Elevation;
		int n = Neighbors[i];
		if (n == 0)
		{
			n = 0;
		}
		if (n >= 0)
		{
			float gradient = (elevation - Terrain[n].Elevation) * NeighborDistInverse[n];
			Gradient[i] = gradient;
		}
	}
}

