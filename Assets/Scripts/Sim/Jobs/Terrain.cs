using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !UpdateTerrainJobDebug
[BurstCompile]
#endif
public struct UpdateTerrainJob : IJobParallelFor {
	public NativeArray<CellTerrain> Terrain;

	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;
	public void Execute(int i)
	{
		Terrain[i] = LastTerrain[i];
	}

}
