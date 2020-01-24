using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

[BurstCompile]
public struct TerrainJob : IJobParallelFor {

	public NativeArray<CellTerrain> Terrain;

	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;

	public void Execute(int i)
	{
		var lastTerrain = LastTerrain[i];
		var nextTerrain = new CellTerrain();

		nextTerrain.Elevation = lastTerrain.Elevation;
		nextTerrain.Roughness = lastTerrain.Roughness;
		nextTerrain.Vegetation = lastTerrain.Vegetation;
		nextTerrain.SoilFertility = lastTerrain.SoilFertility;

		Terrain[i] = nextTerrain;
	}
}
