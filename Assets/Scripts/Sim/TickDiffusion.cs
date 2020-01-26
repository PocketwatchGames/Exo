using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

public struct CellDiffusion {
	public float Water;
}


[BurstCompile]
public struct DiffusionJob : IJobParallelFor {

	public NativeArray<CellDiffusion> CellDiffusion;
	[ReadOnly] public float WaterDiffuseSpeed;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<CellState> Last;
	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;
	[ReadOnly] public NativeArray<CellDependent> LastDependent;

	public void Execute(int i)
	{
		int n = Neighbors[i];
		if (n >= 0)
		{
			int index = i / 6;
			var curTerrain = LastTerrain[index];
			var curDependent = LastDependent[index];
			float water = 0;
			if (curDependent.WaterDepth > 0)
			{
				float waterElevation = curTerrain.Elevation + curDependent.WaterDepth;
				float nWaterElevation = LastTerrain[n].Elevation + LastDependent[n].WaterDepth;
				if (waterElevation > nWaterElevation)
				{
					water = (waterElevation - nWaterElevation) * WaterDiffuseSpeed;
				}
			}
		}
	}
}
[BurstCompile]
public struct DiffusionLimitJob : IJobParallelFor {

	public NativeArray<CellDiffusion> DiffusionLimit;
	[ReadOnly] public NativeArray<CellDiffusion> CellDiffusion;
	[ReadOnly] public NativeArray<CellState> Last;
	[ReadOnly] public NativeArray<CellDependent> LastDependent;

	public void Execute(int i)
	{
		float totalWater = 0;
		float waterLimit;
		for (int j = 0; j < 6; j++)
		{
			int index = i * 6 + j;
			totalWater += CellDiffusion[index].Water;
		}
		float waterDepth = LastDependent[i].WaterDepth;
		if (totalWater > waterDepth && waterDepth > 0)
		{
			waterLimit = waterDepth / totalWater;
		}
		else
		{
			waterLimit = 1;
		}
		DiffusionLimit[i] = new CellDiffusion { Water = waterLimit };
	}
}
