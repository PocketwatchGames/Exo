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
	public float Cloud;
	public float Humidity;
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
			var curCell = Last[index];
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
			float cloud = math.max(0, (curDependent.CloudCoverage - LastDependent[n].CloudCoverage) * WaterDiffuseSpeed);
			float humidity = math.max(0, (curDependent.RelativeHumidity - LastDependent[n].RelativeHumidity) * WaterDiffuseSpeed);

			CellDiffusion[i] = new CellDiffusion { Water = water, Cloud = cloud, Humidity = humidity };
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
		float totalCloud = 0;
		float totalHumidity = 0;
		float waterLimit, cloudLimit, humidityLimit;
		for (int j = 0; j < 6; j++)
		{
			int index = i * 6 + j;
			totalWater += CellDiffusion[index].Water;
			totalCloud += CellDiffusion[index].Cloud;
			totalHumidity += CellDiffusion[index].Humidity;
		}
		float waterDepth = LastDependent[i].WaterDepth;
		float cloud = LastDependent[i].CloudCoverage;
		float humidity = LastDependent[i].RelativeHumidity;
		if (totalWater > waterDepth && waterDepth > 0)
		{
			waterLimit = waterDepth / totalWater;
		}
		else
		{
			waterLimit = 1;
		}
		if (totalCloud > cloud && cloud > 0)
		{
			cloudLimit = cloud / totalCloud;
		}
		else
		{
			cloudLimit = 1;
		}
		if (totalHumidity > humidity && humidity > 0)
		{
			humidityLimit = humidity / totalHumidity;
		}
		else
		{
			humidityLimit = 1;
		}
		DiffusionLimit[i] = new CellDiffusion { Water = waterLimit, Cloud = cloudLimit, Humidity = humidityLimit };
	}
}
