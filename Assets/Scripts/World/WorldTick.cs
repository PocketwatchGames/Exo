using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using Unity.Burst;

public class WorldTick {

	[BurstCompile]
	public struct TickFlowJob : IJobParallelFor {

		[ReadOnly] public NativeArray<int> Neighbors;
		[ReadOnly] public NativeArray<SimStateCell> Last;
		public NativeArray<float> WaterFlow;
		public float WaterDiffuseSpeed;

		public void Execute(int i)
		{
			int index = i / 6;
			var curCell = Last[index];
			float f = 0;
			if (curCell.WaterDepth > 0)
			{
				float waterElevation = curCell.Elevation + curCell.WaterDepth;
				int n = Neighbors[i];
				if (n >= 0)
				{
					float nWaterElevation = Last[n].Elevation + Last[n].WaterDepth;
					if (waterElevation > nWaterElevation)
					{
						f = (waterElevation - nWaterElevation) * WaterDiffuseSpeed;
					}
				}
			}
			WaterFlow[i] = f;
		}
	}
	[BurstCompile]
	public struct TickFlowLimitJob : IJobParallelFor {

		public NativeArray<float> WaterFlowLimit;
		[ReadOnly] public NativeArray<float> WaterFlow;
		[ReadOnly] public NativeArray<SimStateCell> Last;

		public void Execute(int i)
		{
			float total = 0;
			for (int j = 0; j < 6; j++)
			{
				total += WaterFlow[i * 6 + j];
			}
			float waterDepth = Last[i].WaterDepth;
			if (total > waterDepth && waterDepth > 0)
			{
				WaterFlowLimit[i] = waterDepth / total;
			}
			else
			{
				WaterFlowLimit[i] = 1;
			}
		}
	}
	[BurstCompile]
	public struct TickCellJob : IJobParallelFor {

		public int Ticks;
		public float Gravity;
		public float SpinAngle;
		public float SpinSpeed;
		public float OrbitSpeed;
		public float TiltAngle;
		[ReadOnly] public NativeArray<float> WaterFlow;
		[ReadOnly] public NativeArray<float> WaterFlowLimit;
		[ReadOnly] public NativeArray<int> Neighbors;
		[ReadOnly] public NativeArray<SimStateCell> Last;

		public NativeArray<SimStateCell> Cells;

		public void Execute(int i)
		{
			var curCell = Last[i];
			var cell = new SimStateCell();
			cell.CloudCoverage = curCell.CloudCoverage;
			cell.RelativeHumidity = curCell.RelativeHumidity;
			cell.WaterDepth = curCell.WaterDepth;
			cell.RelativeHumidity = curCell.RelativeHumidity;
			cell.CloudElevation = curCell.CloudElevation;
			cell.Elevation = curCell.Elevation;
			cell.Roughness = curCell.Roughness;
			cell.Ice = curCell.Ice;
			cell.Vegetation = curCell.Vegetation;
			cell.SoilFertility = curCell.SoilFertility;


			float evap = math.min(curCell.WaterDepth, 1f) * math.clamp(1.0f - curCell.RelativeHumidity / 100, 0, 1);
			cell.WaterDepth -= evap;
			cell.RelativeHumidity += evap;
			float hToC = curCell.RelativeHumidity / 100;
			cell.RelativeHumidity -= hToC;
			cell.CloudCoverage += hToC;
			if (cell.CloudCoverage > 1000)
			{
				cell.WaterDepth += curCell.CloudCoverage;
				cell.CloudCoverage = 0;
			}
			for (int j = 0; j < 6; j++)
			{
				int n = Neighbors[i * 6 + j];
				if (n >= 0)
				{
					cell.CloudCoverage -= curCell.CloudCoverage * 0.01f;
					cell.CloudCoverage += Last[n].CloudCoverage * 0.01f;
					cell.RelativeHumidity += math.min(math.max(Last[n].RelativeHumidity, curCell.RelativeHumidity), Last[n].RelativeHumidity - curCell.RelativeHumidity) * 0.1f;

					cell.WaterDepth -= WaterFlow[i * 6 + j] * WaterFlowLimit[i];
					for (int k=0;k<6;k++)
					{
						int nToMeIndex = n * 6 + k;
						if (Neighbors[nToMeIndex] == i)
						{
							cell.WaterDepth += WaterFlow[nToMeIndex] * WaterFlowLimit[n];
							break;
						}
					}
				}
			}

			Cells[i] = cell;

		}
	}
}
