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
	public struct TickCellJob : IJobParallelFor {

		public int Ticks;
		public float Gravity;
		public float SpinAngle;
		public float SpinSpeed;
		public float OrbitSpeed;
		public float TiltAngle;
		[ReadOnly] public NativeArray<int> Neighbors;
		[ReadOnly] public NativeArray<SimStateCell> Last;

		public NativeArray<SimStateCell> Cells;

		public void Execute(int i)
		{
			var curCell = Last[i];
			var cell = new SimStateCell();
			cell.CloudCoverage = curCell.CloudCoverage;
			cell.RelativeHumidity = curCell.RelativeHumidity;
			cell.WaterElevation = curCell.WaterElevation;
			cell.RelativeHumidity = curCell.RelativeHumidity;
			cell.CloudElevation = curCell.CloudElevation;
			cell.Elevation = curCell.Elevation;
			cell.Ice = curCell.Ice;
			cell.Vegetation = curCell.Vegetation;


			float depth = math.max(0, curCell.WaterElevation - curCell.Elevation);
			float evap = math.min(depth, 1f) * math.clamp(1.0f - curCell.RelativeHumidity/100,0,1);
			cell.WaterElevation = curCell.WaterElevation - evap;
			cell.RelativeHumidity = curCell.RelativeHumidity + evap;
			float hToC = curCell.RelativeHumidity / 100;
			cell.RelativeHumidity -= hToC;
			cell.CloudCoverage += hToC;
			if (cell.CloudCoverage > 1000)
			{
				cell.WaterElevation += curCell.CloudCoverage;
				cell.CloudCoverage = 0;
			}
			for (int j=0;j<6;j++)
			{
				int n = Neighbors[i * 6 + j];
				if (n >= 0)
				{
					cell.CloudCoverage += math.min(math.max(Last[n].CloudCoverage, curCell.CloudCoverage), Last[n].CloudCoverage - curCell.CloudCoverage) * 0.1f;
					cell.RelativeHumidity += math.min(math.max(Last[n].RelativeHumidity, curCell.RelativeHumidity), Last[n].RelativeHumidity - curCell.RelativeHumidity) * 0.1f;
				}
			}

			Cells[i] = cell;

		}
	}
}
