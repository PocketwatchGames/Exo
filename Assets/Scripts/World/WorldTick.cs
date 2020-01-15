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
		public NativeArray<SimStateCell> Cells;

		public void Execute(int i)
		{
			var curCell = Cells[i];
			var cell = new SimStateCell();
			float depth = math.max(0, curCell.WaterElevation - curCell.Elevation);
			float evap = math.min(depth, 1f) * math.clamp(1.0f - curCell.RelativeHumidity/100,0,1);
			cell.WaterElevation = curCell.WaterElevation - evap;
			cell.RelativeHumidity = curCell.RelativeHumidity + evap;
			cell.CloudCoverage = curCell.CloudCoverage;
			cell.Elevation = curCell.Elevation;
			cell.Ice = curCell.Ice;
			cell.Vegetation = curCell.Vegetation;
			float hToC = curCell.RelativeHumidity / 100;
			cell.RelativeHumidity -= hToC;
			cell.CloudCoverage += hToC;
			float rain = math.clamp(cell.CloudElevation * cell.CloudCoverage / 1000, 0, cell.CloudCoverage);
			cell.CloudElevation = curCell.CloudElevation + 10;
			cell.WaterElevation += rain;
			cell.CloudCoverage -= rain;

			Cells[i] = cell;

		}
	}
}
