using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

public class WorldSim {

	public static void Tick(ref SimState state, ref SimState nextState)
	{
		var tickJob = new TickJob();
		tickJob.Ticks = state.Ticks + 1;
		tickJob.Cells = new NativeArray<SimStateCell>(state.Cells, Allocator.TempJob);
		var jobHandle = tickJob.Schedule(state.Count, 1000);

		jobHandle.Complete();
		nextState.Ticks = tickJob.Ticks;
		for (int i = 0; i < state.Count; i++)
		{
			nextState.Cells[i] = tickJob.Cells[i];
		}

		tickJob.Cells.Dispose();

	}

	struct TickJob : IJobParallelFor {

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
			float evap = math.min(depth, 1f);
			cell.WaterElevation = curCell.WaterElevation - evap;
			cell.RelativeHumidity = curCell.RelativeHumidity + evap;
			cell.CloudCoverage = curCell.CloudCoverage;
			cell.CloudElevation = curCell.CloudElevation + curCell.CloudCoverage;
			cell.Elevation = curCell.Elevation;
			cell.Ice = curCell.Ice;
			cell.Vegetation = curCell.Vegetation;
			if (curCell.RelativeHumidity >= 100)
			{
				cell.RelativeHumidity -= 100.0f;
				cell.CloudCoverage += 100.0f;
			}
			if (curCell.CloudElevation >= 10000)
			{
				cell.CloudElevation -= 5000;
				cell.WaterElevation += curCell.CloudCoverage;
				cell.CloudCoverage -= curCell.CloudCoverage;
			}

			Cells[i] = cell;

		}
	}
}
