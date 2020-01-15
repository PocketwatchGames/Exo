using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using Unity.Entities;
using Unity.Burst;

public class WorldSimSystem : JobComponentSystem {


	[BurstCompile]
	struct WorldSimJob : IJobForEach<CellComponent> {
		public void Execute(ref CellComponent cell)
		{
			var curCell = cell;

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

		}
	}

	protected override JobHandle OnUpdate(JobHandle inputDeps)
	{
		var job = new WorldSimJob() { };
		return job.Schedule(this, inputDeps);
	}

}
