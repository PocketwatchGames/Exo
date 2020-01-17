using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using Unity.Burst;

public static class WorldTick {

	public static void Tick(ref SimState state, ref SimState nextState, int ticksToAdvance, ref StaticState staticState, ref WorldData worldData)
	{
		JobHandle lastJobHandle = default(JobHandle);

		var diffusion = new NativeArray<CellDiffusion>(state.Cells.Length * 6, Allocator.TempJob);
		var diffusionLimit = new NativeArray<CellDiffusion>(state.Cells.Length, Allocator.TempJob);
		var cells = new NativeArray<SimStateCell>[2];
		cells[0] = new NativeArray<SimStateCell>(state.Cells, Allocator.TempJob);
		cells[1] = new NativeArray<SimStateCell>(state.Cells.Length, Allocator.TempJob);
		for (int i = 0; i < ticksToAdvance; i++)
		{
			state.PlanetState.Ticks++;

			// TODO: update
			state.PlanetState.Position = state.PlanetState.Position;
			state.PlanetState.Rotation = state.PlanetState.Rotation;

			int lastStateIndex = i % 2;
			var diffusionJob = new DiffusionJob();
			diffusionJob.Last = cells[lastStateIndex];
			diffusionJob.CellDiffusion = diffusion;
			diffusionJob.Neighbors = staticState.Neighbors;
			diffusionJob.WaterDiffuseSpeed = worldData.WaterDiffuseSpeed;
			lastJobHandle = diffusionJob.Schedule(state.Cells.Length * 6, 100, lastJobHandle);

			var diffusionLimitJob = new DiffusionLimitJob();
			diffusionLimitJob.Last = cells[lastStateIndex];
			diffusionLimitJob.CellDiffusion = diffusion;
			diffusionLimitJob.DiffusionLimit = diffusionLimit;
			lastJobHandle = diffusionLimitJob.Schedule(state.Cells.Length, 100, lastJobHandle);


			var tickJob = new TickCellJob();
			tickJob.Cells = cells[(i + 1) % 2];
			tickJob.Last = cells[lastStateIndex];
			tickJob.worldData = worldData;
			tickJob.staticState = staticState;
			tickJob.Diffusion = diffusion;
			tickJob.DiffusionLimit = diffusionLimit;
			tickJob.PlanetState = state.PlanetState;
			lastJobHandle = tickJob.Schedule(state.Cells.Length, 100, lastJobHandle);
		}
		lastJobHandle.Complete();

		int outputBuffer = ticksToAdvance % 2;

		nextState.PlanetState = state.PlanetState;

		for (int i = 0; i < nextState.Cells.Length; i++)
		{
			nextState.Cells[i] = cells[outputBuffer][i];
		}


		for (int i = 0; i < 2; i++)
		{
			cells[i].Dispose();
		}
		diffusion.Dispose();
		diffusionLimit.Dispose();

	}
}
