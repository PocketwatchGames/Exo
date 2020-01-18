using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using Unity.Burst;
using UnityEngine;

public static class WorldTick {

	public static void Tick(ref SimState state, ref SimState nextState, int ticksToAdvance, ref StaticState staticState, ref WorldData worldData)
	{
		JobHandle lastJobHandle = default(JobHandle);

		var diffusion = new NativeArray<CellDiffusion>(state.Cells.Length * 6, Allocator.TempJob);
		var diffusionLimit = new NativeArray<CellDiffusion>(state.Cells.Length, Allocator.TempJob);
		var cells = new NativeArray<SimCell>[2];
		cells[0] = new NativeArray<SimCell>(state.Cells, Allocator.TempJob);
		cells[1] = new NativeArray<SimCell>(state.Cells.Length, Allocator.TempJob);
		for (int i = 0; i < ticksToAdvance; i++)
		{
			state.PlanetState.Ticks++;

			// TODO: update
			float distanceToSun = math.length(state.PlanetState.Position);
			float angleToSun = math.atan2(state.PlanetState.Position.z, state.PlanetState.Position.x);
			angleToSun += state.PlanetState.OrbitSpeed;
			state.PlanetState.Position = new float3(math.cos(angleToSun), 0, math.sin(angleToSun)) * distanceToSun;
			state.PlanetState.Rotation = new float3(state.PlanetState.Rotation.x, Mathf.Repeat(state.PlanetState.Rotation.y + state.PlanetState.SpinSpeed * worldData.SecondsPerTick, math.PI * 2), 0);

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
			tickJob.DisplayCells = state.DisplayCells;
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
