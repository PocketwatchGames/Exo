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
		NativeArray<JobHandle> tickJobs = new NativeArray<JobHandle>(2, Allocator.TempJob);

		var diffusion = new NativeArray<CellDiffusion>(state.Cells.Length * 6, Allocator.TempJob);
		var diffusionLimit = new NativeArray<CellDiffusion>(state.Cells.Length, Allocator.TempJob);
		var cells = new NativeArray<SimCell>[2];
		cells[0] = new NativeArray<SimCell>(state.Cells, Allocator.TempJob);
		cells[1] = new NativeArray<SimCell>(state.Cells.Length, Allocator.TempJob);
		var wind = new NativeArray<SimWind>[2];
		wind[0] = new NativeArray<SimWind>(state.Wind, Allocator.TempJob);
		wind[1] = new NativeArray<SimWind>(state.Wind.Length, Allocator.TempJob);
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
			var diffusionJobHandle = diffusionJob.Schedule(state.Cells.Length * 6, 100, lastJobHandle);

			var diffusionLimitJob = new DiffusionLimitJob();
			diffusionLimitJob.Last = cells[lastStateIndex];
			diffusionLimitJob.CellDiffusion = diffusion;
			diffusionLimitJob.DiffusionLimit = diffusionLimit;
			var diffusionLimitJobHandle = diffusionLimitJob.Schedule(state.Cells.Length, 100, diffusionJobHandle);

			var cellJob = new TickCellJob();
			cellJob.Cells = cells[(i + 1) % 2];
			cellJob.DisplayCells = nextState.DisplayCells;
			cellJob.Last = cells[lastStateIndex];
			cellJob.LastWind = wind[lastStateIndex];
			cellJob.worldData = worldData;
			cellJob.staticState = staticState;
			cellJob.Diffusion = diffusion;
			cellJob.DiffusionLimit = diffusionLimit;
			cellJob.PlanetState = state.PlanetState;
			var tickJobHandle = cellJob.Schedule(state.Cells.Length, 100, diffusionLimitJobHandle);

			var windJob = new WindJob();
			windJob.Wind = wind[(i + 1) % 2];
			windJob.Last = cells[lastStateIndex];
			windJob.worldData = worldData;
			windJob.staticState = staticState;
			windJob.PlanetState = state.PlanetState;
			var windJobHandle = windJob.Schedule(state.Cells.Length, 100, lastJobHandle);

			tickJobs[0] = tickJobHandle;
			tickJobs[1] = windJobHandle;
			lastJobHandle = JobHandle.CombineDependencies(tickJobs);
		}
		lastJobHandle.Complete();

		int outputBuffer = ticksToAdvance % 2;

		nextState.PlanetState = state.PlanetState;

		for (int i = 0; i < nextState.Cells.Length; i++)
		{
			nextState.Cells[i] = cells[outputBuffer][i];
			nextState.DisplayPlanet.AtmosphericMass += nextState.Cells[i].AirMass;
			nextState.DisplayPlanet.CloudCoverage += nextState.Cells[i].CloudCoverage;
			nextState.DisplayPlanet.CloudMass += nextState.Cells[i].CloudMass;
			nextState.DisplayPlanet.EnergyDelta += nextState.DisplayCells[i].EnergyDelta;
			nextState.DisplayPlanet.EnergyEvapotranspiration += nextState.DisplayCells[i].EnergyEvapotranspiration;
			nextState.DisplayPlanet.EnergyIncoming += nextState.DisplayCells[i].EnergyIncoming;
			nextState.DisplayPlanet.EnergyLand += nextState.Cells[i].GroundEnergy;
			nextState.DisplayPlanet.EnergyLowerAir += nextState.Cells[i].AirEnergy;
			nextState.DisplayPlanet.EnergyOceanConduction += nextState.DisplayCells[i].EnergyOceanConduction;
			nextState.DisplayPlanet.EnergyShallowWater += nextState.Cells[i].WaterEnergy;
			nextState.DisplayPlanet.EnergySolarAbsorbedAtmosphere += nextState.DisplayCells[i].EnergySolarAbsorbedAtmosphere;
			nextState.DisplayPlanet.EnergySolarAbsorbedCloud += nextState.DisplayCells[i].EnergySolarAbsorbedCloud;
			nextState.DisplayPlanet.EnergySolarAbsorbedOcean += nextState.DisplayCells[i].EnergySolarAbsorbedOcean;
			nextState.DisplayPlanet.EnergySolarAbsorbedSurface += nextState.DisplayCells[i].EnergySolarAbsorbedSurface;
			nextState.DisplayPlanet.EnergySolarReflectedAtmosphere += nextState.DisplayCells[i].EnergySolarReflectedAtmosphere;
			nextState.DisplayPlanet.EnergySolarReflectedCloud += nextState.DisplayCells[i].EnergySolarReflectedCloud;
			nextState.DisplayPlanet.EnergySolarReflectedSurface += nextState.DisplayCells[i].EnergySolarReflectedSurface;
			nextState.DisplayPlanet.EnergySurfaceConduction += nextState.DisplayCells[i].EnergySurfaceConduction;
			nextState.DisplayPlanet.EnergyThermalAbsorbedAtmosphere += nextState.DisplayCells[i].EnergySolarAbsorbedAtmosphere;
			nextState.DisplayPlanet.EnergyThermalBackRadiation += nextState.DisplayCells[i].EnergyThermalBackRadiation;
			nextState.DisplayPlanet.EnergyThermalOceanRadiation += nextState.DisplayCells[i].EnergyThermalOceanRadiation;
			nextState.DisplayPlanet.EnergyThermalOutAtmosphere += nextState.DisplayCells[i].EnergyThermalOutAtmosphere;
			nextState.DisplayPlanet.EnergyThermalOutAtmosphericWindow += nextState.DisplayCells[i].EnergyThermalOutAtmosphericWindow;
			nextState.DisplayPlanet.EnergyThermalSurfaceRadiation += nextState.DisplayCells[i].EnergyThermalSurfaceRadiation;
			nextState.DisplayPlanet.Evaporation += nextState.DisplayCells[i].Evaporation;
			nextState.DisplayPlanet.IceMass += nextState.Cells[i].IceMass;
			nextState.DisplayPlanet.OceanCoverage += math.saturate(nextState.Cells[i].WaterDepth / nextState.Cells[i].Roughness);
			nextState.DisplayPlanet.OceanVolume += nextState.Cells[i].WaterAndIceDepth;
			nextState.DisplayPlanet.Rainfall += nextState.DisplayCells[i].Rainfall;
			nextState.DisplayPlanet.SeaLevel += nextState.Cells[i].Elevation + nextState.Cells[i].WaterAndIceDepth;
			nextState.DisplayPlanet.Temperature += nextState.Cells[i].AirTemperature;
			nextState.DisplayPlanet.WaterVapor += nextState.Cells[i].AirWaterMass;

			nextState.Wind[i] = wind[outputBuffer][i];
		}


		for (int i = 0; i < 2; i++)
		{
			cells[i].Dispose();
			wind[i].Dispose();
		}
		diffusion.Dispose();
		diffusionLimit.Dispose();
		tickJobs.Dispose();

	}
}
