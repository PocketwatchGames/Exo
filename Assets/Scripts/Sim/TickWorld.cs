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

		var diffusion = new NativeArray<CellDiffusion>(state.CellDependents.Length * 6, Allocator.TempJob);
		var diffusionLimit = new NativeArray<CellDiffusion>(state.CellDependents.Length, Allocator.TempJob);
		var cells = new NativeArray<CellState>[2];
		cells[0] = new NativeArray<CellState>(state.CellStates, Allocator.TempJob);
		cells[1] = new NativeArray<CellState>(state.CellStates.Length, Allocator.TempJob);
		var dependents = new NativeArray<CellDependent>[2];
		dependents[0] = new NativeArray<CellDependent>(state.CellDependents, Allocator.TempJob);
		dependents[1] = new NativeArray<CellDependent>(state.CellDependents.Length, Allocator.TempJob);
		var terrain = new NativeArray<CellTerrain>[2];
		terrain[0] = new NativeArray<CellTerrain>(state.CellTerrains, Allocator.TempJob);
		terrain[1] = new NativeArray<CellTerrain>(state.CellTerrains.Length, Allocator.TempJob);
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
			diffusionJob.LastDependent = dependents[lastStateIndex];
			diffusionJob.LastTerrain = terrain[lastStateIndex];
			diffusionJob.LastTerrain = terrain[lastStateIndex];
			diffusionJob.CellDiffusion = diffusion;
			diffusionJob.Neighbors = staticState.Neighbors;
			diffusionJob.WaterDiffuseSpeed = worldData.WaterDiffuseSpeed;
			var diffusionJobHandle = diffusionJob.Schedule(state.CellDependents.Length * 6, 100, lastJobHandle);

			var diffusionLimitJob = new DiffusionLimitJob();
			diffusionLimitJob.Last = cells[lastStateIndex];
			diffusionLimitJob.LastDependent = dependents[lastStateIndex];
			diffusionLimitJob.CellDiffusion = diffusion;
			diffusionLimitJob.DiffusionLimit = diffusionLimit;
			var diffusionLimitJobHandle = diffusionLimitJob.Schedule(state.CellDependents.Length, 100, diffusionJobHandle);

			var cellJob = new TickCellJob();
			cellJob.Cells = cells[(i + 1) % 2];
			cellJob.DisplayCells = nextState.CellDisplays;
			cellJob.Last = cells[lastStateIndex];
			cellJob.LastDependent = dependents[lastStateIndex];
			cellJob.LastTerrain = terrain[lastStateIndex];
			cellJob.worldData = worldData;
			cellJob.staticState = staticState;
			cellJob.Diffusion = diffusion;
			cellJob.DiffusionLimit = diffusionLimit;
			cellJob.PlanetState = state.PlanetState;
			var tickJobHandle = cellJob.Schedule(state.CellDependents.Length, 100, diffusionLimitJobHandle);

			var dependentJob = new DependentJob();
			dependentJob.NextDependent = dependents[(i + 1) % 2];
			dependentJob.LastCell = cells[(i + 1) % 2];
			dependentJob.LastDependent = dependents[lastStateIndex];
			dependentJob.LastTerrain = terrain[lastStateIndex];
			dependentJob.worldData = worldData;
			dependentJob.staticState = staticState;
			dependentJob.PlanetState = state.PlanetState;
			var dependentJobHandle = dependentJob.Schedule(state.CellDependents.Length, 100, tickJobHandle);

			var terrainJob = new TerrainJob();
			terrainJob.LastTerrain = terrain[lastStateIndex];
			terrainJob.Terrain = terrain[(i + 1) % 2];
			var terrainJobHandle = terrainJob.Schedule(state.CellTerrains.Length, 100, lastJobHandle);

			tickJobs[0] = terrainJobHandle;
			tickJobs[1] = dependentJobHandle;
			lastJobHandle = JobHandle.CombineDependencies(tickJobs);
		}
		lastJobHandle.Complete();

		int outputBuffer = ticksToAdvance % 2;

		nextState.PlanetState = state.PlanetState;

		for (int i = 0; i < nextState.CellDependents.Length; i++)
		{
			nextState.DisplayPlanet.CloudCoverage += nextState.CellDependents[i].CloudCoverage;
			nextState.DisplayPlanet.CloudMass += nextState.CellStates[i].CloudMass;
			nextState.DisplayPlanet.EnergyDelta += nextState.CellDisplays[i].EnergyDelta;
			nextState.DisplayPlanet.EnergyEvapotranspiration += nextState.CellDisplays[i].EnergyEvapotranspiration;
			nextState.DisplayPlanet.EnergyIncoming += nextState.CellDisplays[i].EnergyIncoming;
			nextState.DisplayPlanet.EnergyLand += nextState.CellStates[i].GroundEnergy;
			nextState.DisplayPlanet.EnergyLowerAir += nextState.CellStates[i].AirEnergy;
			nextState.DisplayPlanet.EnergyOceanConduction += nextState.CellDisplays[i].EnergyOceanConduction;
			nextState.DisplayPlanet.EnergyShallowWater += nextState.CellStates[i].WaterEnergy;
			nextState.DisplayPlanet.EnergySolarAbsorbedAtmosphere += nextState.CellDisplays[i].EnergySolarAbsorbedAtmosphere;
			nextState.DisplayPlanet.EnergySolarAbsorbedCloud += nextState.CellDisplays[i].EnergySolarAbsorbedCloud;
			nextState.DisplayPlanet.EnergySolarAbsorbedOcean += nextState.CellDisplays[i].EnergySolarAbsorbedOcean;
			nextState.DisplayPlanet.EnergySolarAbsorbedSurface += nextState.CellDisplays[i].EnergySolarAbsorbedSurface;
			nextState.DisplayPlanet.EnergySolarReflectedAtmosphere += nextState.CellDisplays[i].EnergySolarReflectedAtmosphere;
			nextState.DisplayPlanet.EnergySolarReflectedCloud += nextState.CellDisplays[i].EnergySolarReflectedCloud;
			nextState.DisplayPlanet.EnergySolarReflectedSurface += nextState.CellDisplays[i].EnergySolarReflectedSurface;
			nextState.DisplayPlanet.EnergySurfaceConduction += nextState.CellDisplays[i].EnergySurfaceConduction;
			nextState.DisplayPlanet.EnergyThermalAbsorbedAtmosphere += nextState.CellDisplays[i].EnergySolarAbsorbedAtmosphere;
			nextState.DisplayPlanet.EnergyThermalBackRadiation += nextState.CellDisplays[i].EnergyThermalBackRadiation;
			nextState.DisplayPlanet.EnergyThermalOceanRadiation += nextState.CellDisplays[i].EnergyThermalOceanRadiation;
			nextState.DisplayPlanet.EnergyThermalOutAtmosphere += nextState.CellDisplays[i].EnergyThermalOutAtmosphere;
			nextState.DisplayPlanet.EnergyThermalSurfaceOutAtmosphericWindow += nextState.CellDisplays[i].EnergyThermalSurfaceOutAtmosphericWindow;
			nextState.DisplayPlanet.EnergyThermalSurfaceRadiation += nextState.CellDisplays[i].EnergyThermalSurfaceRadiation;
			nextState.DisplayPlanet.Evaporation += nextState.CellDisplays[i].Evaporation;
			nextState.DisplayPlanet.IceMass += nextState.CellStates[i].IceMass;
			nextState.DisplayPlanet.OceanCoverage += math.saturate(nextState.CellDependents[i].WaterDepth / nextState.CellTerrains[i].Roughness);
			nextState.DisplayPlanet.OceanVolume += nextState.CellDependents[i].WaterAndIceDepth;
			nextState.DisplayPlanet.Rainfall += nextState.CellDisplays[i].Rainfall;
			nextState.DisplayPlanet.SeaLevel += nextState.CellTerrains[i].Elevation + nextState.CellDependents[i].WaterAndIceDepth;
			nextState.DisplayPlanet.Temperature += nextState.CellDependents[i].AirTemperature;
			nextState.DisplayPlanet.WaterVapor += nextState.CellStates[i].AirWaterMass;

			nextState.CellStates[i] = cells[outputBuffer][i];
			nextState.CellDependents[i] = dependents[outputBuffer][i];
			nextState.CellTerrains[i] = terrain[outputBuffer][i];
		}


		for (int i = 0; i < 2; i++)
		{
			cells[i].Dispose();
			dependents[i].Dispose();
			terrain[i].Dispose();
		}
		diffusion.Dispose();
		diffusionLimit.Dispose();
		tickJobs.Dispose();

	}
}
