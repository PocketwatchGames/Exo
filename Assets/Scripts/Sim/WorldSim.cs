﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using Unity.Burst;
using UnityEngine;
using Unity.Profiling;

[Serializable]
public class WorldSim {

	public JobHelper SimJob;
	public JobHelper NeighborJob;


	private int _cellCount;

	private NativeArray<JobHandle> applyAdvectionAirJobHandles;
	private NativeArray<JobHandle> applyAdvectionWaterJobHandles;
	private NativeArray<JobHandle> energyJobHandles;
	private List<NativeList<JobHandle>> jobHandleDependencies = new List<NativeList<JobHandle>>();
	private List<NativeArray<float>> tempArrays = new List<NativeArray<float>>();

	public WorldSim(int cellCount, ref WorldData worldData)
	{
		_cellCount = cellCount;

		SimJob = new JobHelper(_cellCount);
		NeighborJob = new JobHelper(_cellCount * 6);

		applyAdvectionAirJobHandles = new NativeArray<JobHandle>(worldData.AirLayers, Allocator.Persistent);
		applyAdvectionWaterJobHandles = new NativeArray<JobHandle>(worldData.WaterLayers, Allocator.Persistent);
		energyJobHandles = new NativeArray<JobHandle>(worldData.LayerCount, Allocator.Persistent);

	}

	public void Dispose(ref WorldData worldData)
	{
		energyJobHandles.Dispose();
		applyAdvectionAirJobHandles.Dispose();
		applyAdvectionWaterJobHandles.Dispose();

		foreach (var d in jobHandleDependencies)
		{
			d.Dispose();
		}
		foreach (var a in tempArrays)
		{
			a.Dispose();
		}
	}


	public JobHandle Tick(
		int ticksToAdvance,
		SimState[] states,
		TempState[] tempStates,
		ref StaticState staticState,
		ref WorldData worldData,
		ref SimSettings settings,
		ref int curStateIndex,
		ref int curTempStateIndex)
	{
		var tickJobHandle = default(JobHandle);
		foreach (var d in jobHandleDependencies)
		{
			d.Dispose();
		}
		foreach (var a in tempArrays)
		{
			a.Dispose();
		}
		jobHandleDependencies.Clear();
		tempArrays.Clear();
		int firstStateIndex = curStateIndex;

		#region Init Time step

		int lastStateIndex = curStateIndex;
		ref var lastState = ref states[lastStateIndex];
		while (curStateIndex == lastStateIndex || curStateIndex == firstStateIndex)
		{
			curStateIndex = (curStateIndex + 1) % states.Length;
		}
		ref var nextState = ref states[curStateIndex];

		curTempStateIndex = (curTempStateIndex + 1) % tempStates.Length;
		ref var tempState = ref tempStates[curTempStateIndex];

		tickJobHandle = tempState.Clear(staticState.Count, ref worldData, tickJobHandle);
		tickJobHandle = TempState.Update(SimJob, ref lastState, ref tempState, ref worldData, ref staticState, tickJobHandle);

		#endregion


		#region Update Planetary Globals

		nextState.PlanetState = lastState.PlanetState;
		nextState.PlanetState.Ticks++;

		// TODO: update
		float distanceToSun = math.length(lastState.PlanetState.Position);
		float angleToSun = math.atan2(lastState.PlanetState.Position.z, lastState.PlanetState.Position.x);
		angleToSun += lastState.PlanetState.OrbitSpeed;
		nextState.PlanetState.Position = new float3(math.cos(angleToSun), 0, math.sin(angleToSun)) * distanceToSun;
		nextState.PlanetState.Rotation = new float3(lastState.PlanetState.Rotation.x, Mathf.Repeat(lastState.PlanetState.Rotation.y + lastState.PlanetState.SpinSpeed * worldData.SecondsPerTick, math.PI * 2), 0);

		float coriolisTerm = 2 * lastState.PlanetState.SpinSpeed;
		#endregion

		DoEnergyCycle(energyJobHandles, tickJobHandle, ref nextState, ref lastState, ref tempState, ref staticState, ref worldData, ref settings);

		DoStateChange(energyJobHandles, tickJobHandle, ref nextState, ref lastState, ref tempState, ref staticState, ref worldData, ref settings);

		var groundWaterJob = UpdateGroundWater(energyJobHandles, tickJobHandle, ref nextState, ref lastState, ref tempState, ref staticState, ref worldData, ref settings);

		#region Air Advection and Diffusion
		// Buoyancy, Updrafts, and mixing occur across air layers and water layers
		// TODO: add an empty air layer on top and bottom so we can calculate up/down diffusion in a single step 
		// Temperature and trace elements diffuse into neighboring horizontal cells based on a diffusion constant
		// Air, Water, Cloud
		#region Update Velocity


		var airTerrainFrictionJobHandle = SimJob.Schedule(new AirTerrainFrictionJob()
		{
			Force = tempState.WindFriction,
			IceCoverage = tempState.IceCoverage,
			WaterCoverage = tempState.WaterCoverage[worldData.SurfaceWaterLayer],
			FloraCoverage = tempState.FloraCoverage,
			Roughness = lastState.Roughness,
			IceFriction = worldData.WindIceFriction,
			TerrainFrictionMin = worldData.WindTerrainFrictionMin,
			TerrainFrictionMax = worldData.WindTerrainFrictionMax,
			FloraFriction = worldData.WindFloraFriction,
			WaterFriction = worldData.WindWaterFriction,
			MaxTerrainRoughness = worldData.MaxTerrainRoughnessForWindFriction
		}, tickJobHandle);

		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			int layerIndex = worldData.AirLayer0 + j;
			energyJobHandles[layerIndex] = SimJob.Schedule(new AccelerationAirJob()
			{
				Velocity = nextState.AirVelocity[j],
				Force = tempState.AirAcceleration[j],

				LastVelocity = lastState.AirVelocity[j],
				Friction = tempState.WindFriction,
				Pressure = tempState.AirPressure[j],
				AirMass = tempState.AirMass[j],
				TemperaturePotential = lastState.AirTemperaturePotential[j],
				NewTemperaturePotential = nextState.AirTemperaturePotential[j],
				VaporMass = lastState.AirVapor[j],
				LayerMiddle = tempState.LayerMiddle[j],
				Neighbors = staticState.Neighbors,
				NeighborDiffInverse = staticState.NeighborDiffInverse,
				Positions = staticState.SphericalPosition,
				PlanetRadius = staticState.PlanetRadius,
				Gravity = lastState.PlanetState.Gravity,
				GravityInverse = 1.0f / lastState.PlanetState.Gravity,
				NewUpTemperaturePotential = nextState.AirTemperaturePotential[j + 1],
				UpTemperaturePotential = lastState.AirTemperaturePotential[j + 1],
				UpHumidity = lastState.AirVapor[j + 1],
				UpAirMass = tempState.AirMass[j + 1],
				UpLayerMiddle = tempState.LayerMiddle[j + 1],
				NewDownTemperaturePotential = nextState.AirTemperaturePotential[j - 1],
				DownTemperaturePotential = lastState.AirTemperaturePotential[j - 1],
				DownHumidity = lastState.AirVapor[j - 1],
				DownAirMass = tempState.AirMass[j - 1],
				DownLayerMiddle = tempState.LayerMiddle[j - 1],
				IsTop = j == worldData.AirLayers - 2,
				IsBottom = j == 1,
				FrictionCoefficient = j == 1 ? 1 : 0,
				SecondsPerTick = worldData.SecondsPerTick

			},
			JobHandle.CombineDependencies(airTerrainFrictionJobHandle, 
			JobHandle.CombineDependencies(energyJobHandles[layerIndex], energyJobHandles[layerIndex - 1],energyJobHandles[layerIndex + 1])));
		}

		#endregion

		// Wind and currents move temperature and trace elements horizontally
		// Air, Water, Cloud
		// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
		#region Advection


		if (settings.MakeAirIncompressible)
		{
			for (int j = 1; j < worldData.AirLayers - 1; j++)
			{
				int layer = worldData.AirLayer0 + j;
				energyJobHandles[layer] = SimJob.Schedule(new ResolveAdvectionConflict()
				{
					NewVelocity = tempState.AirVelocityConflictFree[j],
					Neighbors = staticState.Neighbors,
					Position = staticState.SphericalPosition,
					Velocity = nextState.AirVelocity[j],
					VelocityAbove = nextState.AirVelocity[j + 1],
					VelocityBelow = nextState.AirVelocity[j - 1],
					LayerHeight = tempState.LayerHeight[j],
					LayerHeightAbove = tempState.LayerHeight[j + 1],
					LayerHeightBelow = tempState.LayerHeight[j - 1],
					Mass = tempState.AirMass[j],
					MassAbove = tempState.AirMass[j + 1],
					MassBelow = tempState.AirMass[j - 1],
					IsBottom = j == 1,
					IsTop = j == worldData.AirLayers - 2,
					SecondsPerTick = worldData.SecondsPerTick,
				}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer - 1], energyJobHandles[layer + 1]));
			}
			for (int j = 1; j < worldData.AirLayers - 1; j++)
			{
				int layer = worldData.AirLayer0 + j;
				energyJobHandles[layer] = SimJob.Schedule(new GetVectorDestCoordsVerticalJob()
				{
					Destination = tempState.DestinationAir[j],
					Neighbors = staticState.Neighbors,
					Position = staticState.SphericalPosition,
					Velocity = tempState.AirVelocityConflictFree[j],
					LayerHeight = tempState.LayerHeight[j],
					Mass = tempState.AirMass[j],
					MassAbove = tempState.AirMass[j + 1],
					MassBelow = tempState.AirMass[j - 1],
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick,
					MaxWindMove = staticState.CellRadius * 0.9f,
				}, energyJobHandles[layer]);
			}
			for (int j = 1; j < worldData.AirLayers - 1; j++)
			{
				int layer = worldData.AirLayer0 + j;
				energyJobHandles[layer] = SimJob.Schedule(new GetDivergenceJob()
				{
					Divergence = tempState.DivergenceAir[j],
					Destination = tempState.DestinationAir[j],
					DestinationAbove = tempState.DestinationAir[j + 1],
					DestinationBelow = tempState.DestinationAir[j - 1],
					Neighbors = staticState.Neighbors,
					Mass = tempState.AirMass[j],
					MassAbove = tempState.AirMass[j + 1],
					MassBelow = tempState.AirMass[j - 1],
					IsBottom = j == 1,
					IsTop = j == worldData.AirLayers - 2,
				}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer - 1], energyJobHandles[layer + 1]));
			}

			// Calculate Pressure gradient field
			var divergenceJobHandle = default(JobHandle);
			for (int i = 1; i < worldData.AirLayers - 1; i++)
			{
				divergenceJobHandle = JobHandle.CombineDependencies(divergenceJobHandle, energyJobHandles[worldData.AirLayer0 + i]);
			}
			for (int a = 0; a < 20; a++)
			{
				for (int i = 1; i < worldData.AirLayers - 1; i++)
				{
					bool isTop = i == worldData.AirLayers - 2;
					bool isBottom = i == 1;
					var dpj = new GetDivergencePressureJob()
					{
						Pressure = tempState.DivergencePressureAir[i],
						Divergence = tempState.DivergenceAir[i],
						PressureAbove = tempState.DivergencePressureAir[i + 1],
						PressureBelow = tempState.DivergencePressureAir[i - 1],
						IsTop = i == worldData.AirLayers - 2,
						IsBottom = i == 1,
						Neighbors = staticState.Neighbors
					};
					divergenceJobHandle = dpj.Schedule(_cellCount, divergenceJobHandle);
				}
			}

			for (int i = 1; i < worldData.AirLayers - 1; i++)
			{
				int layer = worldData.AirLayer0 + i;
				energyJobHandles[worldData.AirLayer0 + i] = SimJob.Schedule(new GetDivergenceFreeFieldJob()
				{
					VelocityOut = nextState.AirVelocity[i],
					VelocityIn = tempState.AirVelocityConflictFree[i],
					Pressure = tempState.DivergencePressureAir[i],
					PressureAbove = tempState.DivergencePressureAir[i + 1],
					PressureBelow = tempState.DivergencePressureAir[i - 1],
					LayerHeight = tempState.LayerHeight[i],
					Neighbors = staticState.Neighbors,
					NeighborTangent = staticState.NeighborTangent,
					Positions = staticState.SphericalPosition,
					Mass = tempState.AirMass[i],
					IsBottom = i == 1,
					IsTop = i == worldData.AirLayers - 2,
					TicksPerSecond = worldData.TicksPerSecond
				}, divergenceJobHandle);
			}
		}


		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			int layer = worldData.AirLayer0 + i;
			// TODO: we need to deal with coriolis for differently:
			// either apply it as a true force (although in my experience this causes velocity to spin)
			// or we deflect it after advecting it
			energyJobHandles[layer] = SimJob.Schedule(new GetVectorDestCoordsVerticalJob()
			{
				Destination = tempState.DestinationAir[i],
				Neighbors = staticState.Neighbors,
				Position = staticState.SphericalPosition,
				Velocity = nextState.AirVelocity[i],
				LayerHeight = tempState.LayerHeight[i],
				Mass = tempState.AirMass[i],
				MassAbove = tempState.AirMass[i + 1],
				MassBelow = tempState.AirMass[i - 1],
				PlanetRadius = staticState.PlanetRadius,
				SecondsPerTick = worldData.SecondsPerTick,
				MaxWindMove = staticState.CellRadius * 0.9f,
			}, energyJobHandles[layer]);
		}

		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			int layer = worldData.AirLayer0 + i;
			energyJobHandles[layer] = SimJob.Schedule(new AdvectionAirJob()
			{
				Delta = tempState.AdvectionAir[i],
				Temperature = nextState.AirTemperaturePotential[i],
				TemperatureAbove = nextState.AirTemperaturePotential[i + 1],
				TemperatureBelow = nextState.AirTemperaturePotential[i - 1],
				AirMass = tempState.AirMass[i],
				AirMassAbove = tempState.AirMass[i + 1],
				AirMassBelow = tempState.AirMass[i - 1],
				Vapor = nextState.AirVapor[i],
				VaporAbove = nextState.AirVapor[i + 1],
				VaporBelow = nextState.AirVapor[i - 1],
				CarbonDioxide = nextState.AirCarbon[i],
				CarbonDioxideAbove = nextState.AirCarbon[i + 1],
				CarbonDioxideBelow = nextState.AirCarbon[i - 1],
				Dust = nextState.Dust[i],
				DustAbove = nextState.Dust[i + 1],
				DustBelow = nextState.Dust[i - 1],
				Velocity = nextState.AirVelocity[i],
				VelocityAbove = nextState.AirVelocity[i + 1],
				VelocityBelow = nextState.AirVelocity[i - 1],
				Neighbors = staticState.Neighbors,
				Destination = tempState.DestinationAir[i],
				DestinationAbove = tempState.DestinationAir[i + 1],
				DestinationBelow = tempState.DestinationAir[i - 1],
				LayerMiddle = tempState.LayerMiddle[i],
				LayerMiddleAbove = tempState.LayerMiddle[i + 1],
				LayerMiddleBelow = tempState.LayerMiddle[i - 1],
				Positions = staticState.SphericalPosition,
				NeighborDistInverse = staticState.NeighborDistInverse,
				CoriolisMultiplier = staticState.CoriolisMultiplier,
				CoriolisTerm = coriolisTerm,
				SecondsPerTick = worldData.SecondsPerTick,
				TicksPerSecond = worldData.TicksPerSecond
			}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer - 1], energyJobHandles[layer + 1]));
		}

		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			int layer = worldData.AirLayer0 + i;
			applyAdvectionAirJobHandles[i] = SimJob.Schedule(new ApplyAdvectionAirJob()
			{
				Advection = tempState.AdvectionAir[i],
				Vapor = nextState.AirVapor[i],
				Dust = nextState.Dust[i],
				CarbonDioxide = nextState.AirCarbon[i],
				Temperature = nextState.AirTemperaturePotential[i],
				AirVelocity = nextState.AirVelocity[i],
			}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer - 1], energyJobHandles[layer + 1]));
		}
		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			energyJobHandles[worldData.AirLayer0 + i] = applyAdvectionAirJobHandles[i];
		}

		#endregion

		// Diffuse from last time step
		// Air, Water, Cloud
		#region Diffusion

		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			int layer = worldData.AirLayer0 + i;
			// TODO: is it a problem that we are using the dependent variables from last frame while referencing our newly calculated next frame values for temperature and such?
			energyJobHandles[layer] = SimJob.Schedule(new DiffusionAirJob()
			{
				Delta = tempState.DiffusionAir[i],

				AirMass = tempState.AirMass[i],
				AirMassAbove = tempState.AirMass[i + 1],
				AirMassBelow = tempState.AirMass[i - 1],
				Temperature = nextState.AirTemperaturePotential[i],
				TemperatureAbove = nextState.AirTemperaturePotential[i + 1],
				TemperatureBelow = nextState.AirTemperaturePotential[i - 1],
				Vapor = nextState.AirVapor[i],
				VaporAbove = nextState.AirVapor[i + 1],
				VaporBelow = nextState.AirVapor[i - 1],
				CarbonDioxide = nextState.AirCarbon[i],
				CarbonDioxideAbove = nextState.AirCarbon[i + 1],
				CarbonDioxideBelow = nextState.AirCarbon[i - 1],
				Dust = nextState.Dust[i],
				DustAbove = nextState.Dust[i + 1],
				DustBelow = nextState.Dust[i - 1],
				Velocity = nextState.AirVelocity[i],
				AirVelocityAbove = nextState.AirVelocity[i + 1],
				AirVelocityBelow = nextState.AirVelocity[i - 1],
				Neighbors = staticState.Neighbors,
				LayerHeight = tempState.LayerHeight[i],
				NeighborDistInverse = staticState.NeighborDistInverse,
				LayerElevationAbove = tempState.LayerElevation[i + 1],
				LayerHeightAbove = tempState.LayerHeight[i + 1],
				LayerElevationBelow = tempState.LayerElevation[i - 1],
				LayerHeightBelow = tempState.LayerHeight[i - 1],
				IsTop = i == worldData.AirLayers - 2,
				IsBottom = i == 1,
				DiffusionCoefficientHorizontal = worldData.AirDiffusionCoefficientHorizontal,
				DiffusionCoefficientVertical = worldData.AirDiffusionCoefficientVertical,
				CellSurfaceArea = staticState.CellSurfaceArea,
				CellCircumference = staticState.CellCircumference
			}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer - 1], energyJobHandles[layer + 1]));
		}

		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			int layer = worldData.AirLayer0 + i;
			applyAdvectionAirJobHandles[i] = SimJob.Schedule(new ApplyAdvectionAirJob()
			{
				Advection = tempState.DiffusionAir[i],
				Vapor = nextState.AirVapor[i],
				Dust = nextState.Dust[i],
				CarbonDioxide = nextState.AirCarbon[i],
				Temperature = nextState.AirTemperaturePotential[i],
				AirVelocity = nextState.AirVelocity[i],
			}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer - 1], energyJobHandles[layer + 1]));
		}

		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			energyJobHandles[worldData.AirLayer0 + i] = applyAdvectionAirJobHandles[i];
		}

		#endregion

		#endregion

		#region Water and Cloud Advection and Diffusion

		// Buoyancy, Updrafts, and mixing occur across air layers and water layers
		// TODO: add an empty air layer on top and bottom so we can calculate up/down diffusion in a single step 
		// Temperature and trace elements diffuse into neighboring horizontal cells based on a diffusion constant
		// Air, Water, Cloud
		#region Update Velocity

		energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new WaterSurfaceFrictionJob()
		{
			Force = tempState.WaterFriction,

			Position = staticState.SphericalPosition,
			Current = lastState.WaterVelocity[worldData.SurfaceWaterLayer],
			AirVelocityUp = lastState.AirVelocity[worldData.SurfaceAirLayer],
			AirVelocityDown = lastState.WaterVelocity[worldData.SurfaceWaterLayer - 1],
			LayerHeight = tempState.WaterLayerHeight[worldData.SurfaceWaterLayer],
			CoriolisMultiplier = staticState.CoriolisMultiplier,
			FrictionCoefficientUp = worldData.WindToWaterCurrentFrictionCoefficient,
			FrictionCoefficientDown = 0, // TODO: do we want to add a frictional force between layers of water?
			CoriolisTerm = coriolisTerm,
			WaterSurfaceFrictionDepth = worldData.WaterSurfaceFrictionDepth,
			SecondsPerTick = worldData.SecondsPerTick
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.SurfaceAirLayerGlobal]));

		// TODO: since I moved acceleration up here, we will need to add a buoyancy job down after the temperature adjustment to incorporate positive buoyancy in the vertical velocity
		// (see usage of "next state" velocity
		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			energyJobHandles[worldData.WaterLayer0 + j] = SimJob.Schedule(new AccelerationWaterJob()
			{
				Velocity = nextState.WaterVelocity[j],

				LastVelocity = lastState.WaterVelocity[j],
				Positions = staticState.SphericalPosition,
				Neighbors = staticState.Neighbors,
				NeighborDiffInverse = staticState.NeighborDiffInverse,
				WaterDensity = tempState.WaterDensity[j],
				WaterPressure = tempState.WaterPressure[j],
				LayerDepth = tempState.WaterLayerDepth[j],
				LayerHeight = tempState.WaterLayerHeight[j],
				UpWaterDensity = tempState.WaterDensity[j + 1],
				UpWaterPressure = tempState.WaterPressure[j + 1],
				UpLayerDepth = tempState.WaterLayerDepth[j + 1],
				UpLayerHeight = tempState.WaterLayerHeight[j + 1],
				DownWaterDensity = tempState.WaterDensity[j - 1],
				DownWaterPressure = tempState.WaterPressure[j - 1],
				DownLayerDepth = tempState.WaterLayerDepth[j - 1],
				DownLayerHeight = tempState.WaterLayerHeight[j - 1],
				SurfaceElevation = tempState.LayerElevation[worldData.SurfaceAirLayer],
				Friction = tempState.WaterFriction,
				FrictionCoefficient = j == worldData.SurfaceWaterLayer ? 1 : 0,
				Gravity = lastState.PlanetState.Gravity,
				PlanetRadius = staticState.PlanetRadius,
				SecondsPerTick = worldData.SecondsPerTick

			}, JobHandle.CombineDependencies( energyJobHandles[worldData.WaterLayer0 + j], energyJobHandles[worldData.SurfaceWaterLayerGlobal]));
		}


		#endregion

		// Wind and currents move temperature and trace elements horizontally
		// Air, Water, Cloud
		// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
		#region Advection


		if (settings.MakeWaterIncompressible)
		{
			for (int j = 1; j < worldData.WaterLayers - 1; j++)
			{
				int layer = worldData.WaterLayer0 + j;
				energyJobHandles[layer] = SimJob.Schedule(new ResolveAdvectionConflict()
				{
					NewVelocity = tempState.WaterVelocityConflictFree[j],
					Neighbors = staticState.Neighbors,
					Position = staticState.SphericalPosition,
					Velocity = nextState.WaterVelocity[j],
					VelocityAbove = nextState.WaterVelocity[j + 1],
					VelocityBelow = nextState.WaterVelocity[j - 1],
					LayerHeight = tempState.WaterLayerHeight[j],
					LayerHeightAbove = tempState.WaterLayerHeight[j + 1],
					LayerHeightBelow = tempState.WaterLayerHeight[j - 1],
					Mass = nextState.WaterMass[j],
					MassAbove = nextState.WaterMass[j + 1],
					MassBelow = nextState.WaterMass[j - 1],
					IsBottom = j == 1,
					IsTop = j == worldData.WaterLayers - 2,
					SecondsPerTick = worldData.SecondsPerTick,
				}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer - 1], energyJobHandles[layer + 1]));
			}

			for (int j = 1; j < worldData.WaterLayers - 1; j++)
			{
				int layer = worldData.WaterLayer0 + j;
				energyJobHandles[layer] = SimJob.Schedule(new GetVectorDestCoordsVerticalJob()
				{
					Destination = tempState.DestinationWater[j],
					Neighbors = staticState.Neighbors,
					Position = staticState.SphericalPosition,
					Velocity = tempState.WaterVelocityConflictFree[j],
					LayerHeight = tempState.WaterLayerHeight[j],
					Mass = nextState.WaterMass[j],
					MassAbove = nextState.WaterMass[j + 1],
					MassBelow = nextState.WaterMass[j - 1],
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick,
					MaxWindMove = staticState.CellRadius * 0.9f,
				}, energyJobHandles[layer]);
			}
			for (int j = 1; j < worldData.WaterLayers - 1; j++)
			{
				int layer = worldData.WaterLayer0 + j;
				energyJobHandles[layer] = SimJob.Schedule(new GetDivergenceJob()
				{
					Divergence = tempState.DivergenceWater[j],
					Destination = tempState.DestinationWater[j],
					DestinationAbove = tempState.DestinationWater[j + 1],
					DestinationBelow = tempState.DestinationWater[j - 1],
					Neighbors = staticState.Neighbors,
					Mass = nextState.WaterMass[j],
					MassAbove = nextState.WaterMass[j + 1],
					MassBelow = nextState.WaterMass[j - 1],
					IsBottom = j == 1,
					IsTop = j == worldData.SurfaceWaterLayer,
				}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer - 1], energyJobHandles[layer + 1]));
			}

			// Calculate Pressure gradient field
			var divergenceJobHandle = default(JobHandle);
			for (int i = 1; i < worldData.WaterLayers - 1; i++)
			{
				divergenceJobHandle = JobHandle.CombineDependencies(divergenceJobHandle, energyJobHandles[worldData.WaterLayer0 + i]);
			}
			for (int a = 0; a < 20; a++)
			{
				for (int i = 1; i < worldData.WaterLayers - 1; i++)
				{
					bool isTop = i == worldData.WaterLayers - 2;
					bool isBottom = i == 1;
					var dpj = new GetDivergencePressureJob()
					{
						Pressure = tempState.DivergencePressureWater[i],
						Divergence = tempState.DivergenceWater[i],
						PressureAbove = tempState.DivergencePressureWater[i + 1],
						PressureBelow = tempState.DivergencePressureWater[i - 1],
						IsTop = i == worldData.SurfaceWaterLayer,
						IsBottom = i == 1,
						Neighbors = staticState.Neighbors
					};
					divergenceJobHandle = dpj.Schedule(_cellCount, divergenceJobHandle);
				}
			}

			for (int i = 1; i < worldData.WaterLayers - 1; i++)
			{
				int layer = worldData.WaterLayer0 + i;
				energyJobHandles[worldData.WaterLayer0 + i] = SimJob.Schedule(new GetDivergenceFreeFieldJob()
				{
					VelocityOut = nextState.WaterVelocity[i],
					VelocityIn = tempState.WaterVelocityConflictFree[i],
					Pressure = tempState.DivergencePressureWater[i],
					PressureAbove = tempState.DivergencePressureWater[i + 1],
					PressureBelow = tempState.DivergencePressureWater[i - 1],
					LayerHeight = tempState.WaterLayerHeight[i],
					Neighbors = staticState.Neighbors,
					NeighborTangent = staticState.NeighborTangent,
					Positions = staticState.SphericalPosition,
					Mass = nextState.WaterMass[i],
					IsBottom = i == 1,
					IsTop = i == worldData.SurfaceWaterLayer,
					TicksPerSecond = worldData.TicksPerSecond
				}, divergenceJobHandle);
			}
		}



		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			int layer = worldData.WaterLayer0 + j;
			energyJobHandles[layer] = SimJob.Schedule(new GetVectorDestCoordsVerticalJob()
			{
				Destination = tempState.DestinationWater[j],
				Neighbors = staticState.Neighbors,
				Position = staticState.SphericalPosition,
				Velocity = nextState.WaterVelocity[j],
				LayerHeight = tempState.WaterLayerHeight[j],
				Mass = nextState.WaterMass[j],
				MassAbove = nextState.WaterMass[j + 1],
				MassBelow = nextState.WaterMass[j - 1],
				PlanetRadius = staticState.PlanetRadius,
				SecondsPerTick = worldData.SecondsPerTick,
				MaxWindMove = staticState.CellRadius * 0.9f,
			}, energyJobHandles[layer]);
		}

		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			int layer = worldData.WaterLayer0 + j;
			energyJobHandles[layer] = SimJob.Schedule(new AdvectionWaterJob()
			{
				Delta = tempState.AdvectionWater[j],
				Destination = tempState.DestinationWater[j],
				DestinationAbove = tempState.DestinationWater[j + 1],
				DestinationBelow = tempState.DestinationWater[j - 1],
				Velocity = nextState.WaterVelocity[j],
				VelocityAbove = nextState.WaterVelocity[j + 1],
				VelocityBelow = nextState.WaterVelocity[j - 1],
				Temperature = nextState.WaterTemperature[j],
				TemperatureAbove = nextState.WaterTemperature[j + 1],
				TemperatureBelow = nextState.WaterTemperature[j - 1],
				Mass = nextState.WaterMass[j],
				MassAbove = nextState.WaterMass[j + 1],
				MassBelow = nextState.WaterMass[j - 1],
				Salt = nextState.SaltMass[j],
				SaltAbove = nextState.SaltMass[j + 1],
				SaltBelow = nextState.SaltMass[j - 1],
				Carbon = nextState.WaterCarbon[j],
				CarbonAbove = nextState.WaterCarbon[j + 1],
				CarbonBelow = nextState.WaterCarbon[j - 1],
				PlanktonMass = nextState.PlanktonMass[j],
				PlanktonGlucose = nextState.PlanktonGlucose[j],
				Positions = staticState.SphericalPosition,
				Neighbors = staticState.Neighbors,
				CoriolisMultiplier = staticState.CoriolisMultiplier,
				CoriolisTerm = coriolisTerm,
				SecondsPerTick = worldData.SecondsPerTick
			}, JobHandle.CombineDependencies(groundWaterJob,
				JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer + 1], energyJobHandles[layer - 1])));
		}

		energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(new GetVectorDestCoordsJob()
		{
			Destination = tempState.DestinationCloud,
			Neighbors = staticState.Neighbors,
			Position = staticState.SphericalPosition,
			Velocity = tempState.CloudVelocity,
			PlanetRadius = staticState.PlanetRadius,
			SecondsPerTick = worldData.SecondsPerTick,
			MaxWindMove = staticState.CellRadius * 0.9f,
		}, energyJobHandles[worldData.CloudLayer]);
		energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(new AdvectionCloudJob()
		{
			Delta = tempState.AdvectionCloud,
			Destination = tempState.DestinationCloud,
			Mass = nextState.CloudMass,
			Temperature = nextState.CloudTemperature,
			DropletMass = nextState.CloudDropletMass,
			Neighbors = staticState.Neighbors,
		}, energyJobHandles[worldData.CloudLayer]);

		energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(new ApplyAdvectionCloudJob()
		{
			Advection = tempState.AdvectionCloud,
			CloudMass = nextState.CloudMass,
			Temperature = nextState.CloudTemperature,
			DropletMass = nextState.CloudDropletMass,
		}, energyJobHandles[worldData.CloudLayer]);

		for (int i = 1; i < worldData.WaterLayers - 1; i++)
		{
			int layer = worldData.WaterLayer0 + i;
			applyAdvectionWaterJobHandles[i] = SimJob.Schedule(new ApplyAdvectionWaterJob()
			{
				Advection = tempState.AdvectionWater[i],
				SaltMass = nextState.SaltMass[i],
				CarbonMass = nextState.WaterCarbon[i],
				PlanktonMass = nextState.PlanktonMass[i],
				PlanktonGlucose = nextState.PlanktonGlucose[i],
				Temperature = nextState.WaterTemperature[i],
				Velocity = nextState.WaterVelocity[i],
				WaterMass = nextState.WaterMass[i]
			}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer + 1], energyJobHandles[layer - 1]));
		}

		for (int i = 1; i < worldData.WaterLayers - 1; i++)
		{
			energyJobHandles[worldData.WaterLayer0 + i] = applyAdvectionWaterJobHandles[i];
		}

		#endregion

		// Diffuse from last time step
		// Air, Water, Cloud
		#region Diffusion

		for (int i = 1; i < worldData.WaterLayers - 1; i++)
		{
			int layer = worldData.WaterLayer0 + i;
			energyJobHandles[layer] = SimJob.Schedule(new DiffusionWaterJob()
			{
				Delta = tempState.DiffusionWater[i],

				Temperature = nextState.WaterTemperature[i],
				TemperatureAbove = nextState.WaterTemperature[i + 1],
				TemperatureBelow = nextState.WaterTemperature[i - 1],
				SaltMass = nextState.SaltMass[i],
				SaltMassAbove = nextState.SaltMass[i + 1],
				SaltMassBelow = nextState.SaltMass[i - 1],
				PlanktonMass = nextState.PlanktonMass[i],
				PlanktonGlucose = nextState.PlanktonGlucose[i],
				CarbonMass = nextState.WaterCarbon[i],
				CarbonMassAbove = nextState.WaterCarbon[i + 1],
				CarbonMassBelow = nextState.WaterCarbon[i - 1],
				Velocity = nextState.WaterVelocity[i],
				VelocityAbove = nextState.WaterVelocity[i + 1],
				VelocityBelow = nextState.WaterVelocity[i - 1],
				WaterMass = nextState.WaterMass[i],
				WaterMassAbove = nextState.WaterMass[i + 1],
				WaterMassBelow = nextState.WaterMass[i - 1],
				LayerHeight = tempState.LayerHeight[i],
				LayerHeightAbove = tempState.LayerHeight[i + 1],
				LayerHeightBelow = tempState.LayerHeight[i - 1],
				NeighborDistInverse = staticState.NeighborDistInverse,
				Neighbors = staticState.Neighbors,
				DiffusionCoefficientHorizontal = worldData.WaterDiffusionCoefficientHorizontal,
				DiffusionCoefficientVertical = worldData.WaterDiffusionCoefficientVertical,
				CellSurfaceArea = staticState.CellSurfaceArea,
				CellCircumference = staticState.CellCircumference
			}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer + 1], energyJobHandles[layer - 1]));
		}

		energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(new DiffusionCloudJob()
		{
			Delta = tempState.DiffusionCloud,

			LastMass = nextState.CloudMass,
			LastTemperature = nextState.CloudTemperature,
			LastDropletMass = nextState.CloudDropletMass,
			Neighbors = staticState.Neighbors,
			DiffusionCoefficient = worldData.CloudDiffusionCoefficient,
		}, energyJobHandles[worldData.CloudLayer]);

		energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(new ApplyAdvectionCloudJob()
		{
			Advection = tempState.DiffusionCloud,
			CloudMass = nextState.CloudMass,
			Temperature = nextState.CloudTemperature,
			DropletMass = nextState.CloudDropletMass,
		}, energyJobHandles[worldData.CloudLayer]);

		for (int i = 1; i < worldData.WaterLayers - 1; i++)
		{
			int layer = worldData.WaterLayer0 + i;
			applyAdvectionWaterJobHandles[i] = SimJob.Schedule(new ApplyAdvectionWaterJob()
			{
				Advection = tempState.DiffusionWater[i],
				SaltMass = nextState.SaltMass[i],
				CarbonMass = nextState.WaterCarbon[i],
				PlanktonMass = nextState.PlanktonMass[i],
				PlanktonGlucose = nextState.PlanktonGlucose[i],
				Temperature = nextState.WaterTemperature[i],
				Velocity = nextState.WaterVelocity[i],
				WaterMass = nextState.WaterMass[i]
			}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer + 1], energyJobHandles[layer - 1]));
		}

		for (int i = 1; i < worldData.WaterLayers - 1; i++)
		{
			energyJobHandles[worldData.WaterLayer0 + i] = applyAdvectionWaterJobHandles[i];
		}

		#endregion

		#endregion


		tickJobHandle = JobHandle.CombineDependencies(
			groundWaterJob,
			JobHandle.CombineDependencies(energyJobHandles)
			);

		// TODO: we really just want to update depths
		#region Update dependent variables

		tickJobHandle = TempState.Update(SimJob, ref nextState, ref tempState, ref worldData, ref staticState, tickJobHandle);

		#endregion


		#region Flow

		if (settings.WaterSurfaceFlowEnabled)
		{
			// TODO: surface elevation is inaccurate now, we should recalculate (and use water surfae, not ice surface)
			tickJobHandle = NeighborJob.Schedule(new UpdateFlowVelocityJob()
			{
				Flow = nextState.FlowWater,

				LastFlow = lastState.FlowWater,
				SurfaceElevation = tempState.LayerElevation[worldData.SurfaceAirLayer],
				WaterDepth = tempState.WaterLayerHeight[worldData.SurfaceWaterLayer],
				NeighborDistInverse = staticState.NeighborDistInverse,
				Neighbors = staticState.Neighbors,
				Gravity = nextState.PlanetState.Gravity,
				SecondsPerTick = worldData.SecondsPerTick,
				Damping = worldData.SurfaceWaterFlowDamping,
				ViscosityInverse = 1.0f - worldData.WaterViscosity,
			}, tickJobHandle);

			tickJobHandle = SimJob.Schedule(new SumOutgoingFlowJob()
			{
				OutgoingFlow = tempState.OutgoingFlowWater,
				Flow = nextState.FlowWater
			}, tickJobHandle);

			tickJobHandle = NeighborJob.Schedule(new LimitOutgoingFlowJob()
			{
				Flow = nextState.FlowWater,
				FlowPercent = tempState.FlowPercentWater,

				OutgoingFlow = tempState.OutgoingFlowWater,
				Neighbors = staticState.Neighbors,
				WaterDepth = tempState.WaterLayerHeight[worldData.SurfaceWaterLayer]
			}, tickJobHandle);

			tickJobHandle = SimJob.Schedule(new ApplyFlowWaterJob()
			{
				Delta = tempState.DiffusionWater[worldData.SurfaceWaterLayer],

				Mass = nextState.WaterMass[worldData.SurfaceWaterLayer],
				Velocity = nextState.WaterVelocity[worldData.SurfaceWaterLayer],
				Carbon = nextState.WaterCarbon[worldData.SurfaceWaterLayer],
				PlanktonMass = nextState.PlanktonMass[worldData.SurfaceWaterLayer],
				PlanktonGlucose = nextState.PlanktonGlucose[worldData.SurfaceWaterLayer],
				Salt = nextState.SaltMass[worldData.SurfaceWaterLayer],
				Temperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				FlowPercent = tempState.FlowPercentWater,
				Positions = staticState.SphericalPosition,
				Neighbors = staticState.Neighbors,
				CoriolisMultiplier = staticState.CoriolisMultiplier,
				CoriolisTerm = coriolisTerm,
				SecondsPerTick = worldData.SecondsPerTick
			}, tickJobHandle);

			tickJobHandle = SimJob.Schedule(new ApplyAdvectionWaterJob()
			{
				WaterMass = nextState.WaterMass[worldData.SurfaceWaterLayer],
				SaltMass = nextState.SaltMass[worldData.SurfaceWaterLayer],
				CarbonMass = nextState.WaterCarbon[worldData.SurfaceWaterLayer],
				PlanktonMass = nextState.PlanktonMass[worldData.SurfaceWaterLayer],
				PlanktonGlucose = nextState.PlanktonGlucose[worldData.SurfaceWaterLayer],
				Temperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				Velocity = nextState.WaterVelocity[worldData.SurfaceWaterLayer],

				Advection = tempState.DiffusionWater[worldData.SurfaceWaterLayer],
			}, tickJobHandle);
		}

		#endregion

		#region Rebalance

		for (int i = worldData.SurfaceWaterLayer; i >= 2; i--)
		{

			float maxDepth;
			float minDepth;
			// TODO: this doesn't handle more than 3 layers
			if (i == worldData.SurfaceWaterLayer)
			{
				maxDepth = worldData.SurfaceWaterDepth;
				minDepth = worldData.SurfaceWaterDepth - 10;
			} else
			{
				maxDepth = worldData.ThermoclineDepth;
				minDepth = worldData.ThermoclineDepth - 10;
			}
			var h = SimJob.Schedule(new RebalanceWaterLayersLimitJob()
			{
				Delta1 = tempState.DiffusionWater[i],
				Delta2 = tempState.DiffusionWater[i - 1],

				Mass1 = nextState.WaterMass[i],
				Salt1 = nextState.SaltMass[i],
				Carbon1 = nextState.WaterCarbon[i],
				PlanktonMass1 = nextState.PlanktonMass[i],
				PlanktonGlucose1 = nextState.PlanktonGlucose[i],
				Temperature1 = nextState.WaterTemperature[i],
				Velocity1 = nextState.WaterVelocity[i],

				Mass2 = nextState.WaterMass[i - 1],
				Salt2 = nextState.SaltMass[i - 1],
				Carbon2 = nextState.WaterCarbon[i - 1],
				Temperature2 = nextState.WaterTemperature[i - 1],
				Velocity2 = nextState.WaterVelocity[i - 1],

				MaxDepth1 = maxDepth,
				MinDepth1 = minDepth,
			}, tickJobHandle);
			var rebalanceWaterJobHandle1 = SimJob.Schedule(new ApplyAdvectionWaterJob()
			{
				WaterMass = nextState.WaterMass[i],
				SaltMass = nextState.SaltMass[i],
				CarbonMass = nextState.WaterCarbon[i],
				PlanktonMass = nextState.PlanktonMass[i],
				PlanktonGlucose = nextState.PlanktonGlucose[i],
				Temperature = nextState.WaterTemperature[i],
				Velocity = nextState.WaterVelocity[i],

				Advection = tempState.DiffusionWater[i],
			}, h);
			var rebalanceWaterJobHandle2 = SimJob.Schedule(new ApplyAdvectionWaterJob()
			{
				WaterMass = nextState.WaterMass[i - 1],
				SaltMass = nextState.SaltMass[i - 1],
				CarbonMass = nextState.WaterCarbon[i - 1],
				PlanktonMass = nextState.PlanktonMass[i - 1],
				PlanktonGlucose = nextState.PlanktonGlucose[i - 1],
				Temperature = nextState.WaterTemperature[i - 1],
				Velocity = nextState.WaterVelocity[i - 1],

				Advection = tempState.DiffusionWater[i - 1],
			}, h);

			tickJobHandle = JobHandle.CombineDependencies(rebalanceWaterJobHandle1, rebalanceWaterJobHandle2);
		}


		#endregion

		#region Update dependent variables

		tickJobHandle = TempState.Update(SimJob, ref nextState, ref tempState, ref worldData, ref staticState, tickJobHandle);

		#endregion



		return tickJobHandle;

	}

	private void DoEnergyCycle(
		NativeArray<JobHandle> energyJobHandles,
		JobHandle lastJobHandle,
		ref SimState nextState,
		ref SimState lastState,
		ref TempState tempState,
		ref StaticState staticState,
		ref WorldData worldData,
		ref SimSettings settings
		)
	{
		#region Init Solar Radiation Per Cell
		var solarInJobHandle = SimJob.Schedule(new SolarRadiationJob()
		{
			SolarRadiation = tempState.SolarRadiation,
			GeothermalRadiation = tempState.GeothermalRadiation,
			DisplaySolarRadiation = tempState.DisplaySolarRadiation,
			AlbedoSlope = tempState.AlbedoSlope,

			SphericalPosition = staticState.SphericalPosition,
			IncomingSolarRadiation = lastState.PlanetState.SolarRadiation * worldData.SecondsPerTick,
			IncomingGeothermalRadiation = lastState.PlanetState.GeothermalHeat * worldData.SecondsPerTick,
			PlanetRotation = quaternion.Euler(lastState.PlanetState.Rotation),
			SunToPlanetDir = math.normalize(lastState.PlanetState.Position),
		}, lastJobHandle);

		#endregion

		// Calculate emissivity per cell
		// TODO: combine cloud emissivity with cloud conduction
		// TODO: we can probably combine this step with the thermal radiation step
		#region Emissivity Per Cell

		JobHandle[] emissivityJobHandles = new JobHandle[worldData.LayerCount];
		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			int layerIndex = worldData.AirLayer0 + j;
			emissivityJobHandles[layerIndex] = SimJob.Schedule(new EmissivityAirJob()
			{
				Emissivity = tempState.Emissivity[layerIndex],
				AirMass = tempState.AirMass[j],
				VaporMass = lastState.AirVapor[j],
				Dust = lastState.Dust[j],
				CarbonDioxide = lastState.AirCarbon[j],
				EmissivityAir = worldData.ThermalEmissivityAir,
				EmissivityWaterVapor = worldData.ThermalEmissivityWaterVapor,
				EmissivityDust = worldData.ThermalEmissivityDust,
				EmissivityCarbonDioxide = worldData.ThermalEmissivityCarbonDioxide,
				EmissivityOxygen = worldData.ThermalEmissivityOxygen
			}, lastJobHandle);
		}

		// we only do thermal radiation upwards for the surface layer of water,
		// for the bottom we rely on conduction with the terrain for heat transfer (although this might lead to an imbalance!)
		{
			emissivityJobHandles[worldData.SurfaceWaterLayer + worldData.WaterLayer0] = SimJob.Schedule(new EmissivityWaterJob()
			{
				Emissivity = tempState.Emissivity[worldData.SurfaceWaterLayerGlobal],
				WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
				SaltMass = lastState.SaltMass[worldData.SurfaceWaterLayer],
				EmissivitySalt = worldData.ThermalEmissivitySalt,
				EmissivityWater = worldData.ThermalEmissivityWater
			}, lastJobHandle);
		}
		emissivityJobHandles[worldData.TerrainLayer] = SimJob.Schedule(new EmissivityTerrainJob()
		{
			Emissivity = tempState.Emissivity[worldData.TerrainLayer],
			SoilFertility = lastState.GroundCarbon,
			EmissivityDirt = worldData.ThermalEmissivityDirt,
			EmissivitySand = worldData.ThermalEmissivitySand,
		}, lastJobHandle);

		#endregion

		// Calculate how much thermal radition is being emitted out of each layer
		#region Thermal Radiation
		JobHandle[] thermalOutJobHandles = new JobHandle[worldData.LayerCount];

		// ICE
		thermalOutJobHandles[worldData.IceLayer] = SimJob.Schedule(new ThermalEnergyRadiatedConstantEmissivityJob()
		{
			ThermalRadiationDelta = tempState.ThermalRadiationDelta[worldData.IceLayer],
			ThermalRadiationTransmittedUp = tempState.ThermalRadiationTransmittedUp[worldData.IceLayer],
			ThermalRadiationTransmittedDown = tempState.ThermalRadiationTransmittedDown[worldData.IceLayer],
			WindowRadiationTransmittedUp = tempState.WindowRadiationTransmittedUp[worldData.IceLayer],
			WindowRadiationTransmittedDown = tempState.WindowRadiationTransmittedDown[worldData.IceLayer],

			PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
			Emissivity = worldData.ThermalEmissivityIce,
			Energy = tempState.IceEnergy,
			Temperature = lastState.IceTemperature,
			SurfaceArea = tempState.IceCoverage,
			SecondsPerTick = worldData.SecondsPerTick
		}, lastJobHandle);


		// FLORA
		thermalOutJobHandles[worldData.FloraLayer] = SimJob.Schedule(new ThermalEnergyRadiatedConstantEmissivityJob()
		{
			ThermalRadiationDelta = tempState.ThermalRadiationDelta[worldData.FloraLayer],
			ThermalRadiationTransmittedUp = tempState.ThermalRadiationTransmittedUp[worldData.FloraLayer],
			ThermalRadiationTransmittedDown = tempState.ThermalRadiationTransmittedDown[worldData.FloraLayer],
			WindowRadiationTransmittedUp = tempState.WindowRadiationTransmittedUp[worldData.FloraLayer],
			WindowRadiationTransmittedDown = tempState.WindowRadiationTransmittedDown[worldData.FloraLayer],

			PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
			Emissivity = worldData.ThermalEmissivityFlora,
			Energy = tempState.FloraEnergy,
			Temperature = lastState.FloraTemperature,
			SurfaceArea = tempState.FloraCoverage,
			SecondsPerTick = worldData.SecondsPerTick
		}, lastJobHandle);


		// TERRAIN
		thermalOutJobHandles[worldData.TerrainLayer] = SimJob.Schedule(new ThermalEnergyRadiatedTerrainJob()
		{
			ThermalRadiationDelta = tempState.ThermalRadiationDelta[worldData.TerrainLayer],
			ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp[worldData.TerrainLayer],
			WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp[worldData.TerrainLayer],

			PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
			Emissivity = tempState.Emissivity[worldData.TerrainLayer],
			Temperature = lastState.GroundTemperature,
			SecondsPerTick = worldData.SecondsPerTick
		}, emissivityJobHandles[worldData.TerrainLayer]);


		// ATMOSPHERE
		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			int layer = worldData.AirLayer0 + j;
			thermalOutJobHandles[layer] = SimJob.Schedule(new ThermalEnergyRadiatedAirJob()
			{
				ThermalRadiationDelta = tempState.ThermalRadiationDelta[layer],
				ThermalRadiationTransmittedUp = tempState.ThermalRadiationTransmittedUp[layer],
				ThermalRadiationTransmittedDown = tempState.ThermalRadiationTransmittedDown[layer],
				WindowRadiationTransmittedUp = tempState.WindowRadiationTransmittedUp[layer],
				WindowRadiationTransmittedDown = tempState.WindowRadiationTransmittedDown[layer],

				PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
				Energy = tempState.AirPotentialEnergy[j],
				Emissivity = tempState.Emissivity[layer],
				TemperaturePotential = lastState.AirTemperaturePotential[j],
				LayerMiddle = tempState.LayerMiddle[j],
				SecondsPerTick = worldData.SecondsPerTick
			}, emissivityJobHandles[layer]);
		}

		// WATER
		// we only do thermal radiation upwards for the surface layer of water,
		// for the bottom we rely on conduction with the terrain for heat transfer (although this might lead to an imbalance!)
		{
			int layer = worldData.SurfaceWaterLayer + worldData.WaterLayer0;
			thermalOutJobHandles[layer] = SimJob.Schedule(new ThermalEnergyRadiatedWaterJob()
			{
				ThermalRadiationDelta = tempState.ThermalRadiationDelta[layer],
				ThermalRadiationTransmittedUp = tempState.ThermalRadiationTransmittedUp[layer],
				WindowRadiationTransmittedUp = tempState.WindowRadiationTransmittedUp[layer],

				PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
				Emissivity = tempState.Emissivity[layer],
				Energy = tempState.WaterPotentialEnergy[worldData.SurfaceWaterLayer],
				TemperatureAbsolute = lastState.WaterTemperature[worldData.SurfaceWaterLayer],
				SurfaceArea = tempState.WaterCoverage[worldData.SurfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			}, emissivityJobHandles[layer]);
		}
		#endregion


		#region absorptivity

		var cloudAlbedoJobHandle = SimJob.Schedule(new CloudAlbedoJob()
		{
			CloudAlbedo = tempState.CloudAlbedo,
			CloudAbsorptivity = tempState.CloudAbsorptivity,
			CloudMass = lastState.CloudMass,
			CloudDropletMass = lastState.CloudDropletMass,
			DewPoint = tempState.DewPoint,
			CloudElevation = tempState.CloudElevation,
			AlbedoSlope = tempState.AlbedoSlope,
			CloudFreezingTemperatureMax = worldData.maxCloudFreezingTemperature,
			CloudFreezingTemperatureMin = worldData.minCloudFreezingTemperature,
			RainDropSizeAlbedoMax = worldData.rainDropSizeAlbedoMax,
			RainDropSizeAlbedoMin = worldData.rainDropSizeAlbedoMin,
			SolarAbsorptivityCloud = worldData.SolarAbsorptivityCloud,
		}, solarInJobHandle);

		var absorptivityAirJobHandles = new JobHandle[worldData.AirLayers];
		for (int j = worldData.AirLayers - 2; j > 0; j--)
		{
			absorptivityAirJobHandles[j] = SimJob.Schedule(new AbsorptivityAirJob()
			{
				AbsorptivitySolar = tempState.AbsorptivitySolar[j],
				AbsorptivityThermal = tempState.AbsorptivityThermal[j],
				AirMass = tempState.AirMass[j],
				VaporMass = lastState.AirVapor[j],
				AirCarbonDioxide = lastState.AirCarbon[j],
				Dust = lastState.Dust[j],
				CloudMass = lastState.CloudMass,
				CloudAlbedo = tempState.CloudAlbedo,
				CloudAbsorptivity = tempState.CloudAbsorptivity,
				CloudElevation = tempState.CloudElevation,
				LayerElevation = tempState.LayerElevation[j],
				LayerHeight = tempState.LayerHeight[j],
				AlbedoAir = worldData.AlbedoAir,
				AlbedoWaterVapor = worldData.AlbedoWaterVapor,
				AlbedoDust = worldData.AlbedoDust,
				SolarAbsorptivityAir = worldData.SolarAbsorptivityAir,
				SolarAbsorptivityWaterVapor = worldData.SolarAbsorptivityWaterVapor,
				SolarAbsorptivityDust = worldData.SolarAbsorptivityDust,
				ThermalAbsorptivityAir = worldData.ThermalAbsorptivityAir,
				ThermalAbsorptivityWaterVapor = worldData.ThermalAbsorptivityWaterVapor,
				ThermalAbsorptivityOxygen = worldData.ThermalAbsorptivityOxygen,
				ThermalAbsorptivityCarbonDioxide = worldData.ThermalAbsorptivityCarbonDioxide,
				ThermalAbsorptivityDust = worldData.ThermalAbsorptivityDust,
				ThermalAbsorptivityCloud = worldData.ThermalAbsorptivityCloud,
			}, cloudAlbedoJobHandle);
		}

		#endregion


		// Follow the solar radiation down from the top of the atmosphere to ther terrain, and absorb some as it passes through each layer
		#region Solar Radiation Absorbed

		// process each vertical layer in order

		// atmosphere
		JobHandle[] solarInJobHandles = new JobHandle[worldData.LayerCount];
		for (int j = worldData.AirLayer0 + worldData.AirLayers - 2; j > worldData.AirLayer0; j--)
		{
			int airLayerIndex = j - worldData.AirLayer0;
			solarInJobHandles[j] = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedAirJob()
			{
				SolarRadiationAbsorbed = tempState.SolarRadiationIn[j],
				SolarRadiationIncoming = tempState.SolarRadiation,
				SolarRadiationReflected = tempState.SolarReflected[j],
				SolarRadiationAbsorbedCloud = tempState.SolarRadiationIn[worldData.CloudLayer],
				SolarRadiationReflectedCloud = tempState.SolarReflected[worldData.CloudLayer],
				AbsorptivitySolar = tempState.AbsorptivitySolar[airLayerIndex],
			}, JobHandle.CombineDependencies(solarInJobHandle, absorptivityAirJobHandles[airLayerIndex]));
		}


		// ice
		solarInJobHandles[worldData.IceLayer] = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedPartialCoverageConstantAlbedoJob()
		{
			SolarRadiationAbsorbed = tempState.SolarRadiationIn[worldData.IceLayer],
			SolarRadiationIncoming = tempState.SolarRadiation,
			SolarRadiationReflected = tempState.SolarReflected[worldData.IceLayer],
			AlbedoSlope = tempState.AlbedoSlope,
			AlbedoMin = WorldData.AlbedoIce,
			AlbedoRange = 1.0f - WorldData.AlbedoIce,
			Coverage = tempState.IceCoverage
		}, solarInJobHandle);

		for (int j = worldData.SurfaceWaterLayer; j >= 1; j--)
		{
			int layerIndex = worldData.WaterLayer0 + j;
			solarInJobHandles[layerIndex] = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedSlopeJob()
			{
				SolarRadiationAbsorbed = tempState.SolarRadiationIn[layerIndex],
				SolarRadiationIncoming = tempState.SolarRadiation,
				SolarRadiationReflected = tempState.SolarReflected[layerIndex],
				Coverage = tempState.WaterCoverage[layerIndex - worldData.WaterLayer0],
				AlbedoSlope = tempState.AlbedoSlope,
				AlbedoMin = WorldData.AlbedoWater,
			}, solarInJobHandle);
		}

		// flora
		solarInJobHandles[worldData.FloraLayer] = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedPartialCoverageConstantAlbedoJob()
		{
			SolarRadiationAbsorbed = tempState.SolarRadiationIn[worldData.FloraLayer],
			SolarRadiationIncoming = tempState.SolarRadiation,
			SolarRadiationReflected = tempState.SolarReflected[worldData.FloraLayer],
			AlbedoSlope = tempState.AlbedoSlope,
			AlbedoMin = WorldData.AlbedoFloraMin,
			AlbedoRange = WorldData.AlbedoFloraRange,
			Coverage = tempState.FloraCoverage
		}, solarInJobHandle);

		// NOTE: we don't bother with solar radiation in lava

		solarInJobHandles[worldData.TerrainLayer] = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedTerrainJob()
		{
			SolarRadiationAbsorbed = tempState.SolarRadiationIn[worldData.TerrainLayer],
			SolarRadiationIncoming = tempState.SolarRadiation,
			SolarRadiationReflected = tempState.SolarReflected[worldData.TerrainLayer],
			worldData = worldData,
			SoilFertility = lastState.GroundCarbon,
		}, solarInJobHandle);
		#endregion


		// Thermal radiation travels upwards, partially reflecting downwards (clouds), partially absorbed, and partially lost to space
		#region Thermal Radiation Absorbed Up
		// Start at bottom water layer and go up, then go back down
		JobHandle[] thermalInUpJobHandles = new JobHandle[worldData.LayerCount];
		JobHandle[] thermalInDownJobHandles = new JobHandle[worldData.LayerCount];

		// transmit up from land
		for (int j = 1; j < worldData.LayerCount; j++)
		{

			if (j > worldData.AirLayer0 && j < worldData.AirLayer0 + worldData.AirLayers - 1)
			{
				int airLayer = j - worldData.AirLayer0;
				int downIndex = airLayer == 1 ? worldData.IceLayer : (j - 1);
				thermalInUpJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedAirJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[j],
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp[j],
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp[j],

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedUp[downIndex],
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedUp[downIndex],
					AbsorptivityThermal = tempState.AbsorptivityThermal[airLayer],
					LayerElevation = tempState.LayerElevation[airLayer],
					LayerHeight = tempState.LayerHeight[airLayer],
					CloudElevation = tempState.CloudElevation,
					FromTop = false,
				}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex], absorptivityAirJobHandles[airLayer]));
			}
			else if (j == worldData.IceLayer)
			{
				int downIndex = worldData.SurfaceWaterLayerGlobal;
				thermalInUpJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[j],
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp[j],
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp[j],

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedUp[downIndex],
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedUp[downIndex],
					Coverage = tempState.IceCoverage,

				}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex]));
			}
			else if (j == worldData.FloraLayer)
			{
				int downIndex = worldData.TerrainLayer;
				thermalInUpJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[j],
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp[j],
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp[j],

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedUp[downIndex],
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedUp[downIndex],
					Coverage = tempState.FloraCoverage,

				}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex], thermalOutJobHandles[downIndex]));
			}
			else if (j > worldData.WaterLayer0 && j < worldData.WaterLayer0 + worldData.WaterLayers - 1)
			{
				int waterLayerIndex = j - worldData.WaterLayer0;
				int downIndex = (waterLayerIndex == 1) ? worldData.FloraLayer : (j - 1);
				thermalInUpJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[j],
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp[j],
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp[j],

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedUp[downIndex],
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedUp[downIndex],
					Coverage = tempState.WaterCoverage[waterLayerIndex],
				}, JobHandle.CombineDependencies(thermalOutJobHandles[j], thermalInUpJobHandles[downIndex]));
			}
		}

		var thermalInUpJobHandlesCombined = lastJobHandle;
		for (int j = 0; j < worldData.LayerCount; j++)
		{
			thermalInUpJobHandlesCombined = JobHandle.CombineDependencies(thermalInUpJobHandlesCombined, thermalInUpJobHandles[j]);
		}
		#endregion

		// Thermal radiation is absorbed travelling downwards, collecting and then eventually hitting the earth (back radiation)
		// TODO: we need to include the top layer of atmosphere here, since we calculate cloud absorption as part of the air layer step
		#region Thermal Radiation Absorbed Down

		// transmit down from top of atmosphere			
		for (int j = worldData.LayerCount - 1; j >= 0; j--)
		{

			if (j == worldData.TerrainLayer)
			{
				// TERRAIN
				int upIndex = worldData.FloraLayer;
				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
				thermalInDownJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedTerrainJob()
				{
					ThermalRadiationAbsorbed = tempState.ThermalRadiationDelta[j],

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedDown[upIndex],
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedDown[upIndex],
				}, thermalInDependenciesHandle);
			}
			else if (j == worldData.FloraLayer)
			{
				// FLORA
				int upIndex = worldData.SurfaceWaterLayerGlobal;
				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
				thermalInDownJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[j],
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedDown[j],
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedDown[j],

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedDown[upIndex],
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedDown[upIndex],
					Coverage = tempState.FloraCoverage,
				}, thermalInDependenciesHandle);
			}
			else if (j == worldData.IceLayer)
			{
				// ICE
				int upIndex = worldData.SurfaceAirLayerGlobal;
				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
				thermalInDownJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[j],
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedDown[j],
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedDown[j],

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedDown[upIndex],
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedDown[upIndex],
					Coverage = tempState.IceCoverage,
				}, thermalInDependenciesHandle);
			}
			else if (j == worldData.SurfaceWaterLayerGlobal)
			{
				// WATER
				int waterLayerIndex = j - worldData.WaterLayer0;
				int upIndex = (waterLayerIndex == worldData.SurfaceWaterLayer) ? worldData.IceLayer : (j + 1);
				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
				thermalInDownJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedPartialCoverageJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[j],
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedDown[j],
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedDown[j],

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedDown[upIndex],
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedDown[upIndex],
					Coverage = tempState.WaterCoverage[waterLayerIndex],
				}, thermalInDependenciesHandle);
			}
			else if (j > worldData.AirLayer0 && j < worldData.AirLayer0 + worldData.AirLayers - 1)
			{
				int airLayer = j - worldData.AirLayer0;
				int upIndex = j + 1;
				var thermalInDependenciesHandle = JobHandle.CombineDependencies(thermalInDownJobHandles[upIndex], thermalInUpJobHandlesCombined);
				thermalInDownJobHandles[j] = SimJob.Schedule(new ThermalEnergyAbsorbedAirJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[j],
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedDown[j],
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedDown[j],

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedDown[upIndex],
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedDown[upIndex],
					AbsorptivityThermal = tempState.AbsorptivityThermal[airLayer],
					LayerElevation = tempState.LayerElevation[airLayer],
					LayerHeight = tempState.LayerHeight[airLayer],
					CloudElevation = tempState.CloudElevation,
					FromTop = true,
				}, thermalInDependenciesHandle);
			}
		}
		#endregion

		// Conduction is calculated for each Surface that might touch another surface
		// Air to Cloud, Air to Ice, Air to Water, Air to Terrain, Ice to Water, Ice to Terrain, Water to Terrain
		#region Conduction

		JobHandle conductionAirIceJobHandle;
		JobHandle conductionAirWaterJobHandle;
		JobHandle conductionAirFloraJobHandle;
		JobHandle conductionAirTerrainJobHandle;
		JobHandle conductionIceWaterJobHandle;
		JobHandle conductionIceFloraJobHandle;
		JobHandle conductionIceTerrainJobHandle;
		JobHandle conductionFloraTerrainJobHandle;
		JobHandle conductionWaterTerrainJobHandle;

		// air to ice
		conductionAirIceJobHandle = SimJob.ScheduleOrMemset(
			settings.ConductionAirIce,
			tempState.ConductionAirIce,
			0,
			new ConductionBJob()
			{
				EnergyDelta = tempState.ConductionAirIce,
				TemperatureA = tempState.SurfaceAirTemperatureAbsolute,
				TemperatureB = lastState.IceTemperature,
				EnergyB = tempState.IceEnergy,
				ConductionCoefficient = WorldData.ConductivityAirIce,
				SurfaceArea = tempState.SurfaceAreaAirIce,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);

		// air to water
		conductionAirWaterJobHandle = SimJob.ScheduleOrMemset(
			settings.ConductionAirWater,
			tempState.ConductionAirWater,
			0,
			new ConductionBJob()
			{
				EnergyDelta = tempState.ConductionAirWater,
				TemperatureA = tempState.SurfaceAirTemperatureAbsolute,
				TemperatureB = lastState.WaterTemperature[worldData.SurfaceWaterLayer],
				EnergyB = tempState.WaterPotentialEnergy[worldData.SurfaceWaterLayer],
				ConductionCoefficient = WorldData.ConductivityAirWater,
				SurfaceArea = tempState.SurfaceAreaAirWater,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);

		// air to flora
		conductionAirFloraJobHandle = SimJob.ScheduleOrMemset(
			settings.ConductionAirFlora,
			tempState.ConductionAirFlora,
			0,
			new ConductionBJob()
			{
				EnergyDelta = tempState.ConductionAirFlora,
				TemperatureA = tempState.SurfaceAirTemperatureAbsolute,
				TemperatureB = lastState.FloraTemperature,
				ConductionCoefficient = WorldData.ConductivityAirFlora,
				SurfaceArea = tempState.SurfaceAreaAirFlora,
				EnergyB = tempState.FloraEnergy,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);

		// air to terrain
		conductionAirTerrainJobHandle = SimJob.ScheduleOrMemset(
			settings.ConductionAirTerrain,
			tempState.ConductionAirTerrain,
			0,
			new ConductionJob()
			{
				EnergyDelta = tempState.ConductionAirTerrain,
				TemperatureA = tempState.SurfaceAirTemperatureAbsolute,
				TemperatureB = lastState.GroundTemperature,
				ConductionCoefficient = WorldData.ConductivityAirTerrain,
				SurfaceArea = tempState.SurfaceAreaAirTerrain,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);

		// ice to water
		conductionIceWaterJobHandle = SimJob.ScheduleOrMemset(
			settings.ConductionIceWater,
			tempState.ConductionIceWater,
			0,
			new ConductionABJob()
			{
				EnergyDelta = tempState.ConductionIceWater,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.WaterTemperature[worldData.SurfaceWaterLayer],
				EnergyA = tempState.IceEnergy,
				EnergyB = tempState.WaterPotentialEnergy[worldData.SurfaceWaterLayer],
				ConductionCoefficient = WorldData.ConductivityIceWater,
				SurfaceArea = tempState.SurfaceAreaIceWater,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);

		// ice to flora
		conductionIceFloraJobHandle = SimJob.ScheduleOrMemset(
			settings.ConductionIceFlora,
			tempState.ConductionIceFlora,
			0,
			new ConductionABJob()
			{
				EnergyDelta = tempState.ConductionIceFlora,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.FloraTemperature,
				EnergyA = tempState.IceEnergy,
				EnergyB = tempState.FloraEnergy,
				ConductionCoefficient = WorldData.ConductivityIceFlora,
				SurfaceArea = tempState.SurfaceAreaIceFlora,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);

		// ice to terrain
		conductionIceTerrainJobHandle = SimJob.ScheduleOrMemset(
			settings.ConductionIceTerrain,
			tempState.ConductionIceTerrain,
			0,
			new ConductionAJob()
			{
				EnergyDelta = tempState.ConductionIceTerrain,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = lastState.GroundTemperature,
				EnergyA = tempState.IceEnergy,
				ConductionCoefficient = WorldData.ConductivityIceTerrain,
				SurfaceArea = tempState.SurfaceAreaIceTerrain,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);

		// flora to terrain
		conductionFloraTerrainJobHandle = SimJob.ScheduleOrMemset(
			settings.ConductionFloraTerrain,
			tempState.ConductionFloraTerrain,
			0,
			new ConductionAJob()
			{
				EnergyDelta = tempState.ConductionFloraTerrain,
				TemperatureA = lastState.FloraTemperature,
				TemperatureB = lastState.GroundTemperature,
				EnergyA = tempState.FloraEnergy,
				ConductionCoefficient = WorldData.ConductivityFloraTerrain,
				SurfaceArea = tempState.SurfaceAreaFloraTerrain,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);


		// water to terrain
		conductionWaterTerrainJobHandle = lastJobHandle;
		for (int i = 1; i < worldData.WaterLayers - 1; i++)
		{
			conductionWaterTerrainJobHandle = JobHandle.CombineDependencies(conductionWaterTerrainJobHandle, 
				SimJob.ScheduleOrMemset(
					settings.ConductionWaterTerrain,
					tempState.ConductionWaterTerrain[i],
					0,
					new ConductionWaterBottomAJob()
					{
						EnergyDelta = tempState.ConductionWaterTerrain[i],
						EnergyDeltaTotal = tempState.ConductionWaterTerrainTotal,
						TemperatureA = lastState.WaterTemperature[i],
						TemperatureB = lastState.GroundTemperature,
						EnergyA = tempState.WaterPotentialEnergy[i],
						ConductionCoefficient = WorldData.ConductivityWaterTerrain,
						SurfaceArea = tempState.SurfaceAreaWaterTerrain,
						Coverage = tempState.WaterCoverage[i],
						CoverageBelow = tempState.WaterCoverage[i - 1],
						SecondsPerTick = worldData.SecondsPerTick
					}, 
					conductionWaterTerrainJobHandle));
		}

		#endregion

		#region Change temperature due to energy flux

		var terrainEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.Persistent)
			{
				solarInJobHandles[worldData.TerrainLayer],
				thermalOutJobHandles[worldData.TerrainLayer],
				thermalInDownJobHandles[worldData.TerrainLayer],
				thermalInUpJobHandles[worldData.TerrainLayer],
				conductionAirTerrainJobHandle,
				conductionIceTerrainJobHandle,
				conductionWaterTerrainJobHandle,
				conductionFloraTerrainJobHandle,
			};
		jobHandleDependencies.Add(terrainEnergyJobHandleDependencies);
		energyJobHandles[worldData.TerrainLayer] = SimJob.Schedule(new EnergyTerrainJob()
		{
			TerrainTemperature = nextState.GroundTemperature,
			LastTemperature = lastState.GroundTemperature,
			SoilFertility = lastState.GroundCarbon,
			SolarRadiationIn = tempState.SolarRadiationIn[worldData.TerrainLayer],
			ThermalRadiationDelta = tempState.ThermalRadiationDelta[worldData.TerrainLayer],
			ConductionEnergyAir = tempState.ConductionAirTerrain,
			ConductionEnergyIce = tempState.ConductionIceTerrain,
			ConductionEnergyFlora = tempState.ConductionFloraTerrain,
			ConductionEnergyWater = tempState.ConductionWaterTerrainTotal,
			GeothermalEnergy = tempState.GeothermalRadiation,
			HeatingDepth = worldData.SoilHeatDepth,
		}, JobHandle.CombineDependencies(terrainEnergyJobHandleDependencies));

		var energyIceJobHandleDependencies = new NativeList<JobHandle>(Allocator.Persistent)
			{
				solarInJobHandles[worldData.IceLayer],
				thermalOutJobHandles[worldData.IceLayer],
				thermalInDownJobHandles[worldData.IceLayer],
				thermalInUpJobHandles[worldData.IceLayer],
				conductionAirIceJobHandle,
				conductionIceWaterJobHandle,
				conductionIceFloraJobHandle,
				conductionIceTerrainJobHandle,
			};
		jobHandleDependencies.Add(energyIceJobHandleDependencies);
		energyJobHandles[worldData.IceLayer] = SimJob.Schedule(new EnergyIceJob()
		{
			Temperature = nextState.IceTemperature,
			LastTemperature = lastState.IceTemperature,
			LastMass = lastState.IceMass,
			SolarRadiationIn = tempState.SolarRadiationIn[worldData.IceLayer],
			ThermalRadiationDelta = tempState.ThermalRadiationDelta[worldData.IceLayer],
			ConductionEnergyAir = tempState.ConductionAirIce,
			ConductionEnergyTerrain = tempState.ConductionIceTerrain,
			ConductionEnergyWater = tempState.ConductionIceWater,
			ConductionEnergyFlora = tempState.ConductionIceFlora,
		}, JobHandle.CombineDependencies(energyIceJobHandleDependencies));

		var energyFloraJobHandleDependencies = new NativeList<JobHandle>(Allocator.Persistent)
			{
				solarInJobHandles[worldData.FloraLayer],
				thermalOutJobHandles[worldData.FloraLayer],
				thermalInDownJobHandles[worldData.FloraLayer],
				thermalInUpJobHandles[worldData.FloraLayer],
				conductionAirFloraJobHandle,
				conductionIceFloraJobHandle,
				conductionFloraTerrainJobHandle,
			};
		jobHandleDependencies.Add(energyFloraJobHandleDependencies);
		energyJobHandles[worldData.FloraLayer] = SimJob.Schedule(new EnergyFloraJob()
		{
			FloraTemperature = nextState.FloraTemperature,

			LastTemperature = lastState.FloraTemperature,
			FloraMass = lastState.FloraMass,
			FloraWater = lastState.FloraWater,
			ThermalRadiationDelta = tempState.ThermalRadiationDelta[worldData.FloraLayer],
			ConductionEnergyAir = tempState.ConductionAirFlora,
			ConductionEnergyTerrain = tempState.ConductionFloraTerrain,
			ConductionEnergyIce = tempState.ConductionIceFlora,
		}, JobHandle.CombineDependencies(energyFloraJobHandleDependencies));

		energyJobHandles[worldData.LavaLayer] = SimJob.Schedule(new EnergyLavaJob()
		{
			LavaTemperature = nextState.LavaTemperature,
			LastTemperature = lastState.LavaTemperature,
			LavaMass = lastState.LavaMass,
			Emissivity = worldData.ThermalEmissivityLava,
			SecondsPerTick = worldData.SecondsPerTick
		});

		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			int layerIndex = worldData.AirLayer0 + j;
			var airDependencies = new NativeList<JobHandle>(Allocator.Persistent)
				{
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
				};
			if (j == 1)
			{
				airDependencies.Add(conductionAirWaterJobHandle);
				airDependencies.Add(conductionAirIceJobHandle);
				airDependencies.Add(conductionAirFloraJobHandle);
				airDependencies.Add(conductionAirTerrainJobHandle);
				energyJobHandles[layerIndex] = SimJob.Schedule(new EnergyAirSurfaceJob()
				{
					AirTemperaturePotential = nextState.AirTemperaturePotential[j],
					LastTemperaturePotential = lastState.AirTemperaturePotential[j],
					LastVapor = lastState.AirVapor[j],
					AirMass = tempState.AirMass[j],
					ConductionEnergyWater = tempState.ConductionAirWater,
					ConductionEnergyIce = tempState.ConductionAirIce,
					ConductionEnergyFlora = tempState.ConductionAirFlora,
					ConductionEnergyTerrain = tempState.ConductionAirTerrain,
					SolarRadiationIn = tempState.SolarRadiationIn[layerIndex],
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[layerIndex],
				}, JobHandle.CombineDependencies(airDependencies));
			}
			else
			{
				energyJobHandles[layerIndex] = SimJob.Schedule(new EnergyAirJob()
				{
					AirTemperaturePotential = nextState.AirTemperaturePotential[j],
					LastTemperaturePotential = lastState.AirTemperaturePotential[j],
					LastVapor = lastState.AirVapor[j],
					AirMass = tempState.AirMass[j],
					SolarRadiationIn = tempState.SolarRadiationIn[layerIndex],
					ThermalRadiationDelta = tempState.ThermalRadiationDelta[layerIndex],
				}, JobHandle.CombineDependencies(airDependencies));

			}
			jobHandleDependencies.Add(airDependencies);

		}

		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			int layerIndex = worldData.WaterLayer0 + j;
			var waterDependencies = new NativeList<JobHandle>(Allocator.Persistent)
				{
					solarInJobHandles[layerIndex],
					thermalOutJobHandles[layerIndex],
					thermalInDownJobHandles[layerIndex],
					thermalInUpJobHandles[layerIndex],
					conductionAirWaterJobHandle,
					conductionIceWaterJobHandle,
					conductionWaterTerrainJobHandle,
				};
			jobHandleDependencies.Add(waterDependencies);

			energyJobHandles[layerIndex] = SimJob.Schedule(new EnergyWaterJob()
			{
				Temperature = nextState.WaterTemperature[j],
				LastMass = lastState.WaterMass[j],
				LastSaltMass = lastState.SaltMass[j],
				LastTemperature = lastState.WaterTemperature[j],
				ThermalRadiationDelta = tempState.ThermalRadiationDelta[layerIndex],
				CoverageUp = tempState.WaterCoverage[j + 1],
				CoverageDown = tempState.WaterCoverage[j - 1],
				ConductionEnergyAir = tempState.ConductionAirWater,
				ConductionEnergyIce = tempState.ConductionIceWater,
				ConductionEnergyTerrain = tempState.ConductionWaterTerrain[j],
			}, JobHandle.CombineDependencies(waterDependencies));
		}
		lastJobHandle = JobHandle.CombineDependencies(energyJobHandles);

		#endregion


	}


	private void DoStateChange(
		NativeArray<JobHandle> energyJobHandles,
		JobHandle lastJobHandle,
		ref SimState nextState,
		ref SimState lastState,
		ref TempState tempState,
		ref StaticState staticState,
		ref WorldData worldData,
		ref SimSettings settings
		)
	{
		#region State Change (Flux)

		// surface water

		if (settings.Condensation)
		{
			for (int j = 1; j < worldData.AirLayers - 1; j++)
			{
				int layerIndex = j + worldData.AirLayer0;
				energyJobHandles[layerIndex] = SimJob.Schedule(new FluxAirCondensationJob()
				{
					LatentHeat = tempState.LatentHeat[layerIndex],
					CondensationCloudMass = tempState.CondensationCloudMass[j],
					CondensationGroundMass = tempState.CondensationGroundMass[j],

					TemperaturePotential = nextState.AirTemperaturePotential[j],
					LastVapor = lastState.AirVapor[j],
					AirMass = tempState.AirMass[j],
					AirPressure = tempState.AirPressure[j],
					CloudElevation = tempState.CloudElevation,
					LayerElevation = tempState.LayerElevation[j],
					LayerHeight = tempState.LayerHeight[j],
					LayerMiddle = tempState.LayerMiddle[j],
				}, energyJobHandles[layerIndex]);
			}
		}
		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			int layerIndex = j + worldData.AirLayer0;
			energyJobHandles[layerIndex] = SimJob.Schedule(new FluxAirDustJob()
			{
				DustUp = tempState.DustUp[j],
				DustDown = tempState.DustDown[j],

				LayerHeight = tempState.LayerHeight[j],
				LastDust = lastState.Dust[j],
				AirVelocity = lastState.AirVelocity[j],
				Positions = staticState.SphericalPosition,
				DustVerticalVelocity = worldData.DustVerticalVelocity,
				SecondsPerTick = worldData.SecondsPerTick
			}, energyJobHandles[layerIndex]);
		}

		if (settings.Evaporation)
		{
			energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new FluxWaterEvaporationJob()
			{
				EvaporatedWaterMass = tempState.EvaporationMassWater,
				LatentHeatWater = tempState.LatentHeat[worldData.SurfaceWaterLayerGlobal],
				LatentHeatAir = tempState.LatentHeat[worldData.SurfaceAirLayerGlobal],

				WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
				IceCoverage = tempState.IceCoverage,
				WaterCoverage = tempState.WaterCoverage[worldData.SurfaceWaterLayer],
				SurfaceWind = lastState.AirVelocity[worldData.SurfaceAirLayer],
				AirMass = tempState.AirMass[worldData.SurfaceAirLayer],
				AirPressure = tempState.AirPressure[worldData.SurfaceAirLayer],
				AirVapor = lastState.AirVapor[worldData.SurfaceAirLayer],
				WaterHeatingDepth = worldData.WaterHeatingDepth,
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.SurfaceAirLayerGlobal]));
		}

		if (settings.Freezing)
		{
			energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new FluxWaterFreezeJob()
			{
				FrozenMass = tempState.FrozenMass,
				FrozenTemperature = tempState.FrozenTemperature,
				LatentHeatWater = tempState.LatentHeat[worldData.SurfaceWaterLayerGlobal],
				SaltPlume = tempState.SaltPlume,

				WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				AirTemperaturePotential = nextState.AirTemperaturePotential[worldData.SurfaceAirLayer],
				WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
				SaltMass = lastState.SaltMass[worldData.SurfaceWaterLayer],
				AirLayerElevation = tempState.LayerElevation[worldData.SurfaceAirLayer],
				WaterHeatingDepth = worldData.WaterHeatingDepth,
				FreezePointReductionPerSalinity = worldData.FreezePointReductionPerSalinity,
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.SurfaceAirLayerGlobal]));
		}

		energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new FluxFloraWaterConsumeJob()
		{
			FloraWaterConsumed = tempState.WaterConsumedByFlora,

			WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
			FloraMass = lastState.FloraMass,
			FloraWater = lastState.FloraWater,
			FloraWaterConsumptionRate = worldData.FloraWaterConsumptionRate
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.SurfaceAirLayerGlobal]));

		if (settings.Plankton)
		{
			energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new FluxPlanktonJob()
			{
				LatentHeatWater = tempState.LatentHeat[worldData.SurfaceWaterLayerGlobal],
				PlanktonMassDelta = tempState.PlanktonMassDelta,
				PlanktonGlucoseDelta = tempState.PlanktonGlucoseDelta,
				PlanktonDeath = tempState.PlanktonDeath,
				WaterCarbonDelta = tempState.WaterCarbonDelta,

				WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
				SaltMass = lastState.SaltMass[worldData.SurfaceWaterLayer],
				SolarRadiation = tempState.SolarRadiationIn[worldData.SurfaceWaterLayerGlobal],
				WaterCarbon = lastState.WaterCarbon[worldData.SurfaceWaterLayer],
				PlanktonMass = lastState.PlanktonMass[worldData.SurfaceWaterLayer],
				PlanktonGlucoseMass = lastState.PlanktonGlucose[worldData.SurfaceWaterLayer],
				PlanktonDensityMax = worldData.PlanktonDensityMax,
				PlanktonEnergyForPhotosynthesis = worldData.PlanktonEnergyForPhotosynthesis,
				PlanktonCarbonDioxideExtractionEfficiency = worldData.PlanktonCarbonDioxideExtractionEfficiency,
				PlanktonPhotosynthesisSpeed = worldData.PlanktonPhotosynthesisSpeed,
				PlanktonRespirationSpeed = worldData.PlanktonRespirationSpeed,
				PlanktonRespirationPerDegree = worldData.PlanktonRespirationPerDegree,
				PlanktonGrowthRate = worldData.PlanktonGrowthRate,
				PlanktonDeathRate = worldData.PlanktonDeathRate,
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.SurfaceAirLayerGlobal]));
		}

		// CLOUD
		if (settings.Precipitation)
		{
			energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(new FluxCloudJob()
			{
				EvaporationMass = tempState.CloudEvaporationMass,
				PrecipitationMass = tempState.PrecipitationMass,
				PrecipitationTemperature = tempState.PrecipitationTemperature,
				DropletDelta = tempState.DropletDelta,

				SurfaceAirTemperaturePotential = nextState.AirTemperaturePotential[worldData.SurfaceAirLayer],
				SurfaceLayerElevation = tempState.LayerElevation[worldData.SurfaceAirLayer],
				SurfaceLayerMiddle = tempState.LayerMiddle[worldData.SurfaceAirLayer],
				SurfaceSaltMass = lastState.SaltMass[worldData.SurfaceWaterLayer],
				LastCloudMass = lastState.CloudMass,
				LastVelocity = tempState.CloudVelocity,
				LastDropletMass = lastState.CloudDropletMass,
				CloudElevation = tempState.CloudElevation,
				DewPoint = tempState.DewPoint,
				AirDensityCloud = tempState.AirDensityCloud,
				Position = staticState.SphericalPosition,
				Gravity = lastState.PlanetState.Gravity,
				RainDropDragCoefficient = worldData.rainDropDragCoefficient,
				RainDropMaxSize = worldData.rainDropMaxSize,
				RainDropMinSize = worldData.rainDropMinSize,
				RainDropGrowthRate = worldData.RainDropGrowthRate,
				SecondsPerTick = worldData.SecondsPerTick,
				CloudDissapationRateDryAir = worldData.CloudDissapationRateDryAir,
				CloudDissapationRateWind = worldData.CloudDissapationRateWind,
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.CloudLayer], energyJobHandles[worldData.SurfaceAirLayerGlobal]));
		}

		if (settings.Flora)
		{
			energyJobHandles[worldData.FloraLayer] = SimJob.Schedule(new FluxFloraJob()
			{
				LatentHeatAir = tempState.LatentHeat[worldData.SurfaceAirLayerGlobal],
				LatentHeatFlora = tempState.LatentHeat[worldData.FloraLayer],
				EvaporatedWaterMass = tempState.FloraRespirationMassVapor,
				SurfaceWaterDelta = tempState.FloraRespirationMassWater,
				FloraMassDelta = tempState.FloraMassDelta,
				FloraWaterDelta = tempState.FloraWaterDelta,
				FloraGlucoseDelta = tempState.FloraGlucoseDelta,
				FloraDeath = tempState.FloraDeath,
				CarbonDioxideDelta = tempState.AirCarbonDelta,
				OxygenDelta = tempState.OxygenDelta,

				SolarRadiationIn = tempState.SolarRadiationIn[worldData.FloraLayer],
				FloraTemperature = nextState.FloraTemperature,
				FloraMass = lastState.FloraMass,
				FloraGlucose = lastState.FloraGlucose,
				FloraWater = lastState.FloraWater,
				FloraCoverage = tempState.FloraCoverage,
				CarbonDioxide = lastState.AirCarbon[worldData.SurfaceAirLayer],
				LayerElevation = tempState.LayerElevation[worldData.SurfaceAirLayer],
				LayerHeight = tempState.LayerHeight[worldData.SurfaceAirLayer],
				SurfaceWind = lastState.AirVelocity[worldData.SurfaceAirLayer],
				AirMass = tempState.AirMass[worldData.SurfaceAirLayer],
				AirTemperaturePotential = lastState.AirTemperaturePotential[worldData.SurfaceAirLayer],
				AirPressure = tempState.AirPressure[worldData.SurfaceAirLayer],
				AirVapor = lastState.AirVapor[worldData.SurfaceAirLayer],
				SoilFertility = lastState.GroundCarbon,
				FloraGrowthRate = worldData.FloraGrowthRate,
				FloraDeathRate = worldData.FloraDeathRate,
				FloraGrowthTemperatureRangeInverse = worldData.FloraGrowthTemperatureRangeInverse,
				FloraEnergyForPhotosynthesis = worldData.FloraEnergyForPhotosynthesis,
				FloraCarbonDioxideExtractionEfficiency = worldData.FloraCarbonDioxideExtractionEfficiency,
				FloraOxygenExtractionEfficiency = worldData.FloraOxygenExtractionEfficiency,
				FloraPhotosynthesisSpeed = worldData.FloraPhotosynthesisSpeed,
				FloraRespirationSpeed = worldData.FloraRespirationSpeed,
				FloraRespirationPerDegree = worldData.FloraRespirationPerDegree,
				OxygenPercent = lastState.PlanetState.Oxygen,
				Gravity = lastState.PlanetState.Gravity
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.FloraLayer], energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.SurfaceAirLayerGlobal]));
		}

		if (settings.IceMelting)
		{
			energyJobHandles[worldData.IceLayer] = SimJob.Schedule(new FluxIceMeltJob()
			{
				LatentHeatAir = tempState.LatentHeat[worldData.SurfaceAirLayerGlobal],
				LatentHeatWater = tempState.LatentHeat[worldData.SurfaceWaterLayerGlobal],
				LatentHeatTerrain = tempState.LatentHeat[worldData.TerrainLayer],
				LatentHeatIce = tempState.LatentHeat[worldData.IceLayer],
				MeltedMass = tempState.IceMeltedMass,

				Temperature = nextState.IceTemperature,
				LastMass = lastState.IceMass,
				IceHeatingDepth = worldData.IceHeatingDepth,
				AirTemperaturePotential = nextState.AirTemperaturePotential[worldData.SurfaceAirLayer],
				WaterIceSurfaceArea = tempState.SurfaceAreaIceWater,
				WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				TerrainTemperature = nextState.GroundTemperature,
				LayerElevation = tempState.LayerElevation[worldData.SurfaceAirLayer],

			}, JobHandle.CombineDependencies(energyJobHandles[worldData.FloraLayer], energyJobHandles[worldData.IceLayer], energyJobHandles[worldData.TerrainLayer]));
		}

		energyJobHandles[worldData.LavaLayer] = SimJob.Schedule(new FluxLavaJob()
		{
		}, energyJobHandles[worldData.LavaLayer]);

		if (settings.SoilRespiration)
		{
			energyJobHandles[worldData.TerrainLayer] = SimJob.Schedule(new FluxTerrainJob()
			{
				SoilRespiration = tempState.SoilRespiration,
				CrystalizedMass = tempState.LavaCrystalizedMass,
				LavaEjected = tempState.LavaEjected,
				DustEjected = tempState.DustEjected,
				CrustDelta = tempState.CrustDelta,
				LatentHeatLava = tempState.LatentHeat[worldData.LavaLayer],

				LavaTemperature = nextState.LavaTemperature,
				LavaMass = lastState.LavaMass,
				CrustDepth = lastState.CrustDepth,
				MagmaMass = lastState.MagmaMass,
				Elevation = lastState.Elevation,
				SoilCarbon = nextState.GroundCarbon,
				WaterCoverage = tempState.WaterCoverage[worldData.SurfaceWaterLayer],
				LavaCrystalizationTemperature = worldData.LavaCrystalizationTemperature,
				CrustEruptionDepth = worldData.CrustDepthForEruption,
				DustPerLavaEjected = worldData.DustPerLavaEjected,
				MagmaPressureCrustReductionSpeed = worldData.MagmaPressureCrustReductionSpeed,
				LavaEruptionSpeed = worldData.LavaEruptionSpeed,
				SecondsPerTick = worldData.SecondsPerTick,
				SoilRespirationSpeed = worldData.SoilRespirationSpeed,
				OxygenPercent = lastState.PlanetState.Oxygen
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.TerrainLayer], energyJobHandles[worldData.LavaLayer]));
		}

		#endregion

		#region Update Mass - Evaporation, Condensation, Melting, Rainfall

		JobHandle waterDependencies = JobHandle.CombineDependencies(energyJobHandles[worldData.IceLayer], energyJobHandles[worldData.TerrainLayer], energyJobHandles[worldData.SurfaceWaterLayerGlobal]);
		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			energyJobHandles[worldData.WaterLayer0 + j] = SimJob.Schedule(new UpdateMassWaterJob()
			{
				WaterMass = nextState.WaterMass[j],
				SaltMass = nextState.SaltMass[j],
				CarbonMass = nextState.WaterCarbon[j],
				WaterTemperature = nextState.WaterTemperature[j],
				SaltPlume = tempState.SaltPlume,
				SaltPlumeTemperature = tempState.FrozenTemperature,
				LastSaltMass = lastState.SaltMass[j],
				LastCarbonMass = lastState.WaterCarbon[j],
				LastWaterMass = lastState.WaterMass[j],
				DownLastWaterMass = lastState.WaterMass[j - 1],
				SoilRespiration = tempState.SoilRespiration,
				WaterCoverage = tempState.WaterCoverage[j],
				WaterCoverageBelow = tempState.WaterCoverage[j - 1],
			}, JobHandle.CombineDependencies( energyJobHandles[worldData.WaterLayer0 + j], waterDependencies));
		}

		var surfaceWaterJobHandle = energyJobHandles[worldData.SurfaceWaterLayerGlobal];
		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			int layerIndex = worldData.AirLayer0 + j;
			surfaceWaterJobHandle = JobHandle.CombineDependencies(surfaceWaterJobHandle, SimJob.Schedule(new UpdateMassCondensationGroundJob()
			{
				SurfaceWaterMass = nextState.WaterMass[worldData.SurfaceWaterLayer],
				SurfaceWaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],

				AirTemperaturePotential = nextState.AirTemperaturePotential[j],
				GroundCondensation = tempState.CondensationGroundMass[j],
				SurfaceSaltMass = lastState.SaltMass[j],
				LayerMiddle = tempState.LayerMiddle[j],
			}, JobHandle.CombineDependencies(surfaceWaterJobHandle, energyJobHandles[layerIndex])));
			energyJobHandles[layerIndex] = surfaceWaterJobHandle;
		}

		energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new UpdateMassWaterSurfaceJob()
		{
			WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
			WaterMass = nextState.WaterMass[worldData.SurfaceWaterLayer],
			SaltMass = nextState.SaltMass[worldData.SurfaceWaterLayer],
			PlanktonMass = nextState.PlanktonMass[worldData.SurfaceWaterLayer],
			PlanktonGlucose = nextState.PlanktonGlucose[worldData.SurfaceWaterLayer],
			CarbonMass = nextState.WaterCarbon[worldData.SurfaceWaterLayer],

			SaltPlume = tempState.SaltPlume,
			Evaporation = tempState.EvaporationMassWater,
			IceMelted = tempState.IceMeltedMass,
			Precipitation = tempState.PrecipitationMass,
			PrecipitationTemperature = tempState.PrecipitationTemperature,
			FloraRespirationWater = tempState.FloraRespirationMassWater,
			FloraTemperature = lastState.FloraTemperature,
			WaterFrozen = tempState.FrozenMass,
			LastPlanktonMass = lastState.PlanktonMass[worldData.SurfaceWaterLayer],
			LastPlanktonGlucose = lastState.PlanktonGlucose[worldData.SurfaceWaterLayer],
			PlanktonMassDelta = tempState.PlanktonMassDelta,
			PlanktonGlucoseDelta = tempState.PlanktonGlucoseDelta,
			WaterCarbonDelta = tempState.WaterCarbonDelta,
		}, JobHandle.CombineDependencies(surfaceWaterJobHandle, energyJobHandles[worldData.CloudLayer]));

		energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(new UpdateMassCloudJob()
		{
			CloudMass = nextState.CloudMass,
			CloudDropletMass = nextState.CloudDropletMass,
			LastCloudMass = lastState.CloudMass,
			LastDropletMass = lastState.CloudDropletMass,
			CloudEvaporation = tempState.CloudEvaporationMass,
			PrecipitationMass = tempState.PrecipitationMass,
			DropletDelta = tempState.DropletDelta,
		}, energyJobHandles[worldData.CloudLayer]);

		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			int layerIndex = worldData.AirLayer0 + j;
			energyJobHandles[layerIndex] = SimJob.Schedule(new UpdateMassAirJob()
			{
				VaporMass = nextState.AirVapor[j],
				DustMass = nextState.Dust[j],
				CarbonDioxideMass = nextState.AirCarbon[j],

				CloudCondensation = tempState.CondensationCloudMass[j],
				GroundCondensation = tempState.CondensationGroundMass[j],
				LastVaporMass = lastState.AirVapor[j],
				LastDustMass = lastState.Dust[j],
				LastCarbonDioxideMass = lastState.AirCarbon[j],
				DustUp = tempState.DustUp[j],
				DustDown = tempState.DustDown[j],
				DustFromAbove = tempState.DustDown[j + 1],
				DustFromBelow = tempState.DustUp[j - 1],
				IsTop = j == worldData.AirLayers - 2,
				IsBottom = j == 1,
			}, JobHandle.CombineDependencies(energyJobHandles[layerIndex], energyJobHandles[layerIndex-1], energyJobHandles[layerIndex + 1]));
		}

		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			int layerIndex = worldData.AirLayer0 + j;
			energyJobHandles[worldData.CloudLayer] = JobHandle.CombineDependencies(energyJobHandles[worldData.CloudLayer], SimJob.Schedule(new UpdateMassCloudCondensationJob()
			{
				CloudMass = nextState.CloudMass,
				CloudDropletMass = nextState.CloudDropletMass,

				CloudEvaporation = tempState.CloudEvaporationMass,
				CloudElevation = tempState.CloudElevation,
				LayerElevation = tempState.LayerElevation[j],
				LayerHeight = tempState.LayerHeight[j],
				CloudCondensation = tempState.CondensationCloudMass[j],
				GroundCondensation = tempState.CondensationGroundMass[j],
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.CloudLayer], energyJobHandles[layerIndex])));
			energyJobHandles[layerIndex] = energyJobHandles[worldData.CloudLayer];
		}
		energyJobHandles[worldData.SurfaceAirLayerGlobal] = SimJob.Schedule(new UpdateMassAirSurfaceJob()
		{
			AirTemperaturePotential = nextState.AirTemperaturePotential[worldData.SurfaceAirLayer],
			VaporMass = nextState.AirVapor[worldData.SurfaceAirLayer],
			DustMass = nextState.Dust[worldData.SurfaceAirLayer],
			CarbonDioxide = nextState.AirCarbon[worldData.SurfaceAirLayer],

			AirMass = tempState.AirMass[worldData.SurfaceAirLayer],
			EvaporationWater = tempState.EvaporationMassWater,
			EvaporationTemperatureWater = lastState.WaterTemperature[worldData.SurfaceWaterLayer],
			EvaporationFlora = tempState.FloraRespirationMassVapor,
			EvaporationTemperatureFlora = lastState.FloraTemperature,
			DustEjected = tempState.DustEjected,
			AirCarbonDelta = tempState.AirCarbonDelta,
			SoilRespiration = tempState.SoilRespiration,
			WaterCoverage = tempState.WaterCoverage[worldData.SurfaceWaterLayer],
			Elevation = lastState.Elevation,
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceAirLayerGlobal], energyJobHandles[worldData.SurfaceWaterLayerGlobal]));

		energyJobHandles[worldData.IceLayer] = SimJob.Schedule(new UpdateMassIceJob()
		{
			IceMass = nextState.IceMass,
			IceTemperature = nextState.IceTemperature,

			LastIceMass = lastState.IceMass,
			IceMelted = tempState.IceMeltedMass,
			WaterFrozen = tempState.FrozenMass,
			WaterTemperature = tempState.FrozenTemperature,
			Precipitation = tempState.PrecipitationMass,
			PrecipitationTemperature = tempState.PrecipitationTemperature,
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.IceLayer], energyJobHandles[worldData.CloudLayer]));

		energyJobHandles[worldData.TerrainLayer] = SimJob.Schedule(new UpdateTerrainJob()
		{
			SoilCarbon = nextState.GroundCarbon,
			Roughness = nextState.Roughness,
			GroundWater = nextState.GroundWater,
			Elevation = nextState.Elevation,
			LavaMass = nextState.LavaMass,
			CrustDepth = nextState.CrustDepth,
			MagmaMass = nextState.MagmaMass,
			LavaTemperature = nextState.LavaTemperature,

			CrustDelta = tempState.CrustDelta,
			LastElevation = lastState.Elevation,
			LastRoughness = lastState.Roughness,
			LastSoilFertility = lastState.GroundCarbon,
			LastGroundWater = lastState.GroundWater,
			GroundWaterConsumed = tempState.WaterConsumedByFlora,
			SoilRespiration = tempState.SoilRespiration,
			FloraDeath = tempState.FloraDeath,
			PlanktonDeath = tempState.PlanktonDeath,
			WaterCoverage = tempState.WaterCoverage[worldData.SurfaceWaterLayer],
			LastCrustDepth = lastState.CrustDepth,
			LastLavaMass = lastState.LavaMass,
			LastMagmaMass = lastState.MagmaMass,
			DustSettled = tempState.DustDown[worldData.SurfaceAirLayer],
			LavaCrystalized = tempState.LavaCrystalizedMass,
			LavaEjected = tempState.LavaEjected,
			MagmaTemperature = worldData.MagmaTemperature,
			LavaToRockMassAdjustment = worldData.LavaToRockMassAdjustment,
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.TerrainLayer], energyJobHandles[worldData.FloraLayer], energyJobHandles[worldData.SurfaceWaterLayerGlobal]));

		energyJobHandles[worldData.FloraLayer] = SimJob.Schedule(new UpdateFloraJob()
		{
			FloraMass = nextState.FloraMass,
			FloraWater = nextState.FloraWater,
			FloraGlucose = nextState.FloraGlucose,

			FloraGlucoseDelta = tempState.FloraGlucoseDelta,
			FloraMassDelta = tempState.FloraMassDelta,
			FloraWaterDelta = tempState.FloraWaterDelta,
			FloraWaterConsumed = tempState.WaterConsumedByFlora,
			LastMass = lastState.FloraMass,
			LastGlucose = lastState.FloraGlucose,
			LastWater = lastState.FloraWater,
		}, energyJobHandles[worldData.FloraLayer]);

		energyJobHandles[worldData.SurfaceAirLayerGlobal] = SimJob.Schedule(new UpdateWaterAirDiffusionJob()
		{
			AirCarbon = nextState.AirCarbon[worldData.SurfaceAirLayer],
			WaterCarbon = nextState.WaterCarbon[worldData.SurfaceWaterLayer],

			AirMass = tempState.AirMass[worldData.SurfaceAirLayer],
			WaterMass = nextState.WaterMass[worldData.SurfaceWaterLayer],
			SaltMass = nextState.SaltMass[worldData.SurfaceWaterLayer],
			WaterDepth = tempState.WaterLayerHeight[worldData.SurfaceWaterLayer],
			WaterAirCarbonDiffusionCoefficient = worldData.WaterAirCarbonDiffusionCoefficient,
			WaterAirCarbonDiffusionDepth = worldData.WaterAirCarbonDiffusionDepth,
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceAirLayerGlobal], energyJobHandles[worldData.SurfaceWaterLayerGlobal]));
		energyJobHandles[worldData.SurfaceWaterLayerGlobal] = energyJobHandles[worldData.SurfaceAirLayerGlobal];

		#endregion


		#region Apply Latent Heat

		energyJobHandles[worldData.IceLayer] = SimJob.Schedule(new ApplyLatentHeatIceJob()
		{
			IceTemperature = nextState.IceTemperature,
			IceMass = nextState.IceMass,
			LatentHeat = tempState.LatentHeat[worldData.IceLayer]
		}, energyJobHandles[worldData.IceLayer]);

		energyJobHandles[worldData.TerrainLayer] = SimJob.Schedule(new ApplyLatentHeatTerrainJob()
		{
			TerrainTemperature = nextState.GroundTemperature,

			LatentHeat = tempState.LatentHeat[worldData.TerrainLayer],
			SoilFertility = nextState.GroundCarbon,
			HeatingDepth = worldData.SoilHeatDepth
		}, energyJobHandles[worldData.TerrainLayer]);

		energyJobHandles[worldData.LavaLayer] = SimJob.Schedule(new ApplyLatentHeatLavaJob()
		{
			LavaTemperature = nextState.LavaTemperature,

			LatentHeat = tempState.LatentHeat[worldData.LavaLayer],
			LavaMass = nextState.LavaMass
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.LavaLayer], energyJobHandles[worldData.TerrainLayer]));

		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			energyJobHandles[worldData.AirLayer0 + i] = SimJob.Schedule(new ApplyLatentHeatAirJob()
			{
				AirTemperaturePotential = nextState.AirTemperaturePotential[i],
				AirMass = tempState.AirMass[i],
				VaporMass = nextState.AirVapor[i],
				LatentHeat = tempState.LatentHeat[worldData.AirLayer0 + i]
			}, energyJobHandles[worldData.AirLayer0 + i]);
		}
		for (int i = 1; i < worldData.WaterLayers - 1; i++)
		{
			energyJobHandles[worldData.WaterLayer0 + i] = SimJob.Schedule(new ApplyLatentHeatWaterJob()
			{
				WaterTemperature = nextState.WaterTemperature[i],
				WaterMass = nextState.WaterMass[i],
				SaltMass = nextState.SaltMass[i],
				LatentHeat = tempState.LatentHeat[worldData.WaterLayer0 + i]
			}, energyJobHandles[worldData.WaterLayer0 + i]);
		}

		#endregion


	}

	private JobHandle UpdateGroundWater(
		NativeArray<JobHandle> energyJobHandles,
		JobHandle lastJobHandle,
		ref SimState nextState,
		ref SimState lastState,
		ref TempState tempState,
		ref StaticState staticState,
		ref WorldData worldData,
		ref SimSettings settings
		)
	{
		if (settings.GroundWater)
		{

			lastJobHandle = SimJob.Schedule(new GroundWaterFlowJob()
			{
				GroundWater = nextState.GroundWater,
				GroundWaterTemperature = nextState.GroundWaterTemperature,

				LastGroundWater = lastState.GroundWater,
				LastGroundWaterTemperature = lastState.GroundWaterTemperature,
				SurfaceElevation = tempState.LayerElevation[worldData.SurfaceAirLayer],
				Neighbors = staticState.Neighbors,
				NeighborDistInverse = staticState.NeighborDistInverse,
				FlowSpeed = worldData.GroundWaterFlowSpeed,
				GroundWaterMaxInverse = 1.0f / worldData.GroundWaterMax,
			}, JobHandle.CombineDependencies(lastJobHandle, energyJobHandles[worldData.TerrainLayer]));

			lastJobHandle = SimJob.Schedule(new GroundWaterDiffusionJob()
			{
				GroundWater = tempState.GroundWaterFlowMass,
				GroundWaterTemperature = tempState.GroundWaterFlowTemperature,

				LastGroundWater = nextState.GroundWater,
				LastGroundWaterTemperature = nextState.GroundWaterTemperature,
				NeighborDist = staticState.NeighborDist,
				NeighborDistInverse = staticState.NeighborDistInverse,
				Neighbors = staticState.Neighbors,
				DiffusionCoefficient = worldData.GroundWaterDiffusionCoefficient
			}, lastJobHandle);

			for (int i = 1; i < worldData.WaterLayers - 1; i++)
			{
				lastJobHandle = SimJob.Schedule(new GroundWaterAbsorptionJob()
				{
					GroundWater = nextState.GroundWater,
					GroundWaterTemperature = nextState.GroundWaterTemperature,
					WaterMass = nextState.WaterMass[i],
					WaterTemperature = nextState.WaterTemperature[i],

					LastGroundWater = tempState.GroundWaterFlowMass,
					LastGroundWaterTemperature = tempState.GroundWaterFlowTemperature,
					SaltMass = nextState.SaltMass[i],
					WaterBelow = nextState.WaterMass[i - 1],
					GroundWaterAbsorptionRate = worldData.GroundWaterAbsorptionRate * worldData.SecondsPerTick,
					GroundWaterMaxInverse = 1.0f / worldData.GroundWaterMax,
					GroundWaterMax = worldData.GroundWaterMax,
					IsTop = i == worldData.SurfaceWaterLayer
				}, JobHandle.CombineDependencies(lastJobHandle, energyJobHandles[worldData.SurfaceWaterLayerGlobal],
					JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0 + i], energyJobHandles[worldData.WaterLayer0 + i - 1])));
				energyJobHandles[worldData.WaterLayer0 + i] = JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0 + i], lastJobHandle);
			}

			lastJobHandle = SimJob.Schedule(new GroundWaterConductionJob()
			{
				GroundWaterTemperature = tempState.GroundWaterFlowTemperature,
				TerrainTemperature = nextState.GroundTemperature,

				GroundWater = nextState.GroundWater,
				LastGroundWaterTemperature = nextState.GroundWaterTemperature,
				SoilFertility = nextState.GroundCarbon,
				GroundWaterConductionCoefficient = WorldData.ConductivityWaterTerrain,
				HeatingDepth = worldData.SoilHeatDepth,
				SecondsPerTick = worldData.SecondsPerTick,
				GroundWaterSurfaceAreaInverse = worldData.SoilHeatDepth / worldData.GroundWaterMaxDepth
			}, JobHandle.CombineDependencies(lastJobHandle, energyJobHandles[worldData.TerrainLayer]));
			energyJobHandles[worldData.TerrainLayer] = JobHandle.CombineDependencies(energyJobHandles[worldData.TerrainLayer], lastJobHandle);
		}
		return lastJobHandle;
	}
}
