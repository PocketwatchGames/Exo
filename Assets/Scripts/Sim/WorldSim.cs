﻿
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
using Unity.Profiling;

[Serializable]
public class WorldSim {

	public JobHelper SimJob;
	public JobHelper NeighborJob;
	public JobHelper AirJob;
	public JobHelper UpperAirJob;
	public JobHelper AirNeighborJob;
	public JobHelper WaterNeighborJob;
	public JobHelper WaterJob;


	private int _cellCount;

	private NativeArray<JobHandle> energyJobHandles;
	private List<NativeList<JobHandle>> jobHandleDependencies = new List<NativeList<JobHandle>>();
	private List<NativeArray<float>> tempArrays = new List<NativeArray<float>>();

	public WorldSim(int cellCount, ref WorldData worldData)
	{
		_cellCount = cellCount;

		SimJob = new JobHelper(_cellCount);
		NeighborJob = new JobHelper(_cellCount * 6);
		WaterJob = new JobHelper(_cellCount * (worldData.WaterLayers - 2));
		AirJob = new JobHelper(_cellCount * (worldData.AirLayers - 2));
		UpperAirJob = new JobHelper(_cellCount * (worldData.AirLayers - 3));
		AirNeighborJob = new JobHelper(_cellCount * (worldData.AirLayers - 2) * StaticState.MaxNeighborsVert);
		WaterNeighborJob = new JobHelper(_cellCount * (worldData.WaterLayers - 2) * StaticState.MaxNeighborsVert);
		energyJobHandles = new NativeArray<JobHandle>(worldData.LayerCount, Allocator.Persistent);

	}

	public void Dispose(ref WorldData worldData)
	{
		energyJobHandles.Dispose();

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
		tickJobHandle = tempState.Update(ref lastState, ref worldData, ref staticState, ref settings, tickJobHandle);

		for (int i=0;i<energyJobHandles.Length;i++)
		{
			energyJobHandles[i] = tickJobHandle;
		}

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

		tickJobHandle = Utils.MemCopy(nextState.Plate, lastState.Plate, tickJobHandle);

		DoEnergyCycle(energyJobHandles, tickJobHandle, ref nextState, ref lastState, ref tempState, ref staticState, ref worldData, ref settings);

		DoStateChange(energyJobHandles, tickJobHandle, ref nextState, ref lastState, ref tempState, ref staticState, ref worldData, ref settings);

		var groundWaterJob = UpdateGroundWater(energyJobHandles, tickJobHandle, ref nextState, ref lastState, ref tempState, ref staticState, ref worldData, ref settings);

		#region Air Advection and Diffusion
		// Buoyancy, Updrafts, and mixing occur across air layers and water layers
		// TODO: add an empty air layer on top and bottom so we can calculate up/down diffusion in a single step 
		// Temperature and trace elements diffuse into neighboring horizontal cells based on a diffusion constant
		// Air, Water, Cloud
		#region Update Velocity


		var airTerrainFrictionJobHandle = SimJob.Schedule(
			JobType.Schedule, 64,
			new AirTerrainFrictionJob()
			{
				Force = tempState.WindFriction,
				IceCoverage = tempState.IceCoverage,
				WaterCoverage = staticState.GetSliceLayer(tempState.WaterCoverage,worldData.SurfaceWaterLayer),
				FloraCoverage = tempState.FloraCoverage,
				Roughness = lastState.Roughness,
				IceFriction = worldData.WindIceFriction,
				TerrainFrictionMin = worldData.WindTerrainFrictionMin,
				TerrainFrictionMax = worldData.WindTerrainFrictionMax,
				FloraFriction = worldData.WindFloraFriction,
				WaterFriction = worldData.WindWaterFriction,
				MaxTerrainRoughness = worldData.MaxTerrainRoughnessForWindFriction
			}, energyJobHandles[worldData.AirLayer0]);

		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
			JobType.Schedule, 64,
			new AccelerationAirJob()
			{
				Velocity = staticState.GetSliceAir(nextState.AirVelocity),
				Force = staticState.GetSliceAir(tempState.AirAcceleration),

				LastVelocity = staticState.GetSliceAir(lastState.AirVelocity),
				Pressure = staticState.GetSliceAir(tempState.AirPressure),
				AirMass = staticState.GetSliceAir(tempState.AirMass),
				TemperaturePotential = staticState.GetSliceAir(lastState.AirTemperaturePotential),
				NewTemperaturePotential = staticState.GetSliceAir(nextState.AirTemperaturePotential),
				VaporMass = staticState.GetSliceAir(lastState.AirVapor),
				LayerMiddle = staticState.GetSliceAir(tempState.AirLayerMiddle),
				Friction = tempState.WindFriction,
				Neighbors = staticState.Neighbors,
				NeighborDiffInverse = staticState.NeighborDiffInverse,
				Positions = staticState.SphericalPosition,
				PlanetRadius = staticState.PlanetRadius,
				Gravity = lastState.PlanetState.Gravity,
				GravityInverse = 1.0f / lastState.PlanetState.Gravity,
				SecondsPerTick = worldData.SecondsPerTick,
				LayerCount = worldData.AirLayers - 2,
				Count = staticState.Count

			},airTerrainFrictionJobHandle);

		#endregion

		// Wind and currents move temperature and trace elements horizontally
		// Air, Water, Cloud
		// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
		#region Advection

		if (settings.AdvectionAir)
		{
			energyJobHandles[worldData.AirLayer0] = AirNeighborJob.Schedule(
				JobType.Schedule, 64,
				new GetVectorDestCoordsVerticalJob()
				{
					Destination = staticState.GetSliceAirNeighbors(tempState.DestinationAir),
					NeighborsVert = staticState.NeighborsVert,
					Position = staticState.SphericalPosition,
					Velocity = nextState.AirVelocity,
					Mass = tempState.AirMass,
					LayerHeight = tempState.AirLayerHeight,
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick,
					CellsPerLayer = staticState.Count,
					LayerCount = worldData.AirLayers
				}, energyJobHandles[worldData.AirLayer0]);

			energyJobHandles[worldData.AirLayer0] = AirNeighborJob.Schedule(
				JobType.Schedule, 64,
				new ResolveAdvectionConflictVert()
				{
					ResolvedDestination = staticState.GetSliceAirNeighbors(tempState.DestinationAirResolved),
					Destination = staticState.GetSliceAirNeighbors(tempState.DestinationAir),
					ReverseNeighborsVert = staticState.GetSliceAirNeighbors(staticState.ReverseNeighborsVert),
					Count = staticState.Count,
				}, energyJobHandles[worldData.AirLayer0]);

			if (settings.MakeAirIncompressible)
			{

				energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
					JobType.Schedule, 64,
					new GetDivergenceJob()
					{
						Divergence = staticState.GetSliceAir(tempState.DivergenceAir),
						Destination = staticState.GetSliceAirNeighbors(tempState.DestinationAirResolved),
					}, energyJobHandles[worldData.AirLayer0]);

				// Calculate Pressure gradient field
				for (int a = 0; a < settings.IncompressibilityIterations; a++)
				{
					var dpj = new GetDivergencePressureJob()
					{
						Pressure = staticState.GetSliceAir(tempState.DivergencePressureAir),
						Divergence = staticState.GetSliceAir(tempState.DivergenceAir),
						NeighborsVert = staticState.GetSliceAirNeighbors(staticState.NeighborsVert),
						Count = staticState.Count
					};
					energyJobHandles[worldData.AirLayer0] = dpj.Schedule(_cellCount * (worldData.AirLayers - 2), energyJobHandles[worldData.AirLayer0]);
				}

				energyJobHandles[worldData.AirLayer0] = AirNeighborJob.Schedule(
					JobType.Schedule, 64,
					new GetDivergenceFreeFieldJob()
					{
						Destination = staticState.GetSliceAirNeighbors(tempState.DestinationAirResolved),
						NeighborsVert = staticState.GetSliceAirNeighbors(staticState.NeighborsVert),
						Pressure = staticState.GetSliceAir(tempState.DivergencePressureAir),
						Count = staticState.Count
					}, energyJobHandles[worldData.AirLayer0]);

				energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
					JobType.Schedule, 64,
					new SumMassLeavingJob()
					{
						MassLeaving = staticState.GetSliceAir(tempState.AirMassLeaving),
						Destination = staticState.GetSliceAirNeighbors(tempState.DestinationAirResolved),
					}, energyJobHandles[worldData.AirLayer0]);

				energyJobHandles[worldData.AirLayer0] = AirNeighborJob.Schedule(
					JobType.Schedule, 64,
					new CapMassLeavingJob()
					{
						Destination = staticState.GetSliceAirNeighbors(tempState.DestinationAirResolved),
						NeighborsVert = staticState.GetSliceAirNeighbors(staticState.NeighborsVert),
						MassLeaving = staticState.GetSliceAir(tempState.AirMassLeaving),
						Mass = staticState.GetSliceAir(tempState.AirMass),
						Count = staticState.Count
					}, energyJobHandles[worldData.AirLayer0]);

				//energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
				//	JobType.Schedule, 64,
				//	new UpdateDivergenceFreeVelocityJob()
				//	{
				//		Velocity = staticState.GetSliceAir(nextState.AirVelocity),
				//		Mass = staticState.GetSliceAir(tempState.AirMass),
				//		LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
				//		DestinationVert = staticState.GetSliceAirNeighbors(tempState.DestinationAirResolved),
				//		NeighborTangent = staticState.NeighborTangent,
				//		Positions = staticState.SphericalPosition,
				//		TicksPerSecond = worldData.TicksPerSecond,
				//		Count = staticState.Count,
				//	}, energyJobHandles[worldData.AirLayer0]);
			}


			energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
				settings.SynchronousOverrides.AdvectionAir, 64,
				new AdvectionAirJob()
				{
					Delta = staticState.GetSliceAir(tempState.AdvectionAir),
					Temperature = nextState.AirTemperaturePotential,
					AirMass = tempState.AirMass,
					Vapor = nextState.AirVapor,
					CarbonDioxide = nextState.AirCarbonDioxide,
					Dust = nextState.AirDust,
					Velocity = nextState.AirVelocity,
					NeighborsVert = staticState.NeighborsVert,
					Destination = tempState.DestinationAirResolved,
					Positions = staticState.SphericalPosition,
					CoriolisMultiplier = staticState.CoriolisMultiplier,
					CoriolisTerm = coriolisTerm,
					SecondsPerTick = worldData.SecondsPerTick,
					CellsPerLayer = staticState.Count
				}, energyJobHandles[worldData.AirLayer0]);

			energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
				JobType.Schedule, 64,
				new ApplyAdvectionAirJob()
				{
					Advection = staticState.GetSliceAir(tempState.AdvectionAir),
					Vapor = staticState.GetSliceAir(nextState.AirVapor),
					Dust = staticState.GetSliceAir(nextState.AirDust),
					CarbonDioxide = staticState.GetSliceAir(nextState.AirCarbonDioxide),
					Oxygen = staticState.GetSliceAir(nextState.AirOxygen),
					Nitrogen = staticState.GetSliceAir(nextState.AirNitrogen),
					Methane = staticState.GetSliceAir(nextState.AirMethane),
					Minerals = staticState.GetSliceAir(nextState.AirMinerals),
					Temperature = staticState.GetSliceAir(nextState.AirTemperaturePotential),
					AirVelocity = staticState.GetSliceAir(nextState.AirVelocity),
			}, energyJobHandles[worldData.AirLayer0]);
		}
		#endregion

		// Diffuse from last time step
		// Air, Water, Cloud
		#region Diffusion

		if (settings.DiffusionAir)
		{
			// TODO: is it a problem that we are using the dependent variables from last frame while referencing our newly calculated next frame values for temperature and such?
			energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
				settings.SynchronousOverrides.DiffusionAir, 64,
				new DiffusionAirJob()
				{
					Delta = staticState.GetSliceAir(tempState.DiffusionAir),

					AirMass = staticState.GetSliceAir(tempState.AirMass),
					Temperature = staticState.GetSliceAir(nextState.AirTemperaturePotential),
					Vapor = staticState.GetSliceAir(nextState.AirVapor),
					CarbonDioxide = staticState.GetSliceAir(nextState.AirCarbonDioxide),
					Dust = staticState.GetSliceAir(nextState.AirDust),
					Velocity = staticState.GetSliceAir(nextState.AirVelocity),
					LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
					Neighbors = staticState.Neighbors,
					NeighborDistInverse = staticState.NeighborDistInverse,
					DiffusionCoefficientHorizontal = worldData.AirDiffusionCoefficientHorizontal,
					DiffusionCoefficientVertical = worldData.AirDiffusionCoefficientVertical,
					CellSurfaceArea = staticState.CellSurfaceArea,
					CellCircumference = staticState.CellCircumference,
					LayerCount = worldData.AirLayers - 2,
					Count = staticState.Count
				}, energyJobHandles[worldData.AirLayer0]);

			energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
				settings.SynchronousOverrides.DiffusionAir, 64,
				new ApplyAdvectionAirJob()
				{
					Advection = staticState.GetSliceAir(tempState.DiffusionAir),
					Vapor = staticState.GetSliceAir(nextState.AirVapor),
					Dust = staticState.GetSliceAir(nextState.AirDust),
					CarbonDioxide = staticState.GetSliceAir(nextState.AirCarbonDioxide),
					Oxygen = staticState.GetSliceAir(nextState.AirOxygen),
					Nitrogen = staticState.GetSliceAir(nextState.AirNitrogen),
					Methane = staticState.GetSliceAir(nextState.AirMethane),
					Minerals = staticState.GetSliceAir(nextState.AirMinerals),
					Temperature = staticState.GetSliceAir(nextState.AirTemperaturePotential),
					AirVelocity = staticState.GetSliceAir(nextState.AirVelocity),
				}, energyJobHandles[worldData.AirLayer0]);
		}


#endregion

#endregion

#region Water and Cloud Advection and Diffusion

// Buoyancy, Updrafts, and mixing occur across air layers and water layers
// TODO: add an empty air layer on top and bottom so we can calculate up/down diffusion in a single step 
// Temperature and trace elements diffuse into neighboring horizontal cells based on a diffusion constant
// Air, Water, Cloud
#region Update Velocity

		energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
			JobType.Schedule, 64,
			new WaterSurfaceFrictionJob()
			{
				Force = tempState.WaterFriction,

				Position = staticState.SphericalPosition,
				Current = staticState.GetSliceLayer(lastState.WaterVelocity,worldData.SurfaceWaterLayer),
				AirVelocityUp = staticState.GetSliceLayer(lastState.AirVelocity,worldData.SurfaceAirLayer),
				AirVelocityDown = staticState.GetSliceLayer(lastState.WaterVelocity,worldData.SurfaceWaterLayer - 1),
				LayerHeight = staticState.GetSliceLayer(tempState.WaterLayerHeight,worldData.SurfaceWaterLayer),
				CoriolisMultiplier = staticState.CoriolisMultiplier,
				FrictionCoefficientUp = worldData.WindToWaterCurrentFrictionCoefficient,
				FrictionCoefficientDown = 0, // TODO: do we want to add a frictional force between layers of water?
				CoriolisTerm = coriolisTerm,
				WaterSurfaceFrictionDepth = worldData.WaterSurfaceFrictionDepth,
				SecondsPerTick = worldData.SecondsPerTick
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0], energyJobHandles[worldData.SurfaceAirLayerGlobal]));

		// TODO: since I moved acceleration up here, we will need to add a buoyancy job down after the temperature adjustment to incorporate positive buoyancy in the vertical velocity
		// (see usage of "next state" velocity
		energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
			JobType.Schedule, 64,
			new AccelerationWaterJob()
			{
				Velocity = staticState.GetSliceWater(nextState.WaterVelocity),

				LastVelocity = staticState.GetSliceWater(lastState.WaterVelocity),
				WaterDensity = staticState.GetSliceWater(tempState.WaterDensity),
				WaterPressure = staticState.GetSliceWater(tempState.WaterPressure),
				LayerDepth = staticState.GetSliceWater(tempState.WaterLayerDepth),
				LayerHeight = staticState.GetSliceWater(tempState.WaterLayerHeight),
				SurfaceElevation = tempState.SurfaceElevation,
				Friction = tempState.WaterFriction,
				Positions = staticState.SphericalPosition,
				Neighbors = staticState.Neighbors,
				NeighborDiffInverse = staticState.NeighborDiffInverse,
				Gravity = lastState.PlanetState.Gravity,
				PlanetRadius = staticState.PlanetRadius,
				SecondsPerTick = worldData.SecondsPerTick,
				Count = staticState.Count,
				LayerCount = worldData.WaterLayers - 2
			
			}, energyJobHandles[worldData.WaterLayer0]);


		#endregion

		// Wind and currents move temperature and trace elements horizontally
		// Air, Water, Cloud
		// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
		#region Advection

		if (settings.AdvectionWater)
		{
			energyJobHandles[worldData.WaterLayer0] = WaterNeighborJob.Schedule(
			settings.SynchronousOverrides.AdvectionWater, 64,
				new GetVectorDestCoordsVerticalJob()
				{
					Destination = staticState.GetSliceWaterNeighbors(tempState.DestinationWater),
					NeighborsVert = staticState.NeighborsVert,
					Position = staticState.SphericalPosition,
					Velocity = nextState.WaterVelocity,
					Mass = nextState.WaterMass,
					LayerHeight = tempState.WaterLayerHeight,
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick,
					CellsPerLayer = staticState.Count,
					LayerCount = worldData.WaterLayers
				}, energyJobHandles[worldData.WaterLayer0]);

			energyJobHandles[worldData.WaterLayer0] = WaterNeighborJob.Schedule(
				settings.SynchronousOverrides.AdvectionWater, 64,
				new ResolveAdvectionConflictVert()
				{
					ResolvedDestination = staticState.GetSliceWaterNeighbors(tempState.DestinationWaterResolved),
					Destination = staticState.GetSliceWaterNeighbors(tempState.DestinationWater),
					ReverseNeighborsVert = staticState.GetSliceWaterNeighbors(staticState.ReverseNeighborsVert),
					Count = staticState.Count
				}, energyJobHandles[worldData.WaterLayer0]);

			if (settings.MakeWaterIncompressible)
			{
				energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
					settings.SynchronousOverrides.AdvectionWater, 64,
					new GetDivergenceJob()
					{
						Divergence = staticState.GetSliceWater(tempState.DivergenceWater),
						Destination = staticState.GetSliceWaterNeighbors(tempState.DestinationWaterResolved),
					}, energyJobHandles[worldData.WaterLayer0]);

				// Calculate Pressure gradient field
				for (int a = 0; a < settings.IncompressibilityIterations; a++)
				{
					var dpj = new GetDivergencePressureJob()
					{
						Pressure = staticState.GetSliceWater(tempState.DivergencePressureWater),
						Divergence = staticState.GetSliceWater(tempState.DivergenceWater),
						NeighborsVert = staticState.GetSliceWaterNeighbors(staticState.NeighborsVert),
						Count = staticState.Count
					};
					if (settings.SynchronousOverrides.AdvectionWater)
					{
						energyJobHandles[worldData.WaterLayer0].Complete();
						dpj.Run(_cellCount * (worldData.WaterLayers - 2));
						energyJobHandles[worldData.WaterLayer0] = default(JobHandle);
					}
					else
					{
						energyJobHandles[worldData.WaterLayer0] = dpj.Schedule(_cellCount * (worldData.WaterLayers - 2), energyJobHandles[worldData.WaterLayer0]);
					}
				}

				energyJobHandles[worldData.WaterLayer0] = WaterNeighborJob.Schedule(
					settings.SynchronousOverrides.AdvectionWater, 64,
					new GetDivergenceFreeFieldJob()
					{
						Destination = staticState.GetSliceWaterNeighbors(tempState.DestinationWaterResolved),
						NeighborsVert = staticState.GetSliceWaterNeighbors(staticState.NeighborsVert),
						Pressure = staticState.GetSliceWater(tempState.DivergencePressureWater),
						Count = staticState.Count
					}, energyJobHandles[worldData.WaterLayer0]);

				energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
					settings.SynchronousOverrides.AdvectionWater, 64,
					new SumMassLeavingJob()
					{
						MassLeaving = staticState.GetSliceWater(tempState.WaterMassLeaving),
						Destination = staticState.GetSliceWaterNeighbors(tempState.DestinationWaterResolved),
					}, energyJobHandles[worldData.WaterLayer0]);

				energyJobHandles[worldData.WaterLayer0] = WaterNeighborJob.Schedule(
					settings.SynchronousOverrides.AdvectionWater, 64,
					new CapMassLeavingJob()
					{
						Destination = staticState.GetSliceWaterNeighbors(tempState.DestinationWaterResolved),
						NeighborsVert = staticState.GetSliceWaterNeighbors(staticState.NeighborsVert),
						MassLeaving = staticState.GetSliceWater(tempState.WaterMassLeaving),
						Mass = staticState.GetSliceWater(nextState.WaterMass),
						Count = staticState.Count
					}, energyJobHandles[worldData.WaterLayer0]);

				//energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
				//	settings.SynchronousOverrides.AdvectionWater, 64,
				//	new UpdateDivergenceFreeVelocityJob()
				//	{
				//		Velocity = staticState.GetSliceWater(nextState.WaterVelocity),
				//		Mass = staticState.GetSliceWater(nextState.WaterMass),
				//		LayerHeight = staticState.GetSliceWater(tempState.WaterLayerHeight),
				//		DestinationVert = staticState.GetSliceWaterNeighbors(tempState.DestinationWaterResolved),
				//		NeighborTangent = staticState.NeighborTangent,
				//		Positions = staticState.SphericalPosition,
				//		TicksPerSecond = worldData.TicksPerSecond,
				//		Count = staticState.Count
				//	}, energyJobHandles[worldData.WaterLayer0]);
			}


			energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
				settings.SynchronousOverrides.AdvectionWater, 64,
				new AdvectionWaterJob()
				{
					Delta = staticState.GetSliceWater(tempState.AdvectionWater),
					Destination = tempState.DestinationWaterResolved,
					Velocity = nextState.WaterVelocity,
					Temperature = nextState.WaterTemperature,
					Mass = nextState.WaterMass,
					Salt = nextState.WaterSaltMass,
					Carbon = nextState.WaterCarbonDioxide,
					Nitrogen = nextState.WaterNitrogen,
					Glucose = nextState.WaterGlucose,
					Mineral = nextState.WaterMinerals,
					Oxygen = nextState.WaterOxygen,
					Positions = staticState.SphericalPosition,
					NeighborsVert = staticState.NeighborsVert,
					CoriolisMultiplier = staticState.CoriolisMultiplier,
					CoriolisTerm = coriolisTerm,
					SecondsPerTick = worldData.SecondsPerTick,
					CellsPerLayer = staticState.Count
				}, JobHandle.CombineDependencies(groundWaterJob, energyJobHandles[worldData.WaterLayer0]));
			energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
				settings.SynchronousOverrides.AdvectionWater, 64,
				new ApplyAdvectionWaterJob()
				{
					SaltMass = staticState.GetSliceWater(nextState.WaterSaltMass),
					CarbonMass = staticState.GetSliceWater(nextState.WaterCarbonDioxide),
					NitrogenMass = staticState.GetSliceWater(nextState.WaterNitrogen),
					GlucoseMass = staticState.GetSliceWater(nextState.WaterGlucose),
					MineralMass = staticState.GetSliceWater(nextState.WaterMinerals),
					OxygenMass = staticState.GetSliceWater(nextState.WaterOxygen),
					Temperature = staticState.GetSliceWater(nextState.WaterTemperature),
					Velocity = staticState.GetSliceWater(nextState.WaterVelocity),
					WaterMass = staticState.GetSliceWater(nextState.WaterMass),
					Advection = staticState.GetSliceWater(tempState.AdvectionWater),
				}, energyJobHandles[worldData.WaterLayer0]);

		}

		if (settings.AdvectionCloud)
		{
			energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(
				settings.SynchronousOverrides.AdvectionCloud, 64,
				new GetVectorDestCoordsJob()
				{
					Destination = tempState.DestinationCloud,
					Neighbors = staticState.Neighbors,
					Position = staticState.SphericalPosition,
					Velocity = tempState.CloudVelocity,
					PlanetRadius = staticState.PlanetRadius,
					SecondsPerTick = worldData.SecondsPerTick,
					MaxWindMove = staticState.CellRadius * 0.9f,
				}, energyJobHandles[worldData.CloudLayer]);

			energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(
				settings.SynchronousOverrides.AdvectionCloud, 64,
				new AdvectionCloudJob()
				{
					Delta = tempState.AdvectionCloud,
					Destination = tempState.DestinationCloud,
					Mass = nextState.CloudMass,
					Temperature = nextState.CloudTemperature,
					DropletMass = nextState.CloudDropletMass,
					Neighbors = staticState.Neighbors,
				}, energyJobHandles[worldData.CloudLayer]);

			energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(
				settings.SynchronousOverrides.AdvectionCloud, 64,
				new ApplyAdvectionCloudJob()
				{
					Advection = tempState.AdvectionCloud,
					CloudMass = nextState.CloudMass,
					Temperature = nextState.CloudTemperature,
					DropletMass = nextState.CloudDropletMass,
				}, energyJobHandles[worldData.CloudLayer]);
		}
		#endregion

		// Diffuse from last time step
		// Air, Water, Cloud
		#region Diffusion

		if (settings.DiffusionWater)
		{
			energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
				settings.SynchronousOverrides.DiffusionWater, 64,
				new DiffusionWaterJob()
				{
					Delta = staticState.GetSliceWater(tempState.DiffusionWater),

					Temperature = staticState.GetSliceWater(nextState.WaterTemperature),
					SaltMass = staticState.GetSliceWater(nextState.WaterSaltMass),
					CarbonMass = staticState.GetSliceWater(nextState.WaterCarbonDioxide),
					OxygenMass = staticState.GetSliceWater(nextState.WaterOxygen),
					NitrogenMass = staticState.GetSliceWater(nextState.WaterNitrogen),
					GlucoseMass = staticState.GetSliceWater(nextState.WaterGlucose),
					MineralsMass = staticState.GetSliceWater(nextState.WaterMinerals),
					
					Velocity = staticState.GetSliceWater(nextState.WaterVelocity),
					WaterMass = staticState.GetSliceWater(nextState.WaterMass),
					LayerHeight = staticState.GetSliceWater(tempState.WaterLayerHeight),
					NeighborDistInverse = staticState.NeighborDistInverse,
					Neighbors = staticState.Neighbors,
					DiffusionCoefficientHorizontal = worldData.WaterDiffusionCoefficientHorizontal,
					DiffusionCoefficientVertical = worldData.WaterDiffusionCoefficientVertical,
					CellSurfaceArea = staticState.CellSurfaceArea,
					CellCircumference = staticState.CellCircumference,
					Count = staticState.Count,
					LayerCount = worldData.WaterLayers - 2
				}, energyJobHandles[worldData.WaterLayer0]);

			energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
				settings.SynchronousOverrides.DiffusionWater, 64,
				new ApplyAdvectionWaterJob()
				{
					Advection = staticState.GetSliceWater(tempState.DiffusionWater),
					SaltMass = staticState.GetSliceWater(nextState.WaterSaltMass),
					CarbonMass = staticState.GetSliceWater(nextState.WaterCarbonDioxide),
					NitrogenMass = staticState.GetSliceWater(nextState.WaterNitrogen),
					GlucoseMass = staticState.GetSliceWater(nextState.WaterGlucose),
					MineralMass = staticState.GetSliceWater(nextState.WaterMinerals),
					OxygenMass = staticState.GetSliceWater(nextState.WaterOxygen),
					Temperature = staticState.GetSliceWater(nextState.WaterTemperature),
					Velocity = staticState.GetSliceWater(nextState.WaterVelocity),
					WaterMass = staticState.GetSliceWater(nextState.WaterMass)
				}, energyJobHandles[worldData.WaterLayer0]);
		}

		if (settings.DiffusionCloud)
		{
			energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(
				settings.SynchronousOverrides.DiffusionCloud, 64,
				new DiffusionCloudJob()
				{
					Delta = tempState.DiffusionCloud,

					LastMass = nextState.CloudMass,
					LastTemperature = nextState.CloudTemperature,
					LastDropletMass = nextState.CloudDropletMass,
					Neighbors = staticState.Neighbors,
					DiffusionCoefficient = worldData.CloudDiffusionCoefficient,
				}, energyJobHandles[worldData.CloudLayer]);

			energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(
				settings.SynchronousOverrides.DiffusionCloud, 64,
				new ApplyAdvectionCloudJob()
				{
					Advection = tempState.DiffusionCloud,
					CloudMass = nextState.CloudMass,
					Temperature = nextState.CloudTemperature,
					DropletMass = nextState.CloudDropletMass,
				}, energyJobHandles[worldData.CloudLayer]);
		}

#endregion

#endregion


		tickJobHandle = JobHandle.CombineDependencies(
			groundWaterJob,
			JobHandle.CombineDependencies(energyJobHandles)
			);

		// TODO: we really just want to update depths
#region Update dependent variables

		tickJobHandle = tempState.Update(ref nextState, ref worldData, ref staticState, ref settings, tickJobHandle);

		#endregion


		#region Flow

		energyJobHandles[worldData.WaterLayer0] = JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0], tickJobHandle);
		if (settings.WaterSurfaceFlowEnabled)
		{
			// TODO: surface elevation is inaccurate now, we should recalculate (and use water surfae, not ice surface)
			energyJobHandles[worldData.WaterLayer0] = NeighborJob.Schedule(
				JobType.Schedule, 64,
				new UpdateFlowVelocityJob()
				{
					Flow = nextState.FlowWater,

					LastFlow = lastState.FlowWater,
					SurfaceElevation = tempState.SurfaceElevation,
					WaterDepth = staticState.GetSliceLayer(tempState.WaterLayerHeight,worldData.SurfaceWaterLayer),
					NeighborDistInverse = staticState.NeighborDistInverse,
					Neighbors = staticState.Neighbors,
					Gravity = nextState.PlanetState.Gravity,
					SecondsPerTick = worldData.SecondsPerTick,
					Damping = worldData.SurfaceWaterFlowDamping,
					ViscosityInverse = 1.0f - worldData.WaterViscosity,
				}, energyJobHandles[worldData.WaterLayer0]);

			energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
				JobType.Schedule, 64,
				new SumOutgoingFlowJob()
				{
					OutgoingFlow = tempState.OutgoingFlowWater,
					Flow = nextState.FlowWater
				}, energyJobHandles[worldData.WaterLayer0]);

			energyJobHandles[worldData.WaterLayer0] = NeighborJob.Schedule(
				JobType.Schedule, 64,
				new LimitOutgoingFlowJob()
				{
					Flow = nextState.FlowWater,
					FlowPercent = tempState.FlowPercentWater,

					OutgoingFlow = tempState.OutgoingFlowWater,
					Neighbors = staticState.Neighbors,
					WaterDepth = staticState.GetSliceLayer(tempState.WaterLayerHeight,worldData.SurfaceWaterLayer)
				}, energyJobHandles[worldData.WaterLayer0]);

			energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
				JobType.Schedule, 64,
				new ApplyFlowWaterJob()
				{
					Delta = staticState.GetSliceLayer(tempState.DiffusionWater,worldData.SurfaceWaterLayer),

					Mass = staticState.GetSliceLayer(nextState.WaterMass,worldData.SurfaceWaterLayer),
					Velocity = staticState.GetSliceLayer(nextState.WaterVelocity,worldData.SurfaceWaterLayer),
					CarbonDioxide = staticState.GetSliceLayer(nextState.WaterCarbonDioxide,worldData.SurfaceWaterLayer),
					Oxygen = staticState.GetSliceLayer(nextState.WaterOxygen, worldData.SurfaceWaterLayer),
					Nitrogen = staticState.GetSliceLayer(nextState.WaterNitrogen, worldData.SurfaceWaterLayer),
					Glucose = staticState.GetSliceLayer(nextState.WaterGlucose, worldData.SurfaceWaterLayer),
					Minerals = staticState.GetSliceLayer(nextState.WaterMinerals, worldData.SurfaceWaterLayer),
					Salt = staticState.GetSliceLayer(nextState.WaterSaltMass,worldData.SurfaceWaterLayer),
					Temperature = staticState.GetSliceLayer(nextState.WaterTemperature,worldData.SurfaceWaterLayer),
					FlowPercent = tempState.FlowPercentWater,
					Positions = staticState.SphericalPosition,
					Neighbors = staticState.Neighbors,
					ReverseNeighbors = staticState.ReverseNeighbors,
					CoriolisMultiplier = staticState.CoriolisMultiplier,
					CoriolisTerm = coriolisTerm,
					SecondsPerTick = worldData.SecondsPerTick
				}, energyJobHandles[worldData.WaterLayer0]);

			energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
				JobType.Schedule, 64,
				new ApplyAdvectionWaterJob()
				{
					WaterMass = staticState.GetSliceLayer(nextState.WaterMass,worldData.SurfaceWaterLayer),
					SaltMass = staticState.GetSliceLayer(nextState.WaterSaltMass,worldData.SurfaceWaterLayer),
					CarbonMass = staticState.GetSliceLayer(nextState.WaterCarbonDioxide,worldData.SurfaceWaterLayer),
					NitrogenMass = staticState.GetSliceLayer(nextState.WaterNitrogen, worldData.SurfaceWaterLayer),
					GlucoseMass = staticState.GetSliceLayer(nextState.WaterGlucose, worldData.SurfaceWaterLayer),
					MineralMass = staticState.GetSliceLayer(nextState.WaterMinerals, worldData.SurfaceWaterLayer),
					OxygenMass = staticState.GetSliceLayer(nextState.WaterOxygen, worldData.SurfaceWaterLayer),
					Temperature = staticState.GetSliceLayer(nextState.WaterTemperature,worldData.SurfaceWaterLayer),
					Velocity = staticState.GetSliceLayer(nextState.WaterVelocity,worldData.SurfaceWaterLayer),

					Advection = staticState.GetSliceLayer(tempState.DiffusionWater,worldData.SurfaceWaterLayer),
				}, energyJobHandles[worldData.WaterLayer0]);
		}

		#endregion

		#region Rebalance

		if (settings.RebalanceWaterLayers)
		{
			for (int i = worldData.SurfaceWaterLayer; i >= 2; i--)
			{

				float maxDepth;
				float minDepth;
				// TODO: this doesn't handle more than 3 layers
				if (i == worldData.SurfaceWaterLayer)
				{
					maxDepth = worldData.SurfaceWaterDepth;
					minDepth = worldData.SurfaceWaterDepth - 10;
				}
				else
				{
					maxDepth = worldData.ThermoclineDepth;
					minDepth = worldData.ThermoclineDepth - 10;
				}
				energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
					JobType.Schedule, 64,
					new RebalanceWaterLayersLimitJob()
					{
						Delta1 = tempState.RebalanceWater1,
						Delta2 = tempState.RebalanceWater2,

						Mass = nextState.WaterMass,
						Salt = nextState.WaterSaltMass,
						Carbon = nextState.WaterCarbonDioxide,
						Nitrogen = nextState.WaterNitrogen,
						Glucose = nextState.WaterGlucose,
						Minerals = nextState.WaterMinerals,
						Oxygen = nextState.WaterOxygen,
						Temperature = nextState.WaterTemperature,
						Velocity = nextState.WaterVelocity,

						Layer1 = i,
						Layer2 = i - 1,

						MaxDepth1 = maxDepth,
						MinDepth1 = minDepth,

						Count = staticState.Count
					}, energyJobHandles[worldData.WaterLayer0]);
				energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
					JobType.Schedule, 64,
					new ApplyAdvectionWaterJob()
					{
						WaterMass = staticState.GetSliceLayer(nextState.WaterMass, i),
						SaltMass = staticState.GetSliceLayer(nextState.WaterSaltMass, i),
						CarbonMass = staticState.GetSliceLayer(nextState.WaterCarbonDioxide, i),
						NitrogenMass = staticState.GetSliceLayer(nextState.WaterNitrogen, i),
						GlucoseMass = staticState.GetSliceLayer(nextState.WaterGlucose, i),
						MineralMass = staticState.GetSliceLayer(nextState.WaterMinerals, i),
						OxygenMass = staticState.GetSliceLayer(nextState.WaterOxygen, i),
						Temperature = staticState.GetSliceLayer(nextState.WaterTemperature, i),
						Velocity = staticState.GetSliceLayer(nextState.WaterVelocity, i),

						Advection = tempState.RebalanceWater1,
					}, energyJobHandles[worldData.WaterLayer0]);
				energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
					JobType.Schedule, 64,
					new ApplyAdvectionWaterJob()
					{
						WaterMass = staticState.GetSliceLayer(nextState.WaterMass, i - 1),
						SaltMass = staticState.GetSliceLayer(nextState.WaterSaltMass, i - 1),
						CarbonMass = staticState.GetSliceLayer(nextState.WaterCarbonDioxide, i - 1),
						NitrogenMass = staticState.GetSliceLayer(nextState.WaterNitrogen, i - 1),
						GlucoseMass = staticState.GetSliceLayer(nextState.WaterGlucose, i - 1),
						MineralMass = staticState.GetSliceLayer(nextState.WaterMinerals, i - 1),
						OxygenMass = staticState.GetSliceLayer(nextState.WaterOxygen, i - 1),
						Temperature = staticState.GetSliceLayer(nextState.WaterTemperature, i - 1),
						Velocity = staticState.GetSliceLayer(nextState.WaterVelocity, i - 1),

						Advection = tempState.RebalanceWater2,
					}, energyJobHandles[worldData.WaterLayer0]);

			}
		}
		tickJobHandle = JobHandle.CombineDependencies(tickJobHandle, energyJobHandles[worldData.WaterLayer0]);
#endregion

#region Update dependent variables

		tickJobHandle = tempState.Update(ref nextState, ref worldData, ref staticState, ref settings, tickJobHandle);

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
		var solarInJobHandle = SimJob.Schedule(
			JobType.Schedule, 64,
			new SolarRadiationJob()
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
				AlbedoSlopePower = worldData.AlbedoSlopePower
			}, lastJobHandle);

#endregion

		// Calculate emissivity per cell
		// TODO: combine cloud emissivity with cloud conduction
		// TODO: we can probably combine this step with the thermal radiation step
#region Emissivity Per Cell

		JobHandle emissivityJobHandle = default(JobHandle);
		emissivityJobHandle = JobHandle.CombineDependencies(emissivityJobHandle, AirJob.Schedule(
			JobType.Schedule, 64,
			new EmissivityAirJob()
			{
				Emissivity = staticState.GetSliceAir(tempState.EmissivityAir),
				AirMass = staticState.GetSliceAir(tempState.AirMass),
				VaporMass = staticState.GetSliceAir(lastState.AirVapor),
				Dust = staticState.GetSliceAir(lastState.AirDust),
				CarbonDioxide = staticState.GetSliceAir(lastState.AirCarbonDioxide),
				EmissivityAir = worldData.ThermalEmissivityAir,
				EmissivityWaterVapor = worldData.ThermalEmissivityWaterVapor,
				EmissivityDust = worldData.ThermalEmissivityDust,
				EmissivityCarbonDioxide = worldData.ThermalEmissivityCarbonDioxide
			}, lastJobHandle));

		// we only do thermal radiation upwards for the surface layer of water,
		// for the bottom we rely on conduction with the terrain for heat transfer (although this might lead to an imbalance!)
		{
			emissivityJobHandle = JobHandle.CombineDependencies(emissivityJobHandle, SimJob.Schedule(
				JobType.Schedule, 64,
				new EmissivityWaterJob()
				{
					Emissivity = staticState.GetSliceLayer(tempState.EmissivityWater,worldData.SurfaceWaterLayer),
					WaterMass = staticState.GetSliceLayer(lastState.WaterMass,worldData.SurfaceWaterLayer),
					SaltMass = staticState.GetSliceLayer(lastState.WaterSaltMass,worldData.SurfaceWaterLayer),
					EmissivitySalt = worldData.ThermalEmissivitySalt,
					EmissivityWater = worldData.ThermalEmissivityWater
				}, lastJobHandle));
		}
		emissivityJobHandle = JobHandle.CombineDependencies(emissivityJobHandle, SimJob.Schedule(
			JobType.Schedule, 64,
			new EmissivityTerrainJob()
			{
				Emissivity = tempState.EmissivityTerrain,
				SoilFertility = tempState.SoilFertility,
				EmissivityDirt = worldData.ThermalEmissivityDirt,
				EmissivitySand = worldData.ThermalEmissivitySand,
			}, lastJobHandle));

#endregion

		// Calculate how much thermal radition is being emitted out of each layer
#region Thermal Radiation
		JobHandle thermalOutJobHandle = default(JobHandle);

		// ICE
		thermalOutJobHandle = JobHandle.CombineDependencies(thermalOutJobHandle, SimJob.Schedule(
			settings.SynchronousOverrides.ThermalRadiation, 64,
			new ThermalEnergyRadiatedConstantEmissivityJob()
			{
				ThermalRadiationEmitted = tempState.ThermalRadiationEmittedIce,

				Emissivity = worldData.ThermalEmissivityIce,
				Energy = tempState.IceEnergy,
				Temperature = lastState.IceTemperature,
				SurfaceArea = tempState.IceCoverage,
				SecondsPerTick = worldData.SecondsPerTick
			}, lastJobHandle));


		// TERRAIN
		thermalOutJobHandle = JobHandle.CombineDependencies(thermalOutJobHandle, SimJob.Schedule(
			settings.SynchronousOverrides.ThermalRadiation, 64,
			new ThermalEnergyRadiatedTerrainJob()
			{
				ThermalRadiationEmitted = tempState.ThermalRadiationEmittedTerrain,
				
				Emissivity = tempState.EmissivityTerrain,
				Temperature = lastState.GroundTemperature,
				SecondsPerTick = worldData.SecondsPerTick
			}, emissivityJobHandle));


		// ATMOSPHERE
		thermalOutJobHandle = JobHandle.CombineDependencies(thermalOutJobHandle, AirJob.Schedule(
			settings.SynchronousOverrides.ThermalRadiation, 64,
			new ThermalEnergyRadiatedAirJob()
			{
				ThermalRadiationEmitted = staticState.GetSliceAir(tempState.ThermalRadiationEmittedAir),

				Energy = staticState.GetSliceAir(tempState.AirPotentialEnergy),
				Emissivity = staticState.GetSliceAir(tempState.EmissivityAir),
				TemperaturePotential = staticState.GetSliceAir(lastState.AirTemperaturePotential),
				LayerMiddle = staticState.GetSliceAir(tempState.AirLayerMiddle),
				SecondsPerTick = worldData.SecondsPerTick
			}, emissivityJobHandle));

		// WATER
		// we only do thermal radiation upwards for the surface layer of water,
		// for the bottom we rely on conduction with the terrain for heat transfer (although this might lead to an imbalance!)
		{
			thermalOutJobHandle = JobHandle.CombineDependencies(thermalOutJobHandle, WaterJob.Schedule(
				settings.SynchronousOverrides.ThermalRadiation, 64,
				new ThermalEnergyRadiatedWaterJob()
				{
					ThermalRadiationEmitted = staticState.GetSliceWater(tempState.ThermalRadiationEmittedWater),

					Emissivity = staticState.GetSliceWater(tempState.EmissivityWater),
					Energy = staticState.GetSliceWater(tempState.WaterPotentialEnergy),
					TemperatureAbsolute = staticState.GetSliceWater(lastState.WaterTemperature),
					SurfaceArea = staticState.GetSliceWater(tempState.WaterCoverage),
					SecondsPerTick = worldData.SecondsPerTick
				}, emissivityJobHandle));
		}
		#endregion


		#region absorptivity

		solarInJobHandle = SimJob.Schedule(
			settings.SynchronousOverrides.Albedo, 64,
			new AlbedoCloudJob()
			{
				Albedo = tempState.AlbedoCloud,
				Absorptivity = tempState.CloudAbsorptivity,
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
				AlbedoIceMin = worldData.AlbedoIceMin,
				AlbedoWaterMin = worldData.AlbedoWaterMin,
				AlbedoIceRange = worldData.AlbedoIceRange,
				AlbedoWaterRange = worldData.AlbedoWaterRange,
			}, solarInJobHandle);

		solarInJobHandle = AirJob.Schedule(
			settings.SynchronousOverrides.AirAbsorptivity, 64,
			new AbsorptivityAirJob()
			{
				AbsorptivitySolar = staticState.GetSliceAir(tempState.AbsorptivitySolar),
				AbsorptivityThermal = staticState.GetSliceAir(tempState.AbsorptivityThermal),
				AirMass = staticState.GetSliceAir(tempState.AirMass),
				VaporMass = staticState.GetSliceAir(lastState.AirVapor),
				AirCarbonDioxide = staticState.GetSliceAir(lastState.AirCarbonDioxide),
				Dust = staticState.GetSliceAir(lastState.AirDust),
				CloudMass = lastState.CloudMass,
				CloudAbsorptivitySolar = tempState.CloudAbsorptivity,
				CloudAlbedo = tempState.AlbedoCloud,
				CloudElevation = tempState.CloudElevation,
				LayerElevation = staticState.GetSliceAir(tempState.AirLayerElevation),
				LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
				Emissivity = staticState.GetSliceAir( tempState.EmissivityAir),
				EmissivityWater = worldData.ThermalEmissivityWater,
				AlbedoAir = worldData.AlbedoAir,
				AlbedoWaterVapor = worldData.AlbedoWaterVapor,
				AlbedoDust = worldData.AlbedoDust,
				SolarAbsorptivityAir = worldData.SolarAbsorptivityAir,
				SolarAbsorptivityWaterVapor = worldData.SolarAbsorptivityWaterVapor,
				SolarAbsorptivityDust = worldData.SolarAbsorptivityDust,
				Count = staticState.Count
			}, JobHandle.CombineDependencies(solarInJobHandle, emissivityJobHandle));

#endregion


		// Follow the solar radiation down from the top of the atmosphere to ther terrain, and absorb some as it passes through each layer
#region Solar Radiation Absorbed

		// process each vertical layer in order

		// atmosphere
		for (int j = worldData.AirLayer0 + worldData.AirLayers - 2; j > worldData.AirLayer0; j--)
		{
			int airLayerIndex = j - worldData.AirLayer0;
			solarInJobHandle = SimJob.Schedule(
				settings.SynchronousOverrides.SolarRadiationAbsorbed, 64,
				new SolarRadiationAbsorbedAirJob()
				{
					SolarRadiationAbsorbed = staticState.GetSliceLayer(tempState.SolarRadiationInAir,airLayerIndex),
					SolarRadiationIncoming = tempState.SolarRadiation,
					SolarRadiationReflected = staticState.GetSliceLayer(tempState.SolarReflectedAir,airLayerIndex),
					SolarRadiationAbsorbedCloud = tempState.SolarRadiationInCloud,
					SolarRadiationReflectedCloud = tempState.SolarReflectedCloud,
					AbsorptivitySolar = staticState.GetSliceLayer(tempState.AbsorptivitySolar,airLayerIndex),
				}, solarInJobHandle);
		}

		// ice
		solarInJobHandle = SimJob.Schedule(
			settings.SynchronousOverrides.SolarRadiationAbsorbed, 64,
			new SolarRadiationAbsorbedPartialCoverageConstantAlbedoJob()
			{
				SolarRadiationAbsorbed = tempState.SolarRadiationInIce,
				SolarRadiationReflected = tempState.SolarReflectedIce,
				SolarRadiationIncoming = tempState.SolarRadiation,
				AlbedoSlope = tempState.AlbedoSlope,
				AlbedoMin = worldData.AlbedoIceMin,
				AlbedoRange = worldData.AlbedoIceRange,
				Coverage = tempState.IceCoverage
			}, solarInJobHandle);

		// water
		// For data on albedo range due to solar zenith angle:
		// https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2004GL021180
		solarInJobHandle = SimJob.Schedule(
			settings.SynchronousOverrides.SolarRadiationAbsorbed, 64,
			new SolarRadiationAbsorbedPartialCoverageConstantAlbedoJob()
			{
				SolarRadiationAbsorbed = tempState.SolarRadiationInWater,
				SolarRadiationReflected = tempState.SolarReflectedWater,
				SolarRadiationIncoming = tempState.SolarRadiation,
				Coverage = staticState.GetSliceLayer(tempState.WaterCoverage,worldData.SurfaceWaterLayer),
				AlbedoSlope = tempState.AlbedoSlope,
				AlbedoMin = worldData.AlbedoWaterMin,
				AlbedoRange = worldData.AlbedoWaterRange,
			}, solarInJobHandle);


		// NOTE: we don't bother with solar radiation in lava

		solarInJobHandle = SimJob.Schedule(
			settings.SynchronousOverrides.SolarRadiationAbsorbed, 64,
			new SolarRadiationAbsorbedTerrainJob()
			{
				SolarRadiationAbsorbed = tempState.SolarRadiationInTerrain,
				SolarRadiationReflected = tempState.SolarReflectedTerrain,
				SolarRadiationIncoming = tempState.SolarRadiation,
				AlbedoSlope = tempState.AlbedoSlope,
				FloraCoverage = tempState.FloraCoverage,
				SoilFertility = tempState.SoilFertility,
				GroundWaterSaturation = tempState.GroundWaterSaturation,
				AlbedoSandMin = worldData.AlbedoSandMin,
				AlbedoSandRange = worldData.AlbedoSandRange,
				AlbedoSoilMin = worldData.AlbedoSoilMin,
				AlbedoSoilRange = worldData.AlbedoSoilRange,
				AlbedoFloraMin = worldData.AlbedoFloraMin,
				AlbedoFloraRange = worldData.AlbedoFloraRange,
				AlbedoReductionGroundWaterSaturation = worldData.AlbedoReductionGroundWaterSaturation,
			}, solarInJobHandle);
#endregion


		// Thermal radiation travels upwards, partially reflecting downwards (clouds), partially absorbed, and partially lost to space
#region Thermal Radiation Absorbed Up

		// transmit up from land
		for (int j = 0; j < worldData.LayerCount; j++)
		{
			if (j == worldData.TerrainLayer)
			{
				// TERRAIN
				thermalOutJobHandle = SimJob.Schedule(
					settings.SynchronousOverrides.ThermalRadiationAbsorbed, 64,
					new ThermalEnergyRadiatedUpTerrainJob()
					{
						ThermalRadiationDelta = tempState.ThermalRadiationDeltaTerrain,
						WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp,
						ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp,

						ThermalRadiationEmitted = tempState.ThermalRadiationEmittedTerrain,
					}, thermalOutJobHandle);
			}
			else if (j > worldData.WaterLayer0 && j <= worldData.SurfaceWaterLayerGlobal)
			{
				int waterLayerIndex = j - worldData.WaterLayer0;
				thermalOutJobHandle = SimJob.Schedule(
					settings.SynchronousOverrides.ThermalRadiationAbsorbed, 64,
					new ThermalEnergyAbsorbedUpPartialCoverageJob()
					{
						ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaWater, waterLayerIndex),
						ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp,
						WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp,

						ThermalRadiationEmitted = staticState.GetSliceLayer(tempState.ThermalRadiationEmittedWater, waterLayerIndex),
						PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
						Coverage = staticState.GetSliceLayer(tempState.WaterCoverage, waterLayerIndex),
					}, thermalOutJobHandle);
			}
			else if (j == worldData.IceLayer)
			{
				thermalOutJobHandle = SimJob.Schedule(
					settings.SynchronousOverrides.ThermalRadiationAbsorbed, 64,
					new ThermalEnergyAbsorbedUpPartialCoverageJob()
					{
						ThermalRadiationDelta = tempState.ThermalRadiationDeltaIce,
						ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp,
						WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp,

						ThermalRadiationEmitted = tempState.ThermalRadiationEmittedIce,
						PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
						Coverage = tempState.IceCoverage,

					}, thermalOutJobHandle);
			}
			else if (j > worldData.AirLayer0 && j < worldData.AirLayer0 + worldData.AirLayers - 1)
			{
				int airLayer = j - worldData.AirLayer0;
				int downIndex = airLayer == 1 ? worldData.IceLayer : (j - 1);

				thermalOutJobHandle = SimJob.Schedule(
					settings.SynchronousOverrides.ThermalRadiationAbsorbed, 64,
					new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaAir, airLayer),
						ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp,
						WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp,

						ThermalRadiationEmitted = staticState.GetSliceLayer(tempState.ThermalRadiationEmittedAir, airLayer),
						AbsorptivityThermal = staticState.GetSliceLayer(tempState.AbsorptivityThermal, airLayer),
						LayerElevation = staticState.GetSliceLayer(tempState.AirLayerElevation, airLayer),
						LayerHeight = staticState.GetSliceLayer(tempState.AirLayerHeight, airLayer),
						CloudElevation = tempState.CloudElevation,
						PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
						FromTop = false,
					}, JobHandle.CombineDependencies(solarInJobHandle, thermalOutJobHandle));
			}
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
				thermalOutJobHandle = SimJob.Schedule(
					settings.SynchronousOverrides.ThermalRadiationAbsorbed, 64,
					new ThermalEnergyAbsorbedDownTerrainJob()
					{
						ThermalRadiationDelta = tempState.ThermalRadiationDeltaTerrain,

						WindowRadiationIncoming = tempState.WindowRadiationTransmittedDown,
						ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedDown,
					}, thermalOutJobHandle);
			}
			else if (j == worldData.IceLayer)
			{
				// ICE
				thermalOutJobHandle = SimJob.Schedule(
					settings.SynchronousOverrides.ThermalRadiationAbsorbed, 64,
					new ThermalEnergyAbsorbedDownPartialCoverageJob()
					{
						ThermalRadiationDelta = tempState.ThermalRadiationDeltaIce,
						ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedDown,
						WindowRadiationTransmitted = tempState.WindowRadiationTransmittedDown,

						Coverage = tempState.IceCoverage,
					}, thermalOutJobHandle);
			}
			else if (j == worldData.SurfaceWaterLayerGlobal)
			{
				// WATER
				int waterLayerIndex = j - worldData.WaterLayer0;
				thermalOutJobHandle = SimJob.Schedule(
					settings.SynchronousOverrides.ThermalRadiationAbsorbed, 64,
					new ThermalEnergyAbsorbedDownPartialCoverageJob()
					{
						ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaWater, waterLayerIndex),
						ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedDown,
						WindowRadiationTransmitted = tempState.WindowRadiationTransmittedDown,

						Coverage = staticState.GetSliceLayer(tempState.WaterCoverage,waterLayerIndex),
					}, thermalOutJobHandle);
			}
			else if (j > worldData.AirLayer0 && j < worldData.AirLayer0 + worldData.AirLayers - 1)
			{
				int airLayer = j - worldData.AirLayer0;
				thermalOutJobHandle = SimJob.Schedule(
					settings.SynchronousOverrides.ThermalRadiationAbsorbed, 64,
					new ThermalEnergyAbsorbedAirJob()
					{
						ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaAir,airLayer),
						ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedDown,
						WindowRadiationTransmitted = tempState.WindowRadiationTransmittedDown,

						ThermalRadiationEmitted = staticState.GetSliceLayer(tempState.ThermalRadiationEmittedAir, airLayer),
						PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
						AbsorptivityThermal = staticState.GetSliceLayer(tempState.AbsorptivityThermal,airLayer),
						LayerElevation = staticState.GetSliceLayer(tempState.AirLayerElevation,airLayer),
						LayerHeight = staticState.GetSliceLayer(tempState.AirLayerHeight,airLayer),
						CloudElevation = tempState.CloudElevation,
						FromTop = true,
					}, JobHandle.CombineDependencies(thermalOutJobHandle, solarInJobHandle));
			}
		}
#endregion

		// Conduction is calculated for each Surface that might touch another surface
		// Air to Cloud, Air to Ice, Air to Water, Air to Terrain, Ice to Water, Ice to Terrain, Water to Terrain
#region Conduction

		JobHandle conductionAirIceJobHandle;
		JobHandle conductionAirWaterJobHandle;
		JobHandle conductionAirTerrainJobHandle;
		JobHandle conductionIceWaterJobHandle;
		JobHandle conductionIceTerrainJobHandle;
		JobHandle conductionWaterTerrainJobHandle;

		// air to ice
		conductionAirIceJobHandle = SimJob.ScheduleOrMemset(
			settings.SynchronousOverrides.Conduction, 64,
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
			settings.SynchronousOverrides.Conduction, 64,
			settings.ConductionAirWater,
			tempState.ConductionAirWater,
			0,
			new ConductionBJob()
			{
				EnergyDelta = tempState.ConductionAirWater,
				TemperatureA = tempState.SurfaceAirTemperatureAbsolute,
				TemperatureB = staticState.GetSliceLayer(lastState.WaterTemperature,worldData.SurfaceWaterLayer),
				EnergyB = staticState.GetSliceLayer(tempState.WaterPotentialEnergy,worldData.SurfaceWaterLayer),
				ConductionCoefficient = WorldData.ConductivityAirWater,
				SurfaceArea = tempState.SurfaceAreaAirWater,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);

		// air to terrain
		conductionAirTerrainJobHandle = SimJob.ScheduleOrMemset(
			settings.SynchronousOverrides.Conduction, 64,
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
			settings.SynchronousOverrides.Conduction, 64,
			settings.ConductionIceWater,
			tempState.ConductionIceWater,
			0,
			new ConductionABJob()
			{
				EnergyDelta = tempState.ConductionIceWater,
				TemperatureA = lastState.IceTemperature,
				TemperatureB = staticState.GetSliceLayer(lastState.WaterTemperature,worldData.SurfaceWaterLayer),
				EnergyA = tempState.IceEnergy,
				EnergyB = staticState.GetSliceLayer(tempState.WaterPotentialEnergy,worldData.SurfaceWaterLayer),
				ConductionCoefficient = WorldData.ConductivityIceWater,
				SurfaceArea = tempState.SurfaceAreaIceWater,
				SecondsPerTick = worldData.SecondsPerTick
			}, 
			lastJobHandle);


		// ice to terrain
		conductionIceTerrainJobHandle = SimJob.ScheduleOrMemset(
			settings.SynchronousOverrides.Conduction, 64,
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


		// water to terrain
		conductionWaterTerrainJobHandle = lastJobHandle;
		if (!settings.ConductionWaterTerrain)
		{
			conductionWaterTerrainJobHandle = Utils.MemsetArray(staticState.Count, lastJobHandle, tempState.ConductionWaterTerrain, 0);
		}
		else
		{
			for (int i = 1; i < worldData.WaterLayers - 1; i++)
			{
				conductionWaterTerrainJobHandle = SimJob.Schedule(
					settings.SynchronousOverrides.Conduction, 64,
					new ConductionWaterBottomAJob()
					{
						EnergyDelta = tempState.ConductionWaterTerrain,
						TemperatureA = staticState.GetSliceLayer(lastState.WaterTemperature, i),
						TemperatureB = lastState.GroundTemperature,
						EnergyA = staticState.GetSliceLayer(tempState.WaterPotentialEnergy, i),
						ConductionCoefficient = WorldData.ConductivityWaterTerrain,
						SurfaceArea = tempState.SurfaceAreaWaterTerrain,
						Coverage = tempState.WaterCoverage,
						SecondsPerTick = worldData.SecondsPerTick,
						LayerIndex = i,
						Count = staticState.Count
					},
					conductionWaterTerrainJobHandle);
			}
		}

#endregion

#region Change temperature due to energy flux

		var terrainEnergyJobHandleDependencies = new NativeList<JobHandle>(Allocator.Persistent)
			{
				solarInJobHandle,
				thermalOutJobHandle,
				conductionAirTerrainJobHandle,
				conductionIceTerrainJobHandle,
				conductionWaterTerrainJobHandle,
			};
		jobHandleDependencies.Add(terrainEnergyJobHandleDependencies);
		energyJobHandles[worldData.TerrainLayer] = SimJob.Schedule(
			settings.SynchronousOverrides.Energy, 64,
			new EnergyTerrainJob()
			{
				TerrainTemperature = nextState.GroundTemperature,
				LastTemperature = lastState.GroundTemperature,
				SpecificHeatTerrain = tempState.SpecificHeatTerrain,
				SolarRadiationIn = tempState.SolarRadiationInTerrain,
				ThermalRadiationDelta = tempState.ThermalRadiationDeltaTerrain,
				ConductionEnergyAir = tempState.ConductionAirTerrain,
				ConductionEnergyIce = tempState.ConductionIceTerrain,
				ConductionEnergyWater = tempState.ConductionWaterTerrain,
				GeothermalEnergy = tempState.GeothermalRadiation,
			}, JobHandle.CombineDependencies(terrainEnergyJobHandleDependencies));

		var energyIceJobHandleDependencies = new NativeList<JobHandle>(Allocator.Persistent)
			{
				solarInJobHandle,
				thermalOutJobHandle,
				conductionAirIceJobHandle,
				conductionIceWaterJobHandle,
				conductionIceTerrainJobHandle,
			};
		jobHandleDependencies.Add(energyIceJobHandleDependencies);
		energyJobHandles[worldData.IceLayer] = SimJob.Schedule(
			settings.SynchronousOverrides.Energy, 64,
			new EnergyIceJob()
			{
				Temperature = nextState.IceTemperature,
				LastTemperature = lastState.IceTemperature,
				LastMass = lastState.IceMass,
				SolarRadiationIn = tempState.SolarRadiationInIce,
				ThermalRadiationDelta = tempState.ThermalRadiationDeltaIce,
				ConductionEnergyAir = tempState.ConductionAirIce,
				ConductionEnergyTerrain = tempState.ConductionIceTerrain,
				ConductionEnergyWater = tempState.ConductionIceWater,
			}, JobHandle.CombineDependencies(energyIceJobHandleDependencies));


		energyJobHandles[worldData.LavaLayer] = SimJob.Schedule(
			settings.SynchronousOverrides.Energy, 64,
			new EnergyLavaJob()
			{
				LavaTemperature = nextState.LavaTemperature,
				LastTemperature = lastState.LavaTemperature,
				LavaMass = lastState.LavaMass,
				Emissivity = worldData.ThermalEmissivityLava,
				SecondsPerTick = worldData.SecondsPerTick
			});

		energyJobHandles[worldData.AirLayer0] = UpperAirJob.Schedule(
			settings.SynchronousOverrides.Energy, 64,
			new EnergyAirJob()
			{
				AirTemperaturePotential = staticState.GetSliceLayers(nextState.AirTemperaturePotential, 2, worldData.AirLayers-3),
				LastTemperaturePotential = staticState.GetSliceLayers(lastState.AirTemperaturePotential, 2, worldData.AirLayers - 3),
				LastVapor = staticState.GetSliceLayers(lastState.AirVapor, 2, worldData.AirLayers - 3),
				AirMass = staticState.GetSliceLayers(tempState.AirMass, 2, worldData.AirLayers - 3),
				SolarRadiationIn = staticState.GetSliceLayers(tempState.SolarRadiationInAir, 2, worldData.AirLayers - 3),
				ThermalRadiationDelta = staticState.GetSliceLayers(tempState.ThermalRadiationDeltaAir, 2, worldData.AirLayers - 3),
				CloudMass = lastState.CloudMass,
				CloudElevation = tempState.CloudElevation,
				LayerElevation = staticState.GetSliceLayers(tempState.AirLayerElevation, 2, worldData.AirLayers - 3),
				LayerHeight = staticState.GetSliceLayers(tempState.AirLayerHeight, 2, worldData.AirLayers - 3),
				Count = staticState.Count
			}, JobHandle.CombineDependencies(solarInJobHandle, thermalOutJobHandle));

		var airSurfaceDependencies = new NativeList<JobHandle>(Allocator.Persistent)
		{
			energyJobHandles[worldData.AirLayer0],
			solarInJobHandle,
			thermalOutJobHandle,
			conductionAirWaterJobHandle,
			conductionAirIceJobHandle,
			conductionAirTerrainJobHandle,
		};
		jobHandleDependencies.Add(airSurfaceDependencies);
		energyJobHandles[worldData.AirLayer0] = SimJob.Schedule(
			settings.SynchronousOverrides.Energy, 64,
			new EnergyAirSurfaceJob()
			{
				AirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential,worldData.SurfaceAirLayer),
				LastTemperaturePotential = staticState.GetSliceLayer(lastState.AirTemperaturePotential, worldData.SurfaceAirLayer),
				LastVapor = staticState.GetSliceLayer(lastState.AirVapor, worldData.SurfaceAirLayer),
				AirMass = staticState.GetSliceLayer(tempState.AirMass, worldData.SurfaceAirLayer),
				SolarRadiationIn = staticState.GetSliceLayer(tempState.SolarRadiationInAir, worldData.SurfaceAirLayer),
				ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaAir, worldData.SurfaceAirLayer),
				ConductionEnergyWater = tempState.ConductionAirWater,
				ConductionEnergyIce = tempState.ConductionAirIce,
				ConductionEnergyTerrain = tempState.ConductionAirTerrain,
				CloudMass = lastState.CloudMass,
				CloudElevation = tempState.CloudElevation,
				LayerElevation = staticState.GetSliceAir(tempState.AirLayerElevation),
				LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
				Count = staticState.Count
			}, JobHandle.CombineDependencies(airSurfaceDependencies));

		var waterDependencies = new NativeList<JobHandle>(Allocator.Persistent)
				{
					solarInJobHandle,
					thermalOutJobHandle,
					conductionAirWaterJobHandle,
					conductionIceWaterJobHandle,
					conductionWaterTerrainJobHandle,
				};
		jobHandleDependencies.Add(waterDependencies);

		energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
			settings.SynchronousOverrides.Energy, 64,
			new EnergyWaterJob()
			{
				Temperature = staticState.GetSliceWater(nextState.WaterTemperature),
				LastMass = lastState.WaterMass,
				LastSaltMass = lastState.WaterSaltMass,
				LastTemperature = lastState.WaterTemperature,
				ThermalRadiationDelta = tempState.ThermalRadiationDeltaWater,
				Coverage = tempState.WaterCoverage,
				ConductionEnergyAir = tempState.ConductionAirWater,
				ConductionEnergyIce = tempState.ConductionIceWater,
				ConductionEnergyTerrain = tempState.ConductionWaterTerrain,
				Count = staticState.Count
			}, JobHandle.CombineDependencies(waterDependencies));

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
			energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
				settings.SynchronousOverrides.FluxCondensation, 64, new FluxCondensationJob()
				{
					LatentHeatCloud = staticState.GetSliceAir(tempState.LatentHeatCloud),
					CondensationCloudMass = staticState.GetSliceAir(tempState.CondensationCloudMass),
					CondensationGroundMass = staticState.GetSliceAir(tempState.CondensationGroundMass),

					TemperaturePotential = staticState.GetSliceAir(nextState.AirTemperaturePotential),
					LastVapor = staticState.GetSliceAir(lastState.AirVapor),
					AirMass = staticState.GetSliceAir(tempState.AirMass),
					AirPressure = staticState.GetSliceAir(tempState.AirPressure),
					LayerElevation = staticState.GetSliceAir(tempState.AirLayerElevation),
					LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
					LayerMiddle = staticState.GetSliceAir(tempState.AirLayerMiddle),
					CloudElevation = tempState.CloudElevation,
				}, energyJobHandles[worldData.AirLayer0]);
		}
		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
			settings.SynchronousOverrides.FluxDust, 64,
			new FluxDustJob()
			{
				DustUp = staticState.GetSliceAir(tempState.DustUp),
				DustDown = staticState.GetSliceAir(tempState.DustDown),

				LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
				LastDust = staticState.GetSliceAir(lastState.AirDust),
				AirVelocity = staticState.GetSliceAir(lastState.AirVelocity),
				DustVerticalVelocity = worldData.DustVerticalVelocity,
				Positions = staticState.SphericalPosition,
				SecondsPerTick = worldData.SecondsPerTick,
				Count = staticState.Count
			}, energyJobHandles[worldData.AirLayer0]);

		if (settings.Evaporation)
		{
			energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
				settings.SynchronousOverrides.FluxEvaporation, 64,
				new FluxEvaporationJob()
				{
					EvaporatedWaterMass = tempState.EvaporationMassWater,
					LatentHeatWater = tempState.LatentHeatWaterSurface,
					LatentHeatAir = staticState.GetSliceLayer(tempState.LatentHeatAir,worldData.SurfaceAirLayer),
					LatentHeatTerrain = tempState.LatentHeatTerrain,

					WaterTemperature = staticState.GetSliceLayer(nextState.WaterTemperature,worldData.SurfaceWaterLayer),
					WaterMass = staticState.GetSliceLayer(lastState.WaterMass,worldData.SurfaceWaterLayer),
					IceCoverage = tempState.IceCoverage,
					WaterCoverage = staticState.GetSliceLayer(tempState.WaterCoverage,worldData.SurfaceWaterLayer),
					SurfaceWind = staticState.GetSliceLayer(lastState.AirVelocity,worldData.SurfaceAirLayer),
					SurfaceElevation = tempState.SurfaceElevation,
					AirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential, worldData.SurfaceAirLayer),
					AirMass = staticState.GetSliceLayer(tempState.AirMass,worldData.SurfaceAirLayer),
					AirPressure = staticState.GetSliceLayer(tempState.AirPressure,worldData.SurfaceAirLayer),
					AirVapor = staticState.GetSliceLayer(lastState.AirVapor,worldData.SurfaceAirLayer),
					WaterHeatingDepth = worldData.WaterHeatingDepth,
					EvaporationLatentHeatFromAir = worldData.EvaporationLatentHeatFromAir,
				}, JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0], energyJobHandles[worldData.AirLayer0]));
		}

		if (settings.Freezing)
		{
			energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
				settings.SynchronousOverrides.FluxFreeze, 64,
				new FluxFreezeJob()
				{
					FrozenMass = tempState.FrozenMass,
					FrozenTemperature = tempState.FrozenTemperature,
					LatentHeatWater = tempState.LatentHeatWaterSurface,
					SaltPlume = tempState.SaltPlume,

					WaterTemperature = staticState.GetSliceLayer(nextState.WaterTemperature,worldData.SurfaceWaterLayer),
					WaterMass = staticState.GetSliceLayer(lastState.WaterMass,worldData.SurfaceWaterLayer),
					SaltMass = staticState.GetSliceLayer(lastState.WaterSaltMass,worldData.SurfaceWaterLayer),
					AirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential,worldData.SurfaceAirLayer),
					AirLayerElevation = staticState.GetSliceLayer(tempState.AirLayerElevation,worldData.SurfaceAirLayer),
					WaterHeatingDepth = worldData.WaterHeatingDepth,
					FreezePointReductionPerSalinity = worldData.FreezePointReductionPerSalinity,
				}, JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0], energyJobHandles[worldData.AirLayer0]));
		}


		// CLOUD
		if (settings.Precipitation)
		{
			energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(
				settings.SynchronousOverrides.FluxCloud, 64,
				new FluxCloudJob()
				{
					EvaporationMass = tempState.CloudEvaporationMass,
					PrecipitationMass = tempState.PrecipitationMass,
					PrecipitationTemperature = tempState.PrecipitationTemperature,
					DropletDelta = tempState.DropletDelta,

					SurfaceAirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential,worldData.SurfaceAirLayer),
					SurfaceLayerElevation = staticState.GetSliceLayer(tempState.AirLayerElevation,worldData.SurfaceAirLayer),
					SurfaceLayerMiddle = staticState.GetSliceLayer(tempState.AirLayerMiddle,worldData.SurfaceAirLayer),
					SurfaceSaltMass = staticState.GetSliceLayer(lastState.WaterSaltMass,worldData.SurfaceWaterLayer),
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
				}, JobHandle.CombineDependencies(energyJobHandles[worldData.CloudLayer], energyJobHandles[worldData.AirLayer0]));
		}


		if (settings.IceMelting)
		{
			energyJobHandles[worldData.IceLayer] = SimJob.Schedule(
				settings.SynchronousOverrides.FluxIceMelt, 64,
				new FluxIceMeltJob()
				{
					LatentHeatAir = staticState.GetSliceLayer(tempState.LatentHeatAir, worldData.SurfaceAirLayer),
					LatentHeatWater = tempState.LatentHeatWaterSurface,
					LatentHeatTerrain = tempState.LatentHeatTerrain,
					LatentHeatIce = tempState.LatentHeatIce,
					MeltedMass = tempState.IceMeltedMass,

					Temperature = nextState.IceTemperature,
					LastMass = lastState.IceMass,
					IceHeatingDepth = worldData.IceHeatingDepth,
					WaterIceSurfaceArea = tempState.SurfaceAreaIceWater,
					WaterTerrainSurfaceArea = tempState.SurfaceAreaIceTerrain,
					WaterTemperature = staticState.GetSliceLayer(nextState.WaterTemperature,worldData.SurfaceWaterLayer),
					TerrainTemperature = nextState.GroundTemperature,
					SurfaceElevation = tempState.SurfaceElevation,
					AirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential, worldData.SurfaceAirLayer),

				}, JobHandle.CombineDependencies(JobHandle.CombineDependencies(energyJobHandles[worldData.FloraLayer], energyJobHandles[worldData.IceLayer]), energyJobHandles[worldData.TerrainLayer], energyJobHandles[worldData.WaterLayer0]));
		}

		energyJobHandles[worldData.LavaLayer] = SimJob.Schedule(
			settings.SynchronousOverrides.FluxLava, 64,
			new FluxLavaJob()
			{
			}, energyJobHandles[worldData.LavaLayer]);

		if (settings.SoilRespiration)
		{
			energyJobHandles[worldData.TerrainLayer] = SimJob.Schedule(
				settings.SynchronousOverrides.FluxTerrain, 64,
				new FluxTerrainJob()
				{
					CrystalizedMass = tempState.LavaCrystalizedMass,
					LavaEjected = tempState.LavaEjected,
					DustEjected = tempState.DustEjected,
					CrustDelta = tempState.CrustDelta,
					LatentHeatLava = tempState.LatentHeatLava,

					LavaTemperature = nextState.LavaTemperature,
					LavaMass = lastState.LavaMass,
					CrustDepth = lastState.CrustDepth,
					MagmaMass = lastState.MagmaMass,
					Elevation = lastState.Elevation,
					WaterCoverage = staticState.GetSliceLayer(tempState.WaterCoverage,worldData.SurfaceWaterLayer),
					LavaCrystalizationTemperature = worldData.LavaCrystalizationTemperature,
					CrustEruptionDepth = worldData.CrustDepthForEruption,
					DustPerLavaEjected = worldData.DustPerLavaEjected,
					MagmaPressureCrustReductionSpeed = worldData.MagmaPressureCrustReductionSpeed,
					LavaEruptionSpeed = worldData.LavaEruptionSpeed,
					SecondsPerTick = worldData.SecondsPerTick,
					SoilRespirationSpeed = worldData.SoilRespirationSpeed,
				}, JobHandle.CombineDependencies(energyJobHandles[worldData.TerrainLayer], energyJobHandles[worldData.LavaLayer]));
		}
#endregion

#region Update Mass - Evaporation, Condensation, Melting, Rainfall

		JobHandle waterDependencies = JobHandle.CombineDependencies(
			energyJobHandles[worldData.IceLayer], 
			energyJobHandles[worldData.TerrainLayer],
			energyJobHandles[worldData.WaterLayer0]);
		energyJobHandles[worldData.WaterLayer0] = WaterJob.Schedule(
			JobType.Schedule, 64,
			new UpdateMassWaterJob()
			{
				WaterMass = staticState.GetSliceWater(nextState.WaterMass),
				SaltMass = staticState.GetSliceWater(nextState.WaterSaltMass),
				CarbonMass = staticState.GetSliceWater(nextState.WaterCarbonDioxide),
				OxygenMass = staticState.GetSliceWater(nextState.WaterOxygen),
				NitrogenMass = staticState.GetSliceWater(nextState.WaterNitrogen),
				GlucoseMass = staticState.GetSliceWater(nextState.WaterGlucose),
				MineralsMass = staticState.GetSliceWater(nextState.WaterMinerals),
				WaterTemperature = staticState.GetSliceWater(nextState.WaterTemperature),
				SaltPlume = tempState.SaltPlume,
				SaltPlumeTemperature = tempState.FrozenTemperature,
				LastSaltMass = lastState.WaterSaltMass,
				LastCarbonMass = lastState.WaterCarbonDioxide,
				LastOxygenMass = lastState.WaterOxygen,
				LastNitrogenMass = lastState.WaterNitrogen,
				LastGlucoseMass = lastState.WaterGlucose,
				LastMineralsMass = lastState.WaterMinerals,
				LastWaterMass = lastState.WaterMass,
				WaterCoverage = tempState.WaterCoverage,
				Count = staticState.Count
			}, JobHandle.CombineDependencies( energyJobHandles[worldData.WaterLayer0], waterDependencies));

		energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
			JobType.Schedule, 64,
			new UpdateMassCondensationGroundJob()
			{
				SurfaceWaterMass = staticState.GetSliceLayer(nextState.WaterMass,worldData.SurfaceWaterLayer),
				SurfaceWaterTemperature = staticState.GetSliceLayer(nextState.WaterTemperature,worldData.SurfaceWaterLayer),

				SurfaceSaltMass = staticState.GetSliceLayer(lastState.WaterSaltMass,worldData.SurfaceWaterLayer),
				AirTemperaturePotential = staticState.GetSliceAir(nextState.AirTemperaturePotential),
				GroundCondensation = staticState.GetSliceAir(tempState.CondensationGroundMass),
				LayerMiddle = staticState.GetSliceAir(tempState.AirLayerMiddle),
				LayerCount = worldData.AirLayers - 2,
				Count = staticState.Count
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0], energyJobHandles[worldData.AirLayer0]));

		energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
			JobType.Schedule, 64,
			new UpdateMassWaterSurfaceJob()
			{
				WaterTemperature = staticState.GetSliceLayer(nextState.WaterTemperature,worldData.SurfaceWaterLayer),
				WaterMass = staticState.GetSliceLayer(nextState.WaterMass,worldData.SurfaceWaterLayer),
				SaltMass = staticState.GetSliceLayer(nextState.WaterSaltMass,worldData.SurfaceWaterLayer),
				CarbonMass = staticState.GetSliceLayer(nextState.WaterCarbonDioxide,worldData.SurfaceWaterLayer),
				OxygenMass = staticState.GetSliceLayer(nextState.WaterOxygen, worldData.SurfaceWaterLayer),
				NitrogenMass = staticState.GetSliceLayer(nextState.WaterNitrogen, worldData.SurfaceWaterLayer),
				GlucoseMass = staticState.GetSliceLayer(nextState.WaterGlucose, worldData.SurfaceWaterLayer),
				MineralsMass = staticState.GetSliceLayer(nextState.WaterMinerals, worldData.SurfaceWaterLayer),

				SaltPlume = tempState.SaltPlume,
				Evaporation = tempState.EvaporationMassWater,
				IceMelted = tempState.IceMeltedMass,
				Precipitation = tempState.PrecipitationMass,
				PrecipitationTemperature = tempState.PrecipitationTemperature,
				FloraRespirationWater = tempState.FloraRespirationMassWater,
				TerrainTemperature = nextState.GroundTemperature,
				WaterFrozen = tempState.FrozenMass,
				WaterCarbonDelta = tempState.WaterCarbonDelta,
				WaterOxygenDelta = tempState.WaterOxygenDelta,
				WaterNitrogenDelta = tempState.WaterNitrogenDelta,
				WaterGlucoseDelta = tempState.WaterGlucoseDelta,
				WaterMineralsDelta = tempState.WaterMineralsDelta
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0], energyJobHandles[worldData.CloudLayer]));

		energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(
			settings.SynchronousOverrides.UpdateMassCloud, 64,
			new UpdateMassCloudJob()
			{
				CloudMass = nextState.CloudMass,
				CloudDropletMass = nextState.CloudDropletMass,

				LastCloudMass = lastState.CloudMass,
				LastDropletMass = lastState.CloudDropletMass,
				CloudEvaporation = tempState.CloudEvaporationMass,
				PrecipitationMass = tempState.PrecipitationMass,
				DropletDelta = tempState.DropletDelta,
				CloudCondensation = staticState.GetSliceAir(tempState.CondensationCloudMass),
				LayerCount = worldData.AirLayers - 2,
				Count = staticState.Count
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.CloudLayer], energyJobHandles[worldData.AirLayer0]));

		JobHandle airDependencies = JobHandle.CombineDependencies(JobHandle.CombineDependencies(
			energyJobHandles[worldData.AirLayer0],
			energyJobHandles[worldData.CloudLayer],
			energyJobHandles[worldData.IceLayer]),
			energyJobHandles[worldData.TerrainLayer],
			energyJobHandles[worldData.WaterLayer0]);
		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
			settings.SynchronousOverrides.UpdateMassAir, 64,
			new UpdateMassAirJob()
			{
				VaporMass = staticState.GetSliceAir(nextState.AirVapor),
				DustMass = staticState.GetSliceAir(nextState.AirDust),
				CarbonDioxideMass = staticState.GetSliceAir(nextState.AirCarbonDioxide),
				TemperaturePotential = staticState.GetSliceAir(nextState.AirTemperaturePotential),
				
				CloudCondensation = staticState.GetSliceAir(tempState.CondensationCloudMass),
				GroundCondensation = staticState.GetSliceAir(tempState.CondensationGroundMass),
				LastVaporMass = staticState.GetSliceAir(lastState.AirVapor),
				LastDustMass = staticState.GetSliceAir(lastState.AirDust),
				LastCarbonDioxideMass = staticState.GetSliceAir(lastState.AirCarbonDioxide),
				AirMass = staticState.GetSliceAir(tempState.AirMass),
				LayerElevation = staticState.GetSliceAir(tempState.AirLayerElevation),
				LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
				CloudElevation = tempState.CloudElevation,
				CloudMass = lastState.CloudMass,
				DustUp = tempState.DustUp,
				DustDown = tempState.DustDown,
				LayerCount = worldData.AirLayers - 2,
				Count = staticState.Count

			}, airDependencies);

		energyJobHandles[worldData.AirLayer0] = SimJob.Schedule(
			settings.SynchronousOverrides.UpdateMassAir, 64,
			new UpdateMassAirSurfaceJob()
			{
				AirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential, worldData.SurfaceAirLayer),
				VaporMass = staticState.GetSliceLayer(nextState.AirVapor,worldData.SurfaceAirLayer),
				DustMass = staticState.GetSliceLayer(nextState.AirDust,worldData.SurfaceAirLayer),
				CarbonDioxide = staticState.GetSliceLayer(nextState.AirCarbonDioxide,worldData.SurfaceAirLayer),
				Oxygen = staticState.GetSliceLayer(nextState.AirOxygen, worldData.SurfaceAirLayer),
				Nitrogen = staticState.GetSliceLayer(nextState.AirNitrogen, worldData.SurfaceAirLayer),
				Methane = staticState.GetSliceLayer(nextState.AirMethane, worldData.SurfaceAirLayer),
				Minerals = staticState.GetSliceLayer(nextState.AirMinerals, worldData.SurfaceAirLayer),

				AirMass = staticState.GetSliceLayer(tempState.AirMass, worldData.SurfaceAirLayer),
				EvaporationWater = tempState.EvaporationMassWater,
				EvaporationTemperatureWater = staticState.GetSliceLayer(nextState.WaterTemperature,worldData.SurfaceWaterLayer),
				EvaporationFlora = tempState.FloraRespirationMassVapor,
				EvaporationTemperatureFlora = nextState.GroundTemperature,
				DustEjected = tempState.DustEjected,
				CarbonDioxideDelta = tempState.AirCarbonDelta,
				OxygenDelta = tempState.AirOxygenDelta,
				NitrogenDelta = tempState.AirNitrogenDelta,
				MethaneDelta = tempState.AirMethaneDelta,
				MineralsDelta = tempState.AirMineralsDelta,
				WaterCoverage = staticState.GetSliceLayer(tempState.WaterCoverage,worldData.SurfaceWaterLayer),
				Elevation = lastState.Elevation,
				CloudElevation = tempState.CloudElevation,
				CloudMass = nextState.CloudMass,
				LayerElevation = staticState.GetSliceAir(tempState.AirLayerElevation),
				LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
				Count = staticState.Count
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.AirLayer0], energyJobHandles[worldData.WaterLayer0], energyJobHandles[worldData.CloudLayer]));

		energyJobHandles[worldData.IceLayer] = SimJob.Schedule(
			JobType.Schedule, 64,
			new UpdateMassIceJob()
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

		energyJobHandles[worldData.TerrainLayer] = SimJob.Schedule(
			JobType.Schedule, 64,
			new UpdateTerrainJob()
			{
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
				LastGroundWater = lastState.GroundWater,
				GroundWaterConsumed = tempState.WaterConsumedByFlora,
				WaterCoverage = staticState.GetSliceLayer(tempState.WaterCoverage,worldData.SurfaceWaterLayer),
				LastCrustDepth = lastState.CrustDepth,
				LastLavaMass = lastState.LavaMass,
				LastMagmaMass = lastState.MagmaMass,
				DustSettled = staticState.GetSliceLayer(tempState.DustDown,worldData.SurfaceAirLayer),
				LavaCrystalized = tempState.LavaCrystalizedMass,
				LavaEjected = tempState.LavaEjected,
				MagmaTemperature = worldData.MagmaTemperature,
				LavaToRockMassAdjustment = worldData.LavaToRockMassAdjustment,
			}, JobHandle.CombineDependencies(JobHandle.CombineDependencies(
				energyJobHandles[worldData.TerrainLayer], energyJobHandles[worldData.LavaLayer]), 
				energyJobHandles[worldData.FloraLayer], energyJobHandles[worldData.WaterLayer0]));

		if (settings.AirWaterCarbonDioxideDiffusion)
		{
			energyJobHandles[worldData.AirLayer0] = SimJob.Schedule(
				JobType.Schedule, 64,
				new UpdateWaterAirDiffusionJob()
				{
					AirCarbon = staticState.GetSliceLayer(nextState.AirCarbonDioxide, worldData.SurfaceAirLayer),
					WaterCarbon = staticState.GetSliceLayer(nextState.WaterCarbonDioxide, worldData.SurfaceWaterLayer),
					AirOxygen = staticState.GetSliceLayer(nextState.AirOxygen, worldData.SurfaceWaterLayer),
					WaterOxygen = staticState.GetSliceLayer(nextState.WaterOxygen, worldData.SurfaceWaterLayer),
									   
					AirMass = staticState.GetSliceLayer(tempState.AirMass, worldData.SurfaceAirLayer),
					WaterMass = staticState.GetSliceLayer(nextState.WaterMass, worldData.SurfaceWaterLayer),
					SaltMass = staticState.GetSliceLayer(nextState.WaterSaltMass, worldData.SurfaceWaterLayer),
					WaterDepth = staticState.GetSliceLayer(tempState.WaterLayerHeight, worldData.SurfaceWaterLayer),
					WaterAirDiffusionCoefficient = worldData.WaterAirCarbonDiffusionCoefficient,
					WaterAirCarbonDepth = worldData.WaterAirCarbonDiffusionDepth,
				}, JobHandle.CombineDependencies(energyJobHandles[worldData.AirLayer0], energyJobHandles[worldData.WaterLayer0]));
			energyJobHandles[worldData.WaterLayer0] = JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0], energyJobHandles[worldData.AirLayer0]);
		}
#endregion


#region Apply Latent Heat

		energyJobHandles[worldData.IceLayer] = SimJob.Schedule(
			JobType.Schedule, 64,
			new ApplyLatentHeatIceJob()
			{
				IceTemperature = nextState.IceTemperature,
				IceMass = nextState.IceMass,
				LatentHeat = tempState.LatentHeatIce
			}, energyJobHandles[worldData.IceLayer]);

		energyJobHandles[worldData.TerrainLayer] = SimJob.Schedule(
			JobType.Schedule, 64,
			new ApplyLatentHeatTerrainJob()
			{
				TerrainTemperature = nextState.GroundTemperature,

				LatentHeat = tempState.LatentHeatTerrain,
				SpecificHeatTerrain = tempState.SpecificHeatTerrain,
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.TerrainLayer], energyJobHandles[worldData.FloraLayer], energyJobHandles[worldData.AirLayer0]));

		energyJobHandles[worldData.LavaLayer] = SimJob.Schedule(
			JobType.Schedule, 64,
			new ApplyLatentHeatLavaJob()
			{
				LavaTemperature = nextState.LavaTemperature,

				LatentHeat = tempState.LatentHeatLava,
				LavaMass = nextState.LavaMass
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.LavaLayer], energyJobHandles[worldData.TerrainLayer]));

		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(
			JobType.Schedule, 64,
			new ApplyLatentHeatAirJob()
			{
				AirTemperaturePotential = staticState.GetSliceAir(nextState.AirTemperaturePotential),
				AirMass = staticState.GetSliceAir(tempState.AirMass),
				VaporMass = staticState.GetSliceAir(nextState.AirVapor),
				CloudMass = nextState.CloudMass,
				CloudElevation = tempState.CloudElevation,
				LayerElevation = staticState.GetSliceAir(tempState.AirLayerElevation),
				LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
				LatentHeatAir = staticState.GetSliceAir(tempState.LatentHeatAir),
				LatentHeatCloud = staticState.GetSliceAir(tempState.LatentHeatCloud),				
				Count = staticState.Count,
				LayerCount = worldData.AirLayers - 2,
			}, energyJobHandles[worldData.AirLayer0]);
		energyJobHandles[worldData.CloudLayer] = JobHandle.CombineDependencies(energyJobHandles[worldData.CloudLayer], energyJobHandles[worldData.AirLayer0]);

		energyJobHandles[worldData.WaterLayer0] = SimJob.Schedule(
			JobType.Schedule, 64,
			new ApplyLatentHeatWaterJob()
			{
				WaterTemperature = staticState.GetSliceLayer(nextState.WaterTemperature,worldData.SurfaceWaterLayer),
				WaterMass = staticState.GetSliceLayer(nextState.WaterMass,worldData.SurfaceWaterLayer),
				SaltMass = staticState.GetSliceLayer(nextState.WaterSaltMass,worldData.SurfaceWaterLayer),
				LatentHeat = tempState.LatentHeatWaterSurface,
			}, energyJobHandles[worldData.WaterLayer0]);

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

			lastJobHandle = SimJob.Schedule(
				JobType.Schedule, 64,
				new GroundWaterFlowJob()
				{
					GroundWater = nextState.GroundWater,
					GroundWaterTemperature = nextState.GroundWaterTemperature,

					LastGroundWater = lastState.GroundWater,
					LastGroundWaterTemperature = lastState.GroundWaterTemperature,
					SurfaceElevation = tempState.SurfaceElevation,
					Neighbors = staticState.Neighbors,
					NeighborDistInverse = staticState.NeighborDistInverse,
					FlowSpeed = worldData.GroundWaterFlowSpeed,
					GroundWaterMaxInverse = 1.0f / worldData.GroundWaterMax,
				}, JobHandle.CombineDependencies(lastJobHandle, energyJobHandles[worldData.TerrainLayer]));

			lastJobHandle = SimJob.Schedule(
				JobType.Schedule, 64,
				new GroundWaterDiffusionJob()
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
				lastJobHandle = SimJob.Schedule(
					JobType.Schedule, 64,
					new GroundWaterAbsorptionJob()
					{
						GroundWater = nextState.GroundWater,
						GroundWaterTemperature = nextState.GroundWaterTemperature,
						WaterMass = staticState.GetSliceLayer(nextState.WaterMass,i),
						WaterTemperature = staticState.GetSliceLayer(nextState.WaterTemperature,i),

						LastGroundWater = tempState.GroundWaterFlowMass,
						LastGroundWaterTemperature = tempState.GroundWaterFlowTemperature,
						SaltMass = staticState.GetSliceLayer(nextState.WaterSaltMass,i),
						WaterCoverageBelow = staticState.GetSliceLayer(tempState.WaterCoverage, i-1),
						GroundWaterAbsorptionRate = worldData.GroundWaterAbsorptionRate * worldData.SecondsPerTick,
						GroundWaterMaxInverse = 1.0f / worldData.GroundWaterMax,
						GroundWaterMax = worldData.GroundWaterMax,
						IsTop = i == worldData.SurfaceWaterLayer,
						Count = staticState.Count
					}, JobHandle.CombineDependencies(lastJobHandle, energyJobHandles[worldData.WaterLayer0]));
				energyJobHandles[worldData.WaterLayer0] = JobHandle.CombineDependencies(energyJobHandles[worldData.WaterLayer0], lastJobHandle);
			}

			lastJobHandle = SimJob.Schedule(
				JobType.Schedule, 64,
				new GroundWaterConductionJob()
				{
					GroundWaterTemperature = tempState.GroundWaterFlowTemperature,
					TerrainTemperature = nextState.GroundTemperature,

					GroundWater = nextState.GroundWater,
					LastGroundWaterTemperature = nextState.GroundWaterTemperature,
					SpecificHeatTerrain = tempState.SpecificHeatTerrain,
					GroundWaterConductionCoefficient = WorldData.ConductivityWaterTerrain,
					SecondsPerTick = worldData.SecondsPerTick,
					GroundWaterSurfaceAreaInverse = worldData.SoilHeatDepth / worldData.GroundWaterMaxDepth
				}, JobHandle.CombineDependencies(lastJobHandle, energyJobHandles[worldData.TerrainLayer]));
			energyJobHandles[worldData.TerrainLayer] = JobHandle.CombineDependencies(energyJobHandles[worldData.TerrainLayer], lastJobHandle);
		}
		return lastJobHandle;
	}
}
