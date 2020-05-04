#define LayerRefactor

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


	private int _cellCount;

	private NativeArray<JobHandle> applyAdvectionWaterJobHandles;
	private NativeArray<JobHandle> energyJobHandles;
	private List<NativeList<JobHandle>> jobHandleDependencies = new List<NativeList<JobHandle>>();
	private List<NativeArray<float>> tempArrays = new List<NativeArray<float>>();

	public WorldSim(int cellCount, ref WorldData worldData)
	{
		_cellCount = cellCount;

		SimJob = new JobHelper(_cellCount);
		NeighborJob = new JobHelper(_cellCount * 6);
		AirJob = new JobHelper(_cellCount * (worldData.AirLayers - 2));
		UpperAirJob = new JobHelper(_cellCount * (worldData.AirLayers - 3));
		AirNeighborJob = new JobHelper(_cellCount * (worldData.AirLayers - 2) * StaticState.MaxNeighborsVert);

		applyAdvectionWaterJobHandles = new NativeArray<JobHandle>(worldData.WaterLayers, Allocator.Persistent);
		energyJobHandles = new NativeArray<JobHandle>(worldData.LayerCount, Allocator.Persistent);

	}

	public void Dispose(ref WorldData worldData)
	{
		energyJobHandles.Dispose();
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
		tickJobHandle = tempState.Update(ref lastState, ref worldData, ref staticState, tickJobHandle);

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
		}, energyJobHandles[worldData.AirLayer0]);

		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(new AccelerationAirJob()
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


		energyJobHandles[worldData.AirLayer0] = AirNeighborJob.Schedule(new GetVectorDestCoordsVerticalJob()
		{
			Destination = staticState.GetSliceAirNeighbors(tempState.DestinationAir),
			NeighborsVert = staticState.NeighborsVert,
			Position = staticState.SphericalPosition,
			Velocity = nextState.AirVelocity,
			Mass = tempState.AirMass,
			PlanetRadius = staticState.PlanetRadius,
			SecondsPerTick = worldData.SecondsPerTick,
			CellsPerLayer = staticState.Count
		}, energyJobHandles[worldData.AirLayer0]);

		energyJobHandles[worldData.AirLayer0] = AirNeighborJob.Schedule(new ResolveAdvectionConflictVert()
		{
			ResolvedDestination = staticState.GetSliceAirNeighbors(tempState.DestinationAirResolved),
			Destination = tempState.DestinationAir,
			ReverseNeighborsVert = staticState.ReverseNeighborsVert,
			CellsPerLayer = staticState.Count
		}, energyJobHandles[worldData.AirLayer0]);

#if !LayerRefactor
		if (settings.MakeAirIncompressible)
		{

			for (int j = 1; j < worldData.AirLayers - 1; j++)
			{
				int layer = worldData.AirLayer0 + j;
				energyJobHandles[layer] = SimJob.Schedule(new GetDivergenceJob()
				{
					Divergence = tempState.DivergenceAir[j],
					Destination = tempState.DestinationAirResolved[j],
					DestinationAbove = tempState.DestinationAirResolved[j + 1],
					DestinationBelow = tempState.DestinationAirResolved[j - 1],
					Neighbors = staticState.Neighbors,
					Mass = tempState.AirMass[j],
					MassAbove = tempState.AirMass[j + 1],
					MassBelow = tempState.AirMass[j - 1],
				}, JobHandle.CombineDependencies(energyJobHandles[layer], energyJobHandles[layer - 1], energyJobHandles[layer + 1]));
			}

			// Calculate Pressure gradient field
			var divergenceJobHandle = default(JobHandle);
			for (int i = 1; i < worldData.AirLayers - 1; i++)
			{
				divergenceJobHandle = JobHandle.CombineDependencies(divergenceJobHandle, energyJobHandles[worldData.AirLayer0 + i]);
			}
			for (int a = 0; a < settings.IncompressibilityIterations; a++)
			{
				for (int i = 1; i < worldData.AirLayers - 1; i++)
				{
					var dpj = new GetDivergencePressureJob()
					{
						Pressure = tempState.DivergencePressureAir[i],
						Divergence = tempState.DivergenceAir[i],
						Neighbors = staticState.Neighbors
					};
					divergenceJobHandle = dpj.Schedule(_cellCount, divergenceJobHandle);
				}
			}

			for (int i = 1; i < worldData.AirLayers - 1; i++)
			{
				int layer = worldData.AirLayer0 + i;
				energyJobHandles[layer] = NeighborJob.Schedule(new GetDivergenceFreeFieldJob()
				{
					Destination = tempState.DestinationAirResolved[i],
					Pressure = tempState.DivergencePressureAir[i],
					Neighbors = staticState.Neighbors,
					Mass = tempState.AirMass[i],
				}, divergenceJobHandle);
			}
			for (int i = 1; i < worldData.AirLayers - 1; i++)
			{
				int layer = worldData.AirLayer0 + i;
				energyJobHandles[layer] = SimJob.Schedule(new SumMassLeavingJob()
				{
					MassLeaving = tempState.AirMassLeaving[i],
					Destination = tempState.DestinationAirResolved[i],
				}, energyJobHandles[layer]);
			}
			for (int i = 1; i < worldData.AirLayers - 1; i++)
			{
				int layer = worldData.AirLayer0 + i;
				energyJobHandles[layer] = NeighborJob.Schedule(new CapMassLeavingJob()
				{
					Destination = tempState.DestinationAirResolved[i],
					MassLeaving = tempState.AirMassLeaving[i],
					Mass = tempState.AirMass[i],
					Neighbors = staticState.Neighbors
				}, energyJobHandles[layer]);
			}

			for (int i = 1; i < worldData.AirLayers - 1; i++)
			{
				int layer = worldData.AirLayer0 + i;
				energyJobHandles[layer] = SimJob.Schedule(new UpdateDivergenceFreeVelocityJob()
				{
					Velocity = nextState.AirVelocity[i],
					Destination = tempState.DestinationAirResolved[i],
					NeighborTangent = staticState.NeighborTangent,
					Mass = tempState.AirMass[i],
					TicksPerSecond = worldData.TicksPerSecond
				}, energyJobHandles[layer]);
			}

		}

#endif

		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(new AdvectionAirJob()
		{
			Delta = staticState.GetSliceAir(tempState.AdvectionAir),
			Temperature = nextState.AirTemperaturePotential,
			AirMass = tempState.AirMass,
			Vapor = nextState.AirVapor,
			CarbonDioxide = nextState.AirCarbon,
			Dust = nextState.Dust,
			Velocity = nextState.AirVelocity,
			NeighborsVert = staticState.NeighborsVert,
			Destination = tempState.DestinationAirResolved,
			Positions = staticState.SphericalPosition,
			CoriolisMultiplier = staticState.CoriolisMultiplier,
			CoriolisTerm = coriolisTerm,
			SecondsPerTick = worldData.SecondsPerTick,
			CellsPerLayer = staticState.Count
		}, energyJobHandles[worldData.AirLayer0]);

		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(new ApplyAdvectionAirJob()
		{
			Advection = tempState.AdvectionAir,
			Vapor = nextState.AirVapor,
			Dust = nextState.Dust,
			CarbonDioxide = nextState.AirCarbon,
			Temperature = nextState.AirTemperaturePotential,
			AirVelocity = nextState.AirVelocity,
		}, energyJobHandles[worldData.AirLayer0]);

		#endregion

		// Diffuse from last time step
		// Air, Water, Cloud
		#region Diffusion

		// TODO: is it a problem that we are using the dependent variables from last frame while referencing our newly calculated next frame values for temperature and such?
		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(new DiffusionAirJob()
		{
			Delta = staticState.GetSliceAir(tempState.DiffusionAir),

			AirMass = staticState.GetSliceAir(tempState.AirMass),
			Temperature = staticState.GetSliceAir(nextState.AirTemperaturePotential),
			Vapor = staticState.GetSliceAir(nextState.AirVapor),
			CarbonDioxide = staticState.GetSliceAir(nextState.AirCarbon),
			Dust = staticState.GetSliceAir(nextState.Dust),
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

		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(new ApplyAdvectionAirJob()
		{
			Advection = staticState.GetSliceAir(tempState.DiffusionAir),
			Vapor = staticState.GetSliceAir(nextState.AirVapor),
			Dust = staticState.GetSliceAir(nextState.Dust),
			CarbonDioxide = staticState.GetSliceAir(nextState.AirCarbon),
			Temperature = staticState.GetSliceAir(nextState.AirTemperaturePotential),
			AirVelocity = staticState.GetSliceAir(nextState.AirVelocity),
		}, energyJobHandles[worldData.AirLayer0]);



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
			AirVelocityUp = staticState.GetSliceLayer(lastState.AirVelocity,worldData.SurfaceAirLayer),
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
				SurfaceElevation = tempState.SurfaceElevation,
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

#if !LayerRefactor
		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			int layer = worldData.WaterLayer0 + j;
			energyJobHandles[layer] = NeighborJob.Schedule(new GetVectorDestCoordsVerticalJob()
			{
				Destination = tempState.DestinationWater[j],
				NeighborsVert = staticState.NeighborsVert,
				Position = staticState.SphericalPosition,
				Velocity = nextState.WaterVelocity[j],
				Mass = nextState.WaterMass[j],
				PlanetRadius = staticState.PlanetRadius,
				SecondsPerTick = worldData.SecondsPerTick,
				CellsPerLayer = staticState.Count
			}, energyJobHandles[layer]);
		}

		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			int layer = worldData.WaterLayer0 + j;
			energyJobHandles[layer] = NeighborJob.Schedule(new ResolveAdvectionConflictVert()
			{
				ResolvedDestination = tempState.DestinationWaterResolved[j],
				Destination = tempState.DestinationWater[j],
				ReverseNeighborsVert = staticState.ReverseNeighbors,
				CellsPerLayer = staticState.Count
			}, energyJobHandles[layer]);
		}

		if (settings.MakeWaterIncompressible)
		{
			for (int j = 1; j < worldData.WaterLayers - 1; j++)
			{
				int layer = worldData.WaterLayer0 + j;
				energyJobHandles[layer] = SimJob.Schedule(new GetDivergenceJob()
				{
					Divergence = tempState.DivergenceWater[j],
					Destination = tempState.DestinationWater[j],
					CellsPerLayer = staticState.Count
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
						NeighborsVert = staticState.NeighborsVert,
						CellsPerLayer = staticState.Count
					};
					divergenceJobHandle = dpj.Schedule(_cellCount, divergenceJobHandle);
				}
			}

			for (int i = 1; i < worldData.WaterLayers - 1; i++)
			{
				int layer = worldData.WaterLayer0 + i;
				energyJobHandles[worldData.WaterLayer0 + i] = NeighborJob.Schedule(new GetDivergenceFreeFieldJob()
				{
					Destination = tempState.DestinationWater[i],
					Pressure = tempState.DivergencePressureWater[i],
					Neighbors = staticState.Neighbors,
					Mass = nextState.WaterMass[i],
					CellsPerLayer = staticState.Count
				}, divergenceJobHandle);
			}

			for (int i = 1; i < worldData.WaterLayers - 1; i++)
			{
				int layer = worldData.WaterLayer0 + i;
				energyJobHandles[layer] = SimJob.Schedule(new SumMassLeavingJob()
				{
					MassLeaving = tempState.WaterMassLeaving[i],
					Destination = tempState.DestinationWaterResolved[i],
					CellsPerLayer = staticState.Count
				}, energyJobHandles[layer]);
			}
			for (int i = 1; i < worldData.WaterLayers - 1; i++)
			{
				int layer = worldData.WaterLayer0 + i;
				energyJobHandles[layer] = NeighborJob.Schedule(new CapMassLeavingJob()
				{
					Destination = tempState.DestinationWaterResolved[i],
					MassLeaving = tempState.WaterMassLeaving[i],
					Mass = nextState.WaterMass[i],
					NeighborsVert = staticState.NeighborsVert,
					CellsPerLayer = staticState.Count
				}, energyJobHandles[layer]);
			}

			for (int i = 1; i < worldData.WaterLayers - 1; i++)
			{
				int layer = worldData.WaterLayer0 + i;
				energyJobHandles[layer] = SimJob.Schedule(new UpdateDivergenceFreeVelocityJob()
				{
					Velocity = nextState.WaterVelocity[i],
					DestinationVert = tempState.DestinationWaterResolved[i],
					NeighborTangent = staticState.NeighborTangent,
					Mass = nextState.WaterMass[i],
					TicksPerSecond = worldData.TicksPerSecond,
					CellsPerLayer = staticState.Count
				}, energyJobHandles[layer]);
			}
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
#endif

		energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(new GetVectorDestCoordsJob()
		{
			Destination = tempState.DestinationCloud,
			Neighbors = staticState.Neighbors,
			Position = staticState.SphericalPosition,
			Velocity = tempState.CloudVelocity,
			PlanetRadius = staticState.PlanetRadius,
			SecondsPerTick = worldData.SecondsPerTick,
			MaxWindMove = staticState.CellRadius * 0.9f,
			CellsPerLayer = staticState.Count,
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
				LayerHeight = tempState.WaterLayerHeight[i],
				LayerHeightAbove = tempState.WaterLayerHeight[i + 1],
				LayerHeightBelow = tempState.WaterLayerHeight[i - 1],
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

		tickJobHandle = tempState.Update(ref nextState, ref worldData, ref staticState, tickJobHandle);

#endregion


#region Flow

		if (settings.WaterSurfaceFlowEnabled)
		{
			// TODO: surface elevation is inaccurate now, we should recalculate (and use water surfae, not ice surface)
			tickJobHandle = NeighborJob.Schedule(new UpdateFlowVelocityJob()
			{
				Flow = nextState.FlowWater,

				LastFlow = lastState.FlowWater,
				SurfaceElevation = tempState.SurfaceElevation,
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

		tickJobHandle = tempState.Update(ref nextState, ref worldData, ref staticState, tickJobHandle);

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

		JobHandle emissivityJobHandle = default(JobHandle);
		emissivityJobHandle = JobHandle.CombineDependencies(emissivityJobHandle, AirJob.Schedule(new EmissivityAirJob()
		{
			Emissivity = staticState.GetSliceLayers(tempState.EmissivityAir, 1, worldData.AirLayers - 2),
			AirMass = staticState.GetSliceAir(tempState.AirMass),
			VaporMass = staticState.GetSliceAir(lastState.AirVapor),
			Dust = staticState.GetSliceAir(lastState.Dust),
			CarbonDioxide = staticState.GetSliceAir(lastState.AirCarbon),
			EmissivityAir = worldData.ThermalEmissivityAir,
			EmissivityWaterVapor = worldData.ThermalEmissivityWaterVapor,
			EmissivityDust = worldData.ThermalEmissivityDust,
			EmissivityCarbonDioxide = worldData.ThermalEmissivityCarbonDioxide,
			EmissivityOxygen = worldData.ThermalEmissivityOxygen
		}, lastJobHandle));

		// we only do thermal radiation upwards for the surface layer of water,
		// for the bottom we rely on conduction with the terrain for heat transfer (although this might lead to an imbalance!)
		{
			emissivityJobHandle = JobHandle.CombineDependencies(emissivityJobHandle, SimJob.Schedule(new EmissivityWaterJob()
			{
				Emissivity = staticState.GetSliceLayer(tempState.EmissivityWater,worldData.SurfaceWaterLayer),
				WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
				SaltMass = lastState.SaltMass[worldData.SurfaceWaterLayer],
				EmissivitySalt = worldData.ThermalEmissivitySalt,
				EmissivityWater = worldData.ThermalEmissivityWater
			}, lastJobHandle));
		}
		emissivityJobHandle = JobHandle.CombineDependencies(emissivityJobHandle, SimJob.Schedule(new EmissivityTerrainJob()
		{
			Emissivity = tempState.EmissivityTerrain,
			SoilFertility = lastState.GroundCarbon,
			EmissivityDirt = worldData.ThermalEmissivityDirt,
			EmissivitySand = worldData.ThermalEmissivitySand,
		}, lastJobHandle));

#endregion

		// Calculate how much thermal radition is being emitted out of each layer
#region Thermal Radiation
		JobHandle thermalOutJobHandle = default(JobHandle);

		// ICE
		thermalOutJobHandle = JobHandle.CombineDependencies(thermalOutJobHandle, SimJob.Schedule(new ThermalEnergyRadiatedConstantEmissivityJob()
		{
			ThermalRadiationEmitted = tempState.ThermalRadiationEmittedIce,

			Emissivity = worldData.ThermalEmissivityIce,
			Energy = tempState.IceEnergy,
			Temperature = lastState.IceTemperature,
			SurfaceArea = tempState.IceCoverage,
			SecondsPerTick = worldData.SecondsPerTick
		}, lastJobHandle));


		// FLORA
		thermalOutJobHandle = JobHandle.CombineDependencies(thermalOutJobHandle, SimJob.Schedule(new ThermalEnergyRadiatedConstantEmissivityJob()
		{
			ThermalRadiationEmitted = tempState.ThermalRadiationEmittedFlora,

			Emissivity = worldData.ThermalEmissivityFlora,
			Energy = tempState.FloraEnergy,
			Temperature = lastState.FloraTemperature,
			SurfaceArea = tempState.FloraCoverage,
			SecondsPerTick = worldData.SecondsPerTick
		}, lastJobHandle));


		// TERRAIN
		thermalOutJobHandle = JobHandle.CombineDependencies(thermalOutJobHandle, SimJob.Schedule(new ThermalEnergyRadiatedTerrainJob()
		{
			ThermalRadiationEmitted = tempState.ThermalRadiationEmittedTerrain,

			PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
			Emissivity = tempState.EmissivityTerrain,
			Temperature = lastState.GroundTemperature,
			SecondsPerTick = worldData.SecondsPerTick
		}, emissivityJobHandle));


		// ATMOSPHERE
		thermalOutJobHandle = JobHandle.CombineDependencies(thermalOutJobHandle, AirJob.Schedule(new ThermalEnergyRadiatedAirJob()
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
			int layer = worldData.SurfaceWaterLayer + worldData.WaterLayer0;
			thermalOutJobHandle = JobHandle.CombineDependencies(thermalOutJobHandle, SimJob.Schedule(new ThermalEnergyRadiatedWaterJob()
			{
				ThermalRadiationEmitted = staticState.GetSliceLayer(tempState.ThermalRadiationEmittedWater, worldData.SurfaceWaterLayer),

				Emissivity = staticState.GetSliceLayer(tempState.EmissivityWater, worldData.SurfaceWaterLayer),
				Energy = tempState.WaterPotentialEnergy[worldData.SurfaceWaterLayer],
				TemperatureAbsolute = lastState.WaterTemperature[worldData.SurfaceWaterLayer],
				SurfaceArea = tempState.WaterCoverage[worldData.SurfaceWaterLayer],
				SecondsPerTick = worldData.SecondsPerTick
			}, emissivityJobHandle));
		}
#endregion


#region absorptivity

		solarInJobHandle = SimJob.Schedule(new CloudAlbedoJob()
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

		solarInJobHandle = AirJob.Schedule(new AbsorptivityAirJob()
		{
			AbsorptivitySolar = staticState.GetSliceAir(tempState.AbsorptivitySolar),
			AbsorptivityThermal = staticState.GetSliceAir(tempState.AbsorptivityThermal),
			AirMass = staticState.GetSliceAir(tempState.AirMass),
			VaporMass = staticState.GetSliceAir(lastState.AirVapor),
			AirCarbonDioxide = staticState.GetSliceAir(lastState.AirCarbon),
			Dust = staticState.GetSliceAir(lastState.Dust),
			CloudMass = lastState.CloudMass,
			CloudAlbedo = tempState.CloudAlbedo,
			CloudAbsorptivity = tempState.CloudAbsorptivity,
			CloudElevation = tempState.CloudElevation,
			LayerElevation = staticState.GetSliceAir(tempState.AirLayerElevation),
			LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
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
			Count = staticState.Count
		}, solarInJobHandle);

#endregion


		// Follow the solar radiation down from the top of the atmosphere to ther terrain, and absorb some as it passes through each layer
#region Solar Radiation Absorbed

		// process each vertical layer in order

		// atmosphere
		for (int j = worldData.AirLayer0 + worldData.AirLayers - 2; j > worldData.AirLayer0; j--)
		{
			int airLayerIndex = j - worldData.AirLayer0;
			solarInJobHandle = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedAirJob()
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
		solarInJobHandle = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedPartialCoverageConstantAlbedoJob()
		{
			SolarRadiationAbsorbed = tempState.SolarRadiationInIce,
			SolarRadiationReflected = tempState.SolarReflectedIce,
			SolarRadiationIncoming = tempState.SolarRadiation,
			AlbedoSlope = tempState.AlbedoSlope,
			AlbedoMin = WorldData.AlbedoIce,
			AlbedoRange = 1.0f - WorldData.AlbedoIce,
			Coverage = tempState.IceCoverage
		}, solarInJobHandle);

		for (int j = worldData.SurfaceWaterLayer; j >= 1; j--)
		{
			int layerIndex = worldData.WaterLayer0 + j;
			solarInJobHandle = solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedSlopeJob()
			{
				SolarRadiationAbsorbed = staticState.GetSliceLayer(tempState.SolarRadiationInWater,j),
				SolarRadiationReflected = staticState.GetSliceLayer(tempState.SolarReflectedWater,j),
				SolarRadiationIncoming = tempState.SolarRadiation,
				Coverage = tempState.WaterCoverage[layerIndex - worldData.WaterLayer0],
				AlbedoSlope = tempState.AlbedoSlope,
				AlbedoMin = WorldData.AlbedoWater,
			}, solarInJobHandle);
		}

		// flora
		solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedPartialCoverageConstantAlbedoJob()
		{
			SolarRadiationAbsorbed = tempState.SolarRadiationInFlora,
			SolarRadiationReflected = tempState.SolarReflectedFlora,
			SolarRadiationIncoming = tempState.SolarRadiation,
			AlbedoSlope = tempState.AlbedoSlope,
			AlbedoMin = WorldData.AlbedoFloraMin,
			AlbedoRange = WorldData.AlbedoFloraRange,
			Coverage = tempState.FloraCoverage
		}, solarInJobHandle);

		// NOTE: we don't bother with solar radiation in lava

		solarInJobHandle = SimJob.Schedule(new SolarRadiationAbsorbedTerrainJob()
		{
			SolarRadiationAbsorbed = tempState.SolarRadiationInTerrain,
			SolarRadiationReflected = tempState.SolarReflectedTerrain,
			SolarRadiationIncoming = tempState.SolarRadiation,
			worldData = worldData,
			SoilFertility = lastState.GroundCarbon,
		}, solarInJobHandle);
#endregion


		// Thermal radiation travels upwards, partially reflecting downwards (clouds), partially absorbed, and partially lost to space
#region Thermal Radiation Absorbed Up

		// transmit up from land
		for (int j = 0; j < worldData.LayerCount; j++)
		{

			if (j > worldData.AirLayer0 && j < worldData.AirLayer0 + worldData.AirLayers - 1)
			{
				int airLayer = j - worldData.AirLayer0;
				int downIndex = airLayer == 1 ? worldData.IceLayer : (j - 1);

				thermalOutJobHandle = SimJob.Schedule(new ThermalEnergyAbsorbedAirJob()
				{
					ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaAir,airLayer),
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp,
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp,

					ThermalRadiationEmitted = tempState.ThermalRadiationEmittedAir,
					PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
					AbsorptivityThermal = staticState.GetSliceLayer(tempState.AbsorptivityThermal,airLayer),
					LayerElevation = staticState.GetSliceLayer(tempState.AirLayerElevation,airLayer),
					LayerHeight = staticState.GetSliceLayer(tempState.AirLayerHeight,airLayer),
					CloudElevation = tempState.CloudElevation,
					FromTop = false,
				}, JobHandle.CombineDependencies(solarInJobHandle, thermalOutJobHandle));
			}
			else if (j == worldData.IceLayer)
			{
				int downIndex = worldData.SurfaceWaterLayerGlobal;
				thermalOutJobHandle = SimJob.Schedule(new ThermalEnergyAbsorbedUpPartialCoverageJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDeltaIce,
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp,
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp,

					ThermalRadiationEmitted = tempState.ThermalRadiationEmittedIce,
					PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
					Coverage = tempState.IceCoverage,

				}, thermalOutJobHandle);
			}
			else if (j == worldData.FloraLayer)
			{
				thermalOutJobHandle = SimJob.Schedule(new ThermalEnergyAbsorbedUpPartialCoverageJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDeltaFlora,
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp,
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp,

					ThermalRadiationEmitted = tempState.ThermalRadiationEmittedFlora,
					PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
					Coverage = tempState.FloraCoverage,

				}, thermalOutJobHandle);
			}
			else if (j > worldData.WaterLayer0 && j < worldData.WaterLayer0 + worldData.WaterLayers - 1)
			{
				int waterLayerIndex = j - worldData.WaterLayer0;
				int downIndex = (waterLayerIndex == 1) ? worldData.FloraLayer : (j - 1);
				thermalOutJobHandle = SimJob.Schedule(new ThermalEnergyAbsorbedUpPartialCoverageJob()
				{
					ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaWater, waterLayerIndex),
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedUp,
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedUp,

					ThermalRadiationEmitted = staticState.GetSliceLayer(tempState.ThermalRadiationEmittedWater, waterLayerIndex),
					PercentRadiationInAtmosphericWindow = worldData.EnergyLostThroughAtmosphereWindow,
					Coverage = tempState.WaterCoverage[waterLayerIndex],
				}, thermalOutJobHandle);
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
				thermalOutJobHandle = SimJob.Schedule(new ThermalEnergyAbsorbedTerrainJob()
				{
					ThermalRadiationAbsorbed = tempState.ThermalRadiationDeltaTerrain,

					WindowRadiationIncoming = tempState.WindowRadiationTransmittedDown,
					ThermalRadiationIncoming = tempState.ThermalRadiationTransmittedDown,
				}, thermalOutJobHandle);
			}
			else if (j == worldData.FloraLayer)
			{
				// FLORA
				thermalOutJobHandle = SimJob.Schedule(new ThermalEnergyAbsorbedDownPartialCoverageJob()
				{
					ThermalRadiationDelta = tempState.ThermalRadiationDeltaFlora,
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedDown,
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedDown,

					Coverage = tempState.FloraCoverage,
				}, thermalOutJobHandle);
			}
			else if (j == worldData.IceLayer)
			{
				// ICE
				thermalOutJobHandle = SimJob.Schedule(new ThermalEnergyAbsorbedDownPartialCoverageJob()
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
				thermalOutJobHandle = SimJob.Schedule(new ThermalEnergyAbsorbedDownPartialCoverageJob()
				{
					ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaWater, waterLayerIndex),
					ThermalRadiationTransmitted = tempState.ThermalRadiationTransmittedDown,
					WindowRadiationTransmitted = tempState.WindowRadiationTransmittedDown,

					Coverage = tempState.WaterCoverage[waterLayerIndex],
				}, thermalOutJobHandle);
			}
			else if (j > worldData.AirLayer0 && j < worldData.AirLayer0 + worldData.AirLayers - 1)
			{
				int airLayer = j - worldData.AirLayer0;
				thermalOutJobHandle = SimJob.Schedule(new ThermalEnergyAbsorbedAirJob()
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
				solarInJobHandle,
				thermalOutJobHandle,
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
			SolarRadiationIn = tempState.SolarRadiationInTerrain,
			ThermalRadiationDelta = tempState.ThermalRadiationDeltaTerrain,
			ConductionEnergyAir = tempState.ConductionAirTerrain,
			ConductionEnergyIce = tempState.ConductionIceTerrain,
			ConductionEnergyFlora = tempState.ConductionFloraTerrain,
			ConductionEnergyWater = tempState.ConductionWaterTerrainTotal,
			GeothermalEnergy = tempState.GeothermalRadiation,
			HeatingDepth = worldData.SoilHeatDepth,
		}, JobHandle.CombineDependencies(terrainEnergyJobHandleDependencies));

		var energyIceJobHandleDependencies = new NativeList<JobHandle>(Allocator.Persistent)
			{
				solarInJobHandle,
				thermalOutJobHandle,
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
			SolarRadiationIn = tempState.SolarRadiationInIce,
			ThermalRadiationDelta = tempState.ThermalRadiationDeltaIce,
			ConductionEnergyAir = tempState.ConductionAirIce,
			ConductionEnergyTerrain = tempState.ConductionIceTerrain,
			ConductionEnergyWater = tempState.ConductionIceWater,
			ConductionEnergyFlora = tempState.ConductionIceFlora,
		}, JobHandle.CombineDependencies(energyIceJobHandleDependencies));

		var energyFloraJobHandleDependencies = new NativeList<JobHandle>(Allocator.Persistent)
			{
				solarInJobHandle,
				thermalOutJobHandle,
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
			ThermalRadiationDelta = tempState.ThermalRadiationDeltaFlora,
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

		energyJobHandles[worldData.AirLayer0] = UpperAirJob.Schedule(new EnergyAirJob()
		{
			AirTemperaturePotential = staticState.GetSliceLayers(nextState.AirTemperaturePotential, 2, worldData.AirLayers-3),
			LastTemperaturePotential = staticState.GetSliceLayers(lastState.AirTemperaturePotential, 2, worldData.AirLayers - 3),
			LastVapor = staticState.GetSliceLayers(lastState.AirVapor, 2, worldData.AirLayers - 3),
			AirMass = staticState.GetSliceLayers(tempState.AirMass, 2, worldData.AirLayers - 3),
			SolarRadiationIn = staticState.GetSliceLayers(tempState.SolarRadiationInAir, 2, worldData.AirLayers - 3),
			ThermalRadiationDelta = staticState.GetSliceLayers(tempState.ThermalRadiationDeltaAir, 2, worldData.AirLayers - 3),
		}, JobHandle.CombineDependencies(solarInJobHandle, thermalOutJobHandle));

		var airSurfaceDependencies = new NativeList<JobHandle>(Allocator.Persistent)
		{
			energyJobHandles[worldData.AirLayer0],
			solarInJobHandle,
			thermalOutJobHandle,
			conductionAirWaterJobHandle,
			conductionAirIceJobHandle,
			conductionAirFloraJobHandle,
			conductionAirTerrainJobHandle,
		};
		jobHandleDependencies.Add(airSurfaceDependencies);
		energyJobHandles[worldData.AirLayer0] = SimJob.Schedule(new EnergyAirSurfaceJob()
		{
			AirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential,worldData.SurfaceAirLayer),
			LastTemperaturePotential = staticState.GetSliceLayer(lastState.AirTemperaturePotential, worldData.SurfaceAirLayer),
			LastVapor = staticState.GetSliceLayer(lastState.AirVapor, worldData.SurfaceAirLayer),
			AirMass = staticState.GetSliceLayer(tempState.AirMass, worldData.SurfaceAirLayer),
			SolarRadiationIn = staticState.GetSliceLayer(tempState.SolarRadiationInAir, worldData.SurfaceAirLayer),
			ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaAir, worldData.SurfaceAirLayer),
			ConductionEnergyWater = tempState.ConductionAirWater,
			ConductionEnergyIce = tempState.ConductionAirIce,
			ConductionEnergyFlora = tempState.ConductionAirFlora,
			ConductionEnergyTerrain = tempState.ConductionAirTerrain,
		}, JobHandle.CombineDependencies(airSurfaceDependencies));

		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			int layerIndex = worldData.WaterLayer0 + j;
			var waterDependencies = new NativeList<JobHandle>(Allocator.Persistent)
				{
					solarInJobHandle,
					thermalOutJobHandle,
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
				ThermalRadiationDelta = staticState.GetSliceLayer(tempState.ThermalRadiationDeltaWater,j),
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
			energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(new FluxAirCondensationJob()
			{
				LatentHeat = staticState.GetSliceAir(tempState.LatentHeatAir),
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
				Count = staticState.Count
			}, energyJobHandles[worldData.AirLayer0]);
		}
		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(new FluxAirDustJob()
		{
			DustUp = staticState.GetSliceAir(tempState.DustUp),
			DustDown = staticState.GetSliceAir(tempState.DustDown),

			LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
			LastDust = staticState.GetSliceAir(lastState.Dust),
			AirVelocity = staticState.GetSliceAir(lastState.AirVelocity),
			DustVerticalVelocity = worldData.DustVerticalVelocity,
			Positions = staticState.SphericalPosition,
			SecondsPerTick = worldData.SecondsPerTick,
			Count = staticState.Count
		}, energyJobHandles[worldData.AirLayer0]);

		if (settings.Evaporation)
		{
			energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new FluxWaterEvaporationJob()
			{
				EvaporatedWaterMass = tempState.EvaporationMassWater,
				LatentHeatWater = tempState.LatentHeatWater,
				LatentHeatAir = staticState.GetSliceLayer(tempState.LatentHeatAir,worldData.SurfaceAirLayer),

				WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
				IceCoverage = tempState.IceCoverage,
				WaterCoverage = tempState.WaterCoverage[worldData.SurfaceWaterLayer],
				SurfaceWind = staticState.GetSliceLayer(lastState.AirVelocity,worldData.SurfaceAirLayer),
				AirMass = staticState.GetSliceLayer(tempState.AirMass,worldData.SurfaceAirLayer),
				AirPressure = staticState.GetSliceLayer(tempState.AirPressure,worldData.SurfaceAirLayer),
				AirVapor = staticState.GetSliceLayer(lastState.AirVapor,worldData.SurfaceAirLayer),
				WaterHeatingDepth = worldData.WaterHeatingDepth,
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.AirLayer0]));
		}

		if (settings.Freezing)
		{
			energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new FluxWaterFreezeJob()
			{
				FrozenMass = tempState.FrozenMass,
				FrozenTemperature = tempState.FrozenTemperature,
				LatentHeatWater = tempState.LatentHeatWater,
				SaltPlume = tempState.SaltPlume,

				WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
				SaltMass = lastState.SaltMass[worldData.SurfaceWaterLayer],
				AirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential,worldData.SurfaceAirLayer),
				AirLayerElevation = staticState.GetSliceLayer(tempState.AirLayerElevation,worldData.SurfaceAirLayer),
				WaterHeatingDepth = worldData.WaterHeatingDepth,
				FreezePointReductionPerSalinity = worldData.FreezePointReductionPerSalinity,
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.AirLayer0]));
		}

		energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new FluxFloraWaterConsumeJob()
		{
			FloraWaterConsumed = tempState.WaterConsumedByFlora,

			WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
			FloraMass = lastState.FloraMass,
			FloraWater = lastState.FloraWater,
			FloraWaterConsumptionRate = worldData.FloraWaterConsumptionRate
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.AirLayer0]));

		if (settings.Plankton)
		{
			energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new FluxPlanktonJob()
			{
				LatentHeatWater = tempState.LatentHeatWater,
				PlanktonMassDelta = tempState.PlanktonMassDelta,
				PlanktonGlucoseDelta = tempState.PlanktonGlucoseDelta,
				PlanktonDeath = tempState.PlanktonDeath,
				WaterCarbonDelta = tempState.WaterCarbonDelta,

				WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				WaterMass = lastState.WaterMass[worldData.SurfaceWaterLayer],
				SaltMass = lastState.SaltMass[worldData.SurfaceWaterLayer],
				SolarRadiation = staticState.GetSliceLayer(tempState.SolarRadiationInWater,worldData.SurfaceWaterLayer),
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
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.AirLayer0]));
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

				SurfaceAirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential,worldData.SurfaceAirLayer),
				SurfaceLayerElevation = staticState.GetSliceLayer(tempState.AirLayerElevation,worldData.SurfaceAirLayer),
				SurfaceLayerMiddle = staticState.GetSliceLayer(tempState.AirLayerMiddle,worldData.SurfaceAirLayer),
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
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.CloudLayer], energyJobHandles[worldData.AirLayer0]));
		}

		if (settings.Flora)
		{
			energyJobHandles[worldData.FloraLayer] = SimJob.Schedule(new FluxFloraJob()
			{
				LatentHeatAir = staticState.GetSliceLayer(tempState.LatentHeatAir, worldData.SurfaceAirLayer),
				LatentHeatFlora = tempState.LatentHeatFlora,
				EvaporatedWaterMass = tempState.FloraRespirationMassVapor,
				SurfaceWaterDelta = tempState.FloraRespirationMassWater,
				FloraMassDelta = tempState.FloraMassDelta,
				FloraWaterDelta = tempState.FloraWaterDelta,
				FloraGlucoseDelta = tempState.FloraGlucoseDelta,
				FloraDeath = tempState.FloraDeath,
				CarbonDioxideDelta = tempState.AirCarbonDelta,
				OxygenDelta = tempState.OxygenDelta,

				SolarRadiationIn = tempState.SolarRadiationInFlora,
				FloraTemperature = nextState.FloraTemperature,
				FloraMass = lastState.FloraMass,
				FloraGlucose = lastState.FloraGlucose,
				FloraWater = lastState.FloraWater,
				FloraCoverage = tempState.FloraCoverage,
				SoilFertility = lastState.GroundCarbon,
				FloraGrowthRate = worldData.FloraGrowthRate,
				FloraDeathRate = worldData.FloraDeathRate,
				CarbonDioxide = staticState.GetSliceLayer(lastState.AirCarbon,worldData.SurfaceAirLayer),
				LayerElevation = staticState.GetSliceLayer(tempState.AirLayerElevation,worldData.SurfaceAirLayer),
				LayerHeight = staticState.GetSliceLayer(tempState.AirLayerHeight,worldData.SurfaceAirLayer),
				SurfaceWind = staticState.GetSliceLayer(lastState.AirVelocity,worldData.SurfaceAirLayer),
				AirMass = staticState.GetSliceLayer(tempState.AirMass,worldData.SurfaceAirLayer),
				AirTemperaturePotential = staticState.GetSliceLayer(lastState.AirTemperaturePotential,worldData.SurfaceAirLayer),
				AirPressure = staticState.GetSliceLayer(tempState.AirPressure,worldData.SurfaceAirLayer),
				AirVapor = staticState.GetSliceLayer(lastState.AirVapor,worldData.SurfaceAirLayer),
				FloraGrowthTemperatureRangeInverse = worldData.FloraGrowthTemperatureRangeInverse,
				FloraEnergyForPhotosynthesis = worldData.FloraEnergyForPhotosynthesis,
				FloraCarbonDioxideExtractionEfficiency = worldData.FloraCarbonDioxideExtractionEfficiency,
				FloraOxygenExtractionEfficiency = worldData.FloraOxygenExtractionEfficiency,
				FloraPhotosynthesisSpeed = worldData.FloraPhotosynthesisSpeed,
				FloraRespirationSpeed = worldData.FloraRespirationSpeed,
				FloraRespirationPerDegree = worldData.FloraRespirationPerDegree,
				OxygenPercent = lastState.PlanetState.Oxygen,
				Gravity = lastState.PlanetState.Gravity
			}, JobHandle.CombineDependencies(energyJobHandles[worldData.FloraLayer], energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.AirLayer0]));
		}

		if (settings.IceMelting)
		{
			energyJobHandles[worldData.IceLayer] = SimJob.Schedule(new FluxIceMeltJob()
			{
				LatentHeatAir = staticState.GetSliceLayer(tempState.LatentHeatAir, worldData.SurfaceAirLayer),
				LatentHeatWater = tempState.LatentHeatWater,
				LatentHeatTerrain = tempState.LatentHeatTerrain,
				LatentHeatIce = tempState.LatentHeatIce,
				MeltedMass = tempState.IceMeltedMass,

				Temperature = nextState.IceTemperature,
				LastMass = lastState.IceMass,
				IceHeatingDepth = worldData.IceHeatingDepth,
				WaterIceSurfaceArea = tempState.SurfaceAreaIceWater,
				WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
				TerrainTemperature = nextState.GroundTemperature,
				SurfaceElevation = tempState.SurfaceElevation,
				AirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential, worldData.SurfaceAirLayer),

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
				LatentHeatLava = tempState.LatentHeatLava,

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

		energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new UpdateMassCondensationGroundJob()
		{
			SurfaceWaterMass = nextState.WaterMass[worldData.SurfaceWaterLayer],
			SurfaceWaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],

			SurfaceSaltMass = lastState.SaltMass[worldData.SurfaceWaterLayer],
			AirTemperaturePotential = staticState.GetSliceAir(nextState.AirTemperaturePotential),
			GroundCondensation = staticState.GetSliceAir(tempState.CondensationGroundMass),
			LayerMiddle = staticState.GetSliceAir(tempState.AirLayerMiddle),
			LayerCount = worldData.AirLayers - 2,
			Count = staticState.Count
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.AirLayer0]));

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
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.CloudLayer]));

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

		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(new UpdateMassAirJob()
		{
			VaporMass = staticState.GetSliceAir(nextState.AirVapor),
			DustMass = staticState.GetSliceAir(nextState.Dust),
			CarbonDioxideMass = staticState.GetSliceAir(nextState.AirCarbon),

			CloudCondensation = staticState.GetSliceAir(tempState.CondensationCloudMass),
			GroundCondensation = staticState.GetSliceAir(tempState.CondensationGroundMass),
			LastVaporMass = staticState.GetSliceAir(lastState.AirVapor),
			LastDustMass = staticState.GetSliceAir(lastState.Dust),
			LastCarbonDioxideMass = staticState.GetSliceAir(lastState.AirCarbon),
			DustUp = tempState.DustUp,
			DustDown = tempState.DustDown,
			LayerCount = worldData.AirLayers - 2,
			Count = staticState.Count

		}, energyJobHandles[worldData.AirLayer0]);

		energyJobHandles[worldData.CloudLayer] = SimJob.Schedule(new UpdateMassCloudCondensationJob()
		{
			CloudMass = nextState.CloudMass,
			CloudDropletMass = nextState.CloudDropletMass,

			CloudEvaporation = tempState.CloudEvaporationMass,
			CloudElevation = tempState.CloudElevation,
			LayerElevation = staticState.GetSliceAir(tempState.AirLayerElevation),
			LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
			CloudCondensation = staticState.GetSliceAir(tempState.CondensationCloudMass),
			GroundCondensation = staticState.GetSliceAir(tempState.CondensationGroundMass),
			LayerCount = worldData.AirLayers - 2,
			Count = staticState.Count
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.CloudLayer], energyJobHandles[worldData.AirLayer0]));
		energyJobHandles[worldData.AirLayer0] = JobHandle.CombineDependencies(energyJobHandles[worldData.AirLayer0], energyJobHandles[worldData.CloudLayer]);

		energyJobHandles[worldData.AirLayer0] = SimJob.Schedule(new UpdateMassAirSurfaceJob()
		{
			AirTemperaturePotential = staticState.GetSliceLayer(nextState.AirTemperaturePotential, worldData.SurfaceAirLayer),
			VaporMass = staticState.GetSliceLayer(nextState.AirVapor,worldData.SurfaceAirLayer),
			DustMass = staticState.GetSliceLayer(nextState.Dust,worldData.SurfaceAirLayer),
			CarbonDioxide = staticState.GetSliceLayer(nextState.AirCarbon,worldData.SurfaceAirLayer),

			AirMass = staticState.GetSliceLayer(tempState.AirMass, worldData.SurfaceAirLayer),
			EvaporationWater = tempState.EvaporationMassWater,
			EvaporationTemperatureWater = lastState.WaterTemperature[worldData.SurfaceWaterLayer],
			EvaporationFlora = tempState.FloraRespirationMassVapor,
			EvaporationTemperatureFlora = lastState.FloraTemperature,
			DustEjected = tempState.DustEjected,
			AirCarbonDelta = tempState.AirCarbonDelta,
			SoilRespiration = tempState.SoilRespiration,
			WaterCoverage = tempState.WaterCoverage[worldData.SurfaceWaterLayer],
			Elevation = lastState.Elevation,
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.AirLayer0], energyJobHandles[worldData.SurfaceWaterLayerGlobal]));

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
			DustSettled = staticState.GetSliceLayer(tempState.DustDown,worldData.SurfaceAirLayer),
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

		energyJobHandles[worldData.AirLayer0] = SimJob.Schedule(new UpdateWaterAirDiffusionJob()
		{
			AirCarbon = staticState.GetSliceLayer(nextState.AirCarbon, worldData.SurfaceAirLayer),
			WaterCarbon = nextState.WaterCarbon[worldData.SurfaceWaterLayer],

			AirMass = staticState.GetSliceLayer(tempState.AirMass, worldData.SurfaceAirLayer),
			WaterMass = nextState.WaterMass[worldData.SurfaceWaterLayer],
			SaltMass = nextState.SaltMass[worldData.SurfaceWaterLayer],
			WaterDepth = tempState.WaterLayerHeight[worldData.SurfaceWaterLayer],
			WaterAirCarbonDiffusionCoefficient = worldData.WaterAirCarbonDiffusionCoefficient,
			WaterAirCarbonDiffusionDepth = worldData.WaterAirCarbonDiffusionDepth,
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.AirLayer0], energyJobHandles[worldData.SurfaceWaterLayerGlobal]));
		energyJobHandles[worldData.SurfaceWaterLayerGlobal] = JobHandle.CombineDependencies(energyJobHandles[worldData.SurfaceWaterLayerGlobal], energyJobHandles[worldData.AirLayer0]);
#endregion


#region Apply Latent Heat

		energyJobHandles[worldData.IceLayer] = SimJob.Schedule(new ApplyLatentHeatIceJob()
		{
			IceTemperature = nextState.IceTemperature,
			IceMass = nextState.IceMass,
			LatentHeat = tempState.LatentHeatIce
		}, energyJobHandles[worldData.IceLayer]);

		energyJobHandles[worldData.TerrainLayer] = SimJob.Schedule(new ApplyLatentHeatTerrainJob()
		{
			TerrainTemperature = nextState.GroundTemperature,

			LatentHeat = tempState.LatentHeatTerrain,
			SoilFertility = nextState.GroundCarbon,
			HeatingDepth = worldData.SoilHeatDepth
		}, energyJobHandles[worldData.TerrainLayer]);

		energyJobHandles[worldData.LavaLayer] = SimJob.Schedule(new ApplyLatentHeatLavaJob()
		{
			LavaTemperature = nextState.LavaTemperature,

			LatentHeat = tempState.LatentHeatLava,
			LavaMass = nextState.LavaMass
		}, JobHandle.CombineDependencies(energyJobHandles[worldData.LavaLayer], energyJobHandles[worldData.TerrainLayer]));

		energyJobHandles[worldData.AirLayer0] = AirJob.Schedule(new ApplyLatentHeatAirJob()
		{
			AirTemperaturePotential = staticState.GetSliceAir(nextState.AirTemperaturePotential),
			AirMass = staticState.GetSliceAir(tempState.AirMass),
			VaporMass = staticState.GetSliceAir(nextState.AirVapor),
			LatentHeat = staticState.GetSliceAir(tempState.LatentHeatAir)
		}, energyJobHandles[worldData.AirLayer0]);

		energyJobHandles[worldData.SurfaceWaterLayerGlobal] = SimJob.Schedule(new ApplyLatentHeatWaterJob()
		{
			WaterTemperature = nextState.WaterTemperature[worldData.SurfaceWaterLayer],
			WaterMass = nextState.WaterMass[worldData.SurfaceWaterLayer],
			SaltMass = nextState.SaltMass[worldData.SurfaceWaterLayer],
			LatentHeat = tempState.LatentHeatWater,
		}, energyJobHandles[worldData.SurfaceWaterLayerGlobal]);

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
				SurfaceElevation = tempState.SurfaceElevation,
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
