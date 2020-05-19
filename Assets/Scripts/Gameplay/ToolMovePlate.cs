using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Jobs;
using Unity.Burst;
using Unity.Collections;

public class ToolMovePlate : GameTool {

	private JobHelper _jobHelper;
	private NativeArray<bool> Closed;

	public override void OnSelected()
	{
		base.OnSelected();
		Closed = new NativeArray<bool>(Gameplay.Sim.StaticState.Count, Allocator.Persistent);
	}
	public override void OnDeselected()
	{
		base.OnDeselected();
		Closed.Dispose();
	}
	public override void OnPointerDown(Vector3 worldPos, int cellIndex)
	{
		base.OnPointerDown(worldPos, cellIndex);
		Gameplay.ActivePlateIndex = Gameplay.Sim.LastSimState.Plate[cellIndex];
		Gameplay.SetActiveCell(cellIndex, false);

	}
	override public void OnDragMove(Vector3 worldPos, int cellIndex, Vector2 direction)
	{
		if (cellIndex != Gameplay.ActiveCellIndex && Gameplay.ActiveCellIndex >= 0 && cellIndex >= 0)
		{
			float3 worldDir = math.normalize(Gameplay.Sim.StaticState.SphericalPosition[cellIndex] - Gameplay.Sim.StaticState.SphericalPosition[Gameplay.ActiveCellIndex]);

			Gameplay.Sim.Edit((ref SimState last, ref SimState next) => { Activate(ref last, ref next, Gameplay.ActivePlateIndex, worldDir, ref Gameplay.Sim.WorldData, ref Gameplay.Sim.LastTempState, ref Gameplay.Sim.SimSettings); });
		}
		Gameplay.SetActiveCell(cellIndex, false);
	}
	public override void OnUpdate(Vector3 worldPos, int cellIndex)
	{
		base.OnUpdate(worldPos, cellIndex);
		Gameplay.SetActiveCell(cellIndex, false);
	}

	private void Activate(ref SimState lastState, ref SimState nextState, short plateIndex, float3 direction, ref WorldData worldData, ref TempState tempState, ref SimSettings settings)
	{
		nextState.CopyFrom(ref lastState);

		if (plateIndex >= 0)
		{
			if (_jobHelper == null)
			{
				_jobHelper = new JobHelper(Gameplay.Sim.CellCount);
			}
			var movePlateJob = _jobHelper.Schedule(settings.SynchronousOverrides.Tools, 1, new MovePlateJob()
			{
				Plates = nextState.Plate,
				Elevation = nextState.Elevation,
				Roughness = nextState.Roughness,
				CrustDepth = nextState.CrustDepth,
				GroundTemperature = nextState.GroundTemperature,
				GroundCarbon = nextState.GroundCarbon,
				GroundWater = nextState.GroundWater,
				GroundWaterTemperature = nextState.GroundWaterTemperature,
				FloraMass = nextState.FloraMass,
				FloraWater = nextState.FloraWater,
				FloraGlucose = nextState.FloraGlucose,

				PlateIndex = plateIndex,
				Direction = direction,
				NeighborDir = Gameplay.Sim.StaticState.NeighborDir,
				Neighbors = Gameplay.Sim.StaticState.Neighbors,
				LastPlates = lastState.Plate,
				LastElevation = lastState.Elevation,
				LastRoughness = lastState.Roughness,
				LastCrustDepth = lastState.CrustDepth,
				LastGroundTemperature = lastState.GroundTemperature,
				LastGroundCarbon = lastState.GroundCarbon,
				LastGroundWater = lastState.GroundWater,
				LastGroundWaterTemperature = lastState.GroundWaterTemperature,
				LastFloraMass = lastState.FloraMass,
				LastFloraWater = lastState.FloraWater,
				LastFloraGlucose = lastState.FloraGlucose,
				LastWaterDepth = tempState.WaterDepthTotal

			});

			var clearJob = Utils.MemsetArray<bool>(Gameplay.Sim.StaticState.Count, default(JobHandle), Closed, false);
			JobHandle.CombineDependencies(movePlateJob, clearJob).Complete();


			for (int i = 0; i < Gameplay.Sim.StaticState.Count; i++)
			{
				bool erode = false;
				int maxNeighbors = StaticState.GetMaxNeighbors(i, Gameplay.Sim.StaticState.Neighbors);
				for (int j = 0; j < maxNeighbors; j++)
				{
					int n = i * StaticState.MaxNeighbors + j;
					int neighbor = Gameplay.Sim.StaticState.Neighbors[n];
					if (nextState.Plate[neighbor] != plateIndex)
					{
						erode = (nextState.Plate[i] == plateIndex) != (nextState.Plate[neighbor] == plateIndex);
						break;
					}
				}

				if (erode)
				{
					for (int j = 0; j < maxNeighbors; j++)
					{
						int n = i * StaticState.MaxNeighbors + j;
						int neighbor = Gameplay.Sim.StaticState.Neighbors[n];
						float diff = nextState.Elevation[i] - nextState.Elevation[neighbor];
						if (diff > 0)
						{
							nextState.Elevation[i] -= diff / 3;
							nextState.Elevation[neighbor] += diff / 3;
						}
					}
				}
			}


			HashSet<int> allPlates = new HashSet<int>();
			for (int i = 0; i < Gameplay.Sim.StaticState.Count; i++)
			{
				allPlates.Add(nextState.Plate[i]);
			}

			HashSet<int> plates = new HashSet<int>();
			for (int i = 0; i < Gameplay.Sim.StaticState.Count; i++)
			{
				Queue<int> openList = new Queue<int>();
				if (!Closed[i])
				{
					short curPlate = nextState.Plate[i];
					short plateToSet = curPlate;
					if (plates.Contains(curPlate))
					{
						// Find new plate
						plateToSet = 0;
						while (allPlates.Contains(plateToSet) && plateToSet < Gameplay.Sim.StaticState.Count)
						{
							plateToSet++;
						}
					}
					plates.Add(plateToSet);

					openList.Enqueue(i);
					Closed[i] = true;
					while (openList.Count > 0)
					{
						int index = openList.Dequeue();
						int plate = nextState.Plate[index];
						nextState.Plate[index] = plateToSet;

						int maxNeighbors = StaticState.GetMaxNeighbors(index, Gameplay.Sim.StaticState.Neighbors);
						for (int j = 0; j < maxNeighbors; j++)
						{
							int n = Gameplay.Sim.StaticState.Neighbors[index * StaticState.MaxNeighbors + j];
							if (nextState.Plate[n] == plate && !Closed[n])
							{
								Closed[n] = true;
								openList.Enqueue(n);
							}
						}
					}
				}

			}

		}
	}

	private struct MovePlateJob : IJobParallelFor {

		public NativeArray<short> Plates;
		public NativeArray<float> Elevation;
		public NativeArray<float> Roughness;
		public NativeArray<float> CrustDepth;
		public NativeArray<float> GroundTemperature;
		public NativeArray<float> GroundCarbon;
		public NativeArray<float> GroundWater;
		public NativeArray<float> GroundWaterTemperature;
		public NativeArray<float> FloraMass;
		public NativeArray<float> FloraGlucose;
		public NativeArray<float> FloraWater;

		[ReadOnly] public short PlateIndex;
		[ReadOnly] public float3 Direction;
		[ReadOnly] public NativeArray<short> LastPlates;
		[ReadOnly] public NativeArray<float> LastElevation;
		[ReadOnly] public NativeArray<float> LastRoughness;
		[ReadOnly] public NativeArray<float> LastCrustDepth;
		[ReadOnly] public NativeArray<float> LastGroundTemperature;
		[ReadOnly] public NativeArray<float> LastGroundCarbon;
		[ReadOnly] public NativeArray<float> LastGroundWater;
		[ReadOnly] public NativeArray<float> LastGroundWaterTemperature;
		[ReadOnly] public NativeArray<float> LastFloraMass;
		[ReadOnly] public NativeArray<float> LastFloraGlucose;
		[ReadOnly] public NativeArray<float> LastFloraWater;
		[ReadOnly] public NativeArray<float> LastWaterDepth;
		[ReadOnly] public NativeArray<int> Neighbors;
		[ReadOnly] public NativeArray<float3> NeighborDir;

		[BurstCompile]
		public void Execute(int i)
		{
			float fromDot = 1;
			int from = -1;
			float bestRearNeighborDot = -1;
			int bestRearNeighbor = -1;
			int maxNeighbors = StaticState.GetMaxNeighbors(i, Neighbors);
			for (int j = 0; j < maxNeighbors; j++)
			{
				int n = i * StaticState.MaxNeighbors + j;
				int neighbor = Neighbors[n];
				float dot = math.dot(Direction, NeighborDir[n]);
				if (dot <= fromDot)
				{
					fromDot = dot;
					from = neighbor;
				}
				if (dot >= bestRearNeighborDot)
				{
					bestRearNeighborDot = dot;
					bestRearNeighbor = Neighbors[n];
				}
			}
			if (LastPlates[i] == PlateIndex)
			{
				if (from >= 0 && LastPlates[from] == PlateIndex)
				{
					// within same plate
					Plates[i] = LastPlates[from];
					Elevation[i] = LastElevation[from];
					GroundTemperature[i] = LastGroundTemperature[from];
					GroundCarbon[i] = LastGroundCarbon[from];
					GroundWater[i] = LastGroundWater[from];
					GroundWaterTemperature[i] = LastGroundWaterTemperature[from];
					FloraMass[i] = LastFloraMass[from];
					FloraGlucose[i] = LastFloraGlucose[from];
					FloraWater[i] = LastFloraWater[from];
				}
				else
				{
					// trailing edge

					// Divergent				
					Elevation[i] = math.min(LastElevation[i], LastElevation[from]) - 100;
					Plates[i] = LastPlates[from];
					GroundCarbon[i] = 0;
					GroundWater[i] = 0;
					FloraMass[i] = 0;
					FloraGlucose[i] = 0;
					FloraWater[i] = 0;
				}
			}
			else if (from >= 0 && LastPlates[from] == PlateIndex)
			{
				// Leading edge
				GroundTemperature[i] = (LastGroundTemperature[from] + LastGroundTemperature[i]) / 2;
				GroundCarbon[i] = LastGroundCarbon[from] + LastGroundCarbon[i];
				GroundWater[i] = LastGroundWater[from] + LastGroundWater[i];
				GroundWaterTemperature[i] = (LastGroundWaterTemperature[from] * LastGroundWater[from] + LastGroundWaterTemperature[i] * LastGroundWaterTemperature[i]) / (GroundWater[i]);
				FloraMass[i] = LastFloraMass[from] + LastFloraMass[i];
				FloraGlucose[i] = LastFloraGlucose[from] + LastFloraGlucose[i];
				FloraWater[i] = LastFloraWater[from] + LastFloraWater[i];

				bool toContinental = LastWaterDepth[i] < LastRoughness[i];
				bool fromContinental = LastWaterDepth[from] < LastRoughness[from];
				if (LastElevation[from] > LastElevation[i])
				{
					Plates[i] = LastPlates[from];
					if (!toContinental)
					{
						// subduction
						Elevation[i] = (LastElevation[from] + LastElevation[i]) / 2 - 10;
					}
					else
					{
						// Convergent continental plates
						// upthrust
						Elevation[i] = LastElevation[from] + 100;
					}
				}
				else
				{
					Plates[i] = LastPlates[i];
					if (!toContinental)
					{
						// Subduction
						Elevation[i] = (LastElevation[from] + LastElevation[i]) / 2 - 10;
					} else
					{
						// Convergent continental plates
						// upthrust
						Elevation[i] = LastElevation[i] + 100;
					}
				}
			}
		}
	}
}
