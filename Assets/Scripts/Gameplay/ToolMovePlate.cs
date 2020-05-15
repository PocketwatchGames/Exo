using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Jobs;
using Unity.Burst;
using Unity.Collections;

public class ToolMovePlate : GameTool {

	private JobHelper _jobHelper;
	private NativeArray<int> AffectedCells;

	public override void OnSelected()
	{
		base.OnSelected();
		AffectedCells = new NativeArray<int>(Gameplay.Sim.StaticState.Count, Allocator.Persistent);
	}
	public override void OnDeselected()
	{
		base.OnDeselected();
		AffectedCells.Dispose();
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

	private void Activate(ref SimState lastState, ref SimState nextState, sbyte plateIndex, float3 direction, ref WorldData worldData, ref TempState tempState, ref SimSettings settings)
	{
		nextState.CopyFrom(ref lastState);

		if (plateIndex >= 0)
		{
			if (_jobHelper == null)
			{
				_jobHelper = new JobHelper(Gameplay.Sim.CellCount);
			}
			_jobHelper.Schedule(settings.SynchronousOverrides.Tools, 1, new MovePlateJob()
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
				LastFloraGlucose = lastState.FloraGlucose

			}).Complete();
		}
	}

	private struct MovePlateJob : IJobParallelFor {

		public NativeArray<sbyte> Plates;
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

		[ReadOnly] public sbyte PlateIndex;
		[ReadOnly] public float3 Direction;
		[ReadOnly] public NativeArray<sbyte> LastPlates;
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
		[ReadOnly] public NativeArray<int> Neighbors;
		[ReadOnly] public NativeArray<float3> NeighborDir;

		[BurstCompile]
		public void Execute(int i)
		{
			float bestNeighborDot = -1;
			int bestNeighbor = -1;
			float bestRearNeighborDot = 1;
			int bestRearNeighbor = -1;
			int maxNeighbors = StaticState.GetMaxNeighbors(i, Neighbors);
			for (int j = 0; j < maxNeighbors; j++)
			{
				int n = i * StaticState.MaxNeighbors + j;
				float dot = math.dot(Direction, NeighborDir[n]);
				if (dot >= bestNeighborDot)
				{
					bestNeighborDot = dot;
					bestNeighbor = Neighbors[n];
				}
				if (dot <= bestRearNeighborDot)
				{
					bestRearNeighborDot = dot;
					bestRearNeighbor = Neighbors[n];
				}
			}
			if (bestNeighbor >= 0 && LastPlates[bestNeighbor] == PlateIndex)
			{
				// within same plate
				if (LastPlates[i] == PlateIndex)
				{
					Plates[i] = LastPlates[bestNeighbor];
					Elevation[i] = LastElevation[bestNeighbor];
					GroundTemperature[i] = LastGroundTemperature[bestNeighbor];
					GroundCarbon[i] = LastGroundCarbon[bestNeighbor];
					GroundWater[i] = LastGroundWater[bestNeighbor];
					GroundWaterTemperature[i] = LastGroundWaterTemperature[bestNeighbor];
					FloraMass[i] = LastFloraMass[bestNeighbor];
					FloraGlucose[i] = LastFloraGlucose[bestNeighbor];
					FloraWater[i] = LastFloraWater[bestNeighbor];
				} else // moving into new plate
				{
					Plates[i] = LastPlates[bestNeighbor];
					Elevation[i] = math.max(LastElevation[bestNeighbor], LastElevation[i]);
					GroundTemperature[i] = (LastGroundTemperature[bestNeighbor] + LastGroundTemperature[i]) / 2;
					GroundCarbon[i] = LastGroundCarbon[bestNeighbor] + LastGroundCarbon[i];
					GroundWater[i] = LastGroundWater[bestNeighbor] + LastGroundWater[i];
					GroundWaterTemperature[i] = (LastGroundWaterTemperature[bestNeighbor] * LastGroundWater[bestNeighbor] + LastGroundWaterTemperature[i] * LastGroundWaterTemperature[i]) / (GroundWater[i]);
					FloraMass[i] = LastFloraMass[bestNeighbor] + LastFloraMass[i];
					FloraGlucose[i] = LastFloraGlucose[bestNeighbor] + LastFloraGlucose[i];
					FloraWater[i] = LastFloraWater[bestNeighbor] + LastFloraWater[i];
				}
			}
			else if (LastPlates[i] == PlateIndex && bestRearNeighbor >= 0 && LastPlates[bestRearNeighbor] != PlateIndex)
			{
				Elevation[i] -= CrustDepth[i];
				Plates[i] = LastPlates[bestRearNeighbor];
				GroundCarbon[i] = 0;
				GroundWater[i] = 0;
				FloraMass[i] = 0;
				FloraGlucose[i] = 0;
				FloraWater[i] = 0;
			}
		}
	}
}
