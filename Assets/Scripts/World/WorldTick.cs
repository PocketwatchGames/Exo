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

	public struct CellDiffusion {
		public float Water;
		public float Cloud;
		public float Humidity;
	}

	public static void Tick(ref SimState state, ref SimState nextState, int ticksToAdvance, StaticState staticState, WorldData worldData)
	{
		JobHandle lastJobHandle = default(JobHandle);

		var diffusion = new NativeArray<CellDiffusion>(state.Count * 6, Allocator.TempJob);
		var diffusionLimit = new NativeArray<CellDiffusion>(state.Count, Allocator.TempJob);
		var cells = new NativeArray<SimStateCell>[2];
		cells[0] = new NativeArray<SimStateCell>(state.Cells, Allocator.TempJob);
		cells[1] = new NativeArray<SimStateCell>(state.Cells.Length, Allocator.TempJob);
		int ticks = state.Ticks;
		for (int i = 0; i < ticksToAdvance; i++)
		{
			ticks++;
			int lastStateIndex = i % 2;
			var diffusionJob = new DiffusionJob();
			diffusionJob.Last = cells[lastStateIndex];
			diffusionJob.CellDiffusion = diffusion;
			diffusionJob.Neighbors = staticState.Neighbors;
			diffusionJob.WaterDiffuseSpeed = worldData.WaterDiffuseSpeed;
			lastJobHandle = diffusionJob.Schedule(state.Count * 6, 100, lastJobHandle);

			var diffusionLimitJob = new DiffusionLimitJob();
			diffusionLimitJob.Last = cells[lastStateIndex];
			diffusionLimitJob.CellDiffusion = diffusion;
			diffusionLimitJob.DiffusionLimit = diffusionLimit;
			lastJobHandle = diffusionLimitJob.Schedule(state.Count, 100, lastJobHandle);


			var tickJob = new TickCellJob();
			tickJob.Cells = cells[(i + 1) % 2];
			tickJob.Last = cells[lastStateIndex];
			tickJob.Neighbors = staticState.Neighbors;
			tickJob.Diffusion = diffusion;
			tickJob.DiffusionLimit = diffusionLimit;
			tickJob.Ticks = ticks;
			lastJobHandle = tickJob.Schedule(state.Count, 100, lastJobHandle);
		}
		lastJobHandle.Complete();

		int outputBuffer = ticksToAdvance % 2;

		nextState.Ticks = ticks;
		nextState.Gravity = state.Gravity;
		nextState.OrbitSpeed = state.OrbitSpeed;
		nextState.SpinAngle = state.SpinAngle;
		nextState.SpinSpeed = state.SpinSpeed;
		nextState.TiltAngle = state.TiltAngle;
		for (int i = 0; i < nextState.Count; i++)
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

	[BurstCompile]
	public struct DiffusionJob : IJobParallelFor {

		[ReadOnly] public NativeArray<int> Neighbors;
		[ReadOnly] public NativeArray<SimStateCell> Last;
		public NativeArray<CellDiffusion> CellDiffusion;
		public float WaterDiffuseSpeed;

		public void Execute(int i)
		{
			int n = Neighbors[i];
			if (n >= 0)
			{
				int index = i / 6;
				var curCell = Last[index];
				float water = 0;
				if (curCell.WaterDepth > 0)
				{
					float waterElevation = curCell.Elevation + curCell.WaterDepth;
					float nWaterElevation = Last[n].Elevation + Last[n].WaterDepth;
					if (waterElevation > nWaterElevation)
					{
						water = (waterElevation - nWaterElevation) * WaterDiffuseSpeed;
					}
				}
				float cloud = math.max(0, (curCell.CloudCoverage - Last[n].CloudCoverage) * WaterDiffuseSpeed);
				float humidity = math.max(0, (curCell.RelativeHumidity - Last[n].RelativeHumidity) * WaterDiffuseSpeed);

				CellDiffusion[i] = new CellDiffusion { Water = water, Cloud = cloud, Humidity = humidity };
			}
		}
	}
	[BurstCompile]
	public struct DiffusionLimitJob : IJobParallelFor {

		public NativeArray<CellDiffusion> DiffusionLimit;
		[ReadOnly] public NativeArray<CellDiffusion> CellDiffusion;
		[ReadOnly] public NativeArray<SimStateCell> Last;

		public void Execute(int i)
		{
			float totalWater = 0;
			float totalCloud = 0;
			float totalHumidity = 0;
			float waterLimit, cloudLimit, humidityLimit;
			for (int j = 0; j < 6; j++)
			{
				int index = i * 6 + j;
				totalWater += CellDiffusion[index].Water;
				totalCloud += CellDiffusion[index].Cloud;
				totalHumidity += CellDiffusion[index].Humidity;
			}
			float waterDepth = Last[i].WaterDepth;
			float cloud = Last[i].CloudCoverage;
			float humidity = Last[i].RelativeHumidity;
			if (totalWater > waterDepth && waterDepth > 0)
			{
				waterLimit = waterDepth / totalWater;
			}
			else
			{
				waterLimit = 1;
			}
			if (totalCloud > cloud && cloud > 0)
			{
				cloudLimit = cloud / totalCloud;
			}
			else
			{
				cloudLimit = 1;
			}
			if (totalHumidity > humidity && humidity > 0)
			{
				humidityLimit = humidity / totalHumidity;
			}
			else
			{
				humidityLimit = 1;
			}
			DiffusionLimit[i] = new CellDiffusion { Water = waterLimit, Cloud = cloudLimit, Humidity = humidityLimit };
		}
	}
	[BurstCompile]
	public struct TickCellJob : IJobParallelFor {

		public int Ticks;
		public float Gravity;
		public float SpinAngle;
		public float SpinSpeed;
		public float OrbitSpeed;
		public float TiltAngle;
		[ReadOnly] public NativeArray<CellDiffusion> Diffusion;
		[ReadOnly] public NativeArray<CellDiffusion> DiffusionLimit;
		[ReadOnly] public NativeArray<int> Neighbors;
		[ReadOnly] public NativeArray<SimStateCell> Last;

		public NativeArray<SimStateCell> Cells;

		public void Execute(int i)
		{
			var curCell = Last[i];
			var cell = new SimStateCell();
			cell.CloudCoverage = curCell.CloudCoverage;
			cell.RelativeHumidity = curCell.RelativeHumidity;
			cell.WaterDepth = curCell.WaterDepth;
			cell.CloudElevation = curCell.CloudElevation;
			cell.Elevation = curCell.Elevation;
			cell.Roughness = curCell.Roughness;
			cell.IceMass = curCell.IceMass;
			cell.Vegetation = curCell.Vegetation;
			cell.SoilFertility = curCell.SoilFertility;
			cell.WaterMass = curCell.WaterMass;
			cell.WaterEnergy = curCell.WaterEnergy;
			cell.SaltMass = curCell.SaltMass;
			cell.GroundEnergy = curCell.GroundEnergy;
			cell.GroundWater = curCell.GroundWater;
			cell.GroundWaterDepth = curCell.GroundWaterDepth;
			cell.AirMass = curCell.AirMass;
			cell.AirEnergy = curCell.AirEnergy;
			cell.CloudMass = curCell.CloudMass;
			cell.CloudDropletMass = curCell.CloudDropletMass;


			cell.AirWaterMass = curCell.AirWaterMass;
			cell.AirTemperature = curCell.AirTemperature;
			cell.AirPressure = curCell.AirPressure;
			cell.WaterTemperature = curCell.WaterTemperature;
			cell.WindVertical = curCell.WindVertical;
			cell.WindSurface = curCell.WindSurface;
			cell.WindTropopause = curCell.WindTropopause;


			float evap = math.min(curCell.WaterDepth, 1f) * math.clamp(1.0f - curCell.RelativeHumidity / 100, 0, 1);
			cell.WaterDepth -= evap;
			cell.RelativeHumidity += evap;
			float hToC = curCell.RelativeHumidity / 100;
			cell.RelativeHumidity -= hToC;
			cell.CloudCoverage += hToC;
			if (cell.CloudCoverage > 1000)
			{
				cell.WaterDepth += curCell.CloudCoverage;
				cell.CloudCoverage = 0;
			}
			for (int j = 0; j < 6; j++)
			{
				int neighborIndex = i * 6 + j;
				int n = Neighbors[neighborIndex];
				if (n >= 0)
				{
					cell.CloudCoverage -= Diffusion[neighborIndex].Cloud * DiffusionLimit[i].Cloud;
					cell.RelativeHumidity -= Diffusion[neighborIndex].Humidity * DiffusionLimit[i].Humidity;
					cell.WaterDepth -= Diffusion[neighborIndex].Water * DiffusionLimit[i].Water;
					for (int k=0;k<6;k++)
					{
						int nToMeIndex = n * 6 + k;
						if (Neighbors[nToMeIndex] == i)
						{
							cell.CloudCoverage += Diffusion[nToMeIndex].Cloud * DiffusionLimit[n].Cloud;
							cell.RelativeHumidity += Diffusion[nToMeIndex].Humidity * DiffusionLimit[n].Humidity;
							cell.WaterDepth += Diffusion[nToMeIndex].Water * DiffusionLimit[n].Water;
							break;
						}
					}
				}
			}

			Cells[i] = cell;

		}
	}
}
