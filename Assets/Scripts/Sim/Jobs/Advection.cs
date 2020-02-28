﻿//#define DISABLE_VERTICAL_AIR_MOVEMENT
#define DISABLE_AIR_ADVECTION
#define DISABLE_WATER_ADVECTION
//#define AdvectionAirJobDebug
//#define AdvectionCloudJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

public struct DiffusionAir {
	public float Temperature;
	public float WaterVapor;
	public float3 Velocity;
}
public struct DiffusionCloud {
	public float Mass;
	public float DropletMass;
	public float3 Velocity;
}
public struct DiffusionWater {
	public float Temperature;
	public float SaltMass;
	public float3 Velocity;
}


#if !AdvectionAirJobDebug
[BurstCompile]
#endif
public struct AdvectionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float3> Wind;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpHumidity;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> DownTemperature;
	[ReadOnly] public NativeArray<float> DownHumidity;
	[ReadOnly] public NativeArray<float> DownAirMass;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	public void Execute(int i)
	{

		float3 velocity = Wind[i];
		float airMass = AirMass[i];
		float vapor = Vapor[i];
		float temperature = Temperature[i];
		float absoluteHumidity = vapor / (vapor + airMass);
		float3 pos = Position[i];

		float leavingCell = -math.length(velocity - math.cross(velocity, pos)) * InverseCellDiameter * SecondsPerTick;

		// TODO: deal with high speeds here
		leavingCell = math.max(leavingCell, -1);

		float newTemperature = leavingCell * Temperature[i];
		float newWaterVapor = leavingCell * Vapor[i];
		float3 newWind = leavingCell * Wind[i];

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				var nVel = Wind[n];
				float speed = math.length(nVel);
				if (speed > 0)
				{
					// TODO: deal with high speeds here
					if (speed * InverseCellDiameter * SecondsPerTick > 1)
					{
						nVel *= InverseCellDiameter * SecondsPerTick / speed;
					}

					// TODO: this needs to account for cells with 5 neighbors
					const float cosAngleBetweenCells = 0.5f;
					float velDotDir = math.max(0, (math.dot(nVel / speed, math.normalize(pos - Position[n])) - cosAngleBetweenCells) / (1.0f - cosAngleBetweenCells));
					velDotDir *= speed * InverseCellDiameter * SecondsPerTick;

					newWind += nVel * velDotDir;
					newWaterVapor += Vapor[n] * velDotDir;
					newTemperature += Temperature[n] * velDotDir;
				}
			}
		}


#if !DISABLE_VERTICAL_AIR_MOVEMENT

		if (!IsTop)
		{
		}
		if (!IsBottom)
		{
		}


#endif

#if !DISABLE_AIR_ADVECTION
		Delta[i] = new DiffusionAir()
		{
			Temperature = newTemperature,
			WaterVapor = newWaterVapor,
			Velocity = newWind,
		};
#endif



	}
}


#if !AdvectionCloudJobDebug
[BurstCompile]
#endif
public struct AdvectionCloudJob : IJobParallelFor {
	public NativeArray<DiffusionCloud> Delta;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> DropletMass;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	public void Execute(int i)
	{

		float3 velocity = Velocity[i];
		float3 pos = Position[i];
		float mass = Mass[i];
		float dropletMass = DropletMass[i];

		float leavingCell = -math.length(velocity - math.cross(velocity, pos)) * InverseCellDiameter * SecondsPerTick;

		// TODO: deal with high speeds here
		leavingCell = math.max(leavingCell, -1);

		float newMass = leavingCell * mass;
		float newDropletMass = leavingCell * dropletMass;
		float3 newVelocity = leavingCell * velocity;

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				var nVel = Velocity[n];
				float speed = math.length(nVel);
				if (speed > 0)
				{
					// TODO: deal with high speeds here
					if (speed * InverseCellDiameter * SecondsPerTick > 1)
					{
						nVel *= InverseCellDiameter * SecondsPerTick / speed;
					}

					// TODO: this needs to account for cells with 5 neighbors
					const float cosAngleBetweenCells = 0.5f;
					float velDotDir = math.max(0, (math.dot(nVel / speed, math.normalize(pos - Position[n])) - cosAngleBetweenCells) / (1.0f - cosAngleBetweenCells));
					velDotDir *= speed * InverseCellDiameter * SecondsPerTick;

					newMass += Mass[n] * velDotDir;
					newDropletMass += DropletMass[n] * velDotDir;
					newVelocity += Velocity[n] * velDotDir;
				}
			}
		}

#if !DISABLE_CLOUD_ADVECTION
		Delta[i] = new DiffusionCloud()
		{
			Mass = newMass,
			DropletMass = newDropletMass,
			Velocity = newVelocity
		};
#endif
	}
}

#if !AdvectionWaterJobDebug
[BurstCompile]
#endif
public struct AdvectionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float3> Current;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float InverseCoordDiff;
	public void Execute(int i)
	{
		float3 velocity = Current[i];
		float waterMass = Mass[i];
		float salt = Salt[i];
		float mass = Mass[i];
		float temperature = Temperature[i];

		float leavingCell = -math.length(velocity) * InverseCellDiameter * SecondsPerTick;
		float newSaltMass = leavingCell * salt;
		float newTemperature = leavingCell * temperature;
		float3 newVelocity = leavingCell * velocity;
		float3 pos = Position[i];

		for (int j = 0; j < 6; j++)
		{
			int n = Neighbors[i * 6 + j];
			if (n >= 0)
			{
				float nMass = Mass[n];
				if (nMass > 0)
				{
					// TODO: this needs to account for cells with 5 neighbors
					const float cosAngleBetweenCells = 0.5f;
					float velDotDir = math.max(0, (math.dot(Current[n], math.normalize(pos - Position[n])) - cosAngleBetweenCells) / (1.0f - cosAngleBetweenCells)) * InverseCellDiameter * SecondsPerTick;

					newSaltMass += Salt[n] * velDotDir;
					newTemperature += Temperature[n] * velDotDir;
					newVelocity += Current[n] * velDotDir;
				}
			}
		}


#if !DISABLE_WATER_ADVECTION
		Delta[i] = new DiffusionWater()
		{
			Temperature = newTemperature,
			SaltMass = newSaltMass,
			Velocity = newVelocity
		};
#endif

	}
}


#if !ApplyAdvectionAirJobDebug
[BurstCompile]
#endif
public struct ApplyAdvectionAirJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Vapor;
	public NativeArray<float3> Wind;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	public void Execute(int i)
	{
		Temperature[i] += Advection[i].Temperature;
		Vapor[i] += Advection[i].WaterVapor;
		Wind[i] += Advection[i].Velocity;
	}
}


#if !ApplyAdvectionWaterJobDebug
[BurstCompile]
#endif
public struct ApplyAdvectionWaterJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> SaltMass;
	public NativeArray<float3> Velocity;
	public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<DiffusionWater> Advection;
	public void Execute(int i)
	{
		SaltMass[i] += Advection[i].SaltMass;
		Temperature[i] += Advection[i].Temperature;
		Velocity[i] += Advection[i].Velocity;
	}
}


#if !ApplyAdvectionCloudJobDebug
[BurstCompile]
#endif
public struct ApplyAdvectionCloudJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> DropletMass;
	public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<DiffusionCloud> Advection;
	public void Execute(int i)
	{
		CloudMass[i] += Advection[i].Mass;
		DropletMass[i] += Advection[i].DropletMass;
		Velocity[i] += Advection[i].Velocity;
	}
}

