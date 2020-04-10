//#define DISABLE_SURFACE_FLOW

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;


[BurstCompile]
public struct UpdateFlowVelocityJob : IJobParallelFor {
	public NativeArray<float> Flow;

	[ReadOnly] public NativeArray<float> LastFlow;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float> WaterDepth;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float Damping;

	public void Execute(int i)
	{
		int nIndex = Neighbors[i];
		if (nIndex < 0)
		{
			return;
		}
		int cellIndex = i / StaticState.MaxNeighbors;
		float waterDepth = WaterDepth[cellIndex];
		float elevationDiff = SurfaceElevation[cellIndex] - SurfaceElevation[nIndex];
		int upwindIndex = elevationDiff > 0 ? cellIndex : nIndex;
		float acceleration = Gravity * elevationDiff * NeighborDistInverse[i];

		float v = LastFlow[i] * Damping + acceleration * SecondsPerTick;
		Flow[i] = v;
	}
}

[BurstCompile]
public struct SumOutgoingFlowJob : IJobParallelFor {
	public NativeArray<float> OutgoingFlow;
	[ReadOnly] public NativeArray<float> Flow;
	public void Execute(int i)
	{
		float outgoingFlow = 0;
		for (int j = 0; j < StaticState.MaxNeighbors; j++)
		{
			float f = Flow[i * StaticState.MaxNeighbors + j];
			outgoingFlow += math.max(0, f);
		}
		OutgoingFlow[i] = outgoingFlow;
	}
}

[BurstCompile]
public struct LimitOutgoingFlowJob : IJobParallelFor {
	public NativeArray<float> Flow;
	public NativeArray<float> FlowPercent;
	[ReadOnly] public NativeArray<float> OutgoingFlow;
	[ReadOnly] public NativeArray<float> WaterDepth;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		int nIndex = Neighbors[i];
		if (nIndex < 0)
		{
			return;
		}
		float flow = Flow[i];
		if (flow > 0) {
			int cellIndex = i / StaticState.MaxNeighbors;
			float waterDepth = WaterDepth[cellIndex];
			if (waterDepth > 0)
			{
				float outgoing = OutgoingFlow[cellIndex];
				if (outgoing > 0)
				{
					Flow[i] = flow * math.min(1, waterDepth / outgoing);
					FlowPercent[i] = Flow[i] / waterDepth;
				}
			}
			else
			{
				Flow[i] = 0;
				FlowPercent[i] = 0;
			}
		}
		else
		{
			float waterDepth = WaterDepth[nIndex];
			if (waterDepth > 0)
			{
				float outgoing = OutgoingFlow[nIndex];
				if (outgoing > 0)
				{
					Flow[i] = flow * math.min(1, waterDepth / outgoing);
					FlowPercent[i] = Flow[i] / waterDepth;
				}
			}
			else
			{
				Flow[i] = 0;
				FlowPercent[i] = 0;
			}
		}
	}
}


[BurstCompile]
public struct ApplyFlowJob : IJobParallelFor {

	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float> Carbon;
	[ReadOnly] public NativeArray<float> PlanktonMass;
	[ReadOnly] public NativeArray<float> PlanktonGlucose;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> FlowPercent;

	public void Execute(int i)
	{
		float mass = 0;
		float salt = 0;
		float carbon = 0;
		float planktonMass = 0;
		float planktonGlucose = 0;
		float3 velocity = 0;
		float temperature = 0;
		float massPercentRemaining = 1;

#if !DISABLE_SURFACE_FLOW

		for (int j = 0; j < StaticState.MaxNeighbors; j++)
		{
			int n = i * StaticState.MaxNeighbors + j;
			int nIndex = Neighbors[n];
			if (nIndex >= 0)
			{
				float flowPercent = -FlowPercent[n];
				if (flowPercent < 0)
				{
					massPercentRemaining += flowPercent;
				}
				else
				{
					// TODO: cache the neighbor's neighbor inverse
					for (int k = 0; k < StaticState.MaxNeighbors; k++)
					{
						int incomingNIndex = nIndex * StaticState.MaxNeighbors + k;
						if (Neighbors[incomingNIndex] == i)
						{
							float massIncoming = Mass[nIndex] * flowPercent;
							float saltIncoming = Salt[nIndex] * flowPercent;
							mass += massIncoming;
							salt += saltIncoming;
							carbon += Carbon[nIndex] * flowPercent;
							planktonMass += PlanktonMass[nIndex] * flowPercent;
							planktonGlucose += PlanktonGlucose[nIndex] * flowPercent;
							velocity += Velocity[nIndex] * (massIncoming + saltIncoming);
							temperature += Temperature[nIndex] * (massIncoming + saltIncoming);

							break;
						}
					}
				}
			}
		}

		massPercentRemaining = math.max(0, massPercentRemaining);
		float massRemaining = Mass[i] * massPercentRemaining;
		float saltRemaining = Salt[i] * massPercentRemaining;
		mass += massRemaining;
		salt += saltRemaining;
		carbon += Carbon[i] * massPercentRemaining;
		planktonMass += PlanktonMass[i] * massPercentRemaining;
		planktonGlucose += PlanktonGlucose[i] * massPercentRemaining;
		velocity += Velocity[i] * (massRemaining + saltRemaining);
		temperature += Temperature[i] * (massRemaining + saltRemaining);

		if (mass == 0)
		{
			temperature = 0;
			velocity = 0;
			planktonMass = 0;
			planktonGlucose = 0;
		} else
		{
			float inverseMass = 1.0f / (mass + salt);
			temperature *= inverseMass;
			velocity *= inverseMass; 
		}

#endif

		Delta[i] = new DiffusionWater()
		{
			WaterMass = mass,
			SaltMass = salt,
			CarbonMass = carbon,
			Plankton = planktonMass,
			PlanktonGlucose = planktonGlucose,
			Temperature = temperature,
			Velocity = velocity
		};
	}
}


