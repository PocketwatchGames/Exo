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
	[ReadOnly] public float ViscosityInverse;

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

		float v = LastFlow[i] * Damping + acceleration * SecondsPerTick * ViscosityInverse;
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
public struct ApplyFlowWaterJob : IJobParallelFor {

	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float> Carbon;
	[ReadOnly] public NativeArray<float> PlanktonMass;
	[ReadOnly] public NativeArray<float> PlanktonGlucose;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float3> Velocity;
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
		}
		else
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


[BurstCompile]
public struct RebalanceWaterLayersLimitJob : IJobParallelFor {

	public NativeArray<DiffusionWater> Delta1;
	public NativeArray<DiffusionWater> Delta2;
	[ReadOnly] public NativeArray<float> Mass1;
	[ReadOnly] public NativeArray<float> Salt1;
	[ReadOnly] public NativeArray<float> Carbon1;
	[ReadOnly] public NativeArray<float> PlanktonMass1;
	[ReadOnly] public NativeArray<float> PlanktonGlucose1;
	[ReadOnly] public NativeArray<float> Temperature1;
	[ReadOnly] public NativeArray<float3> Velocity1;
	[ReadOnly] public NativeArray<float> Mass2;
	[ReadOnly] public NativeArray<float> Salt2;
	[ReadOnly] public NativeArray<float> Carbon2;
	[ReadOnly] public NativeArray<float> Temperature2;
	[ReadOnly] public NativeArray<float3> Velocity2;
	[ReadOnly] public float MaxDepth1;
	[ReadOnly] public float MinDepth1;

	public void Execute(int i)
	{
		float mass1 = Mass1[i];
		float salt1 = Salt1[i];
		float mass2 = Mass2[i];
		float salt2 = Salt2[i];
		float temperature1 = Temperature1[i];
		float temperature2 = Temperature2[i];
		float3 velocity1 = Velocity1[i];
		float3 velocity2 = Velocity2[i];

		float density1 = Atmosphere.GetWaterDensity(Atmosphere.GetWaterSalinity(mass1, salt1), temperature1);
		float density2 = Atmosphere.GetWaterDensity(Atmosphere.GetWaterSalinity(mass2, salt2), temperature2);
		float depth1 = (mass1 + salt1) / density1;
		float depth2 = (mass2 + salt2) / density2;

		if (depth1 > MaxDepth1)
		{
			float move = 1.0f - MaxDepth1 / depth1;
			float moveMass = mass1 * move;
			float moveSalt = salt1 * move;
			float moveCarbon = Carbon1[i] * move;
			Delta1[i] = new DiffusionWater()
			{
				WaterMass = mass1 - moveMass,
				SaltMass = salt1 - moveSalt,
				CarbonMass = Carbon1[i] - moveCarbon,
				Plankton = PlanktonMass1[i],
				PlanktonGlucose = PlanktonGlucose1[i],
				Temperature = temperature1,
				Velocity = velocity1
			};
			float inverseTotalMass = 1.0f / (moveMass + moveSalt + mass2 + salt2);
			Delta2[i] = new DiffusionWater()
			{
				WaterMass = mass2 + moveMass,
				SaltMass = salt2 + moveSalt,
				CarbonMass = Carbon2[i] + moveCarbon,
				Plankton = 0,
				PlanktonGlucose = 0,
				Temperature = (temperature2 * (mass2 + salt2) + temperature1 * (moveMass + moveSalt)) * inverseTotalMass,
				Velocity = (velocity2 * (mass2 + salt2) + velocity1 * (moveMass + moveSalt)) * inverseTotalMass,
			};
		}
		else if (depth2 > 0 && depth1 < MinDepth1)
		{
			float move = math.min(1.0f, math.min(depth2, MinDepth1 - depth1) / depth2);
			float moveMass = mass2 * move;
			float moveSalt = salt2 * move;
			float moveCarbon = Carbon2[i] * move;
			float inverseTotalMass = 1.0f / (moveMass + moveSalt + mass1 + salt1);
			Delta1[i] = new DiffusionWater()
			{
				WaterMass = mass1 + moveMass,
				SaltMass = salt1 + moveSalt,
				CarbonMass = Carbon1[i] + moveCarbon,
				Plankton = PlanktonMass1[i],
				PlanktonGlucose = PlanktonGlucose1[i],
				Temperature = (temperature1 * (mass1 + salt1) + temperature2 * (moveMass + moveSalt)) * inverseTotalMass,
				Velocity = (velocity1 * (mass1 + salt1) + velocity2 * (moveMass + moveSalt)) * inverseTotalMass,
			};
			float anyMassRemaining = move < 1 ? 1 : 0;
			Delta2[i] = new DiffusionWater()
			{
				WaterMass = mass2 - moveMass,
				SaltMass = salt2 - moveSalt,
				CarbonMass = Carbon2[i] - moveCarbon,
				Plankton = 0,
				PlanktonGlucose = 0,
				Temperature = temperature2 * anyMassRemaining,
				Velocity = velocity2 * anyMassRemaining
			};
		}
		else
		{
			Delta1[i] = new DiffusionWater()
			{
				WaterMass = mass1,
				SaltMass = salt1,
				CarbonMass = Carbon1[i],
				Plankton = PlanktonMass1[i],
				PlanktonGlucose = PlanktonGlucose1[i],
				Temperature = temperature1,
				Velocity = velocity1
			};
			Delta2[i] = new DiffusionWater()
			{
				WaterMass = mass2,
				SaltMass = salt2,
				CarbonMass = Carbon2[i],
				Plankton = 0,
				PlanktonGlucose = 0,
				Temperature = temperature2,
				Velocity = velocity2
			};
		}
	}
}

