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
	[ReadOnly] public NativeSlice<float> WaterDepth;
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
	[ReadOnly] public NativeSlice<float> WaterDepth;
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

	public NativeSlice<DiffusionWater> Delta;
	[ReadOnly] public NativeSlice<float> Mass;
	[ReadOnly] public NativeSlice<float> Salt;
	[ReadOnly] public NativeSlice<float> Carbon;
	[ReadOnly] public NativeSlice<float> PlanktonMass;
	[ReadOnly] public NativeSlice<float> PlanktonGlucose;
	[ReadOnly] public NativeSlice<float> Temperature;
	[ReadOnly] public NativeSlice<float3> Velocity;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> FlowPercent;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float SecondsPerTick;

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
							temperature += Temperature[nIndex] * (massIncoming + saltIncoming);

							// TODO: this is increasing speed, is that right???  Shouldnt it only rotate?
							var deflectedVelocity = Velocity[nIndex] + math.cross(Positions[nIndex], Velocity[nIndex]) * CoriolisMultiplier[nIndex] * CoriolisTerm * SecondsPerTick;

							// TODO: turn velocity along great circle, instead of just erasing the vertical component as we are doing here
							var deflectedVertical = math.dot(deflectedVelocity, Positions[nIndex]);
							deflectedVelocity -= Positions[nIndex] * deflectedVertical;
							deflectedVelocity -= Positions[i] * math.dot(Positions[i], deflectedVelocity);
							deflectedVelocity += deflectedVertical * Positions[i];


							velocity += deflectedVelocity * (massIncoming + saltIncoming);


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
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float> Carbon;
	[ReadOnly] public NativeArray<float> PlanktonMass;
	[ReadOnly] public NativeArray<float> PlanktonGlucose;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public float MaxDepth1;
	[ReadOnly] public float MinDepth1;
	[ReadOnly] public int Layer1;
	[ReadOnly] public int Layer2;
	[ReadOnly] public int Count;

	public void Execute(int i)
	{
		int index1 = Layer1 * Count + i;
		int index2 = Layer2 * Count + i;
		float mass1 = Mass[index1];
		float salt1 = Salt[index1];
		float mass2 = Mass[index2];
		float salt2 = Salt[index2];
		float temperature1 = Temperature[index1];
		float temperature2 = Temperature[index2];
		float3 velocity1 = Velocity[index1];
		float3 velocity2 = Velocity[index2];

		float density1 = Atmosphere.GetWaterDensity(Atmosphere.GetWaterSalinity(mass1, salt1), temperature1);
		float density2 = Atmosphere.GetWaterDensity(Atmosphere.GetWaterSalinity(mass2, salt2), temperature2);
		float depth1 = (mass1 + salt1) / density1;
		float depth2 = (mass2 + salt2) / density2;

		if (depth1 > MaxDepth1)
		{
			float move = 1.0f - MaxDepth1 / depth1;
			float moveMass = mass1 * move;
			float moveSalt = salt1 * move;
			float moveCarbon = Carbon[index1] * move;
			Delta1[i] = new DiffusionWater()
			{
				WaterMass = mass1 - moveMass,
				SaltMass = salt1 - moveSalt,
				CarbonMass = Carbon[index1] - moveCarbon,
				Plankton = PlanktonMass[index1],
				PlanktonGlucose = PlanktonGlucose[index1],
				Temperature = temperature1,
				Velocity = velocity1
			};
			float inverseTotalMass = 1.0f / (moveMass + moveSalt + mass2 + salt2);
			Delta2[i] = new DiffusionWater()
			{
				WaterMass = mass2 + moveMass,
				SaltMass = salt2 + moveSalt,
				CarbonMass = Carbon[index2] + moveCarbon,
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
			float moveCarbon = Carbon[index2] * move;
			float inverseTotalMass = 1.0f / (moveMass + moveSalt + mass1 + salt1);
			Delta1[i] = new DiffusionWater()
			{
				WaterMass = mass1 + moveMass,
				SaltMass = salt1 + moveSalt,
				CarbonMass = Carbon[index1] + moveCarbon,
				Plankton = PlanktonMass[index1],
				PlanktonGlucose = PlanktonGlucose[index1],
				Temperature = (temperature1 * (mass1 + salt1) + temperature2 * (moveMass + moveSalt)) * inverseTotalMass,
				Velocity = (velocity1 * (mass1 + salt1) + velocity2 * (moveMass + moveSalt)) * inverseTotalMass,
			};
			float anyMassRemaining = move < 1 ? 1 : 0;
			Delta2[i] = new DiffusionWater()
			{
				WaterMass = mass2 - moveMass,
				SaltMass = salt2 - moveSalt,
				CarbonMass = Carbon[index2] - moveCarbon,
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
				CarbonMass = Carbon[index1],
				Plankton = PlanktonMass[index1],
				PlanktonGlucose = PlanktonGlucose[index1],
				Temperature = temperature1,
				Velocity = velocity1
			};
			Delta2[i] = new DiffusionWater()
			{
				WaterMass = mass2,
				SaltMass = salt2,
				CarbonMass = Carbon[index2],
				Plankton = 0,
				PlanktonGlucose = 0,
				Temperature = temperature2,
				Velocity = velocity2
			};
		}
	}
}

