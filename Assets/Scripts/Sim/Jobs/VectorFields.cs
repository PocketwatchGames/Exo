#define DISABLE_VERTICAL_DIVERGENCE

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

public struct BarycentricValue {
	public int indexA;
	public int indexB;
	public int indexC;
	public float valueA;
	public float valueB;
	public float valueC;
}


[BurstCompile]
public struct GetVectorDestCoordsJob : IJobParallelFor {
	public NativeArray<BarycentricValue> Destination;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float MaxWindMove;
	public void Execute(int i)
	{
		float3 position = Position[i];
		float3 pos = position * PlanetRadius;

		float3 velocity = Velocity[i];

		float3 move = velocity * SecondsPerTick;
		float windMoveHorizontalSq = math.lengthsq(move);

		// TODO: this limit is bad, we should just iterate to the next neighbor!
		if (windMoveHorizontalSq > MaxWindMove * MaxWindMove)
		{
			move = move / math.sqrt(windMoveHorizontalSq) * MaxWindMove;
		}

		// TODO: move around arc/circumference instead of along a tangent
		float3 movePos = pos + move;
		for (int j = 0; j < 6; j++)
		{
			int indexB = Neighbors[i * 6 + j];
			if (indexB >= 0)
			{
				int indexC = Neighbors[i * 6 + (j + 1) % 6];
				if (indexC < 0)
				{
					indexC = Neighbors[i * 6];
				}

				float a;
				float b;
				float c;
				if (Utils.GetBarycentricIntersection(movePos, pos, Position[indexB] * PlanetRadius, Position[indexC] * PlanetRadius, out a, out b, out c))
				{
					Destination[i] = new BarycentricValue
					{
						indexA = i,
						indexB = indexB,
						indexC = indexC,
						valueA = a,
						valueB = b,
						valueC = c,
					};
					return;
				}
			}
		}

		// TODO: this means the velocity is too high and has skipped over our neighbors!
		Destination[i] = new BarycentricValue
		{
			indexA = i,
			indexB = -1,
			indexC = -1,
			valueA = 1,
			valueB = 0,
			valueC = 0,
		};
	}
}

[BurstCompile]
public struct ResolveAdvectionConflict : IJobParallelFor {
	public NativeArray<float> ResolvedDestination;
	[ReadOnly] public NativeArray<float> Destination;
	[ReadOnly] public NativeArray<int> ReverseNeighbors;
	public void Execute(int i)
	{
		int n = ReverseNeighbors[i];
		if (n >= 0)
		{
			ResolvedDestination[i] = Destination[i] - Destination[n];
		}
	}
}

[BurstCompile]
public struct GetVectorDestCoordsVerticalJob : IJobParallelFor {
	public NativeArray<float> Destination;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float MaxWindMove;
	public void Execute(int i)
	{
		int cellIndex = i / StaticState.MaxNeighbors;
		float moveValue = 0;
		int indexB = Neighbors[i];

		if (indexB >= 0)
		{
			int indexC = Neighbors[cellIndex * StaticState.MaxNeighbors + (i + 1) % StaticState.MaxNeighbors];
			if (indexC < 0)
			{
				indexC = Neighbors[cellIndex * StaticState.MaxNeighbors];
			}

			float3 m;
			var pos = Position[cellIndex] * PlanetRadius;
			if (Utils.GetBarycentricIntersection(
				pos + Velocity[cellIndex] * SecondsPerTick, 
				pos, 
				Position[indexB] * PlanetRadius, 
				Position[indexC] * PlanetRadius, 
				out m.x, out m.y, out m.z))
			{
				float mass = Mass[cellIndex];
				if (mass == 0)
				{
					moveValue = 0;
				} else
				{
					moveValue = m.y * mass;
				}
			} else
			{
				indexC = Neighbors[cellIndex * StaticState.MaxNeighbors + (i + StaticState.MaxNeighbors - 1) % StaticState.MaxNeighbors];
				if (indexC < 0)
				{
					indexC = Neighbors[cellIndex * StaticState.MaxNeighbors + (i + StaticState.MaxNeighbors - 2) % StaticState.MaxNeighbors];
				}

				if (Utils.GetBarycentricIntersection(
					pos + Velocity[cellIndex] * SecondsPerTick,
					pos,
					Position[indexB] * PlanetRadius,
					Position[indexC] * PlanetRadius,
					out m.x, out m.y, out m.z))
				{
					float mass = Mass[cellIndex];
					if (mass == 0)
					{
						moveValue = 0;
					}
					else
					{
						moveValue = m.y * mass;
					}
				}

			}


		}

		// TODO: this means the velocity is too high and has skipped over our neighbors!
		Destination[i] = moveValue;
	
	}
}

//[BurstCompile]
//public struct GetVectorDestCoordsVerticalJob : IJobParallelFor {
//	public NativeArray<float> Destination;
//	[ReadOnly] public NativeArray<float3> Velocity;
//	[ReadOnly] public NativeArray<int> Neighbors;
//	[ReadOnly] public NativeArray<float3> Position;
//	[ReadOnly] public NativeArray<float> LayerHeight;
//	[ReadOnly] public NativeArray<float> Mass;
//	[ReadOnly] public NativeArray<float> MassAbove;
//	[ReadOnly] public NativeArray<float> MassBelow;
//	[ReadOnly] public float SecondsPerTick;
//	[ReadOnly] public float PlanetRadius;
//	[ReadOnly] public float MaxWindMove;
//	public void Execute(int i)
//	{
//		float mass = Mass[i];
//		if (mass > 0)
//		{
//			float inverseMass = 1.0f / mass;
//			float3 position = Position[i];
//			float3 pos = position * PlanetRadius;

//			float3 velocity = Velocity[i];

//			float valueVertical;
//#if !DISABLE_VERTICAL_DIVERGENCE
//			// extract vertical velocity
//			float layerHeight = LayerHeight[i];
//			if (layerHeight > 0)
//			{
//				var velVertical = math.dot(velocity, position);
//				valueVertical = math.clamp(velVertical / layerHeight  /* * SecondsPerTick */, -1, 1);
//			}
//			else
//#endif
//			{
//				valueVertical = 0;
//			}

//			// TODO: deal with high wind speeds appropriately and remove this section
//			float3 move = velocity * SecondsPerTick;
//			float windMoveHorizontalSq = math.lengthsq(move);
//			if (windMoveHorizontalSq > MaxWindMove * MaxWindMove)
//			{
//				move = move / math.sqrt(windMoveHorizontalSq) * MaxWindMove;
//			}


//			float3 movePos = pos + move;
//			for (int j = 0; j < 6; j++)
//			{
//				int indexB = Neighbors[i * 6 + j];
//				if (indexB >= 0)
//				{
//					int indexC = Neighbors[i * 6 + (j + 1) % 6];
//					if (indexC < 0)
//					{
//						indexC = Neighbors[i * 6];
//					}

//					float3 m;
//					if (Utils.GetBarycentricIntersection(movePos, pos, Position[indexB] * PlanetRadius, Position[indexC] * PlanetRadius, out m.x, out m.y, out m.z))
//					{
//#if !DISABLE_VERTICAL_DIVERGENCE
//						// Make sure we aren't moving too much mass into a cell that can't handle it
//						if (valueVertical > 0)
//						{
//							valueVertical -= math.max(0, valueVertical - MassAbove[i] * inverseMass);
//						}
//						else
//						{
//							valueVertical += math.max(0, -valueVertical - MassBelow[i] * inverseMass);
//						}
//						float valueVerticalComplement = 1.0f - math.abs(valueVertical);
//						m *= valueVerticalComplement;
//#endif
//						if (indexB >= 0)
//						{
//							float excessMass = math.max(0, m.y - Mass[indexB] * inverseMass);
//							m.x += excessMass;
//							m.y -= excessMass;
//						}
//						if (indexC >= 0)
//						{
//							float excessMass = math.max(0, m.z - Mass[indexC] * inverseMass);
//							m.x += excessMass;
//							m.z -= excessMass;
//						}

//						Destination[i] = new BarycentricValueVertical
//						{
//							indexA = i,
//							indexB = indexB,
//							indexC = indexC,
//							valueA = m.x,
//							valueB = m.y,
//							valueC = m.z,
//							moveVertical = valueVertical,
//						};
//						return;
//					}
//				}
//			}
//		}

//		// TODO: this means the velocity is too high and has skipped over our neighbors!
//		Destination[i] = new BarycentricValueVertical
//		{
//			indexA = i,
//			indexB = -1,
//			indexC = -1,
//			valueA = 1,
//			valueB = 0,
//			valueC = 0,
//			moveVertical = 0,
//		};
//	}
//}

[BurstCompile]
public struct GetDivergenceJob : IJobParallelFor {
	public NativeArray<float> Divergence;
	[ReadOnly] public NativeArray<float> Destination;
	[ReadOnly] public NativeArray<float> DestinationAbove;
	[ReadOnly] public NativeArray<float> DestinationBelow;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> MassAbove;
	[ReadOnly] public NativeArray<float> MassBelow;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		float value = 0;
		for (int n = 0; n < StaticState.MaxNeighbors; n++)
		{
			int nIndex = i * StaticState.MaxNeighbors + n;
			value -= Destination[nIndex];
		}

		Divergence[i] = value;
	}
}

[BurstCompile]
public struct GetDivergencePressureJob : IJobFor {
	public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> Divergence;
	[ReadOnly] public NativeArray<int> Neighbors;

	public void Execute(int i)
	{
		float pressure = Divergence[i];
		int neighborCount = StaticState.GetMaxNeighbors(i, Neighbors);
		for (int k = 0; k < neighborCount; k++)
		{
			pressure += Pressure[Neighbors[i * StaticState.MaxNeighbors + k]];
		}
		Pressure[i] = pressure / neighborCount;

	}
}


[BurstCompile]
public struct GetDivergenceFreeFieldJob : IJobParallelFor {
	public NativeArray<float> Destination;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		int cellIndex = i / StaticState.MaxNeighbors;
		float mass = Mass[cellIndex];
		int nIndex = Neighbors[i];
		if (mass > 0 && nIndex >= 0)
		{
			float pressureGradient = Pressure[cellIndex] - Pressure[nIndex];
			Destination[i] += pressureGradient;
		}
	}
}
[BurstCompile]
public struct SumMassLeavingJob : IJobParallelFor {
	public NativeArray<float> MassLeaving;
	[ReadOnly] public NativeArray<float> Destination;
	public void Execute(int i)
	{
		float massLeaving = 0;
		for (int j = 0; j < StaticState.MaxNeighbors; j++)
		{
			massLeaving += math.max(0, Destination[i]);
		}
		MassLeaving[i] = massLeaving;
	}
}
[BurstCompile]
public struct CapMassLeavingJob : IJobParallelFor {
	public NativeArray<float> Destination;
	[ReadOnly] public NativeArray<float> MassLeaving;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		int cellIndex = i / StaticState.MaxNeighbors;
		float dest = Destination[i];
		if (dest < 0)
		{
			cellIndex = Neighbors[i];
		}
		float mass = Mass[cellIndex];
		float massLeaving = math.max(mass, MassLeaving[cellIndex]);
		dest *= mass / massLeaving;
		Destination[i] = dest;
	}
}

//[BurstCompile]
public struct UpdateDivergenceFreeVelocityJob : IJobParallelFor {
	public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float> Destination;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float3> NeighborTangent;
	[ReadOnly] public float TicksPerSecond;
	public void Execute(int i)
	{
		float3 vel = 0;
		if (Mass[i] > 0)
		{
			for (int j = 0; j < StaticState.MaxNeighbors; j++)
			{
				vel += -NeighborTangent[i * StaticState.MaxNeighbors + j] * math.max(0, Destination[i * StaticState.MaxNeighbors + j]);
			}
			vel *= TicksPerSecond / Mass[i];
		}
		Velocity[i] = vel;
	}
}


