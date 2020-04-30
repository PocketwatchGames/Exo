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
public struct BarycentricValueVertical {
	public int indexA;
	public int indexB;
	public int indexC;
	public float valueA;
	public float valueB;
	public float valueC;
	public float moveVertical;
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
	public NativeArray<float3> NewVelocity;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> VelocityAbove;
	[ReadOnly] public NativeArray<float3> VelocityBelow;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerHeightAbove;
	[ReadOnly] public NativeArray<float> LayerHeightBelow;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> MassAbove;
	[ReadOnly] public NativeArray<float> MassBelow;
	[ReadOnly] public bool IsBottom;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float mass = Mass[i];
		float3 velocity = Velocity[i];
#if !DISABLE_VERTICAL_DIVERGENCE
		if (mass > 0)
		{
			float3 position = Position[i];


			// extract vertical velocity
			float layerHeight = LayerHeight[i];
			float valueVertical;
			if (layerHeight > 0)
			{
				var velVertical = math.dot(velocity, position);
				valueVertical = math.clamp(velVertical / layerHeight  /* * SecondsPerTick */, -1, 1);
				float conflict = 0;
				if (valueVertical > 0)
				{
					if (IsTop || LayerHeightAbove[i] == 0)
					{
						conflict = 1;
					}
					else
					{
						var velAboveVertical = math.dot(VelocityAbove[i], position);
						if (velAboveVertical < 0)
						{
							float valueAboveVertical = math.min(1, -velAboveVertical / LayerHeightAbove[i]);
							conflict = math.min(valueAboveVertical * MassAbove[i] / (mass * valueVertical), 1);
						}
					}
					velocity -= position * conflict * velVertical;
					valueVertical -= conflict * valueVertical;
				}
				else if (valueVertical < 0)
				{
					if (IsBottom || LayerHeightBelow[i] == 0)
					{
						conflict = 1;
					}
					else
					{
						var velBelowVertical = math.dot(VelocityBelow[i], position);
						if (velBelowVertical > 0)
						{
							float valueBelowVertical = math.min(1, velBelowVertical / LayerHeightBelow[i]);
							conflict = math.min(valueBelowVertical * MassBelow[i] / (mass * -valueVertical), 1);
						}
					}
					velocity -= position * conflict * velVertical;
					valueVertical -= conflict * valueVertical;
				}
			}
		}
#endif
		NewVelocity[i] = velocity;
	}
}

[BurstCompile]
public struct GetVectorDestCoordsVerticalJob : IJobParallelFor {
	public NativeArray<BarycentricValueVertical> Destination;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> MassAbove;
	[ReadOnly] public NativeArray<float> MassBelow;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float MaxWindMove;
	public void Execute(int i)
	{
		float mass = Mass[i];
		if (mass > 0)
		{
			float inverseMass = 1.0f / mass;
			float3 position = Position[i];
			float3 pos = position * PlanetRadius;

			float3 velocity = Velocity[i];

			float valueVertical;
#if !DISABLE_VERTICAL_DIVERGENCE
			// extract vertical velocity
			float layerHeight = LayerHeight[i];
			if (layerHeight > 0)
			{
				var velVertical = math.dot(velocity, position);
				valueVertical = math.clamp(velVertical / layerHeight  /* * SecondsPerTick */, -1, 1);
			}
			else
#endif
			{
				valueVertical = 0;
			}

			// TODO: deal with high wind speeds appropriately and remove this section
			float3 move = velocity * SecondsPerTick;
			float windMoveHorizontalSq = math.lengthsq(move);
			if (windMoveHorizontalSq > MaxWindMove * MaxWindMove)
			{
				move = move / math.sqrt(windMoveHorizontalSq) * MaxWindMove;
			}


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

					float3 m;
					if (Utils.GetBarycentricIntersection(movePos, pos, Position[indexB] * PlanetRadius, Position[indexC] * PlanetRadius, out m.x, out m.y, out m.z))
					{
#if !DISABLE_VERTICAL_DIVERGENCE
						// Make sure we aren't moving too much mass into a cell that can't handle it
						if (valueVertical > 0)
						{
							valueVertical -= math.max(0, valueVertical - MassAbove[i] * inverseMass);
						}
						else
						{
							valueVertical += math.max(0, -valueVertical - MassBelow[i] * inverseMass);
						}
						float valueVerticalComplement = 1.0f - math.abs(valueVertical);
						m *= valueVerticalComplement;
#endif
						if (indexB >= 0)
						{
							float excessMass = math.max(0, m.y - Mass[indexB] * inverseMass);
							m.x += excessMass;
							m.y -= excessMass;
						}
						if (indexC >= 0)
						{
							float excessMass = math.max(0, m.z - Mass[indexC] * inverseMass);
							m.x += excessMass;
							m.z -= excessMass;
						}

						Destination[i] = new BarycentricValueVertical
						{
							indexA = i,
							indexB = indexB,
							indexC = indexC,
							valueA = m.x,
							valueB = m.y,
							valueC = m.z,
							moveVertical = valueVertical,
						};
						return;
					}
				}
			}
		}

		// TODO: this means the velocity is too high and has skipped over our neighbors!
		Destination[i] = new BarycentricValueVertical
		{
			indexA = i,
			indexB = -1,
			indexC = -1,
			valueA = 1,
			valueB = 0,
			valueC = 0,
			moveVertical = 0,
		};
	}
}

[BurstCompile]
public struct GetDivergenceJob : IJobParallelFor {
	public NativeArray<float> Divergence;
	[ReadOnly] public NativeArray<BarycentricValueVertical> Destination;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationAbove;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationBelow;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> MassAbove;
	[ReadOnly] public NativeArray<float> MassBelow;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	public void Execute(int i)
	{
		float value = 0;
		float mass = Mass[i];
		int destIndexA = Destination[i].indexA;
		if (destIndexA == i || destIndexA < 0)
		{
			value += Destination[i].valueA * mass;
		}
		int destIndexB = Destination[i].indexB;
		if (destIndexB == i || destIndexB < 0)
		{
			value += Destination[i].valueB * mass;
		}
		int destIndexC = Destination[i].indexC;
		if (destIndexC == i || destIndexC < 0)
		{
			value += Destination[i].valueC * mass;
		}
#if !DISABLE_VERTICAL_DIVERGENCE
		if (DestinationAbove[i].moveVertical < 0)
		{
			value += -DestinationAbove[i].moveVertical * MassAbove[i];
		}
		if (DestinationBelow[i].moveVertical > 0)
		{
			value += DestinationBelow[i].moveVertical * MassBelow[i];
		}
#endif

		for (int j = 0; j < 6; j++)
		{
			int n = Neighbors[i * 6 + j];
			if (n >= 0)
			{
				if (Destination[n].indexA == i)
				{
					value += Destination[n].valueA * Mass[n];
				}
				else if (Destination[n].indexB == i)
				{
					value += Destination[n].valueB * Mass[n];
				}
				else if (Destination[n].indexC == i)
				{
					value += Destination[n].valueC * Mass[n];
				}
			}
		}
		Divergence[i] = value - mass;
	}
}

[BurstCompile]
public struct GetDivergencePressureJob : IJobFor {
	public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> PressureAbove;
	[ReadOnly] public NativeArray<float> PressureBelow;
	[ReadOnly] public NativeArray<float> Divergence;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;

	public void Execute(int i)
	{
		float pressure = Divergence[i];
		int neighborCount = StaticState.GetMaxNeighbors(i, Neighbors);
		for (int k = 0; k < neighborCount; k++)
		{
			pressure += Pressure[Neighbors[i * StaticState.MaxNeighbors + k]];
		}
#if !DISABLE_VERTICAL_DIVERGENCE
		if (!IsTop)
		{
			pressure += PressureAbove[i];
			neighborCount++;
		}
		if (!IsBottom)
		{
			pressure += PressureBelow[i];
			neighborCount++;
		}
#endif
		Pressure[i] = pressure / neighborCount;

	}
}


[BurstCompile]
public struct GetDivergenceFreeFieldJob : IJobParallelFor {
	public NativeArray<float3> VelocityOut;
	[ReadOnly] public NativeArray<float3> VelocityIn;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> PressureAbove;
	[ReadOnly] public NativeArray<float> PressureBelow;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> NeighborTangent;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	[ReadOnly] public float TicksPerSecond;
	public void Execute(int i)
	{
		int neighborCount = StaticState.GetMaxNeighbors(i, Neighbors);
		float pressure = Pressure[i];
		float3 pressureGradient = 0;
		float3 up = Positions[i];
		for (int j=0;j< neighborCount; j++)
		{
			int nIndex = i * StaticState.MaxNeighbors + j;
			float3 diff = NeighborTangent[nIndex];
			pressureGradient += (pressure - Pressure[Neighbors[nIndex]]) * diff;
		}
#if !DISABLE_VERTICAL_DIVERGENCE
		float height = LayerHeight[i];
		if (!IsTop)
		{
			pressureGradient += (pressure - PressureAbove[i]) * up * height;
		}
		if (!IsBottom)
		{
			pressureGradient += (PressureBelow[i] - pressure) * up * height;
		}
#endif
		VelocityOut[i] = VelocityIn[i] - pressureGradient / Mass[i] * TicksPerSecond;
	}
}

