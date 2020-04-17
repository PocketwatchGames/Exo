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
	public void Execute(int i)
	{
		float3 position = Position[i];
		float3 pos = position * PlanetRadius;

		float3 velocity = Velocity[i];

		float3 move = velocity * SecondsPerTick;
		float windMoveHorizontalSq = math.lengthsq(move);
		const float maxWindMove = 200000;
		if (windMoveHorizontalSq > maxWindMove * maxWindMove)
		{
			move = move / math.sqrt(windMoveHorizontalSq) * maxWindMove;
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
public struct GetVectorDestCoordsVerticalJob : IJobParallelFor {
	public NativeArray<BarycentricValueVertical> Destination;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PlanetRadius;
	public void Execute(int i)
	{
		float3 position = Position[i];
		float3 pos = position * PlanetRadius;

		float3 velocity = Velocity[i];

		// extract vertical velocity
		float layerHeight = LayerHeight[i];
		float valueVertical;
		float valueVerticalComplement;
		if (layerHeight > 0)
		{
			var velVertical = math.dot(velocity, position);
			valueVertical = math.clamp(velVertical / LayerHeight[i], -1, 1);
			valueVerticalComplement = 1.0f - math.abs(valueVertical);
		}
		else
		{
			valueVertical = 0;
			valueVerticalComplement = 1;
		}

		// TODO: deal with high wind speeds appropriately and remove this section
		float3 move = velocity * SecondsPerTick;
		float windMoveHorizontalSq = math.lengthsq(move);
		const float maxWindMove = 200000;
		if (windMoveHorizontalSq > maxWindMove * maxWindMove)
		{
			move = move / math.sqrt(windMoveHorizontalSq) * maxWindMove;
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

				float a;
				float b;
				float c;
				if (Utils.GetBarycentricIntersection(movePos, pos, Position[indexB] * PlanetRadius, Position[indexC] * PlanetRadius, out a, out b, out c))
				{
					Destination[i] = new BarycentricValueVertical
					{
						indexA = i,
						indexB = indexB,
						indexC = indexC,
						valueA = a * valueVerticalComplement,
						valueB = b * valueVerticalComplement,
						valueC = c * valueVerticalComplement,
						moveVertical = valueVertical,
					};
					return;
				}
			}
		}

		// TODO: this means the velocity is too high and has skipped over our neighbors!
		Destination[i] = new BarycentricValueVertical
		{
			indexA = i,
			indexB = -1,
			indexC = -1,
			valueA = valueVerticalComplement,
			valueB = 0,
			valueC = 0,
			moveVertical = valueVertical,
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
		value -= math.abs(Destination[i].moveVertical) * mass;
		if (DestinationAbove[i].moveVertical < 0)
		{
			value += -DestinationAbove[i].moveVertical * MassAbove[i];
		}
		if (DestinationBelow[i].moveVertical > 0)
		{
			value += DestinationBelow[i].moveVertical * MassBelow[i];
		}

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
		float pressure = Pressure[i];
		int neighborCount = StaticState.GetMaxNeighbors(i, Neighbors);
		for (int k = 0; k < neighborCount; k++)
		{
			pressure += Pressure[Neighbors[i * StaticState.MaxNeighbors + k]];
		}
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
		Pressure[i] = pressure / neighborCount;

	}
}


[BurstCompile]
public struct GetDivergenceFreeFieldJob : IJobParallelFor {
	public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> PressureAbove;
	[ReadOnly] public NativeArray<float> PressureBelow;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> NeighborTangent;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	[ReadOnly] public float SecondsPerTickInverse;
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
		float height = LayerHeight[i];
		if (!IsTop)
		{
			pressureGradient += (pressure - PressureAbove[i]) * up * height;
		}
		if (!IsBottom)
		{
			pressureGradient += (PressureBelow[i] - pressure) * up * height;
		}
		Velocity[i] -= pressureGradient / AirMass[i] * SecondsPerTickInverse;
	}
}