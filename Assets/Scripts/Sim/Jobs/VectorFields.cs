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
public struct GetVectorDestCoordsCoriolisJob : IJobParallelFor {
	public NativeArray<BarycentricValue> Destination;
	public NativeArray<float3> VelocityDeflected;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float CoriolisTerm;
	public void Execute(int i)
	{
		float3 position = Position[i];
		float3 pos = position * PlanetRadius;

		float3 velocity = Velocity[i];
		// save out deflected velocity since this is what we will advect
		float3 deflectedVelocity = Velocity[i] + math.cross(position, Velocity[i]) * CoriolisMultiplier[i] * CoriolisTerm * SecondsPerTick;
		VelocityDeflected[i] = deflectedVelocity;

		// TODO: deal with high wind speeds appropriately and remove this section
		float3 move = deflectedVelocity * SecondsPerTick;
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
	public NativeArray<float3> VelocityDeflected;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float CoriolisTerm;
	public void Execute(int i)
	{
		float3 position = Position[i];
		float3 pos = position * PlanetRadius;

		float3 velocity = Velocity[i];
		// save out deflected velocity since this is what we will advect
		float3 deflectedVelocity = Velocity[i] + math.cross(position, Velocity[i]) * CoriolisMultiplier[i] * CoriolisTerm * SecondsPerTick;
		VelocityDeflected[i] = deflectedVelocity;

		// extract vertical velocity
		float layerHeight = LayerHeight[i];
		float valueVertical;
		float valueVerticalComplement;
		if (layerHeight > 0)
		{
			var velVertical = math.dot(deflectedVelocity, position);
			valueVertical = math.clamp(velVertical / LayerHeight[i], -1, 1);
			valueVerticalComplement = 1.0f - math.abs(valueVertical);
		}
		else
		{
			valueVertical = 0;
			valueVerticalComplement = 1;
		}

		// TODO: deal with high wind speeds appropriately and remove this section
		float3 move = deflectedVelocity * SecondsPerTick;
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
public struct GetDivergenceFreeFieldJob : IJobParallelFor {
	public NativeArray<BarycentricValueVertical> Destination;
	[ReadOnly] public NativeArray<float> Divergence;
	[ReadOnly] public NativeArray<float> DivergenceAbove;
	[ReadOnly] public NativeArray<float> DivergenceBelow;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	public void Execute(int i)
	{
		float flux = 0;
		float localDivergence = Divergence[i];
		if (!IsTop)
		{
			flux += math.max(0, localDivergence - DivergenceAbove[i]) / 2;
		}
		if (!IsBottom)
		{
			flux -= math.max(0, localDivergence - DivergenceBelow[i]) / 2;
		}
		var oldDest = Destination[i];
		float valueA = oldDest.valueA;
		float valueB = oldDest.valueB;
		float valueC = oldDest.valueC;
		float m = 1.0f - math.abs(oldDest.moveVertical);
		if (m > 0)
		{
			valueA /= m;
			valueB /= m;
			valueC /= m;
		}

		float newMoveVertical = oldDest.moveVertical + math.clamp(flux / AirMass[i], -1, 1);
		float newMoveVertComplement = 1.0f - math.abs(newMoveVertical);
		Destination[i] = new BarycentricValueVertical()
		{
			indexA = oldDest.indexA,
			indexB = oldDest.indexB,
			indexC = oldDest.indexC,
			valueA = valueA * newMoveVertComplement,
			valueB = valueB * newMoveVertComplement,
			valueC = valueC * newMoveVertComplement,
			moveVertical = newMoveVertical,
		};

	}
}