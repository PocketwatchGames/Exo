//#define DISABLE_VERTICAL_AIR_MOVEMENT
//#define DISABLE_AIR_ADVECTION
//#define DISABLE_WATER_ADVECTION
//#define DISABLE_CLOUD_ADVECTION
//#define AdvectionAirJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

public struct DiffusionAir {
	public float Temperature;
	public float WaterVapor;
	public float Dust;
	public float3 Velocity;
}
public struct DiffusionCloud {
	public float Mass;
	public float Temperature;
	public float DropletMass;
	public float3 Velocity;
}
public struct DiffusionWater {
	public float Temperature;
	public float SaltMass;
	public float3 Velocity;
}

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
	public int indexVertical;
	public float valueA;
	public float valueB;
	public float valueC;
	public float valueVertical;
}


#if !GetVectorDestCoordsJobDebug
[BurstCompile]
#endif
public struct GetVectorDestCoordsJob : IJobParallelFor {
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

#if !GetVectorDestCoordsJobDebug
[BurstCompile]
#endif
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
		int indexVertical;
		float valueVertical;
		float valueVerticalComplement;
		if (layerHeight > 0)
		{
			var velVertical = math.dot(deflectedVelocity, position);
			indexVertical = (velVertical > 0) ? 1 : -1;
			valueVertical = math.min(1, math.abs(velVertical) / LayerHeight[i]);
			valueVerticalComplement = 1.0f - valueVertical;
		} else
		{
			indexVertical = 0;
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
						indexVertical = indexVertical,
						valueA = a * valueVerticalComplement,
						valueB = b * valueVerticalComplement,
						valueC = c * valueVerticalComplement,
						valueVertical = valueVertical,
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
			indexVertical = indexVertical,
			valueA = valueVerticalComplement,
			valueB = 0,
			valueC = 0,
			valueVertical = valueVertical,
		};
	}
}

#if !EnergyAirJobDebug
[BurstCompile]
#endif
public struct UpdateAirVelocityJob : IJobParallelFor {
	public NativeArray<float3> AirVelocity;
	public NativeArray<float> AirMovementVertical;
	[ReadOnly] public NativeArray<float> WindFriction;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float WindFrictionMultiplier;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float3 wind = AirVelocity[i];
		wind *= (1.0f - WindFriction[i] * WindFrictionMultiplier);
		wind += PressureGradientForce[i] * SecondsPerTick;
		AirVelocity[i] = wind;

		AirMovementVertical[i] = math.dot(Position[i], AirVelocity[i]) * SecondsPerTick;
	}
}


#if !AdvectionAirJobDebug
[BurstCompile]
#endif
public struct AdvectionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> TemperatureAbove;
	[ReadOnly] public NativeArray<float> TemperatureBelow;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float> VaporAbove;
	[ReadOnly] public NativeArray<float> VaporBelow;
	[ReadOnly] public NativeArray<float> Dust;
	[ReadOnly] public NativeArray<float> DustAbove;
	[ReadOnly] public NativeArray<float> DustBelow;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> VelocityAbove;
	[ReadOnly] public NativeArray<float3> VelocityBelow;
	[ReadOnly] public NativeArray<float3> VelocityDeflected;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<BarycentricValueVertical> Destination;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationAbove;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationBelow;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		float newTemperature;
		float newWaterVapor;
		float newDust;
		float3 newVelocity;

#if DISABLE_AIR_ADVECTION
		newTemperature = Temperature[i];
		newWaterVapor = Vapor[i];
		newVelocity = Velocity[i];
		newDust = Dust[i];
#else

		newTemperature = 0;
		newWaterVapor = 0;
		newVelocity = 0;
		newDust = 0;

		// TODO: remove this when we have incompressibility
		float totalValue = 0;

		if (Destination[i].indexA == i)
		{
			float v = Destination[i].valueA;
			totalValue += v;
			newTemperature += Temperature[i] * v;
			newWaterVapor += Vapor[i] * v;
			newDust += Dust[i] * v;
			newVelocity += Velocity[i] * v;
		}

		// TODO: subtract vertical motion first before applying horizontal motion (right now it adds up to more than 1)
		for (int j = 0; j < 6; j++)
		{
			int n = Neighbors[i * 6 + j];
			if (n >= 0)
			{

				float incoming = 0;
				if (Destination[n].indexA == i)
				{
					incoming = Destination[n].valueA;
				}
				else if (Destination[n].indexB == i)
				{
					incoming = Destination[n].valueB;
				}
				else if (Destination[n].indexC == i)
				{
					incoming = Destination[n].valueC;
				}
				totalValue += incoming;
				newTemperature += Temperature[n] * incoming;
				newWaterVapor += Vapor[n] * incoming;
				newDust += Dust[n] * incoming;

				// TODO: this is temp
				// need to deal with centrifugal force/gravity so that as air moves horizontally, it can fly into the air or get pulled to earth
				newVelocity += incoming * (VelocityDeflected[n] - Positions[i] * math.dot(Positions[i], VelocityDeflected[n]));

//				newVelocity += DeflectedVelocity[n] * incoming;
			}
		}



#if !DISABLE_VERTICAL_AIR_MOVEMENT

		// from top coming down
		{
			var vertMove = DestinationAbove[i];
			if (vertMove.indexVertical == -1)
			{
				totalValue += vertMove.valueVertical;
				newTemperature += TemperatureAbove[i] * vertMove.valueVertical;
				newWaterVapor += VaporAbove[i] * vertMove.valueVertical;
				newDust += DustAbove[i] * vertMove.valueVertical;
				newVelocity += VelocityAbove[i] * vertMove.valueVertical;
			}
		}

		// from bottom going up
		{
			var vertMove = DestinationBelow[i];
			if (vertMove.indexVertical == 1)
			{
				totalValue += vertMove.valueVertical;
				newTemperature += TemperatureBelow[i] * vertMove.valueVertical;
				newWaterVapor += VaporBelow[i] * vertMove.valueVertical;
				newDust += DustBelow[i] * vertMove.valueVertical;
				newVelocity += VelocityBelow[i] * vertMove.valueVertical;
			}
		}

#endif
		if (totalValue > 0)
		{
			newTemperature /= totalValue;
			newVelocity /= totalValue;
		}
		else
		{
			// TODO: remove once we have incompressibility
			newTemperature = Temperature[i];
			newVelocity = Velocity[i];
		}

#endif

		Delta[i] = new DiffusionAir()
		{
			Temperature = newTemperature,
			WaterVapor = newWaterVapor,
			Dust = newDust,
			Velocity = newVelocity,
		};



	}
}


#if !AdvectionCloudJobDebug
[BurstCompile]
#endif
public struct AdvectionCloudJob : IJobParallelFor {
	public NativeArray<DiffusionCloud> Delta;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> DropletMass;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> VelocityDeflected;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<BarycentricValue> Destination;
	public void Execute(int i)
	{
		float newMass;
		float newTemperature;
		float newDropletMass;
		float3 newVelocity;

#if DISABLE_CLOUD_ADVECTION
		newMass = Mass[i];
		newTemperature = Temperature[i];
		newDropletMass = DropletMass[i];
		newVelocity = Velocity[i];

#else
		newMass = 0;
		newTemperature = 0;
		newDropletMass = 0;
		newVelocity = 0;
		float totalMass = 0;

		if (Destination[i].indexA == i)
		{
			float v = Destination[i].valueA;
			float m = Mass[i] * v;
			totalMass += m;
			newMass += m;
			newTemperature += Temperature[i] * v * m;
			newDropletMass += DropletMass[i] * v;
			newVelocity += Velocity[i] * v;
		}

		for (int j = 0; j < 6; j++)
		{
			int n = Neighbors[i * 6 + j];
			if (n >= 0)
			{
				float incoming = 0;
				if (Destination[n].indexA == i)
				{
					incoming = Destination[n].valueA;
				}
				else if (Destination[n].indexB == i)
				{
					incoming = Destination[n].valueB;
				}
				else if (Destination[n].indexC == i)
				{
					incoming = Destination[n].valueC;
				}
				float m = Mass[n] * incoming;
				newMass += m;
				totalMass += totalMass;
				newTemperature += Temperature[n] * incoming * m;
				newDropletMass += DropletMass[n] * incoming;
				newVelocity += VelocityDeflected[n] * incoming;
			}
		}

		if (totalMass > 0)
		{
			newTemperature /= totalMass;
		}

#endif
		Delta[i] = new DiffusionCloud()
		{
			Mass = newMass,
			Temperature = newTemperature,
			DropletMass = newDropletMass,
			Velocity = newVelocity
		};
	}
}

#if !AdvectionWaterJobDebug
[BurstCompile]
#endif
public struct AdvectionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<BarycentricValueVertical> Destination;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationAbove;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationBelow;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> VelocityAbove;
	[ReadOnly] public NativeArray<float3> VelocityBelow;
	[ReadOnly] public NativeArray<float3> VelocityDeflected;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> MassAbove;
	[ReadOnly] public NativeArray<float> MassBelow;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> TemperatureAbove;
	[ReadOnly] public NativeArray<float> TemperatureBelow;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float> SaltAbove;
	[ReadOnly] public NativeArray<float> SaltBelow;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		float waterMass = Mass[i];

		if (waterMass == 0)
		{
			return;
		}

		float newMass;
		float newSaltMass;
		float newTemperature;
		float3 newVelocity;

#if DISABLE_WATER_ADVECTION
		newMass = waterMass;
		newSaltMass = Salt[i];
		newTemperature = Temperature[i];
		newVelocity = Velocity[i];
#else
		newMass = 0;
		newSaltMass = 0;
		newTemperature = 0;
		newVelocity = 0;


		float valueRemaining = 0;
		int destIndexA = Destination[i].indexA;
		if (destIndexA == i || destIndexA < 0 || Mass[destIndexA] == 0)
		{
			valueRemaining += Destination[i].valueA;
		}
		int destIndexB = Destination[i].indexB;
		if (destIndexB == i || destIndexB < 0 || Mass[destIndexB] == 0)
		{
			valueRemaining += Destination[i].valueB;
		}
		int destIndexC = Destination[i].indexC;
		if (destIndexC == i || destIndexC < 0 || Mass[destIndexC] == 0)
		{
			valueRemaining += Destination[i].valueC;
		}
		if (Destination[i].indexVertical > 0)
		{
			if (MassAbove[i] == 0)
			{
				valueRemaining += Destination[i].valueVertical;
			}
		} else
		{
			if (MassBelow[i] == 0)
			{
				valueRemaining += Destination[i].valueVertical;
			}
		}
		newMass = Mass[i] * valueRemaining;
		newSaltMass += Salt[i] * valueRemaining;
		newTemperature += Temperature[i] * newMass;
		newVelocity += Velocity[i] * newMass;


		for (int j = 0; j < 6; j++)
		{
			int n = Neighbors[i * 6 + j];
			if (n >= 0)
			{
				float nMass = Mass[n];
				if (nMass > 0)
				{
					float incoming = 0;
					if (Destination[n].indexA == i)
					{
						incoming = Destination[n].valueA;
					}
					else if (Destination[n].indexB == i)
					{
						incoming = Destination[n].valueB;
					}
					else if (Destination[n].indexC == i)
					{
						incoming = Destination[n].valueC;
					}
					float massIncoming = nMass * incoming;
					newMass += massIncoming;
					newSaltMass += Salt[n] * incoming;
					newTemperature += Temperature[n] * massIncoming;

					// TODO: this is temp
					// need to deal with centrifugal force/gravity so that as air moves horizontally, it can fly into the air or get pulled to earth
					newVelocity += massIncoming * (VelocityDeflected[n] - Positions[i] * math.dot(Positions[i], VelocityDeflected[n]));

//					newVelocity += VelocityDeflected[n] * massIncoming;
				}
			}
		}


		// from top coming down
		{
			var vertMove = DestinationAbove[i];
			if (vertMove.indexVertical == -1)
			{
				float massIncoming = MassAbove[i] * vertMove.valueVertical;
				newMass += massIncoming;
				newTemperature += TemperatureAbove[i] * massIncoming;
				newVelocity += VelocityAbove[i] * massIncoming;
				newSaltMass += SaltAbove[i] * vertMove.valueVertical;
			}
		}

		// from bottom going up
		{
			var vertMove = DestinationBelow[i];
			if (vertMove.indexVertical == 1)
			{
				float massIncoming = MassBelow[i] * vertMove.valueVertical;
				newMass += massIncoming;
				newTemperature += TemperatureBelow[i] * massIncoming;
				newVelocity += VelocityBelow[i] * massIncoming;
				newSaltMass += SaltBelow[i] * vertMove.valueVertical;
			}
		}


		// TODO: remove this when we have incompressibility??
		// OR NOT? water is compressible!
		if (newMass > 0)
		{
			newTemperature /= newMass;
			newVelocity /= newMass;
		}
		else
		{
			// TODO: remove once we have incompressibility
			newTemperature = Temperature[i];
			newVelocity = Velocity[i];
		}

#endif


		Delta[i] = new DiffusionWater()
		{
			Temperature = newTemperature,
			SaltMass = newSaltMass,
			Velocity = newVelocity
		};

	}
}


#if !ApplyAdvectionAirJobDebug
[BurstCompile]
#endif
public struct ApplyAdvectionAirJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Vapor;
	public NativeArray<float> Dust;
	public NativeArray<float3> AirVelocity;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	public void Execute(int i)
	{
		Temperature[i] = Advection[i].Temperature;
		Vapor[i] = Advection[i].WaterVapor;
		Dust[i] = Advection[i].Dust;
		AirVelocity[i] = Advection[i].Velocity;
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
		SaltMass[i] = Advection[i].SaltMass;
		Temperature[i] = Advection[i].Temperature;
		Velocity[i] = Advection[i].Velocity;
	}
}


#if !ApplyAdvectionCloudJobDebug
[BurstCompile]
#endif
public struct ApplyAdvectionCloudJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> DropletMass;
	public NativeArray<float> Temperature;
	public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<DiffusionCloud> Advection;
	public void Execute(int i)
	{
		CloudMass[i] = Advection[i].Mass;
		Temperature[i] = Advection[i].Temperature;
		DropletMass[i] = Advection[i].DropletMass;
		Velocity[i] = Advection[i].Velocity;
	}
}


