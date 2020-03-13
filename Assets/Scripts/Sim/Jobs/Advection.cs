﻿//#define DISABLE_VERTICAL_AIR_MOVEMENT
//#define DISABLE_AIR_ADVECTION
//#define DISABLE_WATER_ADVECTION
//#define DISABLE_CLOUD_ADVECTION
//#define AdvectionAirJobDebug
//#define AdvectionCloudJobDebug
//#define GetVectorDestCoordsJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

public struct DiffusionAir {
	public float Temperature;
	public float WaterVapor;
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


#if !GetVectorDestCoordsJobDebug
[BurstCompile]
#endif
public struct GetVectorDestCoordsJob : IJobParallelFor {
	public NativeArray<BarycentricValue> Destination;
	public NativeArray<float3> DeflectedVelocity;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public float CoriolisTerm;
	public void Execute(int i)
	{
		float3 position = Position[i];
		float3 velocity = Velocity[i];
		float3 pos = position * PlanetRadius;

		float3 deflectedVelocity = velocity + math.cross(position, velocity) * CoriolisMultiplier[i] * CoriolisTerm * SecondsPerTick;
		float3 velVertical = math.dot(deflectedVelocity, position) * position;
		float3 velHorizontal = deflectedVelocity - velVertical;

		float3 moveHorizontal = velHorizontal * SecondsPerTick;

		// TODO: deal with high wind speeds appropriately and remove this section
		float windMoveHorizontal = math.length(moveHorizontal);

		const float maxWindMove = 200000;
		if (windMoveHorizontal > maxWindMove)
		{
			moveHorizontal = moveHorizontal / windMoveHorizontal * maxWindMove;
		}


		float3 movePos = pos + moveHorizontal;
		DeflectedVelocity[i] = deflectedVelocity;

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
						valueC = c
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
			valueC = 0
		};
	}
}

#if !EnergyAirJobDebug
[BurstCompile]
#endif
public struct UpdateAirVelocityJob : IJobParallelFor {
	public NativeArray<float3> AirVelocity;
	[ReadOnly] public NativeArray<float> WindFriction;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public float WindFrictionMultiplier;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float3 wind = AirVelocity[i];
		wind *= (1.0f - WindFriction[i] * WindFrictionMultiplier);
		wind += PressureGradientForce[i] * SecondsPerTick;
		AirVelocity[i] = wind;
	}
}


#if !AdvectionAirJobDebug
[BurstCompile]
#endif
public struct AdvectionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> DeflectedVelocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<BarycentricValue> Destination;
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
		float newTemperature;
		float newWaterVapor;
		float3 newVelocity;

#if DISABLE_AIR_ADVECTION
		newTemperature = Temperature[i];
		newWaterVapor = Vapor[i];
		newVelocity = Velocity[i];
#else

		newTemperature = 0;
		newWaterVapor = 0;
		newVelocity = 0;

		// TODO: remove this when we have incompressibility
		float totalValue = 0;

		if (Destination[i].indexA == i)
		{
			float v = Destination[i].valueA;
			totalValue += v;
			newTemperature += Temperature[i] * v;
			newWaterVapor += Vapor[i] * v;
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
				totalValue += incoming;
				newTemperature += Temperature[n] * incoming;
				newWaterVapor += Vapor[n] * incoming;
				newVelocity += DeflectedVelocity[n] * incoming;
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
	[ReadOnly] public NativeArray<float3> DeflectedVelocity;
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
				newVelocity += DeflectedVelocity[n] * incoming;
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
	[ReadOnly] public NativeArray<BarycentricValue> Destination;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> DeflectedVelocity;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
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

		// TODO: remove this when we have incompressibility
		float totalMass = 0;


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
		totalMass = Mass[i] * valueRemaining;
		newMass = totalMass;
		newSaltMass += Salt[i] * valueRemaining;
		newTemperature += Temperature[i] * totalMass;
		newVelocity += Velocity[i] * totalMass;


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
				float massIncoming = Mass[n] * incoming;
				newMass += massIncoming;
				newSaltMass += Salt[n] * incoming;
				newTemperature += Temperature[n] * massIncoming;
				newVelocity += DeflectedVelocity[n] * massIncoming;
				totalMass += massIncoming;
			}
		}


		if (totalMass > 0)
		{
			newTemperature /= totalMass;
			newVelocity /= totalMass;
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
	public NativeArray<float3> AirVelocity;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	public void Execute(int i)
	{
		Temperature[i] = Advection[i].Temperature;
		Vapor[i] = Advection[i].WaterVapor;
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


