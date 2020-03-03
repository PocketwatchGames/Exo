//#define DISABLE_VERTICAL_AIR_MOVEMENT
//#define DISABLE_AIR_ADVECTION
#define DISABLE_WATER_ADVECTION
//#define DISABLE_CLOUD_ADVECTION
#define AdvectionAirJobDebug
#define AdvectionCloudJobDebug
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
				newVelocity += Velocity[n] * incoming;
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
		//newTemperature = Temperature[i];
		//newVelocity = Velocity[i];
		//newWaterVapor = Vapor[i];

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
	[ReadOnly] public NativeArray<float> DropletMass;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<BarycentricValue> Destination;
	public void Execute(int i)
	{
		float newMass;
		float newDropletMass;
		float3 newVelocity;

#if DISABLE_CLOUD_ADVECTION
		newMass = Mass[i];
		newDropletMass = DropletMass[i];
		newVelocity = Velocity[i];

#else
		newMass = 0;
		newDropletMass = 0;
		newVelocity = 0;

		if (Destination[i].indexA == i)
		{
			float v = Destination[i].valueA;
			newMass += Mass[i] * v;
			newDropletMass += DropletMass[i] * v;
			newVelocity += Velocity[i] * v;

		}

		if (newMass < 0)
		{
			Debug.Break();
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
				newMass += Mass[n] * incoming;
				newDropletMass += DropletMass[n] * incoming;
				newVelocity += Velocity[n] * incoming;
			}
		}

#endif
		Delta[i] = new DiffusionCloud()
		{
			Mass = newMass,
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
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		float waterMass = Mass[i];
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
		float newMass = 0;
		float newSaltMass = 0;
		float newTemperature = 0;
		float3 newVelocity = 0;

		if (Destination[i].indexA == i)
		{
			float v = Destination[i].valueA;
			newTemperature += Temperature[i] * v;
			newMass += Mass[i] * v;
			newSaltMass += Salt[i] * v;
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
				newTemperature += Temperature[n] * incoming;
				newMass += Mass[n] * incoming;
				newSaltMass += Salt[n] * incoming;
				newVelocity += Velocity[n] * incoming;
			}
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
	public NativeArray<float3> Wind;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	public void Execute(int i)
	{
		Temperature[i] = Advection[i].Temperature;
		Vapor[i] = Advection[i].WaterVapor;
		Wind[i] = Advection[i].Velocity;
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
	public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<DiffusionCloud> Advection;
	public void Execute(int i)
	{
		CloudMass[i] = Advection[i].Mass;
		DropletMass[i] = Advection[i].DropletMass;
		Velocity[i] = Advection[i].Velocity;
	}
}


