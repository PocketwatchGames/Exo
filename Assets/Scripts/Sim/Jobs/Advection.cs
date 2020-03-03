//#define DISABLE_VERTICAL_AIR_MOVEMENT
//#define DISABLE_AIR_ADVECTION
#define DISABLE_WATER_ADVECTION
//#define DISABLE_CLOUD_ADVECTION
//#define AdvectionAirJobDebug
//#define AdvectionCloudJobDebug
//#define GetVectorDestCoordsJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

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


#if !AdvectionAirJobDebug
[BurstCompile]
#endif
public struct AdvectionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float3> Wind;
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

		float3 velocity = Wind[i];
		float vapor = Vapor[i];
		float temperature = Temperature[i];

		float newTemperature = 0;
		float newWaterVapor = 0;
		float3 newWind = 0;

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				if (Destination[i].indexA == i)
				{
					float v = Destination[i].valueA;
					newTemperature += temperature * v;
					newWaterVapor += vapor * v;
					newWind += velocity * v;
				}

				float incoming = 0;
				if (Destination[n].indexA == i)
				{
					incoming = Destination[i].valueA;
				}
				else if (Destination[n].indexB == i)
				{
					incoming = Destination[i].valueB;
				}
				else if (Destination[n].indexC == i)
				{
					incoming = Destination[i].valueC;
				}
				newTemperature += temperature * incoming;
				newWaterVapor += vapor * incoming;
				newWind += velocity * incoming;
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

#if !DISABLE_AIR_ADVECTION
		Delta[i] = new DiffusionAir()
		{
			Temperature = newTemperature,
			WaterVapor = newWaterVapor,
			Velocity = newWind,
		};
#endif



	}
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
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PlanetRadius;
	public void Execute(int i)
	{

		float3 velocity = Velocity[i];
		float3 pos = Position[i] * PlanetRadius;
		float3 velVertical = velocity * pos;
		float3 velHorizontal = velocity - velVertical;
		float3 move = velHorizontal * SecondsPerTick;

		float3 moveHorizontal = math.cross(math.cross(pos, move), pos);
		float3 movePos = pos + moveHorizontal;

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
		float3 velocity = Velocity[i];
		float mass = Mass[i];
		float dropletMass = DropletMass[i];

		float newMass = 0;
		float newDropletMass = 0;
		float3 newVelocity = 0;

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				if (Destination[i].indexA == i)
				{
					float v = Destination[i].valueA;
					newMass += mass * v;
					newDropletMass += dropletMass * v;
					newVelocity += velocity * v;
				}

				float incoming = 0;
				if (Destination[n].indexA == i)
				{
					incoming = Destination[i].valueA;
				} else if (Destination[n].indexB == i)
				{
					incoming = Destination[i].valueB;
				} else if (Destination[n].indexC == i)
				{
					incoming = Destination[i].valueC;
				}
				newMass += Mass[n] * incoming;
				newDropletMass += DropletMass[n] * incoming;
				newVelocity += Velocity[n] * incoming;
			}
		}

#if !DISABLE_CLOUD_ADVECTION
		Delta[i] = new DiffusionCloud()
		{
			Mass = newMass,
			DropletMass = newDropletMass,
			Velocity = newVelocity
		};
#endif
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
		float salt = Salt[i];
		float mass = Mass[i];
		float temperature = Temperature[i];

		float newMass = 0;
		float newSaltMass = 0;
		float newTemperature = 0;
		float3 newVelocity = 0;

		for (int j = 0; j < 6; j++)
		{
			int n = Neighbors[i * 6 + j];
			if (n >= 0)
			{
				if (Destination[i].indexA == i)
				{
					float v = Destination[i].valueA;
					newTemperature += Temperature[i] * v;
					newMass += Mass[i] * v;
					newSaltMass += Salt[i] * v;
					newVelocity += Velocity[i] * v;
				}

				float incoming = 0;
				if (Destination[n].indexA == i)
				{
					incoming = Destination[i].valueA;
				}
				else if (Destination[n].indexB == i)
				{
					incoming = Destination[i].valueB;
				}
				else if (Destination[n].indexC == i)
				{
					incoming = Destination[i].valueC;
				}
				newTemperature += Temperature[i] * incoming;
				newMass += Mass[i] * incoming;
				newSaltMass += Salt[i] * incoming;
				newVelocity += Velocity[i] * incoming;
			}
		}


#if !DISABLE_WATER_ADVECTION
		Delta[i] = new DiffusionWater()
		{
			Temperature = newTemperature,
			SaltMass = newSaltMass,
			Velocity = newVelocity
		};
#endif

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
		Temperature[i] += Advection[i].Temperature;
		Vapor[i] += Advection[i].WaterVapor;
		Wind[i] += Advection[i].Velocity;
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
		SaltMass[i] += Advection[i].SaltMass;
		Temperature[i] += Advection[i].Temperature;
		Velocity[i] += Advection[i].Velocity;
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
		CloudMass[i] += Advection[i].Mass;
		DropletMass[i] += Advection[i].DropletMass;
		Velocity[i] += Advection[i].Velocity;
	}
}


