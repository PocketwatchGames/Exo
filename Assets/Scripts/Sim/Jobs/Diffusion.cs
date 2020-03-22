//#define DISABLE_VERTICAL_AIR_MOVEMENT
//#define DISABLE_AIR_DIFFUSION
//#define DISABLE_CLOUD_DIFFUSION
//#define DISABLE_WATER_DIFFUSION

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

//[BurstCompile]
public struct DiffusionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastDust;
	[ReadOnly] public NativeArray<float3> LastVelocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpVapor;
	[ReadOnly] public NativeArray<float> UpDust;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float3> UpAirVelocity;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> DownTemperature;
	[ReadOnly] public NativeArray<float> DownVapor;
	[ReadOnly] public NativeArray<float> DownDust;
	[ReadOnly] public NativeArray<float> DownAirMass;
	[ReadOnly] public NativeArray<float3> DownAirVelocity;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	[ReadOnly] public float DiffusionCoefficientHorizontal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	[ReadOnly] public float CellSurfaceArea;
	[ReadOnly] public float CellCircumference;
	public void Execute(int i)
	{
		float layerHeight = LayerHeight[i];
		float mass = AirMass[i];
		float temperature = LastTemperature[i];
		float vapor = LastVapor[i];
		float dust = LastDust[i];
		float3 velocity = LastVelocity[i];

		float inverseMass = 1.0f / mass;
		float vaporPercent = vapor * inverseMass;
		float dustPercent = dust * inverseMass;

		float neighborTemperature = temperature;
		float neighborDust = dust;
		float neighborVapor = vapor;
		float3 neighborVelocity = velocity;

#if !DISABLE_AIR_DIFFUSION
		if (mass > 0)
		{
			for (int j = 0; j < 6; j++)
			{
				int nIndex = i * 6 + j;
				int n = Neighbors[nIndex];
				if (n >= 0)
				{
					float nMass = AirMass[n];
					float saToVolume = math.min(layerHeight, LayerHeight[n]) * CellCircumference / (6 * LayerHeight[n] * CellSurfaceArea);
					float inverseTotalMass = 1.0f / (nMass + mass);
					float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientHorizontal;
					neighborTemperature += (LastTemperature[n] - temperature) * diffusion;
					neighborVelocity += (LastVelocity[n] - velocity) * diffusion;
					neighborVapor += (LastVapor[n] / nMass - vaporPercent) * diffusion * mass;
					neighborDust += (LastDust[n] / nMass - dustPercent) * diffusion * mass;
				}
			}
		}


#if !DISABLE_VERTICAL_AIR_MOVEMENT

		// NOTE: we don't diffuse velocity vertically
		if (!IsTop)
		{
			float nMass = UpAirMass[i];
			float saToVolume = 1.0f / UpLayerHeight[i];
			float inverseTotalMass = 1.0f / (nMass + mass);
			float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
			neighborTemperature += (UpTemperature[i] - temperature) * diffusion;
			neighborVapor += (UpVapor[i] / nMass - vaporPercent) * diffusion * nMass;
			neighborDust += (UpDust[i] / nMass - dustPercent) * diffusion * nMass;
		}
		if (!IsBottom)
		{
			float nMass = DownAirMass[i];
			float saToVolume = 1.0f / DownLayerHeight[i];
			float inverseTotalMass = 1.0f / (nMass + mass);
			float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
			neighborTemperature += (DownTemperature[i] - temperature) * diffusion;
			neighborVapor += (DownVapor[i] / nMass - vaporPercent) * diffusion * nMass;
			neighborDust += (DownDust[i] / nMass - dustPercent) * diffusion * nMass;
		}

#endif
#endif

		Delta[i] = new DiffusionAir()
		{
			Temperature = neighborTemperature,
			Velocity = neighborVelocity,
			WaterVapor = neighborVapor,
			Dust = neighborDust,
		};

	}
}


[BurstCompile]
public struct DiffusionCloudJob : IJobParallelFor {
	public NativeArray<DiffusionCloud> Delta;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastDropletMass;
	[ReadOnly] public NativeArray<float3> LastVelocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public float DiffusionCoefficient;
	public void Execute(int i)
	{
		float3 velocity = LastVelocity[i];
		float mass = LastMass[i];
		float temperature = LastTemperature[i];
		float dropletMass = LastDropletMass[i];

		float newMass = 0;
		float newTemperature = 0;
		float newDropletMass = 0;
		float3 newVelocity = float3.zero;

#if !DISABLE_CLOUD_DIFFUSION
		if (mass > 0)
		{
			for (int j = 0; j < 6; j++)
			{
				int neighborIndex = i * 6 + j;
				int n = Neighbors[neighborIndex];
				if (n >= 0)
				{
					float nMass = LastMass[n];


					float diffusion = nMass / (nMass + mass);
					newMass += nMass - mass;
					newTemperature += (LastTemperature[n] - temperature) * diffusion;
					newDropletMass += (LastDropletMass[n] - dropletMass) * diffusion;
					newVelocity += (LastVelocity[n] - velocity) * diffusion;
				}
			}
		}
#endif

		Delta[i] = new DiffusionCloud()
		{
			Mass = newMass * DiffusionCoefficient + mass,
			Temperature = newTemperature * DiffusionCoefficient + temperature,
			DropletMass = newDropletMass * DiffusionCoefficient + LastDropletMass[i],
			Velocity = newVelocity * DiffusionCoefficient + velocity,
		};

	}
}

//[BurstCompile]
public struct DiffusionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastSalt;
	[ReadOnly] public NativeArray<float3> LastVelocity;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> UpMass;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpSalt;
	[ReadOnly] public NativeArray<float3> UpCurrent;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> DownMass;
	[ReadOnly] public NativeArray<float> DownTemperature;
	[ReadOnly] public NativeArray<float> DownSalt;
	[ReadOnly] public NativeArray<float3> DownCurrent;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	[ReadOnly] public float DiffusionCoefficientHorizontal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	[ReadOnly] public float CellSurfaceArea;
	[ReadOnly] public float CellCircumference;
	public void Execute(int i)
	{
		float layerHeight = LayerHeight[i];
		float mass = LastMass[i];
		float temperature = LastTemperature[i];
		float saltMass = LastSalt[i];
		float3 velocity = LastVelocity[i];
		float saltPercent = saltMass / mass;

		float neighborTemperature = LastTemperature[i];
		float neighborSaltMass = LastSalt[i];
		float3 neighborVelocity = LastVelocity[i];


#if !DISABLE_WATER_DIFFUSION
		if (layerHeight > 0)
		{

			for (int j = 0; j < 6; j++)
			{
				int nIndex = i * 6 + j;
				int n = Neighbors[nIndex];
				if (n >= 0)
				{
					float nMass = LastMass[n];
					if (nMass > 0)
					{
						float saToVolume = math.min(layerHeight, LayerHeight[n]) * CellCircumference / (6 * LayerHeight[n] * CellSurfaceArea);
						float inverseTotalMass = 1.0f / (nMass + mass);
						float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientHorizontal;
						neighborTemperature += (LastTemperature[n] - temperature) * diffusion;
						neighborVelocity += (LastVelocity[n] - velocity) * diffusion;
						neighborSaltMass += (LastSalt[n] / nMass - saltPercent) * diffusion * mass;
					}
				}
			}


			//// NOTE: we don't diffuse velocity vertically
			//// NOTE: we are ignoring adibatic processes in the water -- at 10km, the total lapse is less than 1.5 degrees celsius
			//float upMass = UpMass[i];
			//if (upMass > 0)
			//{
			//	float saToVolume = 1.0f / UpLayerHeight[i];
			//	float inverseTotalMass = 1.0f / (upMass + mass);
			//	float diffusion = saToVolume * upMass * inverseTotalMass * DiffusionCoefficientVertical;
			//	neighborTemperature += (UpTemperature[i] - temperature) * diffusion;
			//	neighborSaltMass += (UpSalt[i] / upMass - saltPercent) * diffusion * upMass;
			//}

			//float downMass = DownMass[i];
			//if (downMass > 0)
			//{
			//	float saToVolume = 1.0f / DownLayerHeight[i];
			//	float inverseTotalMass = 1.0f / (downMass + mass);
			//	float diffusion = saToVolume * downMass * inverseTotalMass * DiffusionCoefficientVertical;
			//	neighborTemperature += (DownTemperature[i] - temperature) * diffusion;
			//	neighborSaltMass += (DownSalt[i] / downMass - saltPercent) * diffusion * downMass;
			//}

		}
#endif

		Delta[i] = new DiffusionWater()
		{
			Temperature = neighborTemperature,
			SaltMass = neighborSaltMass,
			Velocity = neighborVelocity,
		};
	}

}

