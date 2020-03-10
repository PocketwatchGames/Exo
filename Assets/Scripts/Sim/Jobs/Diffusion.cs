//#define DISABLE_VERTICAL_AIR_MOVEMENT
#define DISABLE_AIR_DIFFUSION
//#define DISABLE_CLOUD_DIFFUSION
#define DISABLE_WATER_DIFFUSION
//#define DiffusionCloudJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

#if !DiffusionAirJobDebug
[BurstCompile]
#endif
public struct DiffusionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float3> LastVelocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpVapor;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float3> UpAirVelocity;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> DownTemperature;
	[ReadOnly] public NativeArray<float> DownVapor;
	[ReadOnly] public NativeArray<float> DownAirMass;
	[ReadOnly] public NativeArray<float3> DownAirVelocity;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	[ReadOnly] public float DiffusionCoefficientHoriztonal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	public void Execute(int i)
	{
		float airMass = AirMass[i];
		float absoluteHumidity = LastVapor[i] / airMass;

		float newTemperature = LastTemperature[i];
		float newHumidity = absoluteHumidity;
		float3 newVelocity = LastVelocity[i];

#if !DISABLE_AIR_DIFFUSION
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float3 nVelocity = LastVelocity[n];

				float neighborMass = AirMass[n];
				float diffusionAmount = neighborMass / (neighborMass + airMass) * DiffusionCoefficientHoriztonal;
				float neighborHumidity = LastVapor[n] / neighborMass;
				newHumidity += (neighborHumidity - absoluteHumidity) * diffusionAmount;
				newVelocity += (nVelocity - LastVelocity[i]) * diffusionAmount;
				newTemperature += (LastTemperature[n] - LastTemperature[i]) * diffusionAmount;

			}
		}



#if !DISABLE_VERTICAL_AIR_MOVEMENT

		// NOTE: we don't diffuse velocity vertically
		if (!IsTop)
		{
			float diffusionAmount = UpAirMass[i] / (UpAirMass[i] + airMass) * DiffusionCoefficientVertical;

			float absoluteHumidityUp = UpVapor[i] /  UpAirMass[i];
			newHumidity += (absoluteHumidityUp - absoluteHumidity) * diffusionAmount;

			float heightDiff = (UpLayerElevation[i] + UpLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
			float potentialTemperatureUp = UpTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
			newTemperature += (potentialTemperatureUp - LastTemperature[i]) * diffusionAmount;
		}
		if (!IsBottom)
		{
			float diffusionAmount = DownAirMass[i] / (DownAirMass[i] + airMass) * DiffusionCoefficientVertical;

			float absoluteHumidityDown = DownVapor[i] / DownAirMass[i];
			newHumidity += (absoluteHumidityDown - absoluteHumidity) * diffusionAmount;

			float heightDiff = (DownLayerElevation[i] + DownLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
			float potentialTemperatureDown = DownTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
			newTemperature += (potentialTemperatureDown - LastTemperature[i]) * diffusionAmount;
		}

		//		float moveToNeutralBuoyancy = (UpTemperature[i] - Temperature[i]) / WorldData.TemperatureLapseRate - heightDiff;
		//		float vertMovement = math.min(MaxVerticalMovement, math.clamp(moveToNeutralBuoyancy + DiffusionCoefficient, 0, 1));

#endif
#endif

		Delta[i] = new DiffusionAir()
		{
			Temperature = newTemperature,
			WaterVapor = newHumidity * airMass,
			Velocity = newVelocity,
		};

	}
}


#if !DiffusionCloudJobDebug
[BurstCompile]
#endif
public struct DiffusionCloudJob : IJobParallelFor {
	public NativeArray<DiffusionCloud> Delta;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastDropletMass;
	[ReadOnly] public NativeArray<float3> LastVelocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public float DiffusionCoefficient;
	public void Execute(int i)
	{
		float3 velocity = LastVelocity[i];
		float mass = LastMass[i];

		float newMass = 0;
		float newDropletMass = 0;
		float3 newVelocity = float3.zero;

#if !DISABLE_CLOUD_DIFFUSION
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float nMass = LastMass[n];
				if (nMass > 0 || mass > 0)
				{
					float diffusionAmount = nMass / (nMass + mass);

					newMass += (nMass - mass) * diffusionAmount;
					newDropletMass += (LastDropletMass[n] - LastDropletMass[i]) * diffusionAmount;
					newVelocity += (LastVelocity[n] - LastVelocity[i]) * diffusionAmount;
				}
			}
		}
#endif

		Delta[i] = new DiffusionCloud()
		{
			Mass = newMass * DiffusionCoefficient + mass,
			DropletMass = newDropletMass * DiffusionCoefficient + LastDropletMass[i],
			Velocity = newVelocity * DiffusionCoefficient + velocity,
		};

	}
}

#if !DiffusionWaterJobDebug
[BurstCompile]
#endif
public struct DiffusionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastSalt;
	[ReadOnly] public NativeArray<float3> LastCurrent;
	[ReadOnly] public NativeArray<float> UpMass;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpSalt;
	[ReadOnly] public NativeArray<float3> UpCurrent;
	[ReadOnly] public NativeArray<float> DownMass;
	[ReadOnly] public NativeArray<float> DownTemperature;
	[ReadOnly] public NativeArray<float> DownSalt;
	[ReadOnly] public NativeArray<float3> DownCurrent;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public float DiffusionCoefficientHoriztonal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	public void Execute(int i)
	{
		float newSaltMass = LastSalt[i];
		float newTemperature = LastTemperature[i];
		float3 newVelocity = LastCurrent[i];
		float waterMass = LastMass[i];

#if !DISABLE_WATER_DIFFUSION
		if (waterMass > 0)
		{
			float salinity = LastSalt[i] / waterMass;
			for (int j = 0; j < 6; j++)
			{
				int n = Neighbors[i * 6 + j];
				if (n >= 0)
				{
					float nMass = LastMass[n];
					if (nMass > 0)
					{
						float diffusionAmount = nMass / (nMass + waterMass) * DiffusionCoefficientHoriztonal;
						float neighborSalinity = LastSalt[n] / nMass;

						newSaltMass += (neighborSalinity - salinity) * waterMass * diffusionAmount;
						newTemperature += (LastTemperature[n] - LastTemperature[i]) * diffusionAmount;
						newVelocity += (LastCurrent[n] - LastCurrent[i]) * diffusionAmount;
					}
				}
			}

			// NOTE: we don't diffuse velocity vertically
			// NOTE: we are ignoring adibatic processes in the water -- at 10km, the total lapse is less than 1.5 degrees celsius
			float upMass = UpMass[i];
			if (upMass > 0)
			{
				// TODO: shouldn't shorter water columns diffuse faster?
				// If so we should divide (again) by total mass (which is proportional to height) of the shared cells...
				// ...clamp it to 1/6th (or so) of the total, and clamp the height to the diffusion range (10 meters?)
				float diffusionAmount = upMass / (upMass + waterMass) * DiffusionCoefficientVertical;
				float neighborSalinity = UpSalt[i] / upMass;
				newSaltMass += (neighborSalinity - salinity) * waterMass * diffusionAmount;
				newTemperature += (UpTemperature[i] - LastTemperature[i]) * diffusionAmount;
			}
			float downMass = DownMass[i];
			if (downMass > 0)
			{
				float diffusionAmount = downMass / (downMass + waterMass) * DiffusionCoefficientVertical;
				float neighborSalinity = DownSalt[i] / downMass;

				newSaltMass += (neighborSalinity - salinity) * waterMass * diffusionAmount;
				newTemperature += (DownTemperature[i] - LastTemperature[i]) * diffusionAmount;
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

