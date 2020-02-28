//#define DISABLE_VERTICAL_AIR_MOVEMENT

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !DiffusionAirJobDebug
[BurstCompile]
#endif
public struct DiffusionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float3> LastWind;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpHumidity;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float3> UpWind;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> DownTemperature;
	[ReadOnly] public NativeArray<float> DownHumidity;
	[ReadOnly] public NativeArray<float> DownAirMass;
	[ReadOnly] public NativeArray<float3> DownWind;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	[ReadOnly] public float DiffusionCoefficientHoriztonal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	public void Execute(int i)
	{
		float airMass = AirMass[i];
		float absoluteHumidity = LastVapor[i] / (LastVapor[i] + airMass);

		float newTemperature = 0;
		float newWaterVapor = 0;
		float3 newWind = 0;

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float3 nWind = LastWind[n];

				float neighborMass = AirMass[n] + LastVapor[n];
				float diffusionAmount = AirMass[n] / (AirMass[n] + airMass) * DiffusionCoefficientHoriztonal;
				float neighborHumidity = LastVapor[n] / neighborMass;
				newWaterVapor += (neighborHumidity - absoluteHumidity) * diffusionAmount * airMass;
				newWind += (nWind - LastWind[i]) * diffusionAmount;
				newTemperature += (LastTemperature[n] - LastTemperature[i]) * diffusionAmount;

			}
		}



#if !DISABLE_VERTICAL_AIR_MOVEMENT

		// NOTE: we don't diffuse velocity vertically
		if (!IsTop)
		{
			float diffusionAmount = UpAirMass[i] / (UpAirMass[i] + airMass) * DiffusionCoefficientVertical;

			float absoluteHumidityUp = UpHumidity[i] / (UpHumidity[i] + UpAirMass[i]);
			newWaterVapor += (absoluteHumidityUp - absoluteHumidity) * airMass * diffusionAmount;

			float heightDiff = (UpLayerElevation[i] + UpLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
			float potentialTemperatureUp = UpTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
			newTemperature += (potentialTemperatureUp - LastTemperature[i]) * diffusionAmount;
		}
		if (!IsBottom)
		{
			float diffusionAmount = DownAirMass[i] / (DownAirMass[i] + airMass) * DiffusionCoefficientVertical;

			float absoluteHumidityDown = DownHumidity[i] / (DownHumidity[i] + DownAirMass[i]);
			newWaterVapor += (absoluteHumidityDown - absoluteHumidity) * diffusionAmount;

			float heightDiff = (DownLayerElevation[i] + DownLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
			float potentialTemperatureDown = DownTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
			newTemperature += (potentialTemperatureDown - LastTemperature[i]) * diffusionAmount;
		}

		//		float moveToNeutralBuoyancy = (UpTemperature[i] - Temperature[i]) / WorldData.TemperatureLapseRate - heightDiff;
		//		float vertMovement = math.min(MaxVerticalMovement, math.clamp(moveToNeutralBuoyancy + DiffusionCoefficient, 0, 1));

#endif

		Delta[i] = new DiffusionAir()
		{
			Temperature = newTemperature,
			WaterVapor = newWaterVapor,
			Velocity = newWind,
		};

	}
}


#if !DiffusionCloudJob
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
#if !DISABLE_CLOUD_ADVECTION
		float newMass = 0;
		float newDropletMass = 0;
		float3 newVelocity = float3.zero;
		float3 velocity = LastVelocity[i];
		float mass = LastMass[i];

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

		Delta[i] = new DiffusionCloud()
		{
			Mass = newMass * DiffusionCoefficient,
			DropletMass = newDropletMass * DiffusionCoefficient,
			Velocity = newVelocity * DiffusionCoefficient,
		};
#endif
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
		float newSaltMass = 0;
		float newTemperature = 0;
		float3 newVelocity = float3.zero;
		float waterMass = LastMass[i];
		float mass = LastMass[i];

		if (waterMass > 0)
		{
			float salinity = LastSalt[i] / (LastSalt[i] + mass);
			for (int j = 0; j < 6; j++)
			{
				int n = Neighbors[i * 6 + j];
				if (n >= 0)
				{
					float nMass = LastMass[n];
					if (nMass > 0)
					{
						float diffusionAmount = nMass / (nMass + mass) * DiffusionCoefficientHoriztonal;
						float neighborSalinity = LastSalt[n] / (nMass + LastSalt[n]);

						newSaltMass += (neighborSalinity - salinity) * mass * diffusionAmount;
						newTemperature += (LastTemperature[n] - LastTemperature[i]) * diffusionAmount;
						newVelocity += (LastCurrent[n] - LastCurrent[i]) * diffusionAmount;
					}
				}
			}

			// NOTE: we don't diffuse velocity vertically
			float upMass = UpMass[i];
			if (upMass > 0)
			{
				// TODO: shouldn't shorter water columns diffuse faster?
				// If so we should divide (again) by total mass (which is proportional to height) of the shared cells...
				// ...clamp it to 1/6th (or so) of the total, and clamp the height to the diffusion range (10 meters?)
				float diffusionAmount = upMass / (upMass + mass)* DiffusionCoefficientVertical;
				float neighborSalinity = UpSalt[i] / (upMass + LastSalt[i]);
				newSaltMass += (neighborSalinity - salinity) * mass * diffusionAmount;
				newTemperature += (UpTemperature[i] - LastTemperature[i]) * diffusionAmount;
			}
			float downMass = DownMass[i];
			if (downMass > 0)
			{
				float diffusionAmount = downMass / (downMass + mass)* DiffusionCoefficientVertical;
				float neighborSalinity = DownSalt[i] / (downMass + DownSalt[i]);

				newSaltMass += (neighborSalinity - salinity) * mass * diffusionAmount;
				newTemperature += (DownTemperature[i] - LastTemperature[i]) * diffusionAmount;
			}

		}

		Delta[i] = new DiffusionWater()
		{
			Temperature = newTemperature,
			SaltMass = newSaltMass,
			Velocity = newVelocity
		};

	}
}

