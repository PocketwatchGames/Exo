
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

public struct DiffusionAir {
	public float Temperature;
	public float Humidity;
	public float2 Velocity;
	public float VelocityVertical;
}
public struct DiffusionCloud {
	public float Mass;
	public float DropletMass;
	public float2 Velocity;
}
public struct DiffusionWater {
	public float Temperature;
	public float Salinity;
	public float2 Velocity;
}


#if !AdvectionAirJobDebug
[BurstCompile]
#endif
public struct AdvectionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float2> Wind;
	[ReadOnly] public NativeArray<float> WindVertical;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float2> Coords;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpHumidity;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> UpWindVertical;
	[ReadOnly] public NativeArray<float> DownTemperature;
	[ReadOnly] public NativeArray<float> DownHumidity;
	[ReadOnly] public NativeArray<float> DownAirMass;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public NativeArray<float> DownWindVertical;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	[ReadOnly] public float DiffusionCoefficientHoriztonal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	public void Execute(int i)
	{
		float gradientTemperature = 0;
		float gradientWaterVapor = 0;
		float2 gradientVelocity = float2.zero;
		float gradientWindVertical = 0;

		float2 velocity = Wind[i];
		float airMass = AirMass[i];
		float vapor = Vapor[i];
		float temperature = Temperature[i];
		float absoluteHumidity = vapor / (vapor + airMass);

#if !DISABLE_AIR_ADVECTION
		float2 coord = Coords[i];
		//TODO: account for different size air columns, similar to water
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float2 coordDiff = coord - Coords[n];
				float2 diff = math.normalize(math.float2(Utils.WrapAngle(coordDiff.x), Utils.WrapAngle(coordDiff.y)));
				// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
				float velDotDir = math.clamp((math.max(0, math.dot(Wind[n], diff)) + math.max(0, math.dot(velocity, -diff))) * InverseCellDiameter * SecondsPerTick, 0, 0.5f);

				float neighborMass = AirMass[n] + Vapor[n];
				float diffusionAmount = AirMass[n] / (AirMass[n] + airMass);
				float neighborHumidity = Vapor[n] / neighborMass;
				gradientWaterVapor += (neighborHumidity - absoluteHumidity) * airMass * diffusionAmount  * math.min(1, DiffusionCoefficientHoriztonal + velDotDir);
				gradientVelocity += (Wind[n] - Wind[i]) * diffusionAmount * math.min(1, DiffusionCoefficientHoriztonal + velDotDir);
				gradientTemperature += (Temperature[n] - Temperature[i]) * diffusionAmount * math.min(1, DiffusionCoefficientHoriztonal + velDotDir);

			}
		}

#endif


#if !DISABLE_VERTICAL_AIR_MOVEMENT

		if (!IsTop)
		{
			float heightDiff = (UpLayerElevation[i] + UpLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
			float combinedWind = math.min(1, math.max(0, WindVertical[i]) - math.min(0, UpWindVertical[i])) * SecondsPerTick / heightDiff;

			float diffusionAmount = UpAirMass[i] / (UpAirMass[i] + airMass);

			float absoluteHumidityUp = UpHumidity[i] / (UpHumidity[i] + UpAirMass[i]);
			gradientWaterVapor += (absoluteHumidityUp - absoluteHumidity) * diffusionAmount * combinedWind;

			float potentialTemperatureUp = UpTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
			gradientTemperature += (potentialTemperatureUp - Temperature[i]) * diffusionAmount * combinedWind;

			gradientWindVertical += (UpWindVertical[i] - WindVertical[i]) * diffusionAmount * combinedWind;

		}
		if (!IsBottom)
		{
			float heightDiff = (DownLayerElevation[i] + DownLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
			float combinedWind = math.min(1, math.max(0, DownWindVertical[i]) - math.min(0, WindVertical[i])) * SecondsPerTick / heightDiff;

			float diffusionAmount = DownAirMass[i] / (DownAirMass[i] + airMass);

			float absoluteHumidityDown = DownHumidity[i] / (DownHumidity[i] + DownAirMass[i]);
			gradientWaterVapor += (absoluteHumidityDown - absoluteHumidity) * diffusionAmount * math.min(1, combinedWind + DiffusionCoefficientVertical);

			float potentialTemperatureDown = DownTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
			gradientTemperature += (potentialTemperatureDown - Temperature[i]) * diffusionAmount * math.min(1, combinedWind + DiffusionCoefficientVertical);

			gradientWindVertical += (DownWindVertical[i] - WindVertical[i]) * diffusionAmount * math.min(1, combinedWind + DiffusionCoefficientVertical);
		}

		//		float moveToNeutralBuoyancy = (UpTemperature[i] - Temperature[i]) / WorldData.TemperatureLapseRate - heightDiff;
		//		float vertMovement = math.min(MaxVerticalMovement, math.clamp(moveToNeutralBuoyancy + DiffusionCoefficient, 0, 1));

		gradientTemperature *= 0.3f;
		gradientWaterVapor *= 0.3f;
		gradientVelocity *= 0.3f;

#endif
		Delta[i] = new DiffusionAir()
		{
			Temperature = gradientTemperature,
			Humidity = gradientWaterVapor,
			Velocity = gradientVelocity,
			VelocityVertical = gradientWindVertical
		};

	}
}


#if !AdvectionCloudJob
[BurstCompile]
#endif
public struct AdvectionCloudJob : IJobParallelFor {
	public NativeArray<DiffusionCloud> Delta;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> DropletMass;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float2> Coords;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float DiffusionCoefficient;
	public void Execute(int i)
	{
#if !DISABLE_CLOUD_ADVECTION
		float gradientMass = 0;
		float gradientDropletMass = 0;
		float2 gradientVelocity = float2.zero;
		float2 velocity = Velocity[i];
		float2 coord = Coords[i];
		float mass = Mass[i];
		float dropletMass = DropletMass[i];

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float nMass = Mass[n];
				float totalMass = nMass + mass;
				if (totalMass > 0)
				{

					float2 coordDiff = coord - Coords[n];
					float2 diff = math.normalize(math.float2(Utils.WrapAngle(coordDiff.x), Utils.WrapAngle(coordDiff.y)));
					// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
					float velDotDir = math.clamp((math.max(0, math.dot(Velocity[n], diff)) + math.max(0, math.dot(velocity, -diff))) * InverseCellDiameter * SecondsPerTick, 0, 0.5f);
					float diffusionAmount = nMass / (nMass + mass);

					// TODO: account for differing masses having smaller or larger effects on cloud properties
					gradientMass += (Mass[n] - mass) * diffusionAmount * math.min(1, DiffusionCoefficient + velDotDir);
					gradientDropletMass += (DropletMass[n] - dropletMass) * diffusionAmount * math.min(1, DiffusionCoefficient + velDotDir);
					gradientVelocity += (Velocity[n] - velocity) * diffusionAmount * math.min(1, DiffusionCoefficient + velDotDir);


				}
			}
		}
		float outgoing = math.length(velocity) * InverseCellDiameter;
		gradientMass -= Mass[i] * outgoing;
		gradientDropletMass -= DropletMass[i] * outgoing;
		gradientVelocity -= velocity * outgoing;

		Delta[i] = new DiffusionCloud()
		{
			Mass = gradientMass,
			DropletMass = gradientDropletMass,
			Velocity = gradientVelocity
		};

#endif
	}
}

#if !AdvectionWaterJobDebug
[BurstCompile]
#endif
public struct AdvectionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float2> Current;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float2> Coords;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float InverseCoordDiff;
	[ReadOnly] public float DiffusionCoefficientHoriztonal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	public void Execute(int i)
	{
		float gradientSalinity = 0;
		float gradientTemperature = 0;
		float2 gradientVelocity = float2.zero;
		float2 coord = Coords[i];
		float waterMass = Mass[i];
		if (waterMass > 0)
		{
			float salt = Salt[i];
			float mass = Mass[i];
			float temperature = Temperature[i];
			float2 velocity = Current[i];
			float salinity = salt / (salt + mass);
			for (int j = 0; j < 6; j++)
			{
				int n = Neighbors[i * 6 + j];
				if (n >= 0)
				{
					float nMass = Mass[n];
					if (nMass > 0)
					{
						float2 coordDiff = coord - Coords[n];
						float2 diff = math.normalize(math.float2(Utils.WrapAngle(coordDiff.x), Utils.WrapAngle(coordDiff.y)));
						// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
						float velDotDir = math.clamp((math.max(0, math.dot(Current[n], diff)) + math.max(0, math.dot(velocity, -diff))) * InverseCellDiameter * SecondsPerTick, 0, 0.5f);
						float diffusionAmount = nMass / (nMass + mass);
						float neighborSalinity = Salt[n] / (nMass + Salt[n]);

						gradientSalinity += (neighborSalinity - salinity) * mass * diffusionAmount * math.min(1, velDotDir + DiffusionCoefficientHoriztonal);
						gradientTemperature += (Temperature[n] - temperature) * diffusionAmount * math.min(1, velDotDir + DiffusionCoefficientHoriztonal);
						gradientVelocity += (Current[n] - velocity) * diffusionAmount * math.min(1, velDotDir + DiffusionCoefficientHoriztonal);
					}
				}
			}
		}

		Delta[i] = new DiffusionWater()
		{
			Temperature = gradientTemperature * DiffusionCoefficientHoriztonal,
			Salinity = gradientSalinity * DiffusionCoefficientHoriztonal,
			Velocity = gradientVelocity * DiffusionCoefficientHoriztonal
		};


	}
}


#if !AdvectionAirJobDebug
[BurstCompile]
#endif
public struct ApplyAdvectionAirJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Vapor;
	public NativeArray<float2> Wind;
	public NativeArray<float> WindVertical;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	public void Execute(int i)
	{
		WindVertical[i] += Advection[i].VelocityVertical;
		Temperature[i] += Advection[i].VelocityVertical;
		Vapor[i] += Advection[i].Humidity;
		Wind[i] += Advection[i].Velocity;
	}
}


#if !AdvectionWaterJobDebug
[BurstCompile]
#endif
public struct ApplyAdvectionWaterJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> SaltMass;
	public NativeArray<float2> Velocity;
	public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<DiffusionWater> Advection;
	public void Execute(int i)
	{
		SaltMass[i] += Advection[i].Salinity;
		Temperature[i] += Advection[i].Temperature;
		Velocity[i] += Advection[i].Velocity;
	}
}


#if !AdvectionCloudJobDebug
[BurstCompile]
#endif
public struct ApplyAdvectionCloudJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> DropletMass;
	public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<DiffusionCloud> Advection;
	public void Execute(int i)
	{
		CloudMass[i] += Advection[i].Mass;
		DropletMass[i] += Advection[i].DropletMass;
		Velocity[i] += Advection[i].Velocity;
	}
}


