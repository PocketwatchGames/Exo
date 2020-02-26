#define DISABLE_VERTICAL_AIR_MOVEMENT
#define AdvectionAirJobDebug

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
	[ReadOnly] public NativeArray<float2> Coords;
	[ReadOnly] public NativeArray<float3> Position;
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
	[ReadOnly] public NativeArray<float> DownTemperature;
	[ReadOnly] public NativeArray<float> DownHumidity;
	[ReadOnly] public NativeArray<float> DownAirMass;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	public void Execute(int i)
	{
		float gradientTemperature = 0;
		float gradientWaterVapor = 0;
		float3 gradientWind = float3.zero;

		float3 velocity = Wind[i];
		float airMass = AirMass[i];
		float vapor = Vapor[i];
		float temperature = Temperature[i];
		float absoluteHumidity = vapor / (vapor + airMass);

#if !DISABLE_AIR_ADVECTION
		float2 coord = Coords[i];
		float3 pos = Position[i];

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float3 nPos = Position[n];
				float3 diff = math.normalize(pos - nPos);
				float3 nVelocity = Wind[n];

				float cosAngleBetweenCells = 0.5f;

				// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
				float velDotDir = math.clamp((math.max(0, (math.dot(nVelocity, diff)- cosAngleBetweenCells) / (1.0f - cosAngleBetweenCells)) + (math.max(0, math.dot(velocity, -diff))) / (1.0f - cosAngleBetweenCells)) * InverseCellDiameter * SecondsPerTick, 0, 0.5f);

				float neighborMass = AirMass[n] + Vapor[n];
				float diffusionAmount = AirMass[n] / (AirMass[n] + airMass);
				float neighborHumidity = Vapor[n] / neighborMass;
				gradientWaterVapor += (neighborHumidity - absoluteHumidity) * airMass * diffusionAmount * velDotDir;
				gradientWind += (nVelocity - velocity) * diffusionAmount * velDotDir;
				gradientTemperature += (Temperature[n] - Temperature[i]) * diffusionAmount * velDotDir;

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

#endif
		Delta[i] = new DiffusionAir()
		{
			Temperature = gradientTemperature,
			WaterVapor = gradientWaterVapor,
			Velocity = gradientWind,
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
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	public void Execute(int i)
	{
#if !DISABLE_CLOUD_ADVECTION
		float gradientMass = 0;
		float gradientDropletMass = 0;
		float3 gradientVelocity = float3.zero;
		float3 velocity = Velocity[i];
		float3 position = Position[i];
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

					float3 diff = position - Position[n];
					// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
					float velDotDir = math.clamp((math.max(0, math.dot(Velocity[n], diff)) + math.max(0, math.dot(velocity, -diff))) * InverseCellDiameter * SecondsPerTick, 0, 1.0f);
					float diffusionAmount = nMass / (nMass + mass);

					gradientMass += (Mass[n] - mass) * diffusionAmount * velDotDir;
					gradientDropletMass += (DropletMass[n] - dropletMass) * diffusionAmount * velDotDir;
					gradientVelocity += (Velocity[n] - velocity) * diffusionAmount * velDotDir;
				}
			}
		}

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
	[ReadOnly] public NativeArray<float3> Current;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float InverseCoordDiff;
	public void Execute(int i)
	{
		float gradientSalinity = 0;
		float gradientTemperature = 0;
		float3 gradientVelocity = float3.zero;
		float3 position = Position[i];
		float waterMass = Mass[i];
		if (waterMass > 0)
		{
			float salt = Salt[i];
			float mass = Mass[i];
			float temperature = Temperature[i];
			float3 velocity = Current[i];
			float salinity = salt / (salt + mass);
			for (int j = 0; j < 6; j++)
			{
				int n = Neighbors[i * 6 + j];
				if (n >= 0)
				{
					float nMass = Mass[n];
					if (nMass > 0)
					{
						float3 diff = position - Position[n];
						// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
						float velDotDir = math.clamp((math.max(0, math.dot(Current[n], diff)) + math.max(0, math.dot(velocity, -diff))) * InverseCellDiameter * SecondsPerTick, 0, 0.5f);
						float diffusionAmount = nMass / (nMass + mass);
						float neighborSalinity = Salt[n] / (nMass + Salt[n]);

						gradientSalinity += (neighborSalinity - salinity) * mass * diffusionAmount * velDotDir;
						gradientTemperature += (Temperature[n] - temperature) * diffusionAmount * velDotDir;
						gradientVelocity += (Current[n] - velocity) * diffusionAmount * velDotDir;
					}
				}
			}
		}

		Delta[i] = new DiffusionWater()
		{
			Temperature = gradientTemperature,
			SaltMass = gradientSalinity,
			Velocity = gradientVelocity
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


