//#define DISABLE_VERTICAL_AIR_MOVEMENT
//#define DISABLE_AIR_DIFFUSION
//#define DISABLE_CLOUD_DIFFUSION
//#define DISABLE_WATER_DIFFUSION

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

[BurstCompile]
public struct DiffusionAirJob : IJobParallelFor {
	public NativeSlice<DiffusionAir> Delta;
	[ReadOnly] public NativeSlice<float> Temperature;
	[ReadOnly] public NativeSlice<float> Vapor;
	[ReadOnly] public NativeSlice<float> CarbonDioxide;
	[ReadOnly] public NativeSlice<float> Dust;
	[ReadOnly] public NativeSlice<float3> Velocity;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeSlice<float> LayerHeight;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	[ReadOnly] public float DiffusionCoefficientHorizontal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	[ReadOnly] public float CellSurfaceArea;
	[ReadOnly] public float CellCircumference;
	[ReadOnly] public int LayerCount;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		float layerHeight = LayerHeight[i];
		float mass = AirMass[i];
		float temperature = Temperature[i];
		float vapor = Vapor[i];
		float co2 = CarbonDioxide[i];
		float dust = Dust[i];
		float3 velocity = Velocity[i];

		float inverseMass = 1.0f / mass;
		float vaporPercent = vapor * inverseMass;
		float co2Percent = co2 * inverseMass;
		float dustPercent = dust * inverseMass;

		float neighborTemperature = temperature;
		float neighborDust = dust;
		float neighborVapor = vapor;
		float neighborCO2 = co2;
		float3 neighborVelocity = velocity;

		int layer = i / Count;
		bool isTop = layer == LayerCount - 1;
		bool isBottom = layer == 0;
		int cellIndex = i - layer * Count;

		if (mass > 0)
		{
			for (int j = 0; j < StaticState.MaxNeighbors; j++)
			{
				int nIndex = cellIndex * 6 + j;
				int n = Neighbors[nIndex];
				if (n >= 0)
				{
					n += layer * Count;
					float nMass = AirMass[n];
					float nInverseMass = 1.0f / nMass;
					float saToVolume = math.min(layerHeight, LayerHeight[n]) * CellCircumference / (6 * LayerHeight[n] * CellSurfaceArea);
					float inverseTotalMass = 1.0f / (nMass + mass);
					float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientHorizontal;
					neighborTemperature += (Temperature[n] - temperature) * diffusion;
					neighborVelocity += (Velocity[n] - velocity) * diffusion;
					neighborVapor += (Vapor[n] * nInverseMass - vaporPercent) * diffusion * math.min(nMass, mass); // TODO: does this actually have conservation of mass?
					neighborCO2 += (CarbonDioxide[n] * nInverseMass - co2Percent) * diffusion * math.min(nMass, mass); // TODO: does this actually have conservation of mass?
					neighborDust += (Dust[n] * nInverseMass - dustPercent) * diffusion * math.min(nMass, mass);
				}
			}
		}



		// NOTE: we don't diffuse velocity vertically
		if (!isTop)
		{
			int n = i + Count;
			float nMass = AirMass[n];
			float nInverseMass = 1.0f / nMass;
			float saToVolume = 1.0f / LayerHeight[n];
			float inverseTotalMass = 1.0f / (nMass + mass);
			float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
			neighborTemperature += (Temperature[n] - temperature) * diffusion;
			neighborVapor += (Vapor[n] * nInverseMass - vaporPercent) * diffusion * math.min(nMass, mass);
			neighborCO2 += (CarbonDioxide[n] * nInverseMass - co2Percent) * diffusion * math.min(nMass, mass);
			neighborDust += (Dust[n] * nInverseMass - dustPercent) * diffusion * math.min(nMass, mass);
		}
		if (!isBottom)
		{
			int n = i - Count;
			float nMass = AirMass[n];
			float nInverseMass = 1.0f / nMass;
			float saToVolume = 1.0f / LayerHeight[n];
			float inverseTotalMass = 1.0f / (nMass + mass);
			float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
			neighborTemperature += (Temperature[n] - temperature) * diffusion;
			neighborVapor += (Vapor[n] * nInverseMass - vaporPercent) * diffusion * math.min(nMass, mass);
			neighborCO2 += (CarbonDioxide[n] * nInverseMass - co2Percent) * diffusion * math.min(nMass, mass);
			neighborDust += (Dust[n] * nInverseMass - dustPercent) * diffusion * math.min(nMass, mass);
		}


		Delta[i] = new DiffusionAir()
		{
			Temperature = neighborTemperature,
			Velocity = neighborVelocity,
			WaterVapor = neighborVapor,
			CarbonDioxide = neighborCO2,
			Dust = neighborDust,
		};

	}
}


[BurstCompile]
public struct DiffusionWaterJob : IJobParallelFor {
	public NativeSlice<DiffusionWater> Delta;
	[ReadOnly] public NativeSlice<float> WaterMass;
	[ReadOnly] public NativeSlice<float> Temperature;
	[ReadOnly] public NativeSlice<float> SaltMass;
	[ReadOnly] public NativeSlice<float> CarbonMass;
	[ReadOnly] public NativeSlice<float> PlanktonMass;
	[ReadOnly] public NativeSlice<float> PlanktonGlucose;
	[ReadOnly] public NativeSlice<float3> Velocity;
	[ReadOnly] public NativeSlice<float> LayerHeight;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	[ReadOnly] public float DiffusionCoefficientHorizontal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	[ReadOnly] public float CellSurfaceArea;
	[ReadOnly] public float CellCircumference;
	[ReadOnly] public int Count;
	[ReadOnly] public int LayerCount;
	public void Execute(int i)
	{
		int indexUp = i + Count;
		int indexDown = i - Count;

		float mass = WaterMass[i];
		float layerHeight = LayerHeight[i];
		float temperature = Temperature[i];
		float saltMass = SaltMass[i];
		float carbonMass = CarbonMass[i];
		float plankton = PlanktonMass[i];
		float glucose = PlanktonGlucose[i];
		float3 velocity = Velocity[i];

		float neighborTemperature = temperature;
		float neighborSaltMass = saltMass;
		float neighborCarbonMass = carbonMass;
		float neighborPlankton = plankton;
		float neighborGlucose = glucose;
		float3 neighborVelocity = velocity;

		int layer = i / Count;
		bool isTop = layer == LayerCount - 1;
		bool isBottom = layer == 0;
		int cellIndex = i - layer * Count;


		if (mass > 0)
		{
			float inverseMass = 1.0f / mass;
			float saltPercent = saltMass * inverseMass;
			float carbonPercent = carbonMass * inverseMass;
			float planktonPercent = plankton * inverseMass;
			float glucosePercent = glucose * inverseMass;


			for (int j = 0; j < StaticState.MaxNeighbors; j++)
			{
				int nIndex = cellIndex * 6 + j;
				int n = Neighbors[nIndex];
				if (n >= 0)
				{
					float nHeight = LayerHeight[n];
					if (nHeight > 0)
					{
						float nMass = WaterMass[n];
						if (nMass > 0)
						{
							float saToVolume = math.min(layerHeight, nHeight) * CellCircumference / (6 * nHeight * CellSurfaceArea);
							float inverseTotalMass = 1.0f / (nMass + mass);
							float inverseNMass = 1.0f / nMass;
							float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientHorizontal;
							neighborTemperature += (Temperature[n] - temperature) * diffusion;
							neighborVelocity += (Velocity[n] - velocity) * diffusion;
							neighborSaltMass += (SaltMass[n] * inverseNMass - saltPercent) * diffusion * math.min(nMass, mass);
							neighborPlankton += (PlanktonMass[n] * inverseNMass - planktonPercent) * diffusion * math.min(nMass, mass);
							neighborGlucose += (PlanktonGlucose[n] * inverseNMass - glucosePercent) * diffusion * math.min(nMass, mass);
							neighborCarbonMass += (CarbonMass[n] * inverseNMass - carbonPercent) * diffusion * math.min(nMass, mass);
						}
					}
				}
			}


			// NOTE: we don't diffuse velocity vertically
			// NOTE: we are ignoring adibatic processes in the water -- at 10km, the total lapse is less than 1.5 degrees celsius
			if (!isTop)
			{
				float nMass = WaterMass[indexUp];
				if (nMass > 0)
				{
					float nLayerHeight = LayerHeight[indexUp];
					if (nLayerHeight > 0)
					{
						float inverseNMass = 1.0f / nMass;
						float saToVolume = 1.0f / nLayerHeight;
						float inverseTotalMass = 1.0f / (nMass + mass);
						float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
						neighborTemperature += (Temperature[indexUp] - temperature) * diffusion;
						neighborSaltMass += (SaltMass[indexUp] * inverseNMass - saltPercent) * diffusion * math.min(nMass, mass);
						neighborCarbonMass += (CarbonMass[indexUp] * inverseNMass - carbonPercent) * diffusion * math.min(nMass, mass);
					}
				}
			}

			if (!isBottom)
			{
				float nMass = WaterMass[indexDown];
				if (nMass > 0)
				{
					float nLayerHeight = LayerHeight[indexDown];
					if (nLayerHeight > 0)
					{
						float inverseNMass = 1.0f / nMass;
						float saToVolume = 1.0f / nLayerHeight;
						float inverseTotalMass = 1.0f / (nMass + mass);
						float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
						neighborTemperature += (Temperature[indexDown] - temperature) * diffusion;
						neighborSaltMass += (SaltMass[indexDown] * inverseNMass - saltPercent) * diffusion * math.min(nMass, mass);
						neighborCarbonMass += (CarbonMass[indexDown] * inverseNMass - carbonPercent) * diffusion * math.min(nMass, mass);
					}
				}
			}

		}

		Delta[i] = new DiffusionWater()
		{
			WaterMass = mass,
			SaltMass = neighborSaltMass,
			CarbonMass = neighborCarbonMass,
			Plankton = neighborPlankton,
			PlanktonGlucose = neighborGlucose,
			Temperature = neighborTemperature,
			Velocity = neighborVelocity,
		};
	}

}


[BurstCompile]
public struct DiffusionCloudJob : IJobParallelFor {
	public NativeArray<DiffusionCloud> Delta;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastDropletMass;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public float DiffusionCoefficient;
	public void Execute(int i)
	{
		float mass = LastMass[i];
		float temperature = LastTemperature[i];
		float dropletMass = LastDropletMass[i];

		float newMass = 0;
		float newTemperature = 0;
		float newDropletMass = 0;

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
				}
			}
		}
#endif

		Delta[i] = new DiffusionCloud()
		{
			Mass = newMass * DiffusionCoefficient + mass,
			Temperature = newTemperature * DiffusionCoefficient + temperature,
			DropletMass = newDropletMass * DiffusionCoefficient + dropletMass,
		};

	}
}
