﻿//#define DISABLE_VERTICAL_AIR_MOVEMENT
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
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> TemperatureAbove;
	[ReadOnly] public NativeArray<float> TemperatureBelow;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float> VaporAbove;
	[ReadOnly] public NativeArray<float> VaporBelow;
	[ReadOnly] public NativeArray<float> CarbonDioxide;
	[ReadOnly] public NativeArray<float> CarbonDioxideAbove;
	[ReadOnly] public NativeArray<float> CarbonDioxideBelow;
	[ReadOnly] public NativeArray<float> Dust;
	[ReadOnly] public NativeArray<float> DustAbove;
	[ReadOnly] public NativeArray<float> DustBelow;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> AirVelocityAbove;
	[ReadOnly] public NativeArray<float3> AirVelocityBelow;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirMassAbove;
	[ReadOnly] public NativeArray<float> AirMassBelow;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerHeightAbove;
	[ReadOnly] public NativeArray<float> LayerHeightBelow;
	[ReadOnly] public NativeArray<float> LayerElevationAbove;
	[ReadOnly] public NativeArray<float> LayerElevationBelow;
	[ReadOnly] public NativeArray<int> Neighbors;
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


#if !DISABLE_VERTICAL_AIR_MOVEMENT

		// NOTE: we don't diffuse velocity vertically
		if (!IsTop)
		{
			float nMass = AirMassAbove[i];
			float nInverseMass = 1.0f / nMass;
			float saToVolume = 1.0f / LayerHeightAbove[i];
			float inverseTotalMass = 1.0f / (nMass + mass);
			float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
			neighborTemperature += (TemperatureAbove[i] - temperature) * diffusion;
			neighborVapor += (VaporAbove[i] * nInverseMass - vaporPercent) * diffusion * math.min(nMass, mass);
			neighborCO2 += (CarbonDioxideAbove[i] * nInverseMass - co2Percent) * diffusion * math.min(nMass, mass);
			neighborDust += (DustAbove[i] * nInverseMass - dustPercent) * diffusion * math.min(nMass, mass);
		}
		if (!IsBottom)
		{
			float nMass = AirMassBelow[i];
			float nInverseMass = 1.0f / nMass;
			float saToVolume = 1.0f / LayerHeightBelow[i];
			float inverseTotalMass = 1.0f / (nMass + mass);
			float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
			neighborTemperature += (TemperatureBelow[i] - temperature) * diffusion;
			neighborVapor += (VaporBelow[i] * nInverseMass - vaporPercent) * diffusion * math.min(nMass, mass);
			neighborCO2 += (CarbonDioxideBelow[i] * nInverseMass - co2Percent) * diffusion * math.min(nMass, mass);
			neighborDust += (DustBelow[i] * nInverseMass - dustPercent) * diffusion * math.min(nMass, mass);
		}

#endif
#endif

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
			DropletMass = newDropletMass * DiffusionCoefficient + dropletMass,
			Velocity = newVelocity * DiffusionCoefficient + velocity,
		};

	}
}

[BurstCompile]
public struct DiffusionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> WaterMassAbove;
	[ReadOnly] public NativeArray<float> WaterMassBelow;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> TemperatureAbove;
	[ReadOnly] public NativeArray<float> TemperatureBelow;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> SaltMassAbove;
	[ReadOnly] public NativeArray<float> SaltMassBelow;
	[ReadOnly] public NativeArray<float> CarbonMass;
	[ReadOnly] public NativeArray<float> CarbonMassAbove;
	[ReadOnly] public NativeArray<float> CarbonMassBelow;
	[ReadOnly] public NativeArray<float> PlanktonMass;
	[ReadOnly] public NativeArray<float> PlanktonGlucose;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> VelocityAbove;
	[ReadOnly] public NativeArray<float3> VelocityBelow;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerHeightAbove;
	[ReadOnly] public NativeArray<float> LayerHeightBelow;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	[ReadOnly] public float DiffusionCoefficientHorizontal;
	[ReadOnly] public float DiffusionCoefficientVertical;
	[ReadOnly] public float CellSurfaceArea;
	[ReadOnly] public float CellCircumference;
	public void Execute(int i)
	{
		float mass = WaterMass[i];
		float neighborTemperature = Temperature[i];
		float neighborSaltMass = SaltMass[i];
		float neighborCarbonMass = CarbonMass[i];
		float neighborPlankton = PlanktonMass[i];
		float neighborGlucose = PlanktonGlucose[i];
		float3 neighborVelocity = Velocity[i];


#if !DISABLE_WATER_DIFFUSION
		if (mass > 0)
		{
			float layerHeight = LayerHeight[i];
			float temperature = Temperature[i];
			float saltMass = SaltMass[i];
			float carbonMass = CarbonMass[i];
			float plankton = PlanktonMass[i];
			float glucose = PlanktonGlucose[i];
			float3 velocity = Velocity[i];
			float saltPercent = saltMass / mass;
			float carbonPercent = carbonMass / mass;
			float planktonPercent = plankton / mass;
			float glucosePercent = glucose / mass;


			for (int j = 0; j < 6; j++)
			{
				int nIndex = i * 6 + j;
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
							float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientHorizontal;
							neighborTemperature += (Temperature[n] - temperature) * diffusion;
							neighborVelocity += (Velocity[n] - velocity) * diffusion;
							neighborSaltMass += (SaltMass[n] / nMass - saltPercent) * diffusion * math.min(nMass, mass);
							neighborPlankton += (PlanktonMass[n] / nMass - planktonPercent) * diffusion * math.min(nMass, mass);
							neighborGlucose += (PlanktonGlucose[n] / nMass - glucosePercent) * diffusion * math.min(nMass, mass);
							neighborCarbonMass += (CarbonMass[n] / nMass - carbonPercent) * diffusion * math.min(nMass, mass);
						}
					}
				}
			}


			// NOTE: we don't diffuse velocity vertically
			// NOTE: we are ignoring adibatic processes in the water -- at 10km, the total lapse is less than 1.5 degrees celsius
			{
				float nMass = WaterMassAbove[i];
				if (nMass > 0)
				{
					float nLayerHeight = LayerHeightAbove[i];
					if (nLayerHeight > 0)
					{
						float saToVolume = 1.0f / nLayerHeight;
						float inverseTotalMass = 1.0f / (nMass + mass);
						float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
						neighborTemperature += (TemperatureAbove[i] - temperature) * diffusion;
						neighborSaltMass += (SaltMassAbove[i] / nMass - saltPercent) * diffusion * math.min(nMass, mass);
						neighborCarbonMass += (CarbonMassAbove[i] / nMass - carbonPercent) * diffusion * math.min(nMass, mass);
					}
				}
			}

			{
				float nMass = WaterMassBelow[i];
				if (nMass > 0)
				{
					float nLayerHeight = LayerHeightBelow[i];
					if (nLayerHeight > 0)
					{
						float saToVolume = 1.0f / nLayerHeight;
						float inverseTotalMass = 1.0f / (nMass + mass);
						float diffusion = saToVolume * nMass * inverseTotalMass * DiffusionCoefficientVertical;
						neighborTemperature += (TemperatureBelow[i] - temperature) * diffusion;
						neighborSaltMass += (SaltMassBelow[i] / nMass - saltPercent) * diffusion * math.min(nMass, mass);
						neighborCarbonMass += (CarbonMassBelow[i] / nMass - carbonPercent) * diffusion * math.min(nMass, mass);
					}
				}
			}

		}
#endif

		Delta[i] = new DiffusionWater()
		{
			Temperature = neighborTemperature,
			SaltMass = neighborSaltMass,
			CarbonMass = neighborCarbonMass,
			Plankton = neighborPlankton,
			PlanktonGlucose = neighborGlucose,
			Velocity = neighborVelocity,
		};
	}

}

