//#define DISABLE_GROUND_WATER_ABSORPTION

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


[BurstCompile]
public struct GroundWaterDiffusionJob : IJobParallelFor {
	public NativeArray<float> GroundWater;
	public NativeArray<float> GroundWaterTemperature;
	[ReadOnly] public NativeArray<float> LastGroundWaterTemperature;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> NeighborDist;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public float DiffusionCoefficient;
	public void Execute(int i)
	{
		float lastGroundWater = LastGroundWater[i];
		float lastGroundWaterTemperature = LastGroundWaterTemperature[i];
		float newGroundWater = lastGroundWater;
		float newTemperature = lastGroundWaterTemperature;

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float neighborGroundWater = LastGroundWater[n];
				newGroundWater += (neighborGroundWater - lastGroundWater) * DiffusionCoefficient * NeighborDistInverse[neighborIndex];
				float diffusionAmount = Atmosphere.GetDiffusionAmount(lastGroundWater, neighborGroundWater, DiffusionCoefficient, NeighborDist[neighborIndex]);
				newTemperature += (LastGroundWaterTemperature[n] - lastGroundWaterTemperature) * diffusionAmount;
			}
		}

		GroundWaterTemperature[i] = newTemperature;
		GroundWater[i] = newGroundWater;

	}
}

[BurstCompile]
public struct GroundWaterFlowJob : IJobParallelFor {
	public NativeArray<float> GroundWater;
	public NativeArray<float> GroundWaterTemperature;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> LastGroundWaterTemperature;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public float FlowSpeed;
	[ReadOnly] public float GroundWaterMaxInverse;
	public void Execute(int i)
	{
		float groundWater = LastGroundWater[i];
		float newGroundWater = groundWater;
		float inFlow = 0;
		float newTemperature = 0;
		float maxFlowPercentPerNeighbor = 0.1f;
		float elevation = SurfaceElevation[i];

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float terrainGradient = (SurfaceElevation[n] - elevation) * NeighborDistInverse[n];
				if (terrainGradient > 0)
				{
					float flowingIntoUs = LastGroundWater[n] * math.min(maxFlowPercentPerNeighbor, terrainGradient * FlowSpeed);

					newTemperature += LastGroundWaterTemperature[n] * flowingIntoUs;
					inFlow += flowingIntoUs;
				} else
				{
					newGroundWater -= groundWater * math.min(maxFlowPercentPerNeighbor, -terrainGradient * FlowSpeed);
				}
			}
		}

		GroundWater[i] = newGroundWater + inFlow;
		GroundWaterTemperature[i] = (LastGroundWaterTemperature[i] * groundWater + newTemperature) / (inFlow + groundWater);
	}
}

[BurstCompile]
public struct GroundWaterConductionJob : IJobParallelFor {
	public NativeArray<float> GroundWaterTemperature;
	public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeArray<float> GroundWater;
	[ReadOnly] public NativeArray<float> LastGroundWaterTemperature;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public float GroundWaterConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float HeatingDepth;
	[ReadOnly] public float GroundWaterSurfaceAreaInverse;
	public void Execute(int i)
	{
		float groundWater = GroundWater[i];
		float newGroundWaterTemperature = LastGroundWaterTemperature[i];
		float newTerrainTemperature = TerrainTemperature[i];

		if (groundWater > 0)
		{
			float specificHeatTerrain = Atmosphere.GetSpecificHeatTerrain(HeatingDepth, SoilFertility[i]);
			float groundWaterTempDiff = newGroundWaterTemperature - TerrainTemperature[i];
			float specificHeatGroundWater = groundWater * WorldData.SpecificHeatWater;
			float conductionDelta = groundWaterTempDiff * GroundWaterConductionCoefficient * SecondsPerTick * GroundWaterSurfaceAreaInverse;

			// TODO: this can conduct heat past a point of equilibrium
			newGroundWaterTemperature -= conductionDelta / specificHeatGroundWater;
			newTerrainTemperature += conductionDelta / specificHeatTerrain;
		}

		GroundWaterTemperature[i] = newGroundWaterTemperature;
		TerrainTemperature[i] = newTerrainTemperature;

	}
}


[BurstCompile]
public struct GroundWaterAbsorptionJob : IJobParallelFor {
	public NativeArray<float> GroundWater;
	public NativeArray<float> GroundWaterTemperature;
	public NativeArray<float> WaterMass;
	public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> LastGroundWaterTemperature;
	[ReadOnly] public NativeArray<float> WaterBelow;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public float GroundWaterAbsorptionRate;
	[ReadOnly] public float GroundWaterMax;
	[ReadOnly] public float GroundWaterMaxInverse;
	[ReadOnly] public bool IsTop;
	public void Execute(int i)
	{
		float waterMass = WaterMass[i];
		float lastGroundWater = LastGroundWater[i];
		float newGroundWater = lastGroundWater;
		float newGroundWaterTemperature = LastGroundWaterTemperature[i];
		float newWaterMass = WaterMass[i];
		float newWaterTemperature = WaterTemperature[i];

#if !DISABLE_GROUND_WATER_ABSORPTION
		// Check if we are the floor layer
		// TODO: if we put all water layers into one large array, we can access any layer at any time (but may sacrifice ability to process different layers in parallel?)
		// TODO: we can probably find a way to shortcut all of this if a large number of our cells are fully saturated
		if (newGroundWater > GroundWaterMax)
		{
			if ((waterMass > 0 || IsTop) && WaterBelow[i] == 0)
			{
				float waterEmerged = newGroundWater - GroundWaterMax;
				newWaterTemperature = (newWaterTemperature * (newWaterMass * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt) + newGroundWaterTemperature * waterEmerged * WorldData.SpecificHeatWater) / ((newWaterMass + waterEmerged) * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt);
				newWaterMass += waterEmerged;
				newGroundWater = GroundWaterMax;
			}
		}
		else if (waterMass > 0 && WaterBelow[i] == 0)
		{
			float groundWaterAbsorbed = math.min(waterMass, math.max(0, (1.0f - lastGroundWater * GroundWaterMaxInverse) * GroundWaterAbsorptionRate));

			if (groundWaterAbsorbed > 0)
			{
				newGroundWaterTemperature = (lastGroundWater * newGroundWaterTemperature + groundWaterAbsorbed * WaterTemperature[i]) / (lastGroundWater + groundWaterAbsorbed);
				newGroundWater += groundWaterAbsorbed;
				newWaterMass -= groundWaterAbsorbed;
			}
		}
#endif

		GroundWater[i] = newGroundWater;
		GroundWaterTemperature[i] = newGroundWaterTemperature;
		WaterMass[i] = newWaterMass;

	}
}

