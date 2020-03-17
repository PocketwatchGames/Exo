using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !GroundWaterDiffusionJobDebug
[BurstCompile]
#endif
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

#if !GroundWaterFlowJobDebug
[BurstCompile]
#endif
public struct GroundWaterFlowJob : IJobParallelFor {
	public NativeArray<float> GroundWater;
	public NativeArray<float> GroundWaterTemperature;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> LastGroundWaterTemperature;
	[ReadOnly] public NativeArray<float> Elevation;
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
		float saturation = 1.0f - LastGroundWater[i] * GroundWaterMaxInverse;
		float maxFlowPercentPerNeighbor = 0.1f;
		float elevation = Elevation[i];

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float terrainGradient = (Elevation[n] - elevation) * NeighborDistInverse[n];
				if (terrainGradient > 0)
				{
					float flowingIntoUs = LastGroundWater[n] * math.min(maxFlowPercentPerNeighbor, terrainGradient * saturation * FlowSpeed);

					newTemperature += LastGroundWaterTemperature[n] * flowingIntoUs;
					inFlow += flowingIntoUs;
				} else
				{
					float destSaturation = 1.0f - LastGroundWater[n] * GroundWaterMaxInverse;
					newGroundWater -= groundWater * math.min(maxFlowPercentPerNeighbor, -terrainGradient * destSaturation * FlowSpeed);
				}
			}
		}

		GroundWater[i] = newGroundWater + inFlow;
		GroundWaterTemperature[i] = (LastGroundWaterTemperature[i] * groundWater + newTemperature) / (inFlow + groundWater);
	}
}

#if !GroundWaterConductionJobDebug
[BurstCompile]
#endif
public struct GroundWaterConductionJob : IJobParallelFor {
	public NativeArray<float> GroundWaterTemperature;
	public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeArray<float> GroundWater;
	[ReadOnly] public NativeArray<float> LastGroundWaterTemperature;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public float GroundWaterConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		float groundWater = GroundWater[i];
		float newGroundWaterTemperature = LastGroundWaterTemperature[i];
		float newTerrainTemperature = TerrainTemperature[i];

		if (groundWater > 0)
		{
			float specificHeatTerrain = Atmosphere.GetSpecificHeatTerrain(HeatingDepth, Terrain[i].SoilFertility);
			float groundWaterTempDiff = newGroundWaterTemperature - TerrainTemperature[i];
			float specificHeatGroundWater = groundWater * WorldData.SpecificHeatWater;
			float conductionDelta = groundWaterTempDiff * GroundWaterConductionCoefficient * SecondsPerTick;

			// TODO: this can conduct heat past a point of equilibrium
			newGroundWaterTemperature -= conductionDelta / specificHeatGroundWater;
			newTerrainTemperature += conductionDelta / specificHeatTerrain;
		}

		GroundWaterTemperature[i] = newGroundWaterTemperature;
		TerrainTemperature[i] = newTerrainTemperature;

	}
}


#if !GroundWaterAbsorptionJobDebug
[BurstCompile]
#endif
public struct GroundWaterAbsorptionJob : IJobParallelFor {
	public NativeArray<float> GroundWater;
	public NativeArray<float> GroundWaterTemperature;
	public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	[ReadOnly] public NativeArray<float> LastGroundWaterTemperature;
	[ReadOnly] public NativeArray<float> WaterBelow;
	[ReadOnly] public NativeArray<float> WaterTemperature;
	[ReadOnly] public float GroundWaterAbsorptionRate;
	[ReadOnly] public float GroundWaterMaxInverse;
	public void Execute(int i)
	{
		float waterMass = WaterMass[i];
		float lastGroundWater = LastGroundWater[i];
		float newGroundWater = lastGroundWater;
		float newGroundWaterTemperature = LastGroundWaterTemperature[i];
		float newWaterMass = WaterMass[i];

		// Check if we are the floor layer
		// TODO: if we put all water layers into one large array, we can access any layer at any time (but may sacrifice ability to process different layers in parallel?)
		// TODO: we can probably find a way to shortcut all of this if a large number of our cells are fully saturated
		if (waterMass > 0 && WaterBelow[i] == 0)
		{
			float groundWaterAbsorbed = math.min(waterMass, math.max(0, (1.0f - lastGroundWater * GroundWaterMaxInverse) * GroundWaterAbsorptionRate));

			if (groundWaterAbsorbed > 0)
			{
				newGroundWaterTemperature = (lastGroundWater * newGroundWaterTemperature + groundWaterAbsorbed * WaterTemperature[i]) / (lastGroundWater + groundWaterAbsorbed);
				newGroundWater += groundWaterAbsorbed;
				waterMass -= groundWaterAbsorbed;
			}
		}

		GroundWater[i] = newGroundWater;
		GroundWaterTemperature[i] = newGroundWaterTemperature;
		WaterMass[i] = newWaterMass;

	}
}

