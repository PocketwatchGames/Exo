//#define UpdateDependentWaterLayerJobDebug
//#define UpdateDependentAirLayerJobDebug
//#define UpdateWaterDepthJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using System.Collections.Generic;

public static class SimJobs {

	public static JobHandle UpdateWaterDepths(JobHelper jobHelper, ref SimState state, ref DependentState dependent, ref WorldData worldData, JobHandle dependencies, List<NativeArray<float>> arraysToDispose)
	{
		var waterMassSum = new NativeArray<float>(state.IceMass, Allocator.TempJob);
		var waterDepthSum = new NativeArray<float>(state.IceMass.Length, Allocator.TempJob);
		arraysToDispose.Add(waterMassSum);
		arraysToDispose.Add(waterDepthSum);
		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			dependencies = jobHelper.Run(new UpdateWaterDepthJob()
			{
				Density = dependent.WaterDensity[j],
				Pressure = dependent.WaterPressure[j],
				LayerDepth = dependent.WaterLayerDepth[j],
				LayerHeight = dependent.WaterLayerHeight[j],
				WaterMassTotal = dependent.WaterMassTotal,
				WaterDepthTotal = dependent.WaterDepthTotal,
				WaterMassSum = waterMassSum,
				WaterDepthSum = waterDepthSum,

				Temperature = state.WaterTemperature[j],
				SaltMass = state.SaltMass[j],
				WaterMass = state.WaterMass[j],
				UpLayerDepth = dependent.WaterLayerDepth[j + 1],
				UpLayerHeight = dependent.WaterLayerHeight[j + 1],
				WaterDensityPerDegree = worldData.WaterDensityPerDegree,
				WaterDensityPerSalinity = worldData.WaterDensityPerSalinity,
				Gravity = state.PlanetState.Gravity,

			}, dependencies);
		}
		return dependencies;
	}

}

#region Update dependents


#if !UpdateDependentStateJobDebug
[BurstCompile]
#endif
public struct UpdateDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceElevation;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> CloudCoverage;
	public NativeArray<float> SurfaceAirTemperature;
	public NativeArray<float> IceEnergy;
	[ReadOnly] public NativeArray<float> WaterDepth;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		float iceMass = IceMass[i];
		SurfaceElevation[i] = Terrain[i].Elevation + WaterDepth[i] + iceMass / WorldData.MassIce;
		VegetationCoverage[i] = math.saturate(Terrain[i].Vegetation * worldData.inverseFullCoverageVegetation);

		IceCoverage[i] = math.saturate(iceMass * worldData.inverseFullCoverageIce);
		IceEnergy[i] = WorldData.SpecificHeatIce * iceMass * IceTemperature[i];

		float cloudMass = CloudMass[i];
		CloudCoverage[i] = math.saturate(cloudMass * worldData.inverseFullCoverageCloud);
	}

}

#if !UpdateSurfaceDependentStateJobDebug
[BurstCompile]
#endif
public struct UpdateSurfaceDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceAirTemperatureAbsolute;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	public void Execute(int i)
	{
		SurfaceAirTemperatureAbsolute[i] = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], SurfaceElevation[i]);
	}

}


#if !UpdateAirPressureJobDebug
[BurstCompile]
#endif
public struct UpdateAirPressureJob : IJobParallelFor {
	public NativeArray<float> Pressure;
	public NativeArray<float> AirMass;

	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> AirTemperature;
	[ReadOnly] public NativeArray<float2> PressureGradient;
	[ReadOnly] public float Gravity;
	[ReadOnly] public int LayerIndex;

	public void Execute(int i)
	{
	}
}

#if !UpdateDependentAirLayerJobDebug
[BurstCompile]
#endif
public struct UpdateDependentAirLayerJob : IJobParallelFor {
	public NativeArray<float> AbsoluteHumidity;
	public NativeArray<float> RelativeHumidity;
	public NativeArray<float> Pressure;
	public NativeArray<float> PotentialEnergy;
	public NativeArray<float> AirMass;
	public NativeArray<float> AirMassTotal;
	public NativeArray<float> AirMassCloud;
	public NativeArray<float> AirVaporCloud;
	public NativeArray<float> AirPressureCloud;
	public NativeArray<float> AirHumidityRelativeCloud;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> DewPoint;
	public NativeArray<int> AirLayerCloud;

	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public float Gravity;
	[ReadOnly] public int LayerIndex;

	public void Execute(int i)
	{
		float layerMiddle = LayerElevation[i] + LayerHeight[i] / 2;
		float vaporMass = VaporMass[i];
		float airTemperaturePotential = AirTemperaturePotential[i];
		float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(airTemperaturePotential, layerMiddle);
		float airMass = Atmosphere.GetAirMass(LayerElevation[i], LayerHeight[i], airTemperaturePotential, Gravity);

		AirMass[i] = airMass;
		Pressure[i] = (AirMassTotal[i] + airMass / 2) * Gravity;
		AbsoluteHumidity[i] = vaporMass / (vaporMass + airMass);
		RelativeHumidity[i] = Atmosphere.GetRelativeHumidity(airMass, vaporMass, airTemperatureAbsolute, Pressure[i]);
		PotentialEnergy[i] = (airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * airTemperatureAbsolute;

		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];
		float dewPoint = Atmosphere.GetDewPoint(RelativeHumidity[i], airTemperatureAbsolute);
		float cloudElevation = Atmosphere.GetElevationAtDewPoint(dewPoint, airTemperaturePotential);
		float cloudBaseElevation = math.max(cloudElevation, SurfaceElevation[i]);
		if (cloudBaseElevation >= layerElevation && cloudBaseElevation < layerElevation + layerHeight)
		{
			float airMassCloud = AirMass[i];
			float vaporMassCloud = VaporMass[i];
			float airTemperatureCloud = airTemperatureAbsolute; // TODO: this should be modified by cloudelevation right?
			AirLayerCloud[i] = LayerIndex;
			AirMassCloud[i] = airMassCloud;
			AirVaporCloud[i] = VaporMass[i];
			AirPressureCloud[i] = (AirMassTotal[i] + airMass * (1.0f - (cloudBaseElevation - layerElevation) / layerHeight)) * Gravity;
			AirHumidityRelativeCloud[i] = Atmosphere.GetRelativeHumidity(airMassCloud, vaporMassCloud, airTemperatureCloud, AirPressureCloud[i]);
		}
		DewPoint[i] = dewPoint;
		CloudElevation[i] = cloudElevation;
		AirMassTotal[i] += airMass;

	}
}

#if !UpdateDependentWaterLayerJobDebug
[BurstCompile]
#endif
public struct UpdateDependentWaterLayerJob : IJobParallelFor {
	public NativeArray<float> WaterCoverage;
	public NativeArray<float> PotentialEnergy;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Temperature;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			WaterCoverage[i] = math.min(1, LayerHeight[i] / math.max(1, Terrain[i].Roughness));
			PotentialEnergy[i] = (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt) * Temperature[i];
		}
		else
		{
			WaterCoverage[i] = 0;
			PotentialEnergy[i] = 0;
		}
	}
}


#if !UpdateWaterDepthJobDebug
[BurstCompile]
#endif
public struct UpdateWaterDepthJob : IJobParallelFor {
	public NativeArray<float> Pressure;
	public NativeArray<float> LayerDepth;
	public NativeArray<float> LayerHeight;
	public NativeArray<float> Density;
	public NativeArray<float> WaterDepthTotal;
	public NativeArray<float> WaterMassTotal;
	public NativeArray<float> WaterDepthSum;
	public NativeArray<float> WaterMassSum;
	[ReadOnly] public NativeArray<float> UpLayerDepth;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public float WaterDensityPerDegree;
	[ReadOnly] public float WaterDensityPerSalinity;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			float waterMass = WaterMass[i];
			float saltMass = SaltMass[i];
			Density[i] = Atmosphere.GetWaterDensity(waterMass, saltMass, Temperature[i], WaterDensityPerSalinity, WaterDensityPerDegree);
			float waterDepth = (waterMass + saltMass) / Density[i];
			LayerDepth[i] = UpLayerDepth[i] + UpLayerHeight[i] + waterDepth;
			LayerHeight[i] = waterDepth;
			WaterDepthSum[i] += waterDepth;
			float layerMass = waterMass + saltMass;
			Pressure[i] = WaterMassTotal[i] + layerMass / 2 * Gravity;
			WaterMassSum[i] += layerMass;
		}
		else
		{
			LayerDepth[i] = 0;
			LayerHeight[i] = 0;
			Pressure[i] = 0;
			Density[i] = 0;
		}
		WaterDepthTotal[i] = WaterDepthSum[i];
		WaterMassTotal[i] = WaterMassSum[i];
	}
}

#endregion

