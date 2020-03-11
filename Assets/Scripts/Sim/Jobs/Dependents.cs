#define UpdateDependentJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using System.Collections.Generic;

public static class SimJobs {

	public static JobHandle UpdateDependentVariables(JobHelper jobHelper, ref SimState state, ref DependentState dependent, ref WorldData worldData, JobHandle dependencies, List<NativeArray<float>> arraysToDispose)
	{
		var waterMassSum = new NativeArray<float>(state.IceMass, Allocator.TempJob);
		var waterDepthSum = new NativeArray<float>(state.IceMass.Length, Allocator.TempJob);
		var airMassTotal = new NativeArray<float>(state.IceMass.Length, Allocator.TempJob);
		arraysToDispose.Add(waterMassSum);
		arraysToDispose.Add(waterDepthSum);
		arraysToDispose.Add(airMassTotal);
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


		for (int j = worldData.WaterLayers - 2; j >= 1; j--)
		{
			dependencies = jobHelper.Run(new UpdateDependentWaterLayerJob()
			{
				WaterCoverage = dependent.WaterCoverage[j],
				PotentialEnergy = dependent.WaterPotentialEnergy[j],

				LayerHeight = dependent.WaterLayerHeight[j],
				SaltMass = state.SaltMass[j],
				WaterMass = state.WaterMass[j],
				Temperature = state.WaterTemperature[j],
				Terrain = state.Terrain,
			}, dependencies);
		}


		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			dependencies = jobHelper.Run(new UpdateAirLayerHeightsJob()
			{
				LayerHeight = dependent.LayerHeight[j],
				LayerElevation = dependent.LayerElevation[j],

				DownLayerElevation = dependent.LayerElevation[j - 1],
				DownLayerHeight = dependent.LayerHeight[j - 1],
				AirMass = dependent.AirMass[j],
				AirTemperaturePotential = state.AirTemperaturePotential[j],
			}, dependencies);
		}


		dependencies = jobHelper.Run(new UpdateDependentStateJob()
		{
			CloudCoverage = dependent.CloudCoverage,
			IceCoverage = dependent.IceCoverage,
			IceEnergy = dependent.IceEnergy,
			SurfaceElevation = dependent.SurfaceElevation,
			VegetationCoverage = dependent.VegetationCoverage,
			WaterDepth = dependent.WaterDepthTotal,
			StratosphereMass = airMassTotal,

			CloudMass = state.CloudMass,
			IceMass = state.IceMass,
			IceTemperature = state.IceTemperature,
			TropopauseElevation = dependent.LayerElevation[worldData.AirLayers - 2],
			TropopauseHeight = dependent.LayerHeight[worldData.AirLayers - 2],
			Terrain = state.Terrain,
			worldData = worldData,
		}, dependencies);

		for (int j = worldData.AirLayers - 2; j > 0; j--)
		{
			var pressureHandle = jobHelper.Run(new UpdateAirPressureJob()
			{
				Pressure = dependent.AirPressure[j],
				AirMassTotal = airMassTotal,

				AirMass = dependent.AirMass[j],
				CloudMass = state.CloudMass,
				VaporMass = state.AirVapor[j],
				LayerElevation = dependent.LayerElevation[j],
				LayerHeight = dependent.LayerHeight[j],
				SurfaceElevation = dependent.SurfaceElevation,
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
			dependencies = JobHandle.CombineDependencies(dependencies, jobHelper.Run(new UpdateDependentAirLayerJob()
			{
				RelativeHumidity = dependent.AirHumidityRelative[j],
				AbsoluteHumidity = dependent.AirHumidityAbsolute[j],
				AirMass = dependent.AirMass[j],
				PotentialEnergy = dependent.AirPotentialEnergy[j],
				DewPoint = dependent.DewPoint,

				Pressure = dependent.AirPressure[j],
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				VaporMass = state.AirVapor[j],
				LayerElevation = dependent.LayerElevation[j],
				LayerHeight = dependent.LayerHeight[j],
				Gravity = state.PlanetState.Gravity,
			}, pressureHandle));
		}

		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			dependencies = jobHelper.Run(new UpdateDependentCloudJob()
			{
				CloudElevation = dependent.CloudElevation,
				AirDensityCloud = dependent.AirDensityCloud,

				AirMass = dependent.AirMass[j],
				VaporMass = state.AirVapor[j],
				SurfaceElevation = dependent.SurfaceElevation,
				DewPoint = dependent.DewPoint,
				Pressure = dependent.AirPressure[j],
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				LayerElevation = dependent.LayerElevation[j],
				LayerHeight = dependent.LayerHeight[j],
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
		}

		dependencies = jobHelper.Run(new UpdateSurfaceDependentStateJob()
		{
			SurfaceAirTemperatureAbsolute = dependent.SurfaceAirTemperatureAbsolute,

			AirTemperaturePotential = state.AirTemperaturePotential[1],
			SurfaceLayerElevation = dependent.LayerElevation[1]
		}, dependencies);




		return dependencies;
	}

}

#region Update dependents

#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateAirLayerHeightsJob : IJobParallelFor {
	public NativeArray<float> LayerElevation;
	public NativeArray<float> LayerHeight;

	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public NativeArray<float> DownLayerHeight;

	public void Execute(int i)
	{
		float airTemperaturePotential = AirTemperaturePotential[i];
		float layerElevation = DownLayerElevation[i] + DownLayerHeight[i];
		LayerElevation[i] = layerElevation;

		float layerHeight = 2000;
		LayerHeight[i] = layerHeight;
	}
}

#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceElevation;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> CloudCoverage;
	public NativeArray<float> IceEnergy;
	public NativeArray<float> StratosphereMass;
	[ReadOnly] public NativeArray<float> WaterDepth;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> TropopauseElevation;
	[ReadOnly] public NativeArray<float> TropopauseHeight;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public WorldData worldData;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		float iceMass = IceMass[i];
		SurfaceElevation[i] = Terrain[i].Elevation + WaterDepth[i] + iceMass / WorldData.MassIce;
		VegetationCoverage[i] = math.saturate(Terrain[i].Vegetation * worldData.inverseFullCoverageVegetation);

		IceCoverage[i] = math.saturate(iceMass * worldData.inverseFullCoverageIce);
		IceEnergy[i] = WorldData.SpecificHeatIce * iceMass * IceTemperature[i];

		float cloudMass = CloudMass[i];
		CloudCoverage[i] = math.saturate(cloudMass * worldData.inverseFullCoverageCloud);

		// TODO: rebalance stratosphere mass
		float tropopauseElevation = TropopauseElevation[i] + TropopauseHeight[i];
		float stratosphereMass = Atmosphere.GetStandardPressureAtElevation(tropopauseElevation, WorldData.StandardTemperature, Gravity) / Gravity;
		StratosphereMass[i] = stratosphereMass;
	}

}

#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateSurfaceDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceAirTemperatureAbsolute;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> SurfaceLayerElevation;
	public void Execute(int i)
	{
		SurfaceAirTemperatureAbsolute[i] = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], SurfaceLayerElevation[i]);
	}

}


#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateAirPressureJob : IJobParallelFor {
	public NativeArray<float> Pressure;
	public NativeArray<float> AirMassTotal;

	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public float Gravity;

	public void Execute(int i)
	{
		float airMass = AirMass[i];
		Pressure[i] = (AirMassTotal[i] + airMass / 2) * Gravity;
		AirMassTotal[i] += airMass;

	}
}


#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateDependentAirLayerJob : IJobParallelFor {
	public NativeArray<float> AbsoluteHumidity;
	public NativeArray<float> RelativeHumidity;
	public NativeArray<float> PotentialEnergy;
	public NativeArray<float> DewPoint;

	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public float Gravity;

	public void Execute(int i)
	{
		float airTemperaturePotential = AirTemperaturePotential[i];
		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];

		float layerMiddle = layerElevation + LayerHeight[i] / 2;
		float vaporMass = VaporMass[i];
		float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(airTemperaturePotential, layerMiddle);

		float airMass = AirMass[i];
		AbsoluteHumidity[i] = vaporMass / (vaporMass + airMass);
		RelativeHumidity[i] = Atmosphere.GetRelativeHumidity(airMass, vaporMass, airTemperatureAbsolute, Pressure[i]);
		PotentialEnergy[i] = (airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * airTemperatureAbsolute;

		float dewPoint = Atmosphere.GetDewPoint(RelativeHumidity[i], airTemperatureAbsolute);
		DewPoint[i] = dewPoint;

	}
}

#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateDependentCloudJob : IJobParallelFor {
	public NativeArray<float> AirDensityCloud;
	public NativeArray<float> CloudElevation;

	[ReadOnly] public NativeArray<float> DewPoint;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public float Gravity;

	public void Execute(int i)
	{
		float airTemperaturePotential = AirTemperaturePotential[i];
		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];

		float cloudElevation = Atmosphere.GetElevationAtDewPoint(DewPoint[i], airTemperaturePotential);
		float cloudBaseElevation = math.max(cloudElevation, SurfaceElevation[i]);
		if (cloudBaseElevation >= layerElevation && cloudBaseElevation < layerElevation + layerHeight)
		{
			float airMassCloud = AirMass[i];
			float vaporMassCloud = VaporMass[i];
			float layerMiddle = layerElevation + LayerHeight[i] / 2;
			float airPressureCloud = Atmosphere.GetPressureAtElevation(cloudBaseElevation, Gravity, Pressure[i], airTemperaturePotential, layerMiddle);
			AirDensityCloud[i] = Atmosphere.GetAirDensity(airPressureCloud, DewPoint[i], airMassCloud, vaporMassCloud);   
		}
		CloudElevation[i] = cloudElevation;

	}
}

#if !UpdateDependentJobDebug
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


#if !UpdateDependentJobDebug
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

