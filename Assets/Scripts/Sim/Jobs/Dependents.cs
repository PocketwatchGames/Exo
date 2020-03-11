//#define UpdateDependentJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using System.Collections.Generic;

public static class SimJobs {

	public static JobHandle UpdateDependentVariables(JobHelper jobHelper, ref SimState state, ref DependentState dependent, ref WorldData worldData, JobHandle dependencies, List<NativeArray<float>> arraysToDispose)
	{
		var waterMassSum = new NativeArray<float>(state.IceMass, Allocator.TempJob);
		var airMassTotal = new NativeArray<float>(state.IceMass.Length, Allocator.TempJob);
		var standardLayerElevation = new NativeArray<float>(state.IceMass.Length, Allocator.TempJob);
		var waterMassTotal = new NativeArray<float>(state.IceMass.Length, Allocator.Persistent);
		arraysToDispose.Add(waterMassSum);
		arraysToDispose.Add(airMassTotal);
		arraysToDispose.Add(standardLayerElevation);
		arraysToDispose.Add(waterMassTotal);
		for (int j = worldData.WaterLayers - 2; j >= 1; j--)
		{
			dependencies = jobHelper.Run(new UpdateWaterDepthJob()
			{
				Density = dependent.WaterDensity[j],
				LayerDepth = dependent.WaterLayerDepth[j],
				LayerHeight = dependent.WaterLayerHeight[j],
				WaterMassSum = waterMassSum,

				Temperature = state.WaterTemperature[j],
				SaltMass = state.SaltMass[j],
				WaterMass = state.WaterMass[j],
				UpLayerDepth = dependent.WaterLayerDepth[j + 1],
				UpLayerHeight = dependent.WaterLayerHeight[j + 1],
				WaterDensityPerDegree = worldData.WaterDensityPerDegree,
				WaterDensityPerSalinity = worldData.WaterDensityPerSalinity,

			}, dependencies);
		}


		for (int j = worldData.WaterLayers - 2; j >= 1; j--)
		{
			dependencies = jobHelper.Run(new UpdateDependentWaterLayerJob()
			{
				WaterCoverage = dependent.WaterCoverage[j],
				PotentialEnergy = dependent.WaterPotentialEnergy[j],
				Pressure = dependent.WaterPressure[j],
				WaterMassTotal = waterMassTotal,

				LayerHeight = dependent.WaterLayerHeight[j],
				SaltMass = state.SaltMass[j],
				WaterMass = state.WaterMass[j],
				Temperature = state.WaterTemperature[j],
				Terrain = state.Terrain,
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
		}

		dependencies = jobHelper.Run(new UpdateDependentStateJob()
		{
			CloudCoverage = dependent.CloudCoverage,
			IceCoverage = dependent.IceCoverage,
			IceEnergy = dependent.IceEnergy,
			SurfaceElevation = dependent.LayerElevation[1],
			VegetationCoverage = dependent.VegetationCoverage,
			WaterDepth = dependent.WaterLayerDepth[1],

			CloudMass = state.CloudMass,
			IceMass = state.IceMass,
			IceTemperature = state.IceTemperature,
			Terrain = state.Terrain,
			worldData = worldData,
		}, dependencies);


		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			float minumumHeight;
			float columnPercent;
			if (j == 1)
			{
				minumumHeight = worldData.BoundaryZoneElevation;
				columnPercent = 0;
			}
			else
			{
				minumumHeight = 0;
				columnPercent = 1.0f / (worldData.AirLayers - 1 - j);
			}
			dependencies = jobHelper.Run(new UpdateAirLayerHeightsJob()
			{
				StandardLayerElevation = standardLayerElevation,
				LayerHeight = dependent.LayerHeight[j],
				UpLayerElevation = dependent.LayerElevation[j + 1],
				AirMass = dependent.AirMass[j],

				LayerElevation = dependent.LayerElevation[j],				
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				TropopauseElevation = worldData.TropopauseElevation,
				MinimumHeight = minumumHeight,
				ColumnPercent = columnPercent,
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
		}

		dependencies = jobHelper.Run(new UpdateStratosphereJob()
		{
			StratosphereMass = airMassTotal,

			TropopauseElevation = dependent.LayerElevation[worldData.AirLayers - 2],
			TropopauseHeight = dependent.LayerHeight[worldData.AirLayers - 2],
			Gravity = state.PlanetState.Gravity
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
				SurfaceElevation = dependent.LayerElevation[1],
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
				SurfaceElevation = dependent.LayerElevation[1],
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
public struct UpdateDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceElevation;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> CloudCoverage;
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


#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateAirLayerHeightsJob : IJobParallelFor {
	public NativeArray<float> AirMass;
	public NativeArray<float> UpLayerElevation;
	public NativeArray<float> LayerHeight;
	public NativeArray<float> StandardLayerElevation;

	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float TropopauseElevation;
	[ReadOnly] public float MinimumHeight;
	[ReadOnly] public float ColumnPercent;

	public void Execute(int i)
	{
		float standardLayerHeight = math.max(MinimumHeight, (TropopauseElevation - StandardLayerElevation[i]) * ColumnPercent);
		float standardLayerElevation = StandardLayerElevation[i];
		float airMass = (Atmosphere.GetStandardPressureAtElevation(standardLayerElevation, WorldData.StandardTemperature, Gravity) - Atmosphere.GetStandardPressureAtElevation(standardLayerElevation + standardLayerHeight, WorldData.StandardTemperature, Gravity)) / Gravity;
		StandardLayerElevation[i] = standardLayerElevation + standardLayerHeight;
		AirMass[i] = airMass;

		float airTemperaturePotential = AirTemperaturePotential[i];
		float layerElevation = LayerElevation[i];
		float layerHeight = standardLayerHeight * AirTemperaturePotential[i] / WorldData.StandardTemperature;
		LayerHeight[i] = layerHeight;
		UpLayerElevation[i] = layerElevation + layerHeight;
	}
}

#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateStratosphereJob : IJobParallelFor {
	public NativeArray<float> StratosphereMass;
	[ReadOnly] public NativeArray<float> TropopauseElevation;
	[ReadOnly] public NativeArray<float> TropopauseHeight;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
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
	public NativeArray<float> Pressure;
	public NativeArray<float> WaterMassTotal;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			WaterCoverage[i] = math.min(1, LayerHeight[i] / math.max(1, Terrain[i].Roughness));
			PotentialEnergy[i] = (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt) * Temperature[i];
			float layerMass = WaterMass[i] + SaltMass[i];
			Pressure[i] = (WaterMassTotal[i] + layerMass / 2) * Gravity;
			WaterMassTotal[i] += layerMass;
		}
		else
		{
			Pressure[i] = 0;
			WaterCoverage[i] = 0;
			PotentialEnergy[i] = 0;
		}
	}
}


#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateWaterDepthJob : IJobParallelFor {
	public NativeArray<float> LayerDepth;
	public NativeArray<float> LayerHeight;
	public NativeArray<float> Density;
	public NativeArray<float> WaterMassSum;
	[ReadOnly] public NativeArray<float> UpLayerDepth;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public float WaterDensityPerDegree;
	[ReadOnly] public float WaterDensityPerSalinity;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			float waterMass = WaterMass[i];
			float saltMass = SaltMass[i];
			Density[i] = Atmosphere.GetWaterDensity(waterMass, saltMass, Temperature[i], WaterDensityPerSalinity, WaterDensityPerDegree);
			float waterDepth = (waterMass + saltMass) / Density[i];
			LayerDepth[i] = UpLayerDepth[i] + waterDepth;
			LayerHeight[i] = waterDepth;
		}
		else
		{
			LayerDepth[i] = UpLayerDepth[i];
			LayerHeight[i] = 0;
			Density[i] = 0;
		}
	}
}

#endregion

