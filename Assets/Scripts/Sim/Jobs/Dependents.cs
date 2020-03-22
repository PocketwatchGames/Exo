
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
		var waterMassTotal = new NativeArray<float>(state.IceMass.Length, Allocator.Persistent);
		arraysToDispose.Add(waterMassSum);
		arraysToDispose.Add(airMassTotal);
		arraysToDispose.Add(waterMassTotal);
		for (int j = worldData.WaterLayers - 2; j >= 1; j--)
		{
			dependencies = jobHelper.Schedule(new UpdateWaterDepthJob()
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
			dependencies = jobHelper.Schedule(new UpdateDependentWaterLayerJob()
			{
				WaterCoverage = dependent.WaterCoverage[j],
				PotentialEnergy = dependent.WaterPotentialEnergy[j],
				Pressure = dependent.WaterPressure[j],
				WaterMassTotal = waterMassTotal,

				LayerHeight = dependent.WaterLayerHeight[j],
				LayerHeightDown = dependent.WaterLayerHeight[j-1],
				SaltMass = state.SaltMass[j],
				WaterMass = state.WaterMass[j],
				Temperature = state.WaterTemperature[j],
				Roughness = state.Roughness,
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
		}

		dependencies = jobHelper.Schedule(new UpdateDependentStateJob()
		{
			IceEnergy = dependent.IceEnergy,
			FloraEnergy = dependent.FloraEnergy,
			LavaEnergy = dependent.LavaEnergy,
			SurfaceElevation = dependent.LayerElevation[1],

			WaterDepth = dependent.WaterLayerDepth[1],
			Elevation = state.Elevation,
			FloraMass = state.FloraMass,
			FloraWater = state.FloraWater,
			FloraTemperature = state.FloraTemperature,
			LavaMass = state.LavaMass,
			LavaTemperature = state.LavaTemperature,
			IceMass = state.IceMass,
			IceTemperature = state.IceTemperature,
		}, dependencies);
		dependencies.Complete();

		var standardLayerElevation = new NativeArray<float>(dependent.LayerElevation[1], Allocator.TempJob);
		arraysToDispose.Add(standardLayerElevation);

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
			dependencies = jobHelper.Schedule(new UpdateAirLayerHeightsJob()
			{
				StandardLayerElevation = standardLayerElevation,
				LayerHeight = dependent.LayerHeight[j],
				UpLayerElevation = dependent.LayerElevation[j + 1],
				AirMass = dependent.AirMass[j],
				LayerMiddle = dependent.LayerMiddle[j],

				LayerElevation = dependent.LayerElevation[j],				
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				TropopauseElevation = worldData.TropopauseElevation,
				MinimumHeight = minumumHeight,
				ColumnPercent = columnPercent,
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
		}

		dependencies = jobHelper.Schedule(new UpdateStratosphereJob()
		{
			StratosphereMass = airMassTotal,

			TropopauseElevation = dependent.LayerElevation[worldData.AirLayers - 2],
			TropopauseHeight = dependent.LayerHeight[worldData.AirLayers - 2],
			Gravity = state.PlanetState.Gravity
		}, dependencies);


		for (int j = worldData.AirLayers - 2; j > 0; j--)
		{
			var pressureHandle = jobHelper.Schedule(new UpdateAirPressureJob()
			{
				Pressure = dependent.AirPressure[j],
				PressureInverse = dependent.AirPressureInverse[j],
				AirMassTotal = airMassTotal,

				AirMass = dependent.AirMass[j],
				CloudMass = state.CloudMass,
				VaporMass = state.AirVapor[j],
				LayerElevation = dependent.LayerElevation[j],
				LayerHeight = dependent.LayerHeight[j],
				SurfaceElevation = dependent.LayerElevation[1],
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
			dependencies = JobHandle.CombineDependencies(dependencies, jobHelper.Schedule(new UpdateDependentAirLayerJob()
			{
				RelativeHumidity = dependent.AirHumidityRelative[j],
				AbsoluteHumidity = dependent.AirHumidityAbsolute[j],
				AirMass = dependent.AirMass[j],
				PotentialEnergy = dependent.AirPotentialEnergy[j],

				Pressure = dependent.AirPressure[j],
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				VaporMass = state.AirVapor[j],
				LayerMiddle = dependent.LayerMiddle[j],
				Gravity = state.PlanetState.Gravity,
			}, pressureHandle));
		}

		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			dependencies = jobHelper.Schedule(new UpdateDependentCloudJob()
			{
				CloudElevation = dependent.CloudElevation,
				AirDensityCloud = dependent.AirDensityCloud,
				DewPoint = dependent.DewPoint,

				AirMass = dependent.AirMass[j],
				VaporMass = state.AirVapor[j],
				SurfaceElevation = dependent.LayerElevation[1],
				RelativeHumidity = dependent.AirHumidityRelative[j],				
				Pressure = dependent.AirPressure[j],
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				LayerElevation = dependent.LayerElevation[j],
				LayerMiddle = dependent.LayerMiddle[j],
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
		}

		dependencies = jobHelper.Schedule(new UpdateSurfaceDependentStateJob()
		{
			SurfaceAirTemperatureAbsolute = dependent.SurfaceAirTemperatureAbsolute,
			SurfaceAreaAirFlora = dependent.SurfaceAreaAirFlora,
			SurfaceAreaAirLava = dependent.SurfaceAreaAirLava,
			SurfaceAreaAirIce = dependent.SurfaceAreaAirIce,
			SurfaceAreaAirTerrain = dependent.SurfaceAreaAirTerrain,
			SurfaceAreaAirWater = dependent.SurfaceAreaAirWater,
			SurfaceAreaIceLava = dependent.SurfaceAreaIceLava,
			SurfaceAreaIceFlora = dependent.SurfaceAreaIceFlora,
			SurfaceAreaIceTerrain = dependent.SurfaceAreaIceTerrain,
			SurfaceAreaIceWater = dependent.SurfaceAreaIceWater,
			SurfaceAreaWaterLava = dependent.SurfaceAreaWaterLava,
			SurfaceAreaWaterFlora = dependent.SurfaceAreaWaterFlora,
			SurfaceAreaWaterTerrain = dependent.SurfaceAreaWaterTerrain,
			SurfaceAreaFloraLava = dependent.SurfaceAreaFloraLava,
			SurfaceAreaFloraTerrain = dependent.SurfaceAreaFloraTerrain,
			SurfaceAreaLavaTerrain = dependent.SurfaceAreaLavaTerrain,
			LavaCoverage = dependent.LavaCoverage,
			FloraCoverage = dependent.FloraCoverage,
			IceCoverage = dependent.IceCoverage,

			WaterCoverage = dependent.WaterCoverage[worldData.WaterLayers - 2],
			FloraMass = state.FloraMass,
			IceMass = state.IceMass,
			LavaMass = state.LavaMass,
			AirTemperaturePotential = state.AirTemperaturePotential[1],
			SurfaceLayerElevation = dependent.LayerElevation[1],
			inverseFullCoverageFloraMass = 1.0f / worldData.FullCoverageFlora,
			inverseFullCoverageIceMass = 1.0f / (worldData.FullCoverageIce * WorldData.MassIce),
			FloraAirSurfaceArea = worldData.FloraAirSurfaceArea,
			Roughness = state.Roughness
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
	public NativeArray<float> IceEnergy;
	public NativeArray<float> FloraEnergy;
	public NativeArray<float> LavaEnergy;
	[ReadOnly] public NativeArray<float> WaterDepth;
	[ReadOnly] public NativeArray<float> FloraMass;
	[ReadOnly] public NativeArray<float> FloraWater;
	[ReadOnly] public NativeArray<float> FloraTemperature;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> LavaMass;
	[ReadOnly] public NativeArray<float> LavaTemperature;
	[ReadOnly] public NativeArray<float> Elevation;
	public void Execute(int i)
	{
		float iceMass = IceMass[i];
		float lavaMass = LavaMass[i];
		SurfaceElevation[i] = Elevation[i] + WaterDepth[i] + iceMass / WorldData.MassIce + lavaMass / WorldData.MassLava;

		IceEnergy[i] = WorldData.SpecificHeatIce * iceMass * IceTemperature[i];
		LavaEnergy[i] = WorldData.SpecificHeatLava * lavaMass * LavaTemperature[i];
		FloraEnergy[i] = (WorldData.SpecificHeatFlora * FloraMass[i] + WorldData.SpecificHeatWater * FloraWater[i]) * FloraTemperature[i];
	}

}


#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateAirLayerHeightsJob : IJobParallelFor {
	public NativeArray<float> AirMass;
	public NativeArray<float> UpLayerElevation;
	public NativeArray<float> LayerHeight;
	public NativeArray<float> LayerMiddle;
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
		LayerMiddle[i] = layerElevation + layerHeight / 2;
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
public struct UpdateAirPressureJob : IJobParallelFor {
	public NativeArray<float> Pressure;
	public NativeArray<float> PressureInverse;
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
		var pressure = (AirMassTotal[i] + airMass / 2) * Gravity;
		Pressure[i] = pressure;
		PressureInverse[i] = 1.0f / pressure;
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

	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerMiddle;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public float Gravity;

	public void Execute(int i)
	{
		float airTemperaturePotential = AirTemperaturePotential[i];
		float vaporMass = VaporMass[i];
		float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(airTemperaturePotential, LayerMiddle[i]);

		float airMass = AirMass[i];
		AbsoluteHumidity[i] = vaporMass / (vaporMass + airMass);
		RelativeHumidity[i] = Atmosphere.GetRelativeHumidity(airMass, vaporMass, airTemperatureAbsolute, Pressure[i]);
		PotentialEnergy[i] = (airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * airTemperatureAbsolute;

	}
}

#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateDependentCloudJob : IJobParallelFor {
	public NativeArray<float> AirDensityCloud;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> DewPoint;

	[ReadOnly] public NativeArray<float> RelativeHumidity;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> LayerMiddle;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public float Gravity;

	public void Execute(int i)
	{
		float airTemperaturePotential = AirTemperaturePotential[i];
		float layerElevation = LayerElevation[i];

		// TODO: this is overwriting such that it only pays attention to the bottom layer of the atmosphere
		float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(airTemperaturePotential, LayerMiddle[i]);
		float dewPoint = Atmosphere.GetDewPoint(RelativeHumidity[i], airTemperatureAbsolute);
		DewPoint[i] = dewPoint;

		float cloudElevation = Atmosphere.GetElevationAtDewPoint(dewPoint, airTemperaturePotential);
		float cloudBaseElevation = math.max(cloudElevation, SurfaceElevation[i]);
		if (cloudBaseElevation >= layerElevation)
		{
			float airMassCloud = AirMass[i];
			float vaporMassCloud = VaporMass[i];
			float airPressureCloud = Atmosphere.GetPressureAtElevation(cloudBaseElevation, Gravity, Pressure[i], airTemperaturePotential, LayerMiddle[i]);
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
	[ReadOnly] public NativeArray<float> LayerHeightDown;
	[ReadOnly] public NativeArray<float> Roughness;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			WaterCoverage[i] = math.min(1, (LayerHeight[i] + LayerHeightDown[i]) / math.max(1, Roughness[i]));
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
			Density[i] = Atmosphere.GetWaterDensity(Atmosphere.GetWaterSalinity(waterMass, saltMass), Temperature[i]);
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


#if !UpdateDependentJobDebug
[BurstCompile]
#endif
public struct UpdateSurfaceDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceAirTemperatureAbsolute;
	public NativeArray<float> SurfaceAreaAirIce;
	public NativeArray<float> SurfaceAreaAirWater;
	public NativeArray<float> SurfaceAreaAirFlora;
	public NativeArray<float> SurfaceAreaAirLava;
	public NativeArray<float> SurfaceAreaAirTerrain;
	public NativeArray<float> SurfaceAreaIceWater;
	public NativeArray<float> SurfaceAreaIceFlora;
	public NativeArray<float> SurfaceAreaIceLava;
	public NativeArray<float> SurfaceAreaIceTerrain;
	public NativeArray<float> SurfaceAreaWaterFlora;
	public NativeArray<float> SurfaceAreaWaterLava;
	public NativeArray<float> SurfaceAreaWaterTerrain;
	public NativeArray<float> SurfaceAreaFloraLava;
	public NativeArray<float> SurfaceAreaFloraTerrain;
	public NativeArray<float> SurfaceAreaLavaTerrain;
	public NativeArray<float> FloraCoverage;
	public NativeArray<float> LavaCoverage;
	public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> FloraMass;
	[ReadOnly] public NativeArray<float> LavaMass;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> SurfaceLayerElevation;
	[ReadOnly] public NativeArray<float> Roughness;
	[ReadOnly] public float inverseFullCoverageIceMass;
	[ReadOnly] public float inverseFullCoverageFloraMass;
	[ReadOnly] public float FloraAirSurfaceArea;
	public void Execute(int i)
	{
		SurfaceAirTemperatureAbsolute[i] = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], SurfaceLayerElevation[i]);

		float floraCoverage = math.saturate(FloraMass[i] * inverseFullCoverageFloraMass);
		FloraCoverage[i] = floraCoverage;

		float lavaCoverage = math.saturate(LavaMass[i] / (Roughness[i] * WorldData.MassLava));
		LavaCoverage[i] = lavaCoverage;

		float iceCoverage = math.saturate(IceMass[i] * inverseFullCoverageIceMass);
		IceCoverage[i] = iceCoverage;

		float waterCoverage = WaterCoverage[i];

		SurfaceAreaAirIce[i] = iceCoverage;
		SurfaceAreaAirWater[i] = math.max(0, waterCoverage - iceCoverage);
		SurfaceAreaAirFlora[i] = math.min(floraCoverage, 1.0f - math.max(waterCoverage, iceCoverage)) * FloraAirSurfaceArea;
		SurfaceAreaAirLava[i] = math.min(lavaCoverage, 1.0f - math.max(waterCoverage, iceCoverage));
		SurfaceAreaAirTerrain[i] = math.max(0, 1.0f - (math.max(lavaCoverage, floraCoverage) + math.max(waterCoverage, iceCoverage)));
		SurfaceAreaIceWater[i] = math.min(waterCoverage, iceCoverage);
		SurfaceAreaIceFlora[i] = math.max(0, (floraCoverage + math.max(0, iceCoverage - waterCoverage)) - 1.0f) * FloraAirSurfaceArea;
		SurfaceAreaIceLava[i] = math.min(math.max(0, iceCoverage - waterCoverage), lavaCoverage);
		SurfaceAreaIceTerrain[i] = math.max(0, iceCoverage - waterCoverage) - math.max(0, math.max(lavaCoverage, floraCoverage) + iceCoverage - 1.0f);
		SurfaceAreaWaterFlora[i] = math.max(0, (floraCoverage + waterCoverage) - 1.0f) * FloraAirSurfaceArea;
		SurfaceAreaWaterLava[i] = math.min(lavaCoverage, waterCoverage);
		SurfaceAreaWaterTerrain[i] = waterCoverage - math.max(0, (floraCoverage + waterCoverage) - 1.0f);
		SurfaceAreaFloraLava[i] = math.min(lavaCoverage, floraCoverage);
		SurfaceAreaFloraTerrain[i] = floraCoverage;
		SurfaceAreaLavaTerrain[i] = lavaCoverage;
	}

}


#endregion

