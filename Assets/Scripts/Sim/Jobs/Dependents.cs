//#define UpdateDependentWaterLayerJobDebug


using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;



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
	public NativeArray<float> WaterDepth;
	[ReadOnly] public NativeArray<float> WaterDepthTotal;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		float iceMass = IceMass[i];
		float waterDepth = WaterDepthTotal[i];
		WaterDepth[i] = waterDepth;
		SurfaceElevation[i] = Terrain[i].Elevation + waterDepth + iceMass / WorldData.MassIce;
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
	[ReadOnly] public float DewPointZero;
	[ReadOnly] public float WaterVaporMassToAirMassAtDewPoint;
	[ReadOnly] public float InverseDewPointTemperatureRange;
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
	[ReadOnly] public float DewPointZero;
	[ReadOnly] public float WaterVaporMassToAirMassAtDewPoint;
	[ReadOnly] public float InverseDewPointTemperatureRange;
	[ReadOnly] public int LayerIndex;

	public void Execute(int i)
	{
		float layerMiddle = LayerElevation[i] + LayerHeight[i] / 2;
		float vaporMass = VaporMass[i];
		float airTemperaturePotential = AirTemperaturePotential[i];
		float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(airTemperaturePotential, layerMiddle);
		float airMass = Atmosphere.GetAirMass(LayerElevation[i], LayerHeight[i], airTemperaturePotential, Gravity);

		AbsoluteHumidity[i] = vaporMass / (vaporMass + airMass);
		RelativeHumidity[i] = Atmosphere.GetRelativeHumidity(airMass, vaporMass, airTemperatureAbsolute, DewPointZero, WaterVaporMassToAirMassAtDewPoint, InverseDewPointTemperatureRange);
		PotentialEnergy[i] = (airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * airTemperatureAbsolute;
		AirMass[i] = airMass;
		Pressure[i] = (AirMassTotal[i] + airMass / 2) * Gravity;

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
			AirHumidityRelativeCloud[i] = Atmosphere.GetRelativeHumidity(airMassCloud, vaporMassCloud, airTemperatureCloud, DewPointZero, WaterVaporMassToAirMassAtDewPoint, InverseDewPointTemperatureRange);
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
	public NativeArray<float> Salinity;
	public NativeArray<float> Density;
	public NativeArray<float> Pressure;
	public NativeArray<float> LayerDepth;
	public NativeArray<float> LayerHeight;
	public NativeArray<float> WaterCoverage;
	public NativeArray<float> PotentialEnergy;
	public NativeArray<float> WaterDepthTotal;
	public NativeArray<float> WaterMassTotal;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
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
			float waterDepth = WaterMass[i] / WorldData.MassWater + SaltMass[i] / WorldData.MassSalt;
			LayerDepth[i] = UpLayerDepth[i] + UpLayerHeight[i] + waterDepth;
			LayerHeight[i] = waterDepth;
			WaterDepthTotal[i] += waterDepth;
			float layerMass = WaterMass[i] + SaltMass[i];
			Pressure[i] = WaterMassTotal[i] + layerMass / 2 * Gravity;
			WaterMassTotal[i] += layerMass;
			WaterCoverage[i] = math.min(1, waterDepth / math.max(1, Terrain[i].Roughness));
			Salinity[i] = SaltMass[i] / (WaterMass[i] + SaltMass[i]);
			Density[i] = Atmosphere.GetWaterDensity(WaterMass[i], SaltMass[i], Temperature[i], WaterDensityPerSalinity, WaterDensityPerDegree);
			PotentialEnergy[i] = (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt) * Temperature[i];
		}
		else
		{
			LayerDepth[i] = 0;
			LayerHeight[i] = 0;
			Pressure[i] = 0;
			WaterCoverage[i] = 0;
			Salinity[i] = 0;
			PotentialEnergy[i] = 0;
			Density[i] = 0;
		}
	}
}

#endregion

