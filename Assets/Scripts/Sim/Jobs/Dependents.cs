﻿using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

public struct WaterSaltMass {
	public float WaterMass;
	public float SaltMass;
}



#region Update dependents

#if !UpdateWaterSaltMassJobDebug
[BurstCompile]
#endif
public struct UpdateWaterSaltMassJob : IJobParallelFor {
	public NativeArray<WaterSaltMass> WaterSaltMass;
	[ReadOnly] public NativeArray<float> WaterLayerMass;
	[ReadOnly] public NativeArray<float> SaltLayerMass;
	public void Execute(int i)
	{
		WaterSaltMass[i] = new WaterSaltMass()
		{
			WaterMass = WaterSaltMass[i].WaterMass + WaterLayerMass[i],
			SaltMass = WaterSaltMass[i].SaltMass + SaltLayerMass[i],
		};
	}

}

#if !UpdateDependentStateJobDebug
[BurstCompile]
#endif
public struct UpdateDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceElevation;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> CloudCoverage;
	public NativeArray<float> WaterDepth;
	public NativeArray<float> SurfaceAirTemperature;
	public NativeArray<float> IceEnergy;
	[ReadOnly] public NativeArray<WaterSaltMass> WaterSaltMass;
	[ReadOnly] public NativeArray<float> LowerAirTemperature;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public WorldData worldData;
	[ReadOnly] public NativeArray<float> lowerAirHeight;
	public void Execute(int i)
	{
		float waterDepth = WaterSaltMass[i].WaterMass / WorldData.MassWater + WaterSaltMass[i].SaltMass / WorldData.MassSalt;
		float iceMass = IceMass[i];
		WaterDepth[i] = waterDepth;
		SurfaceElevation[i] = Terrain[i].Elevation + waterDepth + iceMass / WorldData.MassIce;
		VegetationCoverage[i] = math.saturate(Terrain[i].Vegetation * worldData.inverseFullCanopyCoverage);

		IceCoverage[i] = math.saturate(iceMass * worldData.inverseFullIceCoverage);
		IceEnergy[i] = WorldData.SpecificHeatIce * iceMass * IceTemperature[i];

		float cloudMass = CloudMass[i];
		CloudCoverage[i] = math.saturate(cloudMass * worldData.inverseCloudMassFullAbsorption);

		SurfaceAirTemperature[i] = LowerAirTemperature[i] - WorldData.TemperatureLapseRate * lowerAirHeight[i] / 2;
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

	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> AirTemperature;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float DewPointZero;
	[ReadOnly] public float WaterVaporMassToAirMassAtDewPoint;
	[ReadOnly] public float InverseDewPointTemperatureRange;
	[ReadOnly] public int LayerIndex;

	public void Execute(int i)
	{
		float layerMiddle = LayerElevation[i] + LayerHeight[i] / 2;
		float vaporMass = VaporMass[i];
		float airTemperature = AirTemperature[i];
		float airMass = Atmosphere.GetAirMass(LayerElevation[i], LayerHeight[i], airTemperature, Gravity);

		AbsoluteHumidity[i] = vaporMass / (vaporMass + airMass);
		RelativeHumidity[i] = Atmosphere.GetRelativeHumidity(airMass, vaporMass, airTemperature, DewPointZero, WaterVaporMassToAirMassAtDewPoint, InverseDewPointTemperatureRange);
		PotentialEnergy[i] = (airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * airTemperature;
		AirMass[i] = airMass;
		Pressure[i] = (AirMassTotal[i] + airMass / 2) * Gravity;

		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];
		float dewPoint = Atmosphere.GetDewPoint(RelativeHumidity[i], AirTemperature[i]);
		float cloudElevation = Atmosphere.GetElevationAtDewPoint(dewPoint, AirTemperature[i], layerElevation + layerHeight / 2);
		if (cloudElevation >= layerElevation && cloudElevation < layerElevation + layerHeight)
		{
			float airMassCloud = AirMass[i];
			float vaporMassCloud = VaporMass[i];
			float airTemperatureCloud = airTemperature;
			AirLayerCloud[i] = LayerIndex;
			AirMassCloud[i] = airMassCloud;
			AirVaporCloud[i] = VaporMass[i];
			AirPressureCloud[i] = (AirMassTotal[i] + airMass * (1.0f - (cloudElevation - layerElevation) / layerHeight)) * Gravity;
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
	public NativeArray<float> WaterCoverage;
	public NativeArray<float> PotentialEnergy;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			float waterDepth = WaterMass[i] / WorldData.MassWater + SaltMass[i] / WorldData.MassSalt;
			WaterCoverage[i] = math.min(1, waterDepth / math.max(1, Terrain[i].Roughness));
			Salinity[i] = SaltMass[i] / (WaterMass[i] + SaltMass[i]);
			PotentialEnergy[i] = (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt) * Temperature[i];
		}
		else
		{
			WaterCoverage[i] = 0;
			Salinity[i] = 0;
			PotentialEnergy[i] = 0;
		}
	}
}

#endregion

