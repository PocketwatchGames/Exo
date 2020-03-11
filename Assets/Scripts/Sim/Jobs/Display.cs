using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;



#if !UpdateDisplayJobDebug
[BurstCompile]
#endif
public struct UpdateDisplayJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbedSurface;
	public NativeArray<float> DisplayPrecipitation;
	public NativeArray<float> DisplayEvaporation;
	public NativeArray<float> Enthalpy;
	[ReadOnly] public NativeArray<float> SolarRadiationInTerrain;
	[ReadOnly] public NativeArray<float> SolarRadiationInIce;
	[ReadOnly] public NativeArray<float> SolarRadiationInWaterSurface;
	[ReadOnly] public NativeArray<float> Precipitation;
	[ReadOnly] public NativeArray<float> Evaporation;
	[ReadOnly] public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		SolarRadiationAbsorbedSurface[i] = SolarRadiationInTerrain[i] + SolarRadiationInIce[i] + SolarRadiationInWaterSurface[i];
		DisplayPrecipitation[i] = Precipitation[i];
		DisplayEvaporation[i] = Evaporation[i];
		Enthalpy[i] = 
			TerrainTemperature[i] * Atmosphere.GetSpecificHeatTerrain(HeatingDepth, Terrain[i].SoilFertility, Terrain[i].Vegetation) 
			+ CloudMass[i] * WorldData.LatentHeatWaterLiquid
			+ IceMass[i] * IceTemperature[i] * WorldData.SpecificHeatIce;
	}
}

#if !InitDisplayJobDebug
[BurstCompile]
#endif
public struct InitDisplayAirLayerJob : IJobParallelFor {
	public NativeArray<float> DisplayPressure;
	public NativeArray<float3> DisplayPressureGradientForce;
	public NativeArray<float> DisplayCondensationGround;
	public NativeArray<float> DisplayCondensationCloud;
	public NativeArray<float> Enthalpy;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CondensationCloud;
	[ReadOnly] public NativeArray<float> CondensationGround;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		DisplayPressure[i] = Atmosphere.GetPressureAtElevation(0, Gravity, AirPressure[i], AirTemperaturePotential[i], LayerElevation[i] + LayerHeight[i] / 2);
		DisplayPressureGradientForce[i] = PressureGradientForce[i];
		DisplayCondensationGround[i] += CondensationCloud[i];
		DisplayCondensationCloud[i] += CondensationGround[i];
		if (AirMass[i] > 0)
		{
			Enthalpy[i] = AirTemperaturePotential[i] * (WorldData.SpecificHeatAtmosphere * AirMass[i] + WorldData.SpecificHeatWaterVapor * VaporMass[i]) + VaporMass[i] * (WorldData.LatentHeatWaterLiquid + WorldData.LatentHeatWaterVapor);
		}
	}
}
#if !InitDisplayJobDebug
[BurstCompile]
#endif
public struct InitDisplayWaterLayerJob : IJobParallelFor {
	public NativeArray<float> Enthalpy;
	public NativeArray<float> Salinity;
	[ReadOnly] public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	public void Execute(int i)
	{
		Salinity[i] = Atmosphere.GetWaterSalinity(WaterMass[i], SaltMass[i]);
		if (WaterMass[i] > 0)
		{
			Enthalpy[i] = WaterTemperature[i] * (WorldData.SpecificHeatWater * WaterMass[i] + WorldData.SpecificHeatSalt * SaltMass[i]) + WaterMass[i] * WorldData.LatentHeatWaterLiquid;
		}
	}
}
