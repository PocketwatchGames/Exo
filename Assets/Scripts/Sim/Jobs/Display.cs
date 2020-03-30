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
	public NativeArray<float> EnthalpyTerrain;
	public NativeArray<float> EnthalpyFlora;
	public NativeArray<float> EnthalpyCloud;
	public NativeArray<float> EnthalpyIce;
	public NativeArray<float> EnthalpyGroundWater;
	[ReadOnly] public NativeArray<float> SolarRadiationInTerrain;
	[ReadOnly] public NativeArray<float> SolarRadiationInIce;
	[ReadOnly] public NativeArray<float> SolarRadiationInWaterSurface;
	[ReadOnly] public NativeArray<float> Precipitation;
	[ReadOnly] public NativeArray<float> EvaporationWater;
	[ReadOnly] public NativeArray<float> EvaporationFlora;
	[ReadOnly] public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> Flora;
	[ReadOnly] public NativeArray<float> FloraWater;
	[ReadOnly] public NativeArray<float> FloraTemperature;
	[ReadOnly] public NativeArray<float> GroundWaterMass;
	[ReadOnly] public NativeArray<float> GroundWaterTemperature;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		SolarRadiationAbsorbedSurface[i] = SolarRadiationInTerrain[i] + SolarRadiationInIce[i] + SolarRadiationInWaterSurface[i];
		DisplayPrecipitation[i] = Precipitation[i];
		DisplayEvaporation[i] = EvaporationWater[i] + EvaporationFlora[i];
		EnthalpyTerrain[i] = TerrainTemperature[i] * Atmosphere.GetSpecificHeatTerrain(HeatingDepth, SoilFertility[i]);
		EnthalpyFlora[i] = FloraTemperature[i] * (Flora[i] * WorldData.SpecificHeatFlora + FloraWater[i] * (WorldData.SpecificHeatWater + WorldData.LatentHeatWaterLiquid));
		EnthalpyCloud[i] = CloudMass[i] * WorldData.LatentHeatWaterLiquid;
		EnthalpyIce[i] = IceMass[i] * IceTemperature[i] * WorldData.SpecificHeatIce;
		EnthalpyGroundWater[i] = GroundWaterMass[i] * (WorldData.LatentHeatWaterLiquid + WorldData.SpecificHeatWater * GroundWaterTemperature[i]);
	}
}

[BurstCompile]
public struct InitDisplayAirLayerJob : IJobParallelFor {
	public NativeArray<float> DisplayPressure;
	public NativeArray<float3> DisplayPressureGradientForce;
	public NativeArray<float> DisplayCondensationGround;
	public NativeArray<float> DisplayCondensationCloud;
	public NativeArray<float> Enthalpy;
	public NativeArray<float> WindVertical;
	public NativeArray<float> DustCoverage;
	public NativeArray<float> CarbonDioxidePercent;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LayerMiddle;
	[ReadOnly] public NativeArray<float> CondensationCloud;
	[ReadOnly] public NativeArray<float> CondensationGround;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> DustMass;
	[ReadOnly] public NativeArray<float> CarbonDioxide;
	[ReadOnly] public NativeArray<BarycentricValueVertical> AdvectionDestination;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		DisplayPressure[i] = Atmosphere.GetPressureAtElevation(0, Gravity, AirPressure[i], AirTemperaturePotential[i], LayerMiddle[i]);
		DisplayPressureGradientForce[i] = PressureGradientForce[i];
		DisplayCondensationGround[i] += CondensationCloud[i];
		DisplayCondensationCloud[i] += CondensationGround[i];
		DustCoverage[i] += DustMass[i];
		CarbonDioxidePercent[i] += CarbonDioxide[i] / AirMass[i];
		if (AirMass[i] > 0)
		{
			Enthalpy[i] = AirTemperaturePotential[i] * (WorldData.SpecificHeatAtmosphere * AirMass[i] + WorldData.SpecificHeatWaterVapor * VaporMass[i]) + VaporMass[i] * (WorldData.LatentHeatWaterLiquid + WorldData.LatentHeatWaterVapor);
		}
		WindVertical[i] = AdvectionDestination[i].moveVertical;
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
