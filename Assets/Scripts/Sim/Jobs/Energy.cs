//#define EnergyAirJobDebug
//#define EnergyWaterSurfaceJobDebug
//#define EnergyIceJobDebug
//#define EnergyTerrainJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !EnergyAirJobDebug
[BurstCompile]
#endif
public struct EnergyAirJob : IJobParallelFor {
	public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastTemperaturePotential;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		float energy =
			SolarRadiationIn[i]
			+ ThermalRadiationDelta[i]
			+ ConductionEnergyIce[i]
			+ ConductionEnergyTerrain[i]
			+ ConductionEnergyWater[i];

		float airMass = AirMass[i];
		float specificHeat = WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWaterVapor * LastVapor[i];
		AirTemperaturePotential[i] = LastTemperaturePotential[i] + energy / specificHeat;
	}
}


#if !EnergyWaterJobDebug
[BurstCompile]
#endif
public struct EnergyWaterJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirVapor;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float3> SurfaceWind;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	public void Execute(int i)
	{
		if (LastMass[i] > 0)
		{
			float energy = 
				- ConductionEnergyAir[i] 
				- ConductionEnergyIce[i] 
				+ SolarRadiationIn[i]
				+ ConductionEnergyTerrain[i] 
				+ ThermalRadiationDelta[i];
			float specificHeat = WorldData.SpecificHeatWater * LastMass[i] + WorldData.SpecificHeatSalt * LastSaltMass[i];
			Temperature[i] = LastTemperature[i] + energy / specificHeat;
		}
		else
		{
			Temperature[i] = 0;
		}
	}
}


#if !EnergyIceJobDebug
[BurstCompile]
#endif
public struct EnergyIceJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;

	public void Execute(int i)
	{
		if (LastMass[i] > 0)
		{
			float energy = 
				-ConductionEnergyAir[i] 
				+ SolarRadiationIn[i] 
				+ ConductionEnergyWater[i] 
				+ ConductionEnergyTerrain[i] 
				+ ThermalRadiationDelta[i];
			Temperature[i] = LastTemperature[i] + energy / (LastMass[i] * WorldData.SpecificHeatIce);
		}
		else
		{
			Temperature[i] = 0;
		}
	}
}

#if !EnergyTerrainJobDebug
[BurstCompile]
#endif
public struct EnergyTerrainJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> GeothermalEnergy;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyAir[i]
			- ConductionEnergyWater[i]
			- ConductionEnergyIce[i];
		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta + GeothermalEnergy[i];
		float specificHeat = Atmosphere.GetSpecificHeatTerrain(HeatingDepth, Terrain[i].SoilFertility, Terrain[i].Vegetation);
		Temperature[i] = LastTemperature[i] + energy / specificHeat;
	}
}

