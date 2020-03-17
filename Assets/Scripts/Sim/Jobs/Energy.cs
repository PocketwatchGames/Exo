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
public struct EnergyAirSurfaceJob : IJobParallelFor {
	public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastTemperaturePotential;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyFlora;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		float energy =
			SolarRadiationIn[i]
			+ ThermalRadiationDelta[i]
			+ ConductionEnergyIce[i]
			+ ConductionEnergyTerrain[i]
			+ ConductionEnergyFlora[i]
			+ ConductionEnergyWater[i];

		float specificHeat = WorldData.SpecificHeatAtmosphere * AirMass[i] + WorldData.SpecificHeatWaterVapor * LastVapor[i];
		AirTemperaturePotential[i] = LastTemperaturePotential[i] + energy / specificHeat;
	}
}
#if !EnergyAirJobDebug
[BurstCompile]
#endif
public struct EnergyAirJob : IJobParallelFor {
	public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastTemperaturePotential;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	public void Execute(int i)
	{
		float energy =
			SolarRadiationIn[i]
			+ ThermalRadiationDelta[i];

		float specificHeat = WorldData.SpecificHeatAtmosphere * AirMass[i] + WorldData.SpecificHeatWaterVapor * LastVapor[i];
		AirTemperaturePotential[i] = LastTemperaturePotential[i] + energy / specificHeat;
	}
}


#if !EnergyWaterJobDebug
[BurstCompile]
#endif
public struct EnergyWaterJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> CoverageUp;
	[ReadOnly] public NativeArray<float> CoverageDown;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyFlora;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		if (LastMass[i] > 0)
		{
			float energy =
				+ SolarRadiationIn[i]
				+ ThermalRadiationDelta[i];

			energy += (1.0f - CoverageDown[i]) * (
				+ConductionEnergyFlora[i]
				+ ConductionEnergyTerrain[i]
				);
			energy += (1.0f - CoverageUp[i]) * (
				- ConductionEnergyAir[i]
				- ConductionEnergyIce[i]
				);

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
	[ReadOnly] public NativeArray<float> ConductionEnergyFlora;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;

	public void Execute(int i)
	{
		if (LastMass[i] > 0)
		{
			float energy =
				//- ConductionEnergyAir[i] 
				//+ ConductionEnergyWater[i]
				//+ ConductionEnergyFlora[i]
				//+ ConductionEnergyTerrain[i]
				+ SolarRadiationIn[i]
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
public struct EnergyFloraJob : IJobParallelFor {
	public NativeArray<float> FloraTemperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> FloraMass;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		if (FloraMass[i] > 0)
		{
			float conductionDelta =
				-ConductionEnergyAir[i]
				- ConductionEnergyWater[i]
				- ConductionEnergyIce[i]
				+ ConductionEnergyTerrain[i];


			float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta;
			FloraTemperature[i] = LastTemperature[i] + energy / (FloraMass[i] * WorldData.SpecificHeatFlora);
		} else
		{
			FloraTemperature[i] = 0;
		}

	}
}



#if !EnergyTerrainJobDebug
[BurstCompile]
#endif
public struct EnergyTerrainJob : IJobParallelFor {
	public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyFlora;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> GeothermalEnergy;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyAir[i]
			- ConductionEnergyWater[i]
			- ConductionEnergyIce[i]
			- ConductionEnergyFlora[i];

		float specificHeatTerrain = Atmosphere.GetSpecificHeatTerrain(HeatingDepth, Terrain[i].SoilFertility);

		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta + GeothermalEnergy[i];
		TerrainTemperature[i] = LastTemperature[i] + energy / specificHeatTerrain;
	}
}

