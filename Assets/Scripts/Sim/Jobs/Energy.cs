
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


[BurstCompile]
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
	[ReadOnly] public NativeArray<float> ConductionEnergyLava;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		float energy =
			SolarRadiationIn[i]
			+ ThermalRadiationDelta[i]
			+ ConductionEnergyIce[i]
			+ ConductionEnergyTerrain[i]
			+ ConductionEnergyLava[i]
			+ ConductionEnergyFlora[i]
			+ ConductionEnergyWater[i];

		float specificHeat = WorldData.SpecificHeatAtmosphere * AirMass[i] + WorldData.SpecificHeatWaterVapor * LastVapor[i];
		AirTemperaturePotential[i] = LastTemperaturePotential[i] + energy / specificHeat;
	}
}
[BurstCompile]
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


[BurstCompile]
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
	[ReadOnly] public NativeArray<float> ConductionEnergyLava;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		if (LastMass[i] > 0)
		{
			float energy =
				+ SolarRadiationIn[i]
				+ ThermalRadiationDelta[i];

			energy += (1.0f - CoverageDown[i]) * (
				+ ConductionEnergyLava[i]
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


[BurstCompile]
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
	[ReadOnly] public NativeArray<float> ConductionEnergyLava;

	public void Execute(int i)
	{
		if (LastMass[i] > 0)
		{
			float energy =
				-ConductionEnergyAir[i]
				+ ConductionEnergyWater[i]
				+ ConductionEnergyFlora[i]
				+ ConductionEnergyLava[i]
				+ ConductionEnergyTerrain[i]
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

[BurstCompile]
public struct EnergyFloraJob : IJobParallelFor {
	public NativeArray<float> FloraTemperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> FloraMass;
	[ReadOnly] public NativeArray<float> FloraWater;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		if (FloraMass[i] > 0)
		{
			float conductionDelta =
				-ConductionEnergyAir[i]
				- ConductionEnergyIce[i]
				+ ConductionEnergyTerrain[i];


			float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta;
			FloraTemperature[i] = LastTemperature[i] + energy / (FloraMass[i] * WorldData.SpecificHeatFlora + FloraWater[i] * WorldData.SpecificHeatWater);
		}
		else
		{
			FloraTemperature[i] = 0;
		}

	}
}


[BurstCompile]
public struct EnergyLavaJob : IJobParallelFor {
	public NativeArray<float> LavaTemperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LavaMass;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		if (LavaMass[i] > 0)
		{
			float conductionDelta =
				-ConductionEnergyAir[i]
				- ConductionEnergyIce[i]
				- ConductionEnergyWater[i]
				+ ConductionEnergyTerrain[i];

			float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta;
			LavaTemperature[i] = LastTemperature[i] + energy / (LavaMass[i] * WorldData.SpecificHeatLava);
		}
		else
		{
			LavaTemperature[i] = 0;
		}

	}
}



[BurstCompile]
public struct EnergyTerrainJob : IJobParallelFor {
	public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyFlora;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyLava;
	[ReadOnly] public NativeArray<float> GeothermalEnergy;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyAir[i]
			- ConductionEnergyWater[i]
			- ConductionEnergyIce[i]
			- ConductionEnergyFlora[i]
			- ConductionEnergyLava[i];

		float specificHeatTerrain = Atmosphere.GetSpecificHeatTerrain(HeatingDepth, SoilFertility[i]);

		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta + GeothermalEnergy[i];
		TerrainTemperature[i] = LastTemperature[i] + energy / specificHeatTerrain;
	}
}

