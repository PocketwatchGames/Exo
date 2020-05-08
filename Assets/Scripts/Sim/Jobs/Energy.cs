
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


[BurstCompile]
public struct EnergyAirSurfaceJob : IJobParallelFor {
	public NativeSlice<float> AirTemperaturePotential;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeSlice<float> LastVapor;
	[ReadOnly] public NativeSlice<float> LastTemperaturePotential;
	[ReadOnly] public NativeSlice<float> ThermalRadiationDelta;
	[ReadOnly] public NativeSlice<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeSlice<float> LayerElevation;
	[ReadOnly] public NativeSlice<float> LayerHeight;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		float energy =
			SolarRadiationIn[i]
			+ ThermalRadiationDelta[i]
			+ ConductionEnergyIce[i]
			+ ConductionEnergyTerrain[i]
			+ ConductionEnergyWater[i];

		int columnIndex = i % Count;
		float cloudMassInLayer = Atmosphere.GetCloudMassInLayer(CloudMass[columnIndex], CloudElevation[columnIndex], LayerElevation[i], LayerHeight[i]);
		float specificHeat = Atmosphere.GetSpecificHeatAir(AirMass[i], LastVapor[i], cloudMassInLayer);
		AirTemperaturePotential[i] = LastTemperaturePotential[i] + energy / specificHeat;
	}
}
[BurstCompile]
public struct EnergyAirJob : IJobParallelFor {
	public NativeSlice<float> AirTemperaturePotential;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeSlice<float> LastVapor;
	[ReadOnly] public NativeSlice<float> LastTemperaturePotential;
	[ReadOnly] public NativeSlice<float> ThermalRadiationDelta;
	[ReadOnly] public NativeSlice<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeSlice<float> LayerElevation;
	[ReadOnly] public NativeSlice<float> LayerHeight;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		float energy =
			SolarRadiationIn[i]
			+ ThermalRadiationDelta[i];

		int columnIndex = i % Count;
		float cloudMassInLayer = Atmosphere.GetCloudMassInLayer(CloudMass[columnIndex], CloudElevation[columnIndex], LayerElevation[i], LayerHeight[i]);
		float specificHeat = Atmosphere.GetSpecificHeatAir(AirMass[i], LastVapor[i], cloudMassInLayer);
		AirTemperaturePotential[i] = LastTemperaturePotential[i] + energy / specificHeat;
	}
}


[BurstCompile]
public struct EnergyWaterJob : IJobParallelFor {
	public NativeSlice<float> Temperature;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		int columnIndex = i % Count;
		int index = i + Count;
		int indexDown = i;
		int indexUp = index + Count;
		if (LastMass[index] > 0)
		{
			float energy = ThermalRadiationDelta[columnIndex];

			energy += (1.0f - Coverage[indexDown]) * (
				+ ConductionEnergyTerrain[columnIndex]
				);
			energy += (1.0f - Coverage[indexUp]) * (
				- ConductionEnergyAir[columnIndex]
				- ConductionEnergyIce[columnIndex]
				);

			float specificHeat = WorldData.SpecificHeatWater * LastMass[index] + WorldData.SpecificHeatSalt * LastSaltMass[index];
			Temperature[i] = LastTemperature[index] + energy / specificHeat;
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
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;

	public void Execute(int i)
	{
		if (LastMass[i] > 0)
		{
			float energy =
				-ConductionEnergyAir[i]
				+ ConductionEnergyWater[i]
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
public struct EnergyLavaJob : IJobParallelFor {
	public NativeArray<float> LavaTemperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LavaMass;
	[ReadOnly] public float Emissivity;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		if (LavaMass[i] > 0)
		{
			float energy = -Atmosphere.GetRadiationRate(LastTemperature[i], Emissivity) * SecondsPerTick;
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
	[ReadOnly] public NativeArray<float> SpecificHeatTerrain;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> GeothermalEnergy;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyAir[i]
			- ConductionEnergyWater[i]
			- ConductionEnergyIce[i];

		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta + GeothermalEnergy[i];
		TerrainTemperature[i] = LastTemperature[i] + energy / SpecificHeatTerrain[i];
	}
}

