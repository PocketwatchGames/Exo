//#define EnergyAirJobDebug
//#define EnergyWaterSurfaceJobDebug
//#define EnergyIceJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !EnergyAirJobDebug
[BurstCompile]
#endif
public struct EnergyAirJob : IJobParallelFor {
	public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> Energy;
	public void Execute(int i)
	{

		// TODO: deal with buoyancy here, but rather than storing a vertical velocity, let's just move to neutral buoyancy every time step

		//float lastWindVertical = WindVertical[i];
		//WindVertical[i] = lastWindVertical;
		//// TODO: this can overshoot
		//WindVertical[i] += Buoyancy[i] /* * SecondsPerTick */;

		float airMass = AirMass[i];
		float specificHeat = WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWaterVapor * Vapor[i];
		AirTemperaturePotential[i] = LastTemperature[i] + Energy[i] / specificHeat;
	}
}


#if !EnergyWaterJobDebug
[BurstCompile]
#endif
public struct EnergyWaterJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Energy;
	public void Execute(int i)
	{
		if (Mass[i] > 0)
		{
			float specificHeat = WorldData.SpecificHeatWater * Mass[i] + WorldData.SpecificHeatSalt * SaltMass[i];
			Temperature[i] = LastTemperature[i] + Energy[i] / specificHeat;
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
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> FluxEnergyIce;
	public void Execute(int i)
	{
		if (Mass[i] > 0)
		{
			Temperature[i] = LastTemperature[i] + FluxEnergyIce[i] / (Mass[i] * WorldData.SpecificHeatIce);
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
	[ReadOnly] public float GeothermalEnergy;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyAir[i]
			- ConductionEnergyWater[i]
			- ConductionEnergyIce[i];
		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta + GeothermalEnergy;
		float specificHeat = Atmosphere.GetSpecificHeatTerrain(HeatingDepth, Terrain[i].SoilFertility, Terrain[i].Vegetation);
		Temperature[i] = LastTemperature[i] + energy / specificHeat;
	}
}

