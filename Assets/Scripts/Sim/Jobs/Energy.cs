//#define EnergyAirJobDebug
//#define EnergyWaterSurfaceJobDebug
//#define DISABLE_CORIOLIS
//#define EnergyIceJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !EnergyAirJobDebug
[BurstCompile]
#endif
public struct EnergyAirJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Vapor;
	public NativeArray<float3> Wind;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float3> LastWind;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> WindFriction;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float WindFrictionMultiplier;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{

		float3 lastWind = LastWind[i];

		float3 wind = lastWind;
		wind *= (1.0f - WindFriction[i] * WindFrictionMultiplier);
		wind += PressureGradientForce[i] * SecondsPerTick;

#if !DISABLE_CORIOLIS
//		var windUp = math.dot(Position[i], lastWind) * Position[i];
		var windRight = math.cross(Position[i], lastWind);
		wind += windRight * math.clamp(CoriolisMultiplier[i] * CoriolisTerm * SecondsPerTick, -1, 1);
#endif

		// TODO: deal with buoyancy here, but rather than storing a vertical velocity, let's just move to neutral buoyancy every time step

		//float lastWindVertical = WindVertical[i];
		//WindVertical[i] = lastWindVertical;
		//// TODO: this can overshoot
		//WindVertical[i] += Buoyancy[i] /* * SecondsPerTick */;


		float airMass = AirMass[i];
		float specificHeat = WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWaterVapor * Vapor[i];
		Temperature[i] += Energy[i] / specificHeat;
		Wind[i] = wind;
	}
}


#if !EnergyWaterJobDebug
[BurstCompile]
#endif
public struct EnergyWaterJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> SaltMass;
	public NativeArray<float3> Velocity;
	public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastWaterMass;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float3> LastVelocity;
	[ReadOnly] public NativeArray<float3> Force;
	[ReadOnly] public NativeArray<float> FrozenMass;
	[ReadOnly] public NativeArray<float> EvaporatedMass;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		Mass[i] = LastWaterMass[i] - FrozenMass[i] - EvaporatedMass[i];
		SaltMass[i] = LastSaltMass[i];
		if (Mass[i] > 0)
		{
			float specificHeat = WorldData.SpecificHeatWater * Mass[i] + WorldData.SpecificHeatSalt * SaltMass[i];
			Temperature[i] = LastTemperature[i] + Energy[i] / specificHeat;
			Velocity[i] = LastVelocity[i] + Force[i] * SecondsPerTick;
		} else
		{
			Temperature[i] = 0;
			Velocity[i] = 0;
		}
	}
}


#if !EnergyIceJobDebug
[BurstCompile]
#endif
public struct EnergyIceJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> FluxEnergyIce;
	[ReadOnly] public float IceHeatingDepth;
	public void Execute(int i)
	{
		if (Mass[i] > 0)
		{
			Temperature[i] += FluxEnergyIce[i] / (Mass[i] * WorldData.SpecificHeatIce);
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
	[ReadOnly] public NativeArray<float> VegetationCoverage;
	[ReadOnly] public float GeothermalEnergy;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyAir[i]
			- ConductionEnergyWater[i]
			- ConductionEnergyIce[i];
		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta + GeothermalEnergy;
		float specificHeat = Atmosphere.GetSpecificHeatTerrain(HeatingDepth, Terrain[i].SoilFertility, VegetationCoverage[i]);
		Temperature[i] = LastTemperature[i] + energy / specificHeat;
	}
}

