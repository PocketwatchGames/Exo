﻿using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

public struct SolarAbsorptivity {
	public float AbsorptivityAirAbove;
	public float ReflectivityAirAbove;
	public float AbsorptivityAirBelow;
	public float ReflectivityAirBelow;
	public float AbsorptivityCloud;
	public float ReflectivityCloud;
}

public struct ThermalAbsorptivity {
	public float AbsorptivityAirAbove;
	public float AbsorptivityAirBelow;
	public float AbsorptivityCloud;
}

#if !SolarRadiationJobDebug
[BurstCompile]
#endif
public struct SolarRadiationJob : IJobParallelFor {
	public NativeArray<float> SolarRadiation;
	public NativeArray<float> GeothermalRadiation;
	public NativeArray<float> DisplaySolarRadiation;
	public NativeArray<float> AlbedoSlope;
	[ReadOnly] public NativeArray<float3> SphericalPosition;
	[ReadOnly] public float3 SunToPlanetDir;
	[ReadOnly] public quaternion PlanetRotation;
	[ReadOnly] public float IncomingSolarRadiation;
	[ReadOnly] public float IncomingGeothermalRadiation;
	public void Execute(int i)
	{
		float sunDotSurface = math.max(0, math.dot(SunToPlanetDir, math.rotate(PlanetRotation, -SphericalPosition[i])));
		AlbedoSlope[i] = math.pow(1.0f - math.max(0, sunDotSurface), 9);
		float r = IncomingSolarRadiation * sunDotSurface;
		GeothermalRadiation[i] = IncomingGeothermalRadiation;
		SolarRadiation[i] = r;
		DisplaySolarRadiation[i] = r;
	}
}



[BurstCompile]
public struct SolarRadiationAbsorbedAirJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationIncoming;
	public NativeArray<float> SolarRadiationReflected;
	public NativeArray<float> SolarRadiationAbsorbedCloud;
	public NativeArray<float> SolarRadiationReflectedCloud;
	[ReadOnly] public NativeArray<SolarAbsorptivity> AbsorptivitySolar;
	public void Execute(int i)
	{
		float incomingRadiation = SolarRadiationIncoming[i];

		float absorbedCloud = 0;
		float reflectedCloud = 0;

		float reflectedAirAbove = incomingRadiation * AbsorptivitySolar[i].ReflectivityAirAbove;
		float absorbedAirAbove = incomingRadiation * AbsorptivitySolar[i].AbsorptivityAirAbove;
		incomingRadiation -= absorbedAirAbove + reflectedAirAbove;

		reflectedCloud = incomingRadiation * AbsorptivitySolar[i].ReflectivityCloud;
		absorbedCloud = incomingRadiation * AbsorptivitySolar[i].AbsorptivityCloud;
		incomingRadiation -= absorbedCloud + reflectedCloud;

		float reflectedAirBelow = incomingRadiation * AbsorptivitySolar[i].ReflectivityAirBelow;
		float absorbedAirBelow = incomingRadiation * AbsorptivitySolar[i].AbsorptivityAirBelow;
		incomingRadiation -= absorbedAirBelow + reflectedAirBelow;

		SolarRadiationAbsorbed[i] = absorbedAirAbove + absorbedAirBelow + absorbedCloud;
		SolarRadiationIncoming[i] = incomingRadiation;
		SolarRadiationReflected[i] = reflectedAirAbove + reflectedAirBelow + reflectedCloud;
		SolarRadiationReflectedCloud[i] = reflectedCloud;
		SolarRadiationAbsorbedCloud[i] = absorbedCloud;
	}
}

#if !SolarRadiationAbsorbedPartialCoverageJobDebug
[BurstCompile]
#endif
public struct SolarRadiationAbsorbedPartialCoverageConstantAlbedoJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationIncoming;
	public NativeArray<float> SolarRadiationReflected;
	[ReadOnly] public NativeArray<float> AlbedoSlope;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public float AlbedoMin;
	[ReadOnly] public float AlbedoRange;
	public void Execute(int i)
	{
		float incoming = SolarRadiationIncoming[i];
		float coverage = Coverage[i];
		float reflected = incoming * coverage * (AlbedoMin + AlbedoRange * AlbedoSlope[i]);
		incoming -= reflected;
		float absorbed = incoming * coverage;
		SolarRadiationAbsorbed[i] = absorbed;
		SolarRadiationIncoming[i] = incoming - absorbed;
		SolarRadiationReflected[i] = reflected;
	}
}

#if !SolarRadiationAbsorbedPartialCoverageJobDebug
[BurstCompile]
#endif
public struct SolarRadiationAbsorbedSlopeJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationIncoming;
	public NativeArray<float> SolarRadiationReflected;
	[ReadOnly] public float AlbedoMin;
	[ReadOnly] public float AlbedoRange;
	[ReadOnly] public NativeArray<float> AlbedoSlope;
	[ReadOnly] public NativeArray<float> Coverage;
	public void Execute(int i)
	{
		float incoming = SolarRadiationIncoming[i];
		float coverage = Coverage[i];
		float reflected = incoming * coverage * (AlbedoMin + AlbedoRange * AlbedoSlope[i]);
		incoming -= reflected;
		float absorbed = incoming * coverage;
		SolarRadiationAbsorbed[i] = absorbed;
		SolarRadiationIncoming[i] = incoming - absorbed;
		SolarRadiationReflected[i] = reflected;
	}
}

#if !SolarRadiationAbsorbedTerrainJobDebug
[BurstCompile]
#endif
public struct SolarRadiationAbsorbedTerrainJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationReflected;
	[ReadOnly] public NativeArray<float> SolarRadiationIncoming;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		float incoming = SolarRadiationIncoming[i];

		float slopeAlbedo = 0;
		float soilReflectivity = Atmosphere.GetAlbedo(WorldData.AlbedoLand - worldData.AlbedoReductionSoilQuality * SoilFertility[i], slopeAlbedo);
		float reflected = incoming * soilReflectivity;
		incoming -= reflected;

		SolarRadiationReflected[i] = reflected;
		SolarRadiationAbsorbed[i] = incoming;
	}
}