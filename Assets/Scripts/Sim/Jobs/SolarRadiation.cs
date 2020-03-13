//#define SolarRadiationAbsorbedAirJobDebug
//#define SolarRadiationAbsorbedTerrainJobDebug

using Unity.Burst;
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
	public NativeArray<float> WaterSlopeAlbedo;
	[ReadOnly] public NativeArray<float3> SphericalPosition;
	[ReadOnly] public float3 SunToPlanetDir;
	[ReadOnly] public quaternion PlanetRotation;
	[ReadOnly] public float IncomingSolarRadiation;
	[ReadOnly] public float IncomingGeothermalRadiation;
	public void Execute(int i)
	{
		float sunDotSurface = math.max(0, math.dot(SunToPlanetDir, math.rotate(PlanetRotation, -SphericalPosition[i])));
		WaterSlopeAlbedo[i] = math.pow(1.0f - math.max(0, sunDotSurface), 9);
		float r = IncomingSolarRadiation * sunDotSurface;
		GeothermalRadiation[i] = IncomingGeothermalRadiation;
		SolarRadiation[i] = r;
		DisplaySolarRadiation[i] = r;
	}
}



#if !SolarRadiationAbsorbedAirJobDebug
[BurstCompile]
#endif
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
		incomingRadiation -= reflectedAirAbove;
		float absorbedAirAbove = incomingRadiation * AbsorptivitySolar[i].AbsorptivityAirAbove;
		incomingRadiation -= absorbedAirAbove;

		reflectedCloud = incomingRadiation * AbsorptivitySolar[i].ReflectivityCloud;
		incomingRadiation -= reflectedCloud;
		absorbedCloud = incomingRadiation * AbsorptivitySolar[i].AbsorptivityCloud;
		incomingRadiation -= absorbedCloud;

		float reflectedAirBelow = incomingRadiation * AbsorptivitySolar[i].ReflectivityAirBelow;
		incomingRadiation -= reflectedAirBelow;
		float absorbedAirBelow = incomingRadiation * AbsorptivitySolar[i].AbsorptivityAirBelow;
		incomingRadiation -= absorbedAirBelow;

		SolarRadiationAbsorbed[i] = absorbedAirAbove + absorbedAirBelow + absorbedCloud;
		SolarRadiationIncoming[i] = incomingRadiation;
		SolarRadiationReflected[i] = reflectedAirAbove + reflectedAirBelow + reflectedCloud;
		SolarRadiationReflectedCloud[i] = reflectedCloud;
		SolarRadiationAbsorbedCloud[i] = absorbedCloud;
	}
}

#if !SolarRadiationAbsorbedIceJobDebug
[BurstCompile]
#endif
public struct SolarRadiationAbsorbedIceJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationIncoming;
	public NativeArray<float> SolarRadiationReflected;
	[ReadOnly] public float AlbedoIce;
	[ReadOnly] public NativeArray<float> IceCoverage;
	public void Execute(int i)
	{
		float incoming = SolarRadiationIncoming[i];
		float iceCoverage = IceCoverage[i];
		float reflected = incoming * (AlbedoIce * iceCoverage);
		incoming -= reflected;
		float absorbed = incoming * iceCoverage;
		SolarRadiationAbsorbed[i] = absorbed;
		SolarRadiationIncoming[i] = incoming - absorbed;
		SolarRadiationReflected[i] = reflected;
	}
}

#if !SolarRadiationAbsorbedWaterJobDebug
[BurstCompile]
#endif
public struct SolarRadiationAbsorbedWaterJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationIncoming;
	public NativeArray<float> SolarRadiationReflected;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> WaterSlopeAlbedo;
	public void Execute(int i)
	{
		float incoming = SolarRadiationIncoming[i];
		float waterCoverage = WaterCoverage[i];
		float reflected = incoming * (WaterSlopeAlbedo[i] * waterCoverage);
		incoming -= reflected;
		float absorbed = incoming * waterCoverage;
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
	[ReadOnly] public NativeArray<float> VegetationCoverage;
	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		float incoming = SolarRadiationIncoming[i];

		float slopeAlbedo = 0;
		float vegetationCoverage = VegetationCoverage[i];
		float soilReflectivity = Atmosphere.GetAlbedo(WorldData.AlbedoLand - worldData.AlbedoReductionSoilQuality * LastTerrain[i].SoilFertility, slopeAlbedo);
		float reflected = incoming * (vegetationCoverage * WorldData.AlbedoFoliage + (1.0f - vegetationCoverage) * soilReflectivity);
		incoming -= reflected;

		SolarRadiationReflected[i] = reflected;
		SolarRadiationAbsorbed[i] = incoming;
	}
}