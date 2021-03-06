﻿using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


[BurstCompile]
public struct EmissivityAirJob : IJobParallelFor {
	public NativeSlice<float> Emissivity;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeSlice<float> VaporMass;
	[ReadOnly] public NativeSlice<float> Dust;
	[ReadOnly] public NativeSlice<float> CarbonDioxide;
	[ReadOnly] public float EmissivityAir;
	[ReadOnly] public float EmissivityCarbonDioxide;
	[ReadOnly] public float EmissivityWaterVapor;
	[ReadOnly] public float EmissivityDust;
	public void Execute(int i)
	{
		// TODO: this should calculate the chance of hitting and interacting with each consituent, then combine to build an albedo/absorptivity profile
		float emissivity =
			((AirMass[i] - CarbonDioxide[i]) * EmissivityAir
			+ CarbonDioxide[i] * EmissivityCarbonDioxide
			+ VaporMass[i] * EmissivityWaterVapor
			+ Dust[i] * EmissivityDust);
		Emissivity[i] = math.exp(-emissivity);
	}
}

[BurstCompile]
public struct EmissivityWaterJob : IJobParallelFor {
	public NativeSlice<float> Emissivity;
	[ReadOnly] public NativeSlice<float> SaltMass;
	[ReadOnly] public NativeSlice<float> WaterMass;
	[ReadOnly] public float EmissivityWater;
	[ReadOnly] public float EmissivitySalt;
	public void Execute(int i)
	{
		Emissivity[i] = math.saturate((EmissivityWater * WaterMass[i] + EmissivitySalt * SaltMass[i]) / (WaterMass[i] + SaltMass[i]));
	}
}

#if !EmissivityTerrainJobDebug
[BurstCompile]
#endif
public struct EmissivityTerrainJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public float EmissivitySand;
	[ReadOnly] public float EmissivityDirt;
	public void Execute(int i)
	{
		Emissivity[i] = math.saturate(math.lerp(EmissivitySand, EmissivityDirt, SoilFertility[i]));
	}
}


[BurstCompile]
public struct AlbedoCloudJob : IJobParallelFor {
	public NativeArray<float> Absorptivity;
	public NativeArray<float> Albedo;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> DewPoint;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> AlbedoSlope;
	[ReadOnly] public float CloudFreezingTemperatureMin;
	[ReadOnly] public float CloudFreezingTemperatureMax;
	[ReadOnly] public float RainDropSizeAlbedoMin;
	[ReadOnly] public float RainDropSizeAlbedoMax;
	[ReadOnly] public float SolarAbsorptivityCloud;
	[ReadOnly] public float AlbedoIceMin;
	[ReadOnly] public float AlbedoIceRange;
	[ReadOnly] public float AlbedoWaterMin;
	[ReadOnly] public float AlbedoWaterRange;
	public void Execute(int i)
	{
		float cloudCollision = math.exp(-SolarAbsorptivityCloud * CloudMass[i]);
		Absorptivity[i] = cloudCollision;

		float cloudIceContent = math.saturate((DewPoint[i] - CloudFreezingTemperatureMin) / (CloudFreezingTemperatureMax - CloudFreezingTemperatureMin));
		float cloudTemperatureAlbedoMin = math.lerp(AlbedoIceMin, AlbedoWaterMin, cloudIceContent);
		float cloudTemperatureAlbedoRange = math.lerp(AlbedoIceRange, AlbedoWaterRange, cloudIceContent);
		float rainDropSizeAlbedo = math.saturate(1.0f - CloudDropletMass[i] / CloudMass[i]) * (RainDropSizeAlbedoMax - RainDropSizeAlbedoMin) + RainDropSizeAlbedoMin;
		float albedo = (cloudTemperatureAlbedoMin + cloudTemperatureAlbedoRange * AlbedoSlope[i]) * rainDropSizeAlbedo;
		Albedo[i] = albedo;
	}
}



[BurstCompile]
public struct AbsorptivityAirJob : IJobParallelFor {
	public NativeSlice<SolarAbsorptivity> AbsorptivitySolar;
	public NativeSlice<ThermalAbsorptivity> AbsorptivityThermal;
	[ReadOnly] public NativeSlice<float> VaporMass;
	[ReadOnly] public NativeSlice<float> Dust;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeSlice<float> AirCarbonDioxide;
	[ReadOnly] public NativeSlice<float> Emissivity;
	[ReadOnly] public NativeSlice<float> LayerElevation;
	[ReadOnly] public NativeSlice<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudAbsorptivitySolar;
	[ReadOnly] public NativeArray<float> CloudAlbedo;
	[ReadOnly] public float AlbedoAir;
	[ReadOnly] public float AlbedoWaterVapor;
	[ReadOnly] public float AlbedoDust;
	[ReadOnly] public float SolarAbsorptivityAir;
	[ReadOnly] public float SolarAbsorptivityWaterVapor;
	[ReadOnly] public float SolarAbsorptivityDust;
	[ReadOnly] public float EmissivityWater;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		int cloudIndex = i % Count;
		float cloudMass = CloudMass[cloudIndex];
		float cloudElevation = CloudElevation[cloudIndex];
		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];

		float belowCloud = math.saturate((cloudElevation - layerElevation) / layerHeight);
		float solarAbsorptivityCloud = 0;
		float albedoCloud = 0;
		float thermalAbsorptivityCloud = 0;
		if (cloudMass > 0)
		{
			if (cloudElevation >= layerElevation && cloudElevation < layerElevation + layerHeight)
			{
				thermalAbsorptivityCloud = 1.0f - math.exp(-EmissivityWater * cloudMass);
				solarAbsorptivityCloud = CloudAbsorptivitySolar[cloudIndex] * (1.0f - CloudAlbedo[cloudIndex]);

				albedoCloud = CloudAbsorptivitySolar[cloudIndex] * CloudAlbedo[cloudIndex];
			}
		}

		// TODO: this should calculate the chance of hitting and interacting with each consituent, then combine to build an albedo/absorptivity profile

		float thermalAbsorptivityAirAbove = 1.0f - Emissivity[i] * (1.0f - belowCloud);
		float thermalAbsorptivityAirBelow = 1.0f - Emissivity[i] * belowCloud;
		float solarAbsorptivityMass = SolarAbsorptivityAir * AirMass[i] + SolarAbsorptivityWaterVapor * VaporMass[i] + SolarAbsorptivityDust * Dust[i];
		float solarAbsorptivityAbove = math.saturate(1.0f - math.exp(-solarAbsorptivityMass * (1.0f - belowCloud)));
		float solarAbsorptivityBelow = math.saturate(1.0f - math.exp(-solarAbsorptivityMass * belowCloud));

		AbsorptivitySolar[i] = new SolarAbsorptivity
		{
			AbsorptivityAirAbove = solarAbsorptivityAbove * (1.0f - AlbedoAir),
			ReflectivityAirAbove = solarAbsorptivityAbove * AlbedoAir,
			AbsorptivityAirBelow = solarAbsorptivityBelow * (1.0f - AlbedoAir),
			ReflectivityAirBelow = solarAbsorptivityBelow * AlbedoAir,
			AbsorptivityCloud = solarAbsorptivityCloud,
			ReflectivityCloud = albedoCloud,			
		};

		AbsorptivityThermal[i] = new ThermalAbsorptivity() {
			AbsorptivityCloud = thermalAbsorptivityCloud,
			AbsorptivityAirAbove = thermalAbsorptivityAirAbove,
			AbsorptivityAirBelow = thermalAbsorptivityAirBelow
		};

	}
}


