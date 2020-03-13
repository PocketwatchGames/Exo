using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !EmissivityAirJobDebug
[BurstCompile]
#endif
public struct EmissivityAirJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	public void Execute(int i)
	{
		Emissivity[i] = (AirMass[i] * WorldData.EmissivityAir + VaporMass[i] * WorldData.EmissivityWaterVapor) / (AirMass[i] + VaporMass[i]);
	}
}

#if !EmissivityWaterJobDebug
[BurstCompile]
#endif
public struct EmissivityWaterJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> WaterMass;
	public void Execute(int i)
	{
		Emissivity[i] = math.saturate((WorldData.EmissivityWater * WaterMass[i] + WorldData.EmissivitySalt * SaltMass[i]) / (WaterMass[i] + SaltMass[i]));
	}
}

#if !EmissivityTerrainJobDebug
[BurstCompile]
#endif
public struct EmissivityTerrainJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> VegetationCoverage;
	public void Execute(int i)
	{
		Emissivity[i] = math.saturate(math.lerp(
			math.lerp(WorldData.EmissivitySand, WorldData.EmissivityDirt, Terrain[i].SoilFertility),
			WorldData.EmissivityVegetation,
			VegetationCoverage[i]));
	}
}


#if !SolarAbsorptivityAirJobDebug
[BurstCompile]
#endif
public struct AbsorptivityAirJob : IJobParallelFor {
	public NativeArray<SolarAbsorptivity> AbsorptivitySolar;
	public NativeArray<ThermalAbsorptivity> AbsorptivityThermal;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> DewPoint;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> WaterSlopeAlbedo;
	[ReadOnly] public float CarbonDioxide;
	[ReadOnly] public float SolarReflectivityAir;
	[ReadOnly] public float SolarReflectivityCloud;
	[ReadOnly] public float SolarAbsorptivityAir;
	[ReadOnly] public float SolarAbsorptivityWaterVapor;
	[ReadOnly] public float SolarAbsorptivityCloud;
	[ReadOnly] public float ThermalAbsorptivityAir;
	[ReadOnly] public float ThermalAbsorptivityWaterVapor;
	[ReadOnly] public float ThermalAbsorptivityWater;
	[ReadOnly] public float ThermalAbsorptivityCarbonDioxide;
	[ReadOnly] public float CloudFreezingTemperatureMin;
	[ReadOnly] public float CloudFreezingTemperatureMax;
	[ReadOnly] public float RainDropSizeAlbedoMin;
	[ReadOnly] public float RainDropSizeAlbedoMax;
	[ReadOnly] public float CloudSlopeAlbedoMax;
	public void Execute(int i)
	{
		float airMass = AirMass[i];
		float waterVaporMass = VaporMass[i];


		float cloudMass = CloudMass[i];
		float cloudElevation = CloudElevation[i];
		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];
		bool isCloudLayer = cloudElevation >= layerElevation && cloudElevation < layerElevation + layerHeight;
		float solarAbsorptivityCloud = 0;
		float solarReflectivityCloud = 0;
		float thermalAbsorptivityCloud = 0;

		if (isCloudLayer && cloudMass > 0)
		{
			float cloudIceContent = math.saturate((DewPoint[i] - CloudFreezingTemperatureMin) / (CloudFreezingTemperatureMax - CloudFreezingTemperatureMin));
			float cloudTemperatureAlbedo = WorldData.AlbedoIce + (WorldData.AlbedoWater - WorldData.AlbedoIce) * cloudIceContent;
			float rainDropSizeAlbedo = math.saturate(1.0f - CloudDropletMass[i] / cloudMass) * (RainDropSizeAlbedoMax - RainDropSizeAlbedoMin) + RainDropSizeAlbedoMin;
			float cloudAlbedo = math.min(1.0f, SolarReflectivityCloud * cloudTemperatureAlbedo * cloudMass * rainDropSizeAlbedo / math.max(CloudSlopeAlbedoMax, 1.0f - WaterSlopeAlbedo[i]));

			solarReflectivityCloud = math.saturate(1.0f - math.exp10(-cloudAlbedo * cloudMass));
			solarAbsorptivityCloud = math.saturate(1.0f - math.exp10(-SolarAbsorptivityCloud * cloudMass));
			thermalAbsorptivityCloud = math.saturate(1.0f - math.exp10(-ThermalAbsorptivityWater * cloudMass));
		}

		AbsorptivitySolar[i] = new SolarAbsorptivity
		{
			AbsorptivityAir = math.saturate(1.0f - math.exp10(-SolarAbsorptivityAir * airMass + SolarAbsorptivityWaterVapor * waterVaporMass)),
			ReflectivityAir = math.saturate(1.0f - math.exp10(-math.saturate(SolarReflectivityAir * (airMass + waterVaporMass)))),
			AbsorptivityCloud = solarAbsorptivityCloud,
			ReflectivityCloud = solarReflectivityCloud
		};

		AbsorptivityThermal[i] = new ThermalAbsorptivity() {
			AbsorptivityCloud = thermalAbsorptivityCloud,
			AbsorptivityAir = math.saturate(1.0f - math.exp10(-AirMass[i] * ((1.0f - CarbonDioxide) * ThermalAbsorptivityAir + CarbonDioxide * ThermalAbsorptivityCarbonDioxide) + VaporMass[i] * ThermalAbsorptivityWaterVapor))
		};

	}
}


