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
	[ReadOnly] public NativeArray<float> Dust;
	[ReadOnly] public float CarbonDioxide;
	[ReadOnly] public float EmissivityAir;
	[ReadOnly] public float EmissivityWaterVapor;
	[ReadOnly] public float EmissivityDust;
	[ReadOnly] public float EmissivityCarbonDioxide;
	public void Execute(int i)
	{
		Emissivity[i] = (AirMass[i] * ((1.0f - CarbonDioxide) * EmissivityAir + CarbonDioxide * EmissivityCarbonDioxide) + VaporMass[i] * EmissivityWaterVapor + Dust[i] * EmissivityDust) / (AirMass[i] + VaporMass[i] + Dust[i]);
	}
}

#if !EmissivityWaterJobDebug
[BurstCompile]
#endif
public struct EmissivityWaterJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> WaterMass;
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


#if !CloudAlbedoJobDebug
[BurstCompile]
#endif
public struct CloudAlbedoJob : IJobParallelFor {
	public NativeArray<float> CloudAlbedo;
	public NativeArray<float> CloudAbsorptivity;
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
	public void Execute(int i)
	{
		float cloudIceContent = math.saturate((DewPoint[i] - CloudFreezingTemperatureMin) / (CloudFreezingTemperatureMax - CloudFreezingTemperatureMin));
		float cloudTemperatureAlbedo = WorldData.AlbedoIce * cloudIceContent + WorldData.AlbedoWater * (1.0f - cloudIceContent);
		float rainDropSizeAlbedo = math.saturate(1.0f - CloudDropletMass[i] / CloudMass[i]) * (RainDropSizeAlbedoMax - RainDropSizeAlbedoMin) + RainDropSizeAlbedoMin;

		float cloudCollision = math.saturate(1.0f - math.exp10(-SolarAbsorptivityCloud * CloudMass[i]));
		float albedo = cloudTemperatureAlbedo * rainDropSizeAlbedo * (1.0f - AlbedoSlope[i]) + AlbedoSlope[i];
		CloudAlbedo[i] = albedo;
		CloudAbsorptivity[i] = cloudCollision;
	}
}




#if !SolarAbsorptivityAirJobDebug
[BurstCompile]
#endif
public struct AbsorptivityAirJob : IJobParallelFor {
	public NativeArray<SolarAbsorptivity> AbsorptivitySolar;
	public NativeArray<ThermalAbsorptivity> AbsorptivityThermal;
	[ReadOnly] public NativeArray<float> EmissivityAir;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> Dust;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudAlbedo;
	[ReadOnly] public NativeArray<float> CloudAbsorptivity;
	[ReadOnly] public float CarbonDioxide;
	[ReadOnly] public float AlbedoAir;
	[ReadOnly] public float AlbedoWaterVapor;
	[ReadOnly] public float AlbedoDust;
	[ReadOnly] public float SolarAbsorptivityAir;
	[ReadOnly] public float SolarAbsorptivityWaterVapor;
	[ReadOnly] public float SolarAbsorptivityDust;
	[ReadOnly] public float ThermalAbsorptivityCloud;
	[ReadOnly] public float ThermalAbsorptivityAir;
	[ReadOnly] public float ThermalAbsorptivityWaterVapor;
	[ReadOnly] public float ThermalAbsorptivityDust;
	[ReadOnly] public float EmissivityCloud;
	public void Execute(int i)
	{
		float cloudMass = CloudMass[i];
		float cloudElevation = CloudElevation[i];
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
				thermalAbsorptivityCloud = EmissivityCloud * math.saturate(1.0f - math.exp10(-ThermalAbsorptivityCloud * cloudMass));
				solarAbsorptivityCloud = CloudAbsorptivity[i] * (1.0f - CloudAlbedo[i]);
				albedoCloud = CloudAbsorptivity[i] * CloudAlbedo[i];
			}
		}

		// TODO: this should calculate the chance of hitting and interacting with each consituent, then combine to build an albedo/absorptivity profile
		float fullAbsorptivity = AirMass[i] * ThermalAbsorptivityAir + VaporMass[i] * ThermalAbsorptivityWaterVapor + Dust[i] * ThermalAbsorptivityDust;
		float thermalAbsorptivityAirAbove = EmissivityAir[i] * math.saturate(1.0f - math.exp10(-fullAbsorptivity * (1.0f - belowCloud)));
		float thermalAbsorptivityAirBelow = EmissivityAir[i] * math.saturate(1.0f - math.exp10(-fullAbsorptivity * belowCloud));
		float solarAbsorptivityMass = SolarAbsorptivityAir * AirMass[i] + SolarAbsorptivityWaterVapor * VaporMass[i] + SolarAbsorptivityDust * Dust[i];
		float solarAbsorptivityAbove = math.saturate(1.0f - math.exp10(-(solarAbsorptivityMass) * (1.0f - belowCloud)));
		float solarAbsorptivityBelow = math.saturate(1.0f - math.exp10(-(solarAbsorptivityMass) * belowCloud));

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


