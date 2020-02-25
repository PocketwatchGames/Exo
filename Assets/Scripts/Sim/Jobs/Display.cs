using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;



#if !UpdateDisplayJobDebug
[BurstCompile]
#endif
public struct UpdateDisplayJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbedSurface;
	public NativeArray<float> DisplayRainfall;
	public NativeArray<float> DisplayEvaporation;
	[ReadOnly] public NativeArray<float> SolarRadiationInTerrain;
	[ReadOnly] public NativeArray<float> SolarRadiationInIce;
	[ReadOnly] public NativeArray<float> SolarRadiationInWaterSurface;
	[ReadOnly] public NativeArray<float> RainfallWater;
	[ReadOnly] public NativeArray<float> Evaporation;
	public void Execute(int i)
	{
		SolarRadiationAbsorbedSurface[i] = SolarRadiationInTerrain[i] + SolarRadiationInIce[i] + SolarRadiationInWaterSurface[i];
		DisplayRainfall[i] = RainfallWater[i];
		DisplayEvaporation[i] = Evaporation[i];
	}
}

#if !InitDisplayJobDebug
[BurstCompile]
#endif
public struct InitDisplayAirLayerJob : IJobParallelFor {
	public NativeArray<float> DisplayPressure;
	public NativeArray<float3> DisplayPressureGradientForce;
	[ReadOnly] public NativeArray<float> AirLayerElevation;
	[ReadOnly] public NativeArray<float> AirLayerHeight;
	[ReadOnly] public NativeArray<float> AirTemperature;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		DisplayPressure[i] = Atmosphere.GetPressureAtElevation(0, Gravity, AirPressure[i], AirTemperature[i], AirLayerElevation[i] + AirLayerHeight[i] / 2);
		DisplayPressureGradientForce[i] = PressureGradientForce[i];
	}
}
