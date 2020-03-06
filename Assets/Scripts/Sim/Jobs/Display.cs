using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;



#if !UpdateDisplayJobDebug
[BurstCompile]
#endif
public struct UpdateDisplayJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbedSurface;
	public NativeArray<float> DisplayPrecipitation;
	public NativeArray<float> DisplayEvaporation;
	[ReadOnly] public NativeArray<float> SolarRadiationInTerrain;
	[ReadOnly] public NativeArray<float> SolarRadiationInIce;
	[ReadOnly] public NativeArray<float> SolarRadiationInWaterSurface;
	[ReadOnly] public NativeArray<float> Precipitation;
	[ReadOnly] public NativeArray<float> Evaporation;
	public void Execute(int i)
	{
		SolarRadiationAbsorbedSurface[i] = SolarRadiationInTerrain[i] + SolarRadiationInIce[i] + SolarRadiationInWaterSurface[i];
		DisplayPrecipitation[i] = Precipitation[i];
		DisplayEvaporation[i] = Evaporation[i];
	}
}

#if !InitDisplayJobDebug
[BurstCompile]
#endif
public struct InitDisplayAirLayerJob : IJobParallelFor {
	public NativeArray<float> DisplayPressure;
	public NativeArray<float3> DisplayPressureGradientForce;
	public NativeArray<float> DisplayCondensationGround;
	public NativeArray<float> DisplayCondensationCloud;
	[ReadOnly] public NativeArray<float> AirTemperaturePotential;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CondensationCloud;
	[ReadOnly] public NativeArray<float> CondensationGround;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		DisplayPressure[i] = Atmosphere.GetPressureAtElevation(0, Gravity, AirPressure[i], AirTemperaturePotential[i], LayerElevation[i] + LayerHeight[i] / 2);
		DisplayPressureGradientForce[i] = PressureGradientForce[i];
		DisplayCondensationGround[i] += CondensationCloud[i];
		DisplayCondensationCloud[i] += CondensationGround[i];
	}
}
