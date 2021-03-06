﻿//#define DISABLE_CONDUCTION

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;



[BurstCompile]
public struct ConductionJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> SurfaceArea;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
#if !DISABLE_CONDUCTION
		// TODO: this can conduct heat past a point of equilibrium
		EnergyDelta[i] = (TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * SurfaceArea[i];
#endif
	}
}

[BurstCompile]
public struct ConductionAJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> SurfaceArea;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
#if !DISABLE_CONDUCTION
		// TODO: this can conduct heat past a point of equilibrium
		EnergyDelta[i] = math.max((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * SurfaceArea[i], -EnergyA[i]);
#endif
	}
}
[BurstCompile]
public struct ConductionBJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeSlice<float> TemperatureB;
	[ReadOnly] public NativeSlice<float> EnergyB;
	[ReadOnly] public NativeArray<float> SurfaceArea;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
#if !DISABLE_CONDUCTION
		// TODO: this can conduct heat past a point of equilibrium
		EnergyDelta[i] = math.min((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * SurfaceArea[i], EnergyB[i]);
#endif
	}
}

[BurstCompile]
public struct ConductionABJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeSlice<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeSlice<float> EnergyB;
	[ReadOnly] public NativeArray<float> SurfaceArea;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
#if !DISABLE_CONDUCTION
		// TODO: this can conduct heat past a point of equilibrium
		EnergyDelta[i] = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * SurfaceArea[i], -EnergyA[i], EnergyB[i]);
#endif
	}
}


[BurstCompile]
public struct ConductionWaterBottomAJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeSlice<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeSlice<float> EnergyA;
	[ReadOnly] public NativeArray<float> SurfaceArea;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public int LayerIndex;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{
		float coverage = Coverage[i + LayerIndex * Count];
		float coverageBelow = Coverage[i + (LayerIndex - 1) * Count];

		// TODO: this can conduct heat past a point of equilibrium
		EnergyDelta[i] += math.max((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * SurfaceArea[i] * coverage * (1.0f - coverageBelow), -EnergyA[i]);
	}
}

[BurstCompile]
public struct ConductionWaterBottomABJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	public NativeArray<float> EnergyDeltaTotal;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> EnergyB;
	[ReadOnly] public NativeArray<float> SurfaceArea;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public NativeArray<float> CoverageBelow;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
#if !DISABLE_CONDUCTION
		// TODO: this can conduct heat past a point of equilibrium
		float delta = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * SurfaceArea[i] * Coverage[i] * (1.0f - CoverageBelow[i]), -EnergyA[i], EnergyB[i]);
		EnergyDelta[i] = delta;
		EnergyDeltaTotal[i] += delta;
#endif
	}
}

