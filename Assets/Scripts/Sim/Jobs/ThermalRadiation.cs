
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

[BurstCompile]
public struct ThermalEnergyRadiatedJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeArray<float> TemperatureAbsolute;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		float maxRadiationPercent = 0.01f;

		// radiate half up and half down
		ThermalRadiationEmitted[i] = math.min(Energy[i] * maxRadiationPercent, Atmosphere.GetRadiationRate(TemperatureAbsolute[i], Emissivity[i]) * SecondsPerTick);
	}
}

[BurstCompile]
public struct ThermalEnergyRadiatedWaterJob : IJobParallelFor {
	public NativeSlice<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeSlice<float> TemperatureAbsolute;
	[ReadOnly] public NativeSlice<float> Energy;
	[ReadOnly] public NativeSlice<float> Emissivity;
	[ReadOnly] public NativeSlice<float> SurfaceArea;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		float maxRadiationPercent = 0.01f;

		// radiate half up
		ThermalRadiationEmitted[i] = math.min(Energy[i] * maxRadiationPercent, Atmosphere.GetRadiationRate(TemperatureAbsolute[i], Emissivity[i]) * SurfaceArea[i] * SecondsPerTick);
	}
}

[BurstCompile]
public struct ThermalEnergyRadiatedAirJob : IJobParallelFor {
	public NativeSlice<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeSlice<float> TemperaturePotential;
	[ReadOnly] public NativeSlice<float> Energy;
	[ReadOnly] public NativeSlice<float> Emissivity;
	[ReadOnly] public NativeSlice<float> LayerMiddle;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float maxRadiationPercent = 0.01f;

		// radiate half up and half down
		ThermalRadiationEmitted[i] = math.min(Energy[i] * maxRadiationPercent, Atmosphere.GetRadiationRate(Atmosphere.GetAbsoluteTemperature(TemperaturePotential[i], LayerMiddle[i]), Emissivity[i]) * SecondsPerTick);
	}
}

[BurstCompile]
public struct ThermalEnergyRadiatedConstantEmissivityJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> SurfaceArea;
	[ReadOnly] public float Emissivity;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		// radiate half up and half down
		ThermalRadiationEmitted[i] = math.min(Energy[i] / 2, Atmosphere.GetRadiationRate(Temperature[i], Emissivity) * SurfaceArea[i] * SecondsPerTick);
	}
}

[BurstCompile]
public struct ThermalEnergyRadiatedTerrainJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		// radiate half up and half down
		ThermalRadiationEmitted[i] = Atmosphere.GetRadiationRate(Temperature[i], Emissivity[i]) * SecondsPerTick;

		// TODO: should we radiate more given that the surface area of flora can be greater than 1?
	}
}

[BurstCompile]
public struct ThermalEnergyAbsorbedAirJob : IJobParallelFor {
	public NativeSlice<float> ThermalRadiationDelta;
	public NativeSlice<float> ThermalRadiationTransmitted;
	public NativeSlice<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeSlice<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeSlice<float> LayerElevation;
	[ReadOnly] public NativeSlice<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeSlice<ThermalAbsorptivity> AbsorptivityThermal;
	[ReadOnly] public bool FromTop;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		float emitted = ThermalRadiationEmitted[i];

		float cloudElevation = CloudElevation[i];
		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];
		float beforeCloud = math.saturate((cloudElevation - layerElevation) / layerHeight);
		float thermalRadiationDelta = ThermalRadiationDelta[i];

		float absorptivityBefore;
		float absorptivityAfter;
		if (!FromTop)
		{
			beforeCloud = 1.0f - beforeCloud;
			absorptivityBefore = AbsorptivityThermal[i].AbsorptivityAirBelow;
			absorptivityAfter = AbsorptivityThermal[i].AbsorptivityAirAbove;
		} else
		{
			absorptivityBefore = AbsorptivityThermal[i].AbsorptivityAirAbove;
			absorptivityAfter = AbsorptivityThermal[i].AbsorptivityAirBelow;
		}

		float emitting = ThermalRadiationEmitted[i];
		float emittingWindowRadiation = emitting * PercentRadiationInAtmosphericWindow;
		float emittingThermalRadiation = emitting * (1.0f - PercentRadiationInAtmosphericWindow);
		thermalRadiationDelta -= emitting;

		float incoming = ThermalRadiationTransmitted[i];
		float absorbedBeforeCloud = incoming * absorptivityBefore;
		incoming -= absorbedBeforeCloud;
		thermalRadiationDelta += absorbedBeforeCloud;

		float incomingAbsorbedByCloud = incoming * AbsorptivityThermal[i].AbsorptivityCloud;
		float emittingAbsorbedByCloud = beforeCloud * emittingThermalRadiation * AbsorptivityThermal[i].AbsorptivityCloud;
		thermalRadiationDelta += incomingAbsorbedByCloud + emittingAbsorbedByCloud;

		incoming -= incomingAbsorbedByCloud;
		emittingThermalRadiation -= emittingAbsorbedByCloud;

		float absorbedAfterCloud = incoming * absorptivityAfter;
		incoming -= absorbedAfterCloud;
		thermalRadiationDelta += absorbedAfterCloud;

		ThermalRadiationDelta[i] = thermalRadiationDelta;
		ThermalRadiationTransmitted[i] = emittingThermalRadiation + incoming;
		WindowRadiationTransmitted[i] += emittingWindowRadiation;
	}
}


[BurstCompile]
public struct ThermalEnergyRadiatedUpTerrainJob : IJobParallelFor {
	public NativeSlice<float> ThermalRadiationDelta;
	public NativeSlice<float> ThermalRadiationTransmitted;
	public NativeSlice<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeSlice<float> ThermalRadiationEmitted;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		float emitted = ThermalRadiationEmitted[i];
		ThermalRadiationDelta[i] = -emitted;
		ThermalRadiationTransmitted[i] = emitted * (1.0f - PercentRadiationInAtmosphericWindow);
		WindowRadiationTransmitted[i] = emitted * PercentRadiationInAtmosphericWindow;
	}
}

[BurstCompile]
public struct ThermalEnergyAbsorbedDownTerrainJob : IJobParallelFor {
	public NativeSlice<float> ThermalRadiationDelta;
	[ReadOnly] public NativeSlice<float> WindowRadiationIncoming;
	[ReadOnly] public NativeSlice<float> ThermalRadiationIncoming;
	public void Execute(int i)
	{
		ThermalRadiationDelta[i] += ThermalRadiationIncoming[i] + WindowRadiationIncoming[i];
	}
}

[BurstCompile]
public struct ThermalEnergyAbsorbedDownPartialCoverageJob : IJobParallelFor {
	public NativeSlice<float> ThermalRadiationDelta;
	public NativeSlice<float> ThermalRadiationTransmitted;
	public NativeSlice<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeSlice<float> Coverage;
	public void Execute(int i)
	{
		float atmosphericWindowAbsorbed = WindowRadiationTransmitted[i] * Coverage[i];
		WindowRadiationTransmitted[i] -= atmosphericWindowAbsorbed;

		float absorbed = ThermalRadiationTransmitted[i] * Coverage[i];
		ThermalRadiationTransmitted[i] -= absorbed;

		ThermalRadiationDelta[i] += absorbed + atmosphericWindowAbsorbed;
	}
}

[BurstCompile]
public struct ThermalEnergyAbsorbedUpPartialCoverageJob : IJobParallelFor {
	public NativeSlice<float> ThermalRadiationDelta;
	public NativeSlice<float> ThermalRadiationTransmitted;
	public NativeSlice<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeSlice<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeSlice<float> Coverage;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		float absorptivity = Coverage[i];

		float emitted = ThermalRadiationEmitted[i];

		float atmosphericWindowAbsorbed = WindowRadiationTransmitted[i] * absorptivity;
		WindowRadiationTransmitted[i] = WindowRadiationTransmitted[i] - atmosphericWindowAbsorbed + emitted * PercentRadiationInAtmosphericWindow;

		float incoming = ThermalRadiationTransmitted[i];
		float absorbed = incoming * absorptivity;
		ThermalRadiationTransmitted[i] = incoming - absorbed + emitted * (1.0f - PercentRadiationInAtmosphericWindow);
		ThermalRadiationDelta[i] += absorbed + atmosphericWindowAbsorbed - emitted;
	}
}

