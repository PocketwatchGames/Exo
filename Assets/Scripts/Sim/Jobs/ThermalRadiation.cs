using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;

#if !ThermalEnergyRadiatedJobDebug
[BurstCompile]
#endif
public struct ThermalEnergyRadiatedJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> ThermalRadiationTransmittedDown;
	public NativeArray<float> WindowRadiationTransmittedUp;
	public NativeArray<float> WindowRadiationTransmittedDown;
	[ReadOnly] public NativeArray<float> TemperatureAbsolute;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		float maxRadiationPercent = 0.01f;

		// radiate half up and half down
		float transmittedUp = math.min(Energy[i] * maxRadiationPercent, Atmosphere.GetRadiationRate(TemperatureAbsolute[i], Emissivity[i]) * SecondsPerTick);
		ThermalRadiationDelta[i] = -2 * transmittedUp;

		float windowTransmittedUp = transmittedUp * PercentRadiationInAtmosphericWindow;
		transmittedUp -= windowTransmittedUp;

		WindowRadiationTransmittedUp[i] = windowTransmittedUp;
		WindowRadiationTransmittedDown[i] = windowTransmittedUp;
		ThermalRadiationTransmittedUp[i] = transmittedUp;
		ThermalRadiationTransmittedDown[i] = transmittedUp;
	}
}

#if !ThermalEnergyRadiatedWaterJobDebug
[BurstCompile]
#endif
public struct ThermalEnergyRadiatedWaterJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> WindowRadiationTransmittedUp;
	[ReadOnly] public NativeArray<float> TemperatureAbsolute;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> SurfaceArea;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		float maxRadiationPercent = 0.01f;

		// radiate half up
		float transmittedUp = math.min(Energy[i] * maxRadiationPercent, Atmosphere.GetRadiationRate(TemperatureAbsolute[i], Emissivity[i]) * SurfaceArea[i] * SecondsPerTick);
		ThermalRadiationDelta[i] = -transmittedUp;

		float windowTransmittedUp = transmittedUp * PercentRadiationInAtmosphericWindow;
		transmittedUp -= windowTransmittedUp;

		WindowRadiationTransmittedUp[i] = windowTransmittedUp;
		ThermalRadiationTransmittedUp[i] = transmittedUp;
	}
}

#if !ThermalEnergyRadiatedAirJobDebug
[BurstCompile]
#endif
public struct ThermalEnergyRadiatedAirJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> ThermalRadiationTransmittedDown;
	public NativeArray<float> WindowRadiationTransmittedUp;
	public NativeArray<float> WindowRadiationTransmittedDown;
	[ReadOnly] public NativeArray<float> TemperaturePotential;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> LayerMiddle;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		float maxRadiationPercent = 0.01f;

		// radiate half up and half down
		float transmittedUp = math.min(Energy[i] * maxRadiationPercent, Atmosphere.GetRadiationRate(Atmosphere.GetAbsoluteTemperature(TemperaturePotential[i], LayerMiddle[i]), Emissivity[i]) * SecondsPerTick);
		ThermalRadiationDelta[i] = -2 * transmittedUp;

		float windowTransmittedUp = transmittedUp * PercentRadiationInAtmosphericWindow;
		transmittedUp -= windowTransmittedUp;

		WindowRadiationTransmittedUp[i] = windowTransmittedUp;
		WindowRadiationTransmittedDown[i] = windowTransmittedUp;
		ThermalRadiationTransmittedUp[i] = transmittedUp;
		ThermalRadiationTransmittedDown[i] = transmittedUp;
	}
}

#if !ThermalEnergyRadiatedConstantEmissivityJobDebug
[BurstCompile]
#endif
public struct ThermalEnergyRadiatedConstantEmissivityJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> ThermalRadiationTransmittedDown;
	public NativeArray<float> WindowRadiationTransmittedUp;
	public NativeArray<float> WindowRadiationTransmittedDown;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> SurfaceArea;
	[ReadOnly] public float Emissivity;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		// radiate half up and half down
		float transmittedUp = math.min(Energy[i] / 2, Atmosphere.GetRadiationRate(Temperature[i], Emissivity) * SurfaceArea[i] * SecondsPerTick);
		ThermalRadiationDelta[i] = -2 * transmittedUp;

		float windowTransmittedUp = transmittedUp * PercentRadiationInAtmosphericWindow;
		transmittedUp -= windowTransmittedUp;

		WindowRadiationTransmittedUp[i] = windowTransmittedUp;
		WindowRadiationTransmittedDown[i] = windowTransmittedUp;
		ThermalRadiationTransmittedUp[i] = transmittedUp;
		ThermalRadiationTransmittedDown[i] = transmittedUp;
	}
}

#if !ThermalEnergyRadiatedTerrainJobDebug
[BurstCompile]
#endif
public struct ThermalEnergyRadiatedTerrainJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		// radiate half up and half down
		float emitted = Atmosphere.GetRadiationRate(Temperature[i], Emissivity[i]) * SecondsPerTick;
		ThermalRadiationDelta[i] = -emitted;

		float emittedOutAtmosphericWindow = emitted * PercentRadiationInAtmosphericWindow;
		emitted -= emittedOutAtmosphericWindow;
		WindowRadiationTransmitted[i] = emittedOutAtmosphericWindow;
		ThermalRadiationTransmitted[i] = emitted;
	}
}

#if !ThermalEnergyAbsorbedAirJobDebug
[BurstCompile]
#endif
public struct ThermalEnergyAbsorbedAirJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeArray<float> WindowRadiationIncoming;
	[ReadOnly] public NativeArray<float> ThermalRadiationIncoming;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<ThermalAbsorptivity> AbsorptivityThermal;
	[ReadOnly] public bool FromTop;
	public void Execute(int i)
	{
		float cloudElevation = CloudElevation[i];
		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];
		float beforeCloud = math.min(1, (cloudElevation - layerElevation) / layerHeight);
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

		float transmitting = ThermalRadiationTransmitted[i];

		float incoming = ThermalRadiationIncoming[i];
		float absorbedBeforeCloud = incoming * absorptivityBefore;
		incoming -= absorbedBeforeCloud;
		thermalRadiationDelta += absorbedBeforeCloud;

		float incomingAbsorbedByCloud = incoming * AbsorptivityThermal[i].AbsorptivityCloud;
		float transmittingAbsorbedByCloud = beforeCloud * transmitting * AbsorptivityThermal[i].AbsorptivityCloud;
		thermalRadiationDelta += incomingAbsorbedByCloud + transmittingAbsorbedByCloud;

		incoming -= incomingAbsorbedByCloud;
		transmitting -= transmittingAbsorbedByCloud;

		float absorbedAfterCloud = incoming * absorptivityAfter;
		incoming -= absorbedAfterCloud;
		thermalRadiationDelta += absorbedAfterCloud;

		ThermalRadiationDelta[i] = thermalRadiationDelta;
		ThermalRadiationTransmitted[i] = transmitting + incoming;
		WindowRadiationTransmitted[i] += WindowRadiationIncoming[i];
	}
}

#if !ThermalEnergyAbsorbedJob
[BurstCompile]
#endif
public struct ThermalEnergyAbsorbedJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeArray<float> WindowRadiationIncoming;
	[ReadOnly] public NativeArray<float> ThermalRadiationIncoming;
	public void Execute(int i)
	{
		float absorptivity = 1;

		float windowIncoming = WindowRadiationIncoming[i];
		float windowAbsorbed = windowIncoming * absorptivity;
		WindowRadiationTransmitted[i] += windowIncoming - windowAbsorbed;

		float incoming = ThermalRadiationIncoming[i];
		float absorbed = incoming * absorptivity;
		ThermalRadiationTransmitted[i] += incoming - absorbed;

		ThermalRadiationDelta[i] += absorbed + windowAbsorbed;
	}
}

#if !ThermalEnergyAbsorbedPartialCoverageJobDebug
[BurstCompile]
#endif
public struct ThermalEnergyAbsorbedPartialCoverageJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeArray<float> WindowRadiationIncoming;
	[ReadOnly] public NativeArray<float> ThermalRadiationIncoming;
	[ReadOnly] public NativeArray<float> Coverage;
	public void Execute(int i)
	{
		float absorptivity = Coverage[i];

		float windowRadiationIncoming = WindowRadiationIncoming[i];
		float atmosphericWindowAbsorbed = windowRadiationIncoming * absorptivity;
		WindowRadiationTransmitted[i] += windowRadiationIncoming - atmosphericWindowAbsorbed;

		float incoming = ThermalRadiationIncoming[i];
		float absorbed = incoming * absorptivity;
		ThermalRadiationTransmitted[i] += incoming - absorbed;
		ThermalRadiationDelta[i] += absorbed + atmosphericWindowAbsorbed;
	}
}

#if !ThermalEnergyAbsorbedTerrainJobDebug
[BurstCompile]
#endif
public struct ThermalEnergyAbsorbedTerrainJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationAbsorbed;
	[ReadOnly] public NativeArray<float> WindowRadiationIncoming;
	[ReadOnly] public NativeArray<float> ThermalRadiationIncoming;
	public void Execute(int i)
	{
		ThermalRadiationAbsorbed[i] += ThermalRadiationIncoming[i] + WindowRadiationIncoming[i];
	}
}

