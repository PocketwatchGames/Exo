//#define DISABLE_CLOUD_DISSAPATION
#define DISABLE_CLOUD_ELEVATION
//#define DISABLE_RAINFALL
//#define DISABLE_CONDENSATION
#define DISABLE_VERTICAL_AIR_MOVEMENT
#define DISABLE_CLOUD_ADVECTION
#define DISABLE_AIR_ADVECTION
//#define DISABLE_EVAPORATION
//#define DISABLE_FREEZE_TOP
//#define DISABLE_FREEZE_BOTTOM
//#define DISABLE_MELTING_TOP
//#define DISABLE_MELTING_BOTTOM

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

public struct DiffusionAir {
	public float Temperature;
	public float Humidity;
	public float2 Velocity;
}
public struct DiffusionAirVertical {
	public float Temperature;
	public float Humidity;
}
public struct DiffusionCloud {
	public float Temperature;
	public float Mass;
	public float DropletMass;
	public float Elevation;
	public float2 Velocity;
}
public struct DiffusionWater {
	public float Temperature;
	public float Salinity;
	public float2 Velocity;
}
public struct WaterSaltMass {
	public float WaterMass;
	public float SaltMass;
}

#if !SolarRadiationJobDebug
[BurstCompile]
#endif
public struct SolarRadiationJob : IJobParallelFor {
	public NativeArray<float> SolarRadiation;
	public NativeArray<float> DisplaySolarRadiation;
	public NativeArray<float> WaterSlopeAlbedo;
	[ReadOnly] public NativeArray<float3> SphericalPosition;
	[ReadOnly] public float3 SunToPlanetDir;
	[ReadOnly] public quaternion PlanetRotation;
	[ReadOnly] public float IncomingSolarRadiation;
	public void Execute(int i)
	{
		float sunDotSurface = math.max(0, math.dot(SunToPlanetDir, math.rotate(PlanetRotation, -SphericalPosition[i])));
		WaterSlopeAlbedo[i] = math.pow(1.0f - math.max(0, sunDotSurface), 9);
		float r = IncomingSolarRadiation * sunDotSurface;
		SolarRadiation[i] = r;
		DisplaySolarRadiation[i] = r;
	}
}

#if !SolarRadiationAbsorbedAirJobDebug
[BurstCompile]
#endif
public struct SolarRadiationAbsorbedAirJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationAbsorbedCloud;
	public NativeArray<float> SolarRadiationIncoming;
	public NativeArray<float> SolarRadiationReflected;
	public NativeArray<float> SolarRadiationReflectedCloud;
	[ReadOnly] public float SolarReflectivityAir;
	[ReadOnly] public float SolarAbsorptivityAir;
	[ReadOnly] public float SolarAbsorptivityWaterVapor;
	[ReadOnly] public float SolarAbsorptivityCloud;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CloudCoverage;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudTemperature;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> WaterSlopeAlbedo;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public int LayerIndex;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		float airMass = AirMass[i];
		float waterVaporMass = VaporMass[i];
		float incomingRadiation = SolarRadiationIncoming[i];


		float cloudMass = CloudMass[i];
		float cloudElevation = CloudElevation[i];
		float layerElevation = LayerElevation[LayerIndex];
		float layerHeight = LayerHeight[LayerIndex];
		bool isCloudLayer = cloudElevation >= layerElevation && cloudElevation < layerElevation + layerHeight;
		float beforeCloud = math.min(1, (cloudElevation - layerElevation) / layerHeight);
		float afterCloud = 1.0f - beforeCloud;

		float reflectivity = math.min(1, SolarReflectivityAir * (airMass + waterVaporMass));
		float absorptivity = SolarAbsorptivityAir * airMass + SolarAbsorptivityWaterVapor * waterVaporMass;
		float absorbed = 0;
		float energyReflectedAtmosphere = 0;
		float absorbedByCloudsIncoming = 0;
		float energyReflectedClouds = 0;

		if (beforeCloud > 0)
		{
			energyReflectedAtmosphere = incomingRadiation * math.max(0, 1.0f - 1.0f / math.exp10(reflectivity * beforeCloud));
			incomingRadiation -= energyReflectedAtmosphere;
			absorbed += incomingRadiation * math.max(0, 1.0f - 1.0f / math.exp10(absorptivity * beforeCloud));
			incomingRadiation -= absorbed;
		}


		if (isCloudLayer)
		{
			float cloudIceContent = math.saturate((CloudTemperature[i] - worldData.minCloudFreezingTemperature) / (worldData.maxCloudFreezingTemperature - worldData.minCloudFreezingTemperature));
			float cloudTemperatureAlbedo = WorldData.AlbedoIce + (WorldData.AlbedoWater - WorldData.AlbedoIce) * cloudIceContent;
			float rainDropSizeAlbedo = math.saturate(1.0f - CloudDropletMass[i] / CloudMass[i]) * (worldData.rainDropSizeAlbedoMax - worldData.rainDropSizeAlbedoMin) + worldData.rainDropSizeAlbedoMin;
			float cloudAlbedo = math.min(1.0f, worldData.SolarReflectivityCloud * cloudTemperatureAlbedo * CloudMass[i] * rainDropSizeAlbedo / math.max(worldData.maxCloudSlopeAlbedo, 1.0f - WaterSlopeAlbedo[i]));

			energyReflectedClouds = incomingRadiation * math.max(0, 1.0f - 1.0f / math.exp10(cloudAlbedo * cloudMass));
			incomingRadiation -= energyReflectedClouds;
			absorbedByCloudsIncoming = incomingRadiation * math.max(0, 1.0f - 1.0f / math.exp10(SolarAbsorptivityCloud * cloudMass));
			incomingRadiation -= absorbedByCloudsIncoming;

			SolarRadiationReflectedCloud[i] = energyReflectedClouds;
			//			SolarRadiationAbsorbedCloud[i] = absorbedByCloudsIncoming;
			absorbed += absorbedByCloudsIncoming;
		}

		// below cloud
		if (afterCloud > 0)
		{
			float reflectedBelow = incomingRadiation * math.max(0, 1.0f - 1.0f / math.exp10(reflectivity * afterCloud));
			incomingRadiation -= reflectedBelow;
			energyReflectedAtmosphere += reflectedBelow;
			float absorbedBelow = incomingRadiation * math.max(0, 1.0f - 1.0f / math.exp10(absorptivity * afterCloud));
			absorbed += absorbedBelow;
			incomingRadiation -= absorbedBelow;
		}

		SolarRadiationAbsorbed[i] = absorbed;
		SolarRadiationIncoming[i] = incomingRadiation;
		SolarRadiationReflected[i] = energyReflectedAtmosphere;
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
		float reflected = incoming * AlbedoIce * iceCoverage;
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
		float reflected = incoming * WaterSlopeAlbedo[i] * waterCoverage;
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
		float reflected = incoming * vegetationCoverage * WorldData.AlbedoFoliage + (1.0f - vegetationCoverage) * soilReflectivity;
		incoming -= reflected;

		SolarRadiationReflected[i] = reflected;
		SolarRadiationAbsorbed[i] = incoming;
	}
}

#if !EmissivityAirJobDebug
[BurstCompile]
#endif
public struct EmissivityAirJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	public void Execute(int i)
	{
		Emissivity[i] = (AirMass[i] * WorldData.EmissivityAir +	VaporMass[i] * WorldData.EmissivityWater) / (AirMass[i] + VaporMass[i]);
	}
}

#if !EmissivityWaterJobDebug
[BurstCompile]
#endif
public struct EmissivityWaterJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> WaterSaltMass;
	[ReadOnly] public NativeArray<float> WaterMass;
	public void Execute(int i)
	{
		Emissivity[i] = math.saturate((WorldData.EmissivityWater * WaterMass[i] + WorldData.EmissivitySalt * WaterSaltMass[i]) / (WaterMass[i] + WaterSaltMass[i]));
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

#if !ThermalEnergyRadiatedJobDebug
[BurstCompile]
#endif
public struct ThermalEnergyRadiatedJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> ThermalRadiationTransmittedDown;
	public NativeArray<float> WindowRadiationTransmittedUp;
	public NativeArray<float> WindowRadiationTransmittedDown;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		float maxRadiationPercent = 0.01f;

		// radiate half up and half down
		float transmittedUp = math.min(Energy[i] * maxRadiationPercent, Atmosphere.GetRadiationRate(Temperature[i], Emissivity[i]) * SecondsPerTick);
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
		float transmittedUp = math.min(Energy[i] / 2,  Atmosphere.GetRadiationRate(Temperature[i], Emissivity) * SurfaceArea[i] * SecondsPerTick);
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
	public NativeArray<float> ThermalRadiationDeltaCloud;
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeArray<float> WindowRadiationIncoming;
	[ReadOnly] public NativeArray<float> ThermalRadiationIncoming;
	[ReadOnly] public NativeArray<float> ThermalRadiationTransmittedCloud;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public float AirAbsorptivity;
	[ReadOnly] public float VaporAbsorptivity;
	[ReadOnly] public float WaterAbsorptivity;
	[ReadOnly] public float CarbonAbsorptivity;
	[ReadOnly] public float CarbonDioxide;
	[ReadOnly] public int LayerIndex;
	[ReadOnly] public bool FromTop;
	public void Execute(int i)
	{
		WindowRadiationTransmitted[i] += WindowRadiationIncoming[i];

		float absorptivity = AirMass[i] * ((1.0f - CarbonDioxide) * AirAbsorptivity + CarbonDioxide * CarbonAbsorptivity) + VaporMass[i] * VaporAbsorptivity;

		float cloudMass = CloudMass[i];
		float cloudElevation = CloudElevation[i];
		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];
		bool isCloudLayer = cloudElevation >= layerElevation && cloudElevation < layerElevation + layerHeight;
		float beforeCloud = math.min(1, (cloudElevation - layerElevation) / layerHeight);
		if (!FromTop)
		{
			beforeCloud = 1.0f - beforeCloud;
		}
		float afterCloud = 1.0f - beforeCloud;

		float transmitting = ThermalRadiationTransmitted[i];

		float incoming = ThermalRadiationIncoming[i];
		if (beforeCloud > 0)
		{
			float absorbedBeforeCloud = incoming * math.max(0, 1.0f - 1.0f / math.exp10(absorptivity * beforeCloud));
			incoming -= absorbedBeforeCloud;
			ThermalRadiationDelta[i] += absorbedBeforeCloud;
		}

		if (cloudMass > 0 && isCloudLayer)
		{
			float absorbance = math.max(0, 1.0f - 1.0f / math.exp10(WaterAbsorptivity * cloudMass));
			float incomingAbsorbedByCloud = incoming * absorbance;
			float transmittingAbsorbedByCloud = beforeCloud * transmitting * absorbance;
			//	ThermalRadiationDeltaCloud[i] += incomingAbsorbedByCloud + transmittingAbsorbedByCloud;
			ThermalRadiationDelta[i] += incomingAbsorbedByCloud + transmittingAbsorbedByCloud;


			incoming += ThermalRadiationTransmittedCloud[i] - incomingAbsorbedByCloud;
			transmittingAbsorbedByCloud -= transmittingAbsorbedByCloud;
		}

		if (afterCloud > 0)
		{
			float absorbedAfterCloud = incoming * math.max(0, 1.0f - 1.0f / math.exp10(absorptivity * afterCloud));
			ThermalRadiationDelta[i] += absorbedAfterCloud;
			incoming -= absorbedAfterCloud;
		}
		ThermalRadiationTransmitted[i] = transmitting + incoming;
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

#if !DiffusionAirJobDebug
[BurstCompile]
#endif
public struct DiffusionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public float DiffusionCoefficient;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Humidity;
	[ReadOnly] public NativeArray<float2> Wind;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		float gradientTemperature = 0;
		float gradientWaterVapor = 0;
		float2 gradientVelocity = float2.zero;
		int neighborCount = 0;
		//TODO: account for different size air columns, similar to water
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				gradientWaterVapor += Humidity[n];
				gradientVelocity += Wind[n];
				gradientTemperature += Temperature[n];
				neighborCount++;
			}
		}
		gradientTemperature -= Temperature[i] * neighborCount;
		gradientWaterVapor -= Humidity[i] * neighborCount;
		gradientVelocity -= Wind[i] * neighborCount;

		Delta[i] = new DiffusionAir()
		{
			Temperature = gradientTemperature * DiffusionCoefficient,
			Humidity = gradientWaterVapor * DiffusionCoefficient,
			Velocity = gradientVelocity * DiffusionCoefficient
		};
	}
}

#if !AirVerticalMovementUpJobDebug
[BurstCompile]
#endif
public struct AirVerticalMovementUpJob : IJobParallelFor {
	public NativeArray<DiffusionAirVertical> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Humidity;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpHumidity;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public float DiffusionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float MaxVerticalMovement;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		float gradientTemperature = 0;
		float gradientWaterVapor = 0;

		float heightDiff = (UpLayerElevation[i] + UpLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);

		//TODO: account for different size air columns, similar to water


		float atmosphereMass = Humidity[i] + AirMass[i];
		float atmosphereMassUp = UpHumidity[i] + UpAirMass[i];
		float absoluteHumidity = Humidity[i] / atmosphereMass;
		float absoluteHumidityUp = UpHumidity[i] / atmosphereMassUp;
		float gradientHumidity = absoluteHumidityUp - absoluteHumidity;
		gradientWaterVapor = gradientHumidity * math.min(atmosphereMass, atmosphereMassUp);
		float potentialTemperatureUp = UpTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
		gradientTemperature = potentialTemperatureUp - Temperature[i];
		float moveToNeutralBouyancy = (UpTemperature[i] - Temperature[i]) / WorldData.TemperatureLapseRate - heightDiff;
		float vertMovement = math.min(MaxVerticalMovement, math.clamp(moveToNeutralBouyancy + DiffusionCoefficient, 0, 1));

		//Debug.Log("I: " + i + " G: " + gradientWaterVapor + " HA: " + Humidity[i] + " HB: " + UpHumidity[i] + " AA: " + atmosphereMass + " AB: " + atmosphereMassUp);

		float bouyancyForce = math.max(0, Gravity * (Temperature[i] / potentialTemperatureUp - 1));
		Delta[i] = new DiffusionAirVertical()
		{
			Temperature = gradientTemperature * vertMovement,
			Humidity = gradientWaterVapor * vertMovement,
		};
	}
}

#if !AirVerticalMovementDownJobDebug
[BurstCompile]
#endif
public struct AirVerticalMovementDownJob : IJobParallelFor {
	public NativeArray<DiffusionAirVertical> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Humidity;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpHumidity;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public float DiffusionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float MaxVerticalMovement;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		float gradientTemperature = 0;
		float gradientWaterVapor = 0;

		float heightDiff = (UpLayerElevation[i] + UpLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);

		//TODO: account for different size air columns, similar to water


		float atmosphereMass = Humidity[i] + AirMass[i];
		float atmosphereMassUp = UpHumidity[i] + UpAirMass[i];
		float absoluteHumidity = Humidity[i] / atmosphereMass;
		float absoluteHumidityUp = UpHumidity[i] / atmosphereMassUp;
		float gradientHumidity = absoluteHumidityUp - absoluteHumidity;
		gradientWaterVapor = gradientHumidity * math.min(atmosphereMass, atmosphereMassUp);
		float potentialTemperatureUp = UpTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
		gradientTemperature = potentialTemperatureUp - Temperature[i];
		float moveToNeutralBouyancy = (UpTemperature[i] - Temperature[i]) / WorldData.TemperatureLapseRate - heightDiff;
		float vertMovement = math.min(MaxVerticalMovement, math.clamp(moveToNeutralBouyancy - DiffusionCoefficient, -1, 0));

		Delta[i] = new DiffusionAirVertical()
		{
			Temperature = gradientTemperature * vertMovement,
			Humidity = gradientWaterVapor * vertMovement,
		};
	}
}

#if !DiffusionCloudJobDebug
[BurstCompile]
#endif
public struct DiffusionCloudJob : IJobParallelFor {
	public NativeArray<DiffusionCloud> Delta;
	[ReadOnly] public float DiffusionCoefficient;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Elevation;
	[ReadOnly] public NativeArray<float> DropletSize;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		float gradientMass = 0;
		float gradientTemperature = 0;
		float gradientElevation = 0;
		float gradientDropletSize = 0;
		float2 gradientVelocity = float2.zero;
		int neighborCount = 0;
		float neighborMassTotal = 0;
		float mass = Mass[i];
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float nMass = Mass[n];
				gradientMass += nMass;
				gradientTemperature += Temperature[n] * nMass;
				gradientElevation += Elevation[n] * nMass;
				gradientDropletSize += DropletSize[n] * nMass;
				gradientVelocity += Velocity[n] * nMass;
				neighborCount++;
			}
		}
		if (neighborMassTotal > 0)
		{
			float massDivisor = 1.0f / neighborMassTotal;
			gradientTemperature = gradientTemperature * massDivisor - Temperature[i];
			gradientElevation = gradientElevation * massDivisor - Elevation[i];
			gradientDropletSize = gradientDropletSize * massDivisor - DropletSize[i];
			gradientVelocity = gradientVelocity * massDivisor - Velocity[i];
		}
		gradientMass -= mass * neighborCount;

		Delta[i] = new DiffusionCloud()
		{
			Mass = gradientMass * DiffusionCoefficient,
			Temperature = gradientTemperature * DiffusionCoefficient,
			Elevation = gradientElevation * DiffusionCoefficient,
			DropletMass = gradientDropletSize * DiffusionCoefficient,
			Velocity = gradientVelocity * DiffusionCoefficient,
		};
	}
}

#if !DiffusionWaterJobDebug
[BurstCompile]
#endif
public struct DiffusionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public float DiffusionCoefficient;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float2> Current;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		float gradientSalinity = 0;
		float gradientTemperature = 0;
		float2 gradientVelocity = float2.zero;
		float neighborMassToDiffuse = 0;
		float waterMass = WaterMass[i];
		if (waterMass > 0) {
			for (int j = 0; j < 6; j++)
			{
				int n = Neighbors[i * 6 + j];
				if (n >= 0)
				{
					float nMass = WaterMass[n];
					if (nMass > 0)
					{
						// TODO: probably want to deal with "salinity" here not salt mass
						float diffuseMass = math.min(nMass, waterMass) / nMass;
						float diffuse = math.min(nMass, waterMass) / waterMass;
						gradientSalinity += Salt[n] * diffuseMass;
						gradientTemperature += Temperature[n] * diffuse;
						gradientVelocity += Current[n] * diffuse;

						neighborMassToDiffuse += diffuse;
					}
				}
			}
			gradientSalinity -= Salt[i] * neighborMassToDiffuse;
			gradientTemperature -= Temperature[i] * neighborMassToDiffuse;
			gradientVelocity -= Current[i] * neighborMassToDiffuse;
		}

		Delta[i] = new DiffusionWater()
		{
			Temperature = gradientTemperature * DiffusionCoefficient,
			Salinity = gradientSalinity * DiffusionCoefficient,
			Velocity = gradientVelocity * DiffusionCoefficient
		};

	}
}

#if !AdvectionAirJobDebug
[BurstCompile]
#endif
public struct AdvectionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float2> Wind;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float2> Coords;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float InverseCoordDiff;
	public void Execute(int i)
	{
#if !DISABLE_AIR_ADVECTION
		float gradientTemperature = 0;
		float gradientWaterVapor = 0;
		float2 gradientVelocity = float2.zero;
		float2 velocity = Wind[i];
		float2 coord = Coords[i];
		//TODO: account for different size air columns, similar to water
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float2 coordDiff = Coords[n] - coord;
				float2 diff = math.normalize(math.float2(Utils.WrapAngle(coordDiff.x), Utils.WrapAngle(coordDiff.y)));
				// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
				float velDotDir = math.max(0, math.dot(Wind[n], diff)) * InverseCellDiameter;
				if (velDotDir < 0)
				{
					gradientWaterVapor -= velDotDir * Vapor[n];
					gradientVelocity -= velDotDir * Wind[n];
					gradientTemperature -= velDotDir * Temperature[n];
				}
			}
		}
		float outgoing = math.length(velocity) * InverseCellDiameter;
		gradientTemperature -= Temperature[i] * outgoing;
		gradientWaterVapor -= Vapor[i] * outgoing;
		gradientVelocity -= Wind[i] * outgoing;

		Delta[i] = new DiffusionAir()
		{
			Temperature = gradientTemperature,
			Humidity = gradientWaterVapor,
			Velocity = gradientVelocity
		};
#endif
	}
}

#if !AdvectionCloudJob
[BurstCompile]
#endif
public struct AdvectionCloudJob : IJobParallelFor {
	public NativeArray<DiffusionCloud> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Elevation;
	[ReadOnly] public NativeArray<float> DropletMass;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float2> Coords;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float InverseCoordDiff;
	public void Execute(int i)
	{
#if !DISABLE_CLOUD_ADVECTION
		float gradientTemperature = 0;
		float gradientMass = 0;
		float gradientElevation = 0;
		float gradientDropletMass = 0;
		float2 gradientVelocity = float2.zero;
		float2 velocity = Velocity[i];
		float2 coord = Coords[i];
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float2 coordDiff = Coords[n] - coord;
				float2 diff = math.normalize(math.float2(Utils.WrapAngle(coordDiff.x), Utils.WrapAngle(coordDiff.y)));
				// TODO: Dot product doesn't actually work given the hexagonal nature of the grid
				float velDotDir = math.max(0, math.dot(Velocity[n], diff)) * InverseCellDiameter;
				if (velDotDir < 0)
				{
					// TODO: account for differing masses having smaller or larger effects on cloud properties
					gradientTemperature -= velDotDir * Temperature[n];
					gradientMass -= velDotDir * Mass[n];
					gradientElevation -= velDotDir * Elevation[n];
					gradientDropletMass -= velDotDir * DropletMass[n];
					gradientVelocity -= velDotDir * Velocity[n];
				}
			}
		}
		float outgoing = math.length(velocity) * InverseCellDiameter;
		gradientTemperature -= Temperature[i] * outgoing;
		gradientMass -= Mass[i] * outgoing;
		gradientElevation -= Elevation[i] * outgoing;
		gradientDropletMass -= DropletMass[i] * outgoing;
		gradientVelocity -= velocity * outgoing;

		Delta[i] = new DiffusionCloud()
		{
			Temperature = gradientTemperature,
			Mass = gradientMass,
			Elevation = gradientElevation,
			DropletMass = gradientDropletMass,
			Velocity = gradientVelocity
		};

#endif
	}
}

#if !AdvectionWaterJobDebug
[BurstCompile]
#endif
public struct AdvectionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float2> Current;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{

	}
}

#if !PressureGradientForceAirJobDebug
[BurstCompile]
#endif
public struct PressureGradientForceAirJob : IJobParallelFor {
	public NativeArray<float2> Delta;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float2> Coords;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float InverseCoordDiff;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		float2 gradientPressure = 0;
		float2 coord = Coords[i];
		float elevation = LayerElevation[i] + LayerHeight[i] / 2;
		float pressure = Pressure[i];

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float2 coordDiff = Coords[n] - coord;
				float2 diff = math.normalize(math.float2(Utils.WrapAngle(coordDiff.x), Utils.WrapAngle(coordDiff.y))); // TODO: cache the normalized vectors in staticstate

				float neighborMidElevation = LayerElevation[n] + LayerHeight[n] / 2;
				float neighborElevationAtPressure = Atmosphere.GetElevationAtPressure(pressure, Temperature[n], Pressure[n], neighborMidElevation, Gravity);
				gradientPressure += diff * (neighborElevationAtPressure - elevation);
				//	Debug.Log("I: " + i + " D: " + diff + " P: " + pressure + " E: " + (neighborElevationAtPressure - elevation) + " NT: " + Temperature[n] + " NP: " + Pressure[n] + " NE: " + neighborMidElevation + " NEAP: " + neighborElevationAtPressure);
			}
		}
		//float density = Atmosphere.GetAirDensity(pressure, Temperature[i], AirMass[i], VaporMass[i]);
		float inverseDensity = WorldData.InverseDensityAir;
		Delta[i] = gradientPressure * Gravity * InverseCellDiameter * inverseDensity;
	}
}

#if !WindFrictionJobDebug
[BurstCompile]
#endif
public struct WindFrictionJob : IJobParallelFor {
	public NativeArray<float> Friction;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> VegetationCoverage;
	[ReadOnly] public float IceFriction;
	[ReadOnly] public float WaterFriction;
	[ReadOnly] public float TerrainFrictionMin;
	[ReadOnly] public float TerrainFrictionMax;
	[ReadOnly] public float VegetationFriction;
	[ReadOnly] public float MaxTerrainRoughness;
	public void Execute(int i)
	{
		float exposedIce = IceCoverage[i];
		float exposedWater = math.max(0, WaterCoverage[i] - exposedIce);
		float exposedVegetation = math.max(0, VegetationCoverage[i] - exposedIce - exposedWater);
		float exposedTerrain = 1.0f - exposedWater - exposedVegetation - exposedIce;
		Friction[i] = exposedIce * IceFriction +
			exposedWater * WaterFriction +
			exposedVegetation * VegetationFriction +
			exposedTerrain * (TerrainFrictionMin + (TerrainFrictionMax - TerrainFrictionMin) * math.saturate(Terrain[i].Roughness / MaxTerrainRoughness));
		//Debug.Log("I: " + i + " F: " + Friction[i]);
	}
}

#if !ConductionCloudAirJobDebug
[BurstCompile]
#endif
public struct ConductionCloudAirJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;

	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudDropletSize;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> ThermalDelta;
	[ReadOnly] public NativeArray<float> SolarIn;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> AirTemperature;
	[ReadOnly] public NativeArray<float> CloudTemperature;
	[ReadOnly] public NativeArray<float> CloudEnergy;
	[ReadOnly] public NativeArray<float> CloudCoverage;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;

	public void Execute(int i)
	{
		float cloudElevation = CloudElevation[i];
		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];
		if (cloudElevation >= layerElevation && cloudElevation < layerElevation + layerHeight)
		{
			float specificHeat = CloudMass[i] * WorldData.SpecificHeatWater;
			float cloudMass = CloudMass[i];
			float cloudTemp = CloudTemperature[i] + (ThermalDelta[i] + SolarIn[i]) / specificHeat;
			float tempDiff = AirTemperature[i] - cloudTemp;

			float delta = math.max(tempDiff * ConductionCoefficient * SecondsPerTick * CloudCoverage[i], -CloudEnergy[i]);
			float energyToEquilibrium = tempDiff * specificHeat;
			if (delta > 0)
			{
				delta = math.min(delta, energyToEquilibrium);
			}
			else
			{
				delta = math.max(delta, energyToEquilibrium);
			}
			EnergyDelta[i] = delta;
		}
	}
}

#if !ConductionAirIceJobDebug
[BurstCompile]
#endif
public struct ConductionAirIceJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyB;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		// TODO: this can conduct heat past a point of equilibrium
		EnergyDelta[i] = math.min((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * Coverage[i], EnergyB[i]);
	}
}

#if !ConductionAirWaterJobDebug
[BurstCompile]
#endif
public struct ConductionAirWaterJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> EnergyB;
	[ReadOnly] public NativeArray<float> CoverageIce;
	[ReadOnly] public NativeArray<float> CoverageWater;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		// TODO: this can conduct heat past a point of equilibrium
		EnergyDelta[i] = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * math.min(1.0f - CoverageIce[i], CoverageWater[i]), -EnergyA[i], EnergyB[i]);
	}
}

#if !ConductionIceWaterJobDebug
[BurstCompile]
#endif
public struct ConductionIceWaterJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> EnergyB;
	[ReadOnly] public NativeArray<float> CoverageA;
	[ReadOnly] public NativeArray<float> CoverageB;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		// TODO: this can conduct heat past a point of equilibrium
		EnergyDelta[i] = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * math.min(CoverageA[i], CoverageB[i]), -EnergyA[i], EnergyB[i]);
	}
}

#if !ConductionIceTerrainJobDebug
[BurstCompile]
#endif
public struct ConductionIceTerrainJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> CoverageIce;
	[ReadOnly] public NativeArray<float> CoverageWater;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		// TODO: this can conduct heat past a point of equilibrium
		EnergyDelta[i] = math.max((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * math.min(CoverageIce[i], 1.0f - CoverageWater[i]), -EnergyA[i]);
	}
}

#if !ConductionAirTerrainJobDebug
[BurstCompile]
#endif
public struct ConductionAirTerrainJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> CoverageIce;
	[ReadOnly] public NativeArray<float> CoverageWater;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		EnergyDelta[i] = (TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * (1.0f - math.max(CoverageIce[i], CoverageWater[i]));
	}
}

#if !ConductionWaterTerrainJobDebug
[BurstCompile]
#endif
public struct ConductionWaterTerrainJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	public NativeArray<float> EnergyDeltaWaterTotal;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public NativeArray<float> CoverageBelow;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		// TODO: this can conduct heat past a point of equilibrium
		float delta = math.max((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * Coverage[i] * (1.0f - CoverageBelow[i]), -EnergyA[i]);
		EnergyDelta[i] = delta;
		EnergyDeltaWaterTotal[i] += delta;
	}
}

#if !ConductionWaterBottomJobDebug
[BurstCompile]
#endif
public struct ConductionWaterBottomJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	public NativeArray<float> EnergyDeltaWaterTotal;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		// TODO: this can conduct heat past a point of equilibrium
		float delta = math.max((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * Coverage[i], -EnergyA[i]);
		EnergyDelta[i] = delta;
		EnergyDeltaWaterTotal[i] += delta;
	}
}


#if !StateChangeJobDebug
[BurstCompile]
#endif
public struct StateChangeJob : IJobParallelFor {
	public NativeArray<float> SurfaceAirTemperature;
	public NativeArray<float> SurfaceAirVapor;
	public NativeArray<float> SurfaceWaterMass;
	public NativeArray<float> SurfaceWaterTemperature;
	public NativeArray<float> IceMass;
	public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> IceMeltedMass;
	[ReadOnly] public NativeArray<float> WaterEvaporatedMass;
	[ReadOnly] public NativeArray<float> WaterFrozenMass;
	[ReadOnly] public NativeArray<float> RainfallWaterMass;
	[ReadOnly] public NativeArray<float> RainfallTemperature;
	[ReadOnly] public NativeArray<float> SurfaceAirMass;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	public void Execute(int i)
	{
		float vaporMass = SurfaceAirVapor[i];
		float airMass = SurfaceAirMass[i];
		float evaporatedMass = WaterEvaporatedMass[i];
		float waterMass = SurfaceWaterMass[i];
		float rainfallMass = RainfallWaterMass[i];
		float waterTemperature = SurfaceWaterTemperature[i];
		float rainfallTemperature = RainfallTemperature[i];
		float airTemperature = SurfaceAirTemperature[i];
		float iceMeltedMass = IceMeltedMass[i];
		float waterFrozenMass = WaterFrozenMass[i];
		float iceMass = IceMass[i];
		float saltMass = SurfaceSaltMass[i];

		float newWaterMass = waterMass + rainfallMass + iceMeltedMass;
		float newVaporMass = vaporMass + evaporatedMass;
		float newAirTemperature = airTemperature;

		float specificHeatAir = WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWater * vaporMass;
		newAirTemperature = (newAirTemperature * (airMass + vaporMass) + evaporatedMass * waterTemperature) / (airMass + vaporMass + evaporatedMass);

		float newWaterTemperature;
		if (newWaterMass > 0)
		{
			newWaterTemperature = (
				waterTemperature * (saltMass + waterMass)
				+ rainfallTemperature * rainfallMass
				+ WorldData.FreezingTemperature * iceMeltedMass)
				/ (waterMass + saltMass + rainfallMass + iceMeltedMass);
		} else
		{
			newWaterTemperature = 0;
		}

		float newIceMass = iceMass + waterFrozenMass;
		float newIceTemperature;
		if (newIceMass > 0)
		{
			newIceTemperature = (IceTemperature[i] * iceMass + WorldData.FreezingTemperature * waterFrozenMass) / newIceMass;
		} else
		{
			newIceTemperature = 0;
		}

		SurfaceWaterTemperature[i] = newWaterTemperature;
		SurfaceAirTemperature[i] = newAirTemperature;
		SurfaceWaterMass[i] = newWaterMass;
		SurfaceAirVapor[i] = newVaporMass;
		IceTemperature[i] = newIceTemperature;
		IceMass[i] = newIceMass;


	}
}

#if !StateChangeAirLayerJobDebug
[BurstCompile]
#endif
public struct StateChangeAirLayerJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudTemperature;
	public NativeArray<float> CloudDropletMass;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> SurfaceWaterMass;
	public NativeArray<float> SurfaceWaterTemperature;
	public NativeArray<float> AirTemperature;
	[ReadOnly] public NativeArray<float> GroundCondensationMass;
	[ReadOnly] public NativeArray<float> CloudCondensationMass;
	[ReadOnly] public NativeArray<float> CloudEvaporationMass;
	[ReadOnly] public NativeArray<float> SurfaceSaltMass;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	// TODO: implement
	//[ReadOnly] public NativeArray<float> AirTemperatureTop;
	//[ReadOnly] public NativeArray<float> AirTemperatureBottom;
	[ReadOnly] public int LayerIndex;
	public void Execute(int i)
	{
		float cloudElevation = CloudElevation[i];

		//float layerElevation = LayerElevation[LayerIndex];
		//float layerHeight = LayerHeight[LayerIndex];
		//bool isCloudLayer = newCloudElevation >= layerElevation && newCloudElevation < layerElevation + layerHeight;

		float cloudMass = CloudMass[i];
		float cloudTemperature = CloudTemperature[i];
		float cloudDropletMass = CloudDropletMass[i];
		float groundCondensationMass = GroundCondensationMass[i];
		float cloudCondensationMass = CloudCondensationMass[i];
		float cloudEvaporationMass = CloudEvaporationMass[i];
		float airTemperature = AirTemperature[i];


		float cloudMassDelta = math.max(-cloudMass, cloudCondensationMass - cloudEvaporationMass);
		float newCloudMass = cloudMass + cloudMassDelta;

		float newDropletSize = 0;
		if (newCloudMass > 0)
		{
			// TODO: use dew point and dew point elevation for condensation values
	//		cloudTemperature = ((cloudMass - cloudEvaporationMass) * cloudTemperature + airTemperature * cloudCondensationMass) / (cloudMass - cloudEvaporationMass + cloudCondensationMass); 
			cloudElevation = ((cloudMass - cloudEvaporationMass) * cloudElevation + (LayerElevation[i] + LayerHeight[i] / 2) * cloudCondensationMass) / (cloudMass - cloudEvaporationMass + cloudCondensationMass);
			newDropletSize = cloudDropletMass * (cloudMass - cloudEvaporationMass) / (cloudMass + cloudCondensationMass);
		}
		else
		{
			cloudTemperature = 0;
			cloudElevation = 0;
			newDropletSize = 0;
		}

		CloudMass[i] = newCloudMass;
		CloudDropletMass[i] = newDropletSize;
		CloudTemperature[i] = cloudTemperature;
		CloudElevation[i] = cloudElevation;

		float newWaterMass = SurfaceWaterMass[i] + groundCondensationMass;
		if (newWaterMass == 0)
		{
			SurfaceWaterMass[i] = 0;
			SurfaceWaterTemperature[i] = 0;
		} else
		{
			SurfaceWaterTemperature[i] = ((SurfaceWaterMass[i] + SurfaceSaltMass[i]) * SurfaceWaterTemperature[i] + groundCondensationMass * AirTemperature[i]) / (newWaterMass + SurfaceSaltMass[i]);
			SurfaceWaterMass[i] = newWaterMass;
		}
	}
}


#if !EnergySurfaceAirJobDebug
[BurstCompile]
#endif
public struct EnergySurfaceAirJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Vapor;
	public NativeArray<float2> Wind;
	public NativeArray<float> CondensationGroundMass;
	public NativeArray<float> CondensationCloudMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float2> LastWind;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	[ReadOnly] public NativeArray<DiffusionAir> Diffusion;
	[ReadOnly] public NativeArray<DiffusionAirVertical> VerticalAirMovement;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyCloud;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> WindFriction;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<float2> PressureGradientForce;
	[ReadOnly] public float WindFrictionMultiplier;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float DewPointZero;
	[ReadOnly] public float WaterVaporMassToAirMassAtDewPoint;
	[ReadOnly] public float InverseDewPointTemperatureRange;
	public void Execute(int i)
	{
		float airMass = AirMass[i];
		float conductionDelta =
			-ConductionEnergyCloud[i]
			+ ConductionEnergyIce[i]
			+ ConductionEnergyTerrain[i]
			+ ConductionEnergyWater[i];
		float specificHeat = WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWaterVapor * LastVapor[i];

		// TODO: account for differing Air masses
		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta;
		float temperature = LastTemperature[i] + Diffusion[i].Temperature - Advection[i].Temperature + energy / specificHeat;
		float vapor = LastVapor[i] + Diffusion[i].Humidity - Advection[i].Humidity;
#if !DISABLE_VERTICAL_AIR_MOVEMENT
		temperature += + VerticalAirMovement[i].Temperature;
		vapor += VerticalAirMovement[i].Humidity;

#endif

		float2 lastWind = LastWind[i];

		float2 wind = lastWind * (1.0f - WindFriction[i] * WindFrictionMultiplier) + Diffusion[i].Velocity - Advection[i].Velocity + PressureGradientForce[i] * SecondsPerTick;
		float2 coriolisForce = CoriolisMultiplier[i] * CoriolisTerm * new float2(wind.y, -wind.x);
		wind += coriolisForce * SecondsPerTick;

		float condensationGroundMass = 0;
		float condensationCloudMass = 0;

#if !DISABLE_CONDENSATION
		var relativeHumidity = Atmosphere.GetRelativeHumidity(airMass, vapor, temperature, DewPointZero, WaterVaporMassToAirMassAtDewPoint, InverseDewPointTemperatureRange);
		if (relativeHumidity > 1.0f)
		{
			float aboveCloud = math.saturate((CloudElevation[i] - LayerElevation[i]) / LayerHeight[i]);
			float vaporToCondense = (1.0f - 1.0f / relativeHumidity) * vapor;
			condensationCloudMass = aboveCloud * vaporToCondense;
			condensationGroundMass = (1.0f - aboveCloud) * vaporToCondense;
			vapor -= vaporToCondense;
		}
#endif

		CondensationGroundMass[i] = condensationGroundMass;
		CondensationCloudMass[i] = condensationCloudMass;

		Temperature[i] = temperature;
		Vapor[i] = vapor;
		Wind[i] = wind;
	}
}

#if !EnergyAirJobDebug
[BurstCompile]
#endif
public struct EnergyAirJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Vapor;
	public NativeArray<float2> Wind;
	public NativeArray<float> CondensationGroundMass;
	public NativeArray<float> CondensationCloudMass;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float2> LastWind;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	[ReadOnly] public NativeArray<DiffusionAir> Diffusion;
	[ReadOnly] public NativeArray<DiffusionAirVertical> VerticalAirMovementUp;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyCloud;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> WindFriction;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<float2> PressureGradientForce;
	[ReadOnly] public float WindFrictionMultiplier;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float DewPointZero;
	[ReadOnly] public float WaterVaporMassToAirMassAtDewPoint;
	[ReadOnly] public float InverseDewPointTemperatureRange;
	public void Execute(int i)
	{
		float airMass = AirMass[i];
		float conductionDelta =	-ConductionEnergyCloud[i];
		float specificHeat = WorldData.SpecificHeatAtmosphere * airMass + WorldData.SpecificHeatWaterVapor * LastVapor[i];

		// TODO: account for differing Air masses
		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta;
		float temperature = LastTemperature[i] + Diffusion[i].Temperature - Advection[i].Temperature /*+ VerticalAirMovementUp[i].Temperature*/ + energy / specificHeat;
		float vapor = LastVapor[i] + Diffusion[i].Humidity - Advection[i].Humidity /*+ VerticalAirMovementUp[i].Humidity*/;
		float2 lastWind = LastWind[i];

		float2 wind = lastWind * (1.0f - WindFriction[i] * WindFrictionMultiplier) + Diffusion[i].Velocity - Advection[i].Velocity + PressureGradientForce[i] * SecondsPerTick;
		float2 coriolisForce = CoriolisMultiplier[i] * CoriolisTerm * new float2(wind.y, -wind.x);
		wind += coriolisForce * SecondsPerTick;

		float condensationGroundMass = 0;
		float condensationCloudMass = 0;

#if !DISABLE_CONDENSATION
		var relativeHumidity = Atmosphere.GetRelativeHumidity(airMass, vapor, temperature, DewPointZero, WaterVaporMassToAirMassAtDewPoint, InverseDewPointTemperatureRange);
		if (relativeHumidity > 1.0f)
		{
			float aboveCloud = math.saturate((CloudElevation[i] - LayerElevation[i]) / LayerHeight[i]);
			float vaporToCondense = (1.0f - 1.0f / relativeHumidity) * vapor;
			condensationCloudMass = aboveCloud * vaporToCondense;
			condensationGroundMass = (1.0f - aboveCloud) * vaporToCondense;
			vapor -= vaporToCondense;
		}
#endif

		CondensationGroundMass[i] = condensationGroundMass;
		CondensationCloudMass[i] = condensationCloudMass;

		Temperature[i] = temperature;
		Vapor[i] = vapor;
		Wind[i] = wind;
	}
}

#if !EnergyWaterJobDebug
[BurstCompile]
#endif
public struct EnergyWaterJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> SaltMass;
	public NativeArray<float2> Velocity;
	public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float2> LastVelocity;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<DiffusionWater> Advection;
	[ReadOnly] public NativeArray<DiffusionWater> Diffusion;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float conductionDelta = ConductionEnergyTerrain[i];
		float specificHeat = WorldData.SpecificHeatWater * LastMass[i] + WorldData.SpecificHeatSalt * LastSaltMass[i];
		float energySources = 0;
		if (specificHeat > 0)
		{
			energySources = (SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta) / specificHeat;
		}
		Temperature[i] = LastTemperature[i] + Diffusion[i].Temperature - Advection[i].Temperature + energySources;
		SaltMass[i] = LastSaltMass[i];
		Mass[i] = LastMass[i];

		float2 pressureGradientForce = float2.zero;
		float2 coriolisForce = float2.zero;

		Velocity[i] = LastVelocity[i] + Diffusion[i].Velocity - Advection[i].Velocity + (pressureGradientForce + coriolisForce) * SecondsPerTick;
	}
}

#if !EnergyWaterJobSurfaceDebug
[BurstCompile]
#endif
public struct EnergyWaterJobSurface : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> SaltMass;
	public NativeArray<float2> Velocity;
	public NativeArray<float> WaterMass;
	public NativeArray<float> EvaporatedWaterMass;
	public NativeArray<float> FrozenMass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float2> LastVelocity;
	[ReadOnly] public NativeArray<float2> LastSurfaceWind;
	[ReadOnly] public NativeArray<DiffusionWater> Advection;
	[ReadOnly] public NativeArray<DiffusionWater> Diffusion;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaTop;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaBottom;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public NativeArray<float> RelativeHumidity;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public float EvaporationRate;
	[ReadOnly] public float EvapTemperatureMin;
	[ReadOnly] public float EvapTemperatureMax;
	[ReadOnly] public float WaterHeatingDepth;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float lastTemperature = LastTemperature[i];
		float lastMassWater = LastMass[i];
		float lastMassSalt = LastSaltMass[i];
		float energyTop = -ConductionEnergyAir[i] - ConductionEnergyIce[i] + ThermalRadiationDeltaTop[i] + SolarRadiationIn[i];
		float energyBottom = ConductionEnergyTerrain[i] + ThermalRadiationDeltaBottom[i];

		float temperature = lastTemperature;
		float saltMass = LastSaltMass[i];
		float2 velocity = LastVelocity[i];
		float waterMass = lastMassWater;

		float evapMass = 0;
		float frozenTopMass = 0;
		float frozenBottomMass = 0;
		if (waterMass > 0)
		{

#if !DISABLE_EVAPORATION
			if (energyTop > 0)
			{
				float evapRate = math.saturate((1.0f - RelativeHumidity[i]) * math.max(0, WaterCoverage[i] - IceCoverage[i]) * EvaporationRate * (LastTemperature[i] - EvapTemperatureMin) / (EvapTemperatureMax - EvapTemperatureMin));
				evapMass = math.clamp(evapRate * energyTop / WorldData.LatentHeatWaterVapor, 0, waterMass);
				if (evapRate > 0)
				{
					float evapEnergy = evapMass * WorldData.LatentHeatWaterVapor / evapRate;
					energyTop -= evapEnergy;
				}
				waterMass -= evapMass;
			}
#endif

#if !DISABLE_FREEZE_TOP
			if (waterMass > 0)
			{
				float heatingMass = math.min(waterMass, WaterHeatingDepth * WorldData.MassWater);
				float specificHeatSaltWater = (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatWater * saltMass) / (waterMass + saltMass);
				float newTempTop = temperature + energyTop / (specificHeatSaltWater * heatingMass);
				if (newTempTop < WorldData.FreezingTemperature)
				{
					float energyToRelease = (WorldData.FreezingTemperature - newTempTop) * specificHeatSaltWater * heatingMass;
					frozenTopMass = math.min(waterMass, energyToRelease / WorldData.LatentHeatWaterLiquid);
					energyTop += frozenTopMass * WorldData.LatentHeatWaterLiquid;
					waterMass -= frozenTopMass;
				}
			}
#endif

#if !DISABLE_FREEZE_BOTTOM
			if (waterMass > 0)
			{
				float heatingMass = math.min(waterMass, WaterHeatingDepth * WorldData.MassWater);
				float specificHeatSaltWater = (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatWater * saltMass) / (waterMass + saltMass);
				float newTempBottom = temperature + energyBottom / (specificHeatSaltWater * heatingMass);
				if (newTempBottom < WorldData.FreezingTemperature)
				{
					float energyToRelease = (WorldData.FreezingTemperature - newTempBottom) * specificHeatSaltWater * heatingMass;
					frozenBottomMass = math.min(waterMass, energyToRelease / WorldData.LatentHeatWaterLiquid);
					energyBottom += frozenBottomMass * WorldData.LatentHeatWaterLiquid;
					waterMass -= frozenBottomMass;
				}
			}
#endif
		}


		EvaporatedWaterMass[i] = evapMass;
		FrozenMass[i] = frozenTopMass + frozenBottomMass;

		SaltMass[i] = saltMass + Diffusion[i].Salinity - Advection[i].Salinity;
		WaterMass[i] = waterMass;
		if (waterMass == 0)
		{
			Velocity[i] = 0;
			Temperature[i] = 0;
		} else
		{
			float specificHeat = WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatWater * saltMass;
			float energySources = (energyTop + energyBottom) / specificHeat;
			Temperature[i] = temperature + Diffusion[i].Temperature - Advection[i].Temperature + energySources;

			float2 pressureGradientForce = float2.zero;
			float2 coriolisForce = float2.zero;

			Velocity[i] = velocity + Diffusion[i].Velocity - Advection[i].Velocity + (pressureGradientForce + coriolisForce) * SecondsPerTick;
		}
	}
}

#if !EnergyCloudJobDebug
[BurstCompile]
#endif
public struct EnergyCloudJob : IJobParallelFor {
	public NativeArray<float> CloudTemperature;
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> DropletMass;
	public NativeArray<float2> Velocity;
	public NativeArray<float> RainfallWaterMass;
	public NativeArray<float> CloudEvaporationMass;
	[ReadOnly] public NativeArray<DiffusionCloud> Advection;
	[ReadOnly] public NativeArray<DiffusionCloud> Diffusion;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> LastCloudTemperature;
	[ReadOnly] public NativeArray<float2> LastVelocity;
	[ReadOnly] public NativeArray<float> LastCloudElevation;
	[ReadOnly] public NativeArray<float> LastDropletMass;
	[ReadOnly] public NativeArray<float> LastCloudMass;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> AirMassCloud;
	[ReadOnly] public NativeArray<float> WaterVaporCloud;
	[ReadOnly] public NativeArray<float> AirTemperatureCloud;
	[ReadOnly] public NativeArray<float> AirPressureCloud;
	[ReadOnly] public NativeArray<float> RelativeHumidityCloud;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float> WindFriction;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<float2> PressureGradientForce;
	[ReadOnly] public float WindFrictionMultiplier;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float RainDropMinSize;
	[ReadOnly] public float RainDropMaxSize;
	[ReadOnly] public float RainDropDragCoefficient;
	[ReadOnly] public float RainfallRate;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float TicksPerSecond;
	[ReadOnly] public float CloudDissapationRateDryAir;
	[ReadOnly] public float CloudDissapationRateWind;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float surfaceElevation = SurfaceElevation[i];
		float airTemperatureAtCloud = AirTemperatureCloud[i];

		float cloudMass = LastCloudMass[i] + Diffusion[i].Mass - Advection[i].Mass;
		float cloudTemperature = LastCloudTemperature[i];
		float dropletMass = LastDropletMass[i];
		float cloudElevation = LastCloudElevation[i];
		float velocityVertical = 0;
		float2 lastVelocity = LastVelocity[i];
		float rainfallWaterMass = 0;
		float cloudEvaporationMass = 0;

		float energySources = 0;
		if (cloudMass > 0)
		{
			energySources = (SolarRadiationIn[i] + ThermalRadiationDelta[i] + ConductionEnergyAir[i]) / (WorldData.SpecificHeatWater * cloudMass);
			cloudTemperature += Diffusion[i].Temperature - Advection[i].Temperature + energySources;
		}
		dropletMass += Diffusion[i].DropletMass - Advection[i].DropletMass;
#if !DISABLE_CLOUD_ELEVATION
		cloudElevation = math.max(surfaceElevation, cloudElevation + Diffusion[i].Elevation - Advection[i].Elevation);
#endif

		float2 velocity = lastVelocity * (1.0f - WindFriction[i] * WindFrictionMultiplier) + Diffusion[i].Velocity - Advection[i].Velocity + PressureGradientForce[i] * SecondsPerTick;
		float2 coriolisForce = CoriolisMultiplier[i] * CoriolisTerm * new float2(velocity.y, -velocity.x);
		velocity += coriolisForce * SecondsPerTick;

#if !DISABLE_CLOUD_ELEVATION
		float airDensityAtElevation = Atmosphere.GetAirDensity(AirPressureCloud[i], airTemperatureAtCloud, AirMassCloud[i], WaterVaporCloud[i]);
		float waterDensityAtElevation = Atmosphere.GetWaterDensityAtElevation(CloudTemperature[i], cloudElevation);
		float rainDropRadius = math.clamp(math.pow(dropletMass / waterDensityAtElevation * 3 / (4 * math.PI), 0.333f), RainDropMinSize, RainDropMaxSize);
		float rainDropVolume = 4 / 3 * math.PI * math.pow(rainDropRadius, 3);
		velocityVertical = math.sqrt(8 * rainDropRadius * waterDensityAtElevation * Gravity / (3 * airDensityAtElevation * RainDropDragCoefficient));

		cloudElevation = math.max(surfaceElevation, cloudElevation + velocityVertical * SecondsPerTick);
#endif

#if !DISABLE_RAINFALL
		if (velocityVertical < 0 && dropletMass > 0)
		{
			// TODO: use this to detemine temperature when it hits the ground (snow or rain)
			float rainDropFallTime = (cloudElevation - surfaceElevation) / -velocityVertical;
			if (rainDropFallTime < SecondsPerTick)
			{
				// TODO: account for drop size variance so we don't just dump it all at once
				rainfallWaterMass = cloudMass * (1.0f - rainDropFallTime / SecondsPerTick);
				cloudMass -= rainfallWaterMass;
				dropletMass *= 1.0f - rainfallWaterMass / cloudMass;
			}
		}
#endif

#if !DISABLE_CLOUD_DISSAPATION
		// dissapation
		if (cloudMass > 0)
		{
			float dissapationSpeed = math.min(1.0f, CloudDissapationRateWind * math.max(0, -velocityVertical) + CloudDissapationRateDryAir) * math.max(0, 1.0f - RelativeHumidityCloud[i]);
			cloudEvaporationMass = cloudMass * dissapationSpeed;
			dropletMass *= 1.0f - dissapationSpeed;
			cloudMass -= cloudEvaporationMass;
		}
#endif





		if (cloudMass <= 0)
		{
			cloudTemperature = 0;
			cloudElevation = surfaceElevation;
			dropletMass = 0;
			velocityVertical = 0;
			velocity = 0;
		}

		CloudTemperature[i] = cloudTemperature;
		DropletMass[i] = dropletMass;
		CloudMass[i] = cloudMass;
		CloudElevation[i] = cloudElevation;
		Velocity[i] = velocity;

		RainfallWaterMass[i] = rainfallWaterMass;
		CloudEvaporationMass[i] = cloudEvaporationMass;
	}
}

#if !EnergyIceJobDebug
[BurstCompile]
#endif
public struct EnergyIceJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Mass;
	public NativeArray<float> MeltedMass;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaBottom;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaTop;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	[ReadOnly] public float IceHeatingDepth;
	public void Execute(int i)
	{
		float energyTop = -ConductionEnergyAir[i] + ThermalRadiationDeltaTop[i] + SolarRadiationIn[i];
		float energyBottom = ConductionEnergyWater[i] + ConductionEnergyTerrain[i] + ThermalRadiationDeltaBottom[i];

		float meltedTopMass = 0;
		float meltedBottomMass = 0;
		float iceMass = LastMass[i];
		float iceTemperature = LastTemperature[i];

#if !DISABLE_MELTING_TOP
		if (iceMass > 0)
		{
			float iceDepth = iceMass / WorldData.MassIce;
			float heatingDepth = math.min(IceHeatingDepth, iceDepth);
			float heatingMass = heatingDepth * WorldData.MassIce;
			float newTempTop = iceTemperature + energyTop / (WorldData.SpecificHeatIce * heatingMass);
			if (newTempTop > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (newTempTop - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * heatingMass;
				meltedTopMass = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				energyTop -= meltedTopMass * WorldData.LatentHeatWaterLiquid;
				iceMass -= meltedTopMass;
			}
		}
#endif

#if !DISABLE_MELTING_BOTTOM
		if (iceMass > 0)
		{
			float iceDepth = iceMass / WorldData.MassIce;
			float heatingDepth = math.min(IceHeatingDepth, iceDepth);
			float heatingMass = heatingDepth * WorldData.MassIce;
			float newTempBottom = iceTemperature + energyBottom / (WorldData.SpecificHeatIce * heatingMass);
			if (newTempBottom > WorldData.FreezingTemperature)
			{
				float energyToAbsorb = (newTempBottom - WorldData.FreezingTemperature) * WorldData.SpecificHeatIce * heatingMass;
				meltedBottomMass = math.min(iceMass, energyToAbsorb / WorldData.LatentHeatWaterLiquid);
				energyBottom -= meltedBottomMass * WorldData.LatentHeatWaterLiquid;
				iceMass -= meltedBottomMass;
			}
		}
#endif

		Mass[i] = iceMass;
		MeltedMass[i] = meltedTopMass + meltedBottomMass;
		Temperature[i] = iceMass > 0 ? iceTemperature + (energyTop + energyBottom) / (iceMass * WorldData.SpecificHeatIce) : 0;
	}
}

#if !EnergyTerrainJobDebug
[BurstCompile]
#endif
public struct EnergyTerrainJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> VegetationCoverage;
	[ReadOnly] public float GeothermalEnergy;
	[ReadOnly] public float HeatingDepth;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyAir[i]
			-ConductionEnergyWater[i]
			-ConductionEnergyIce[i];
		float energy = SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta + GeothermalEnergy;
		float specificHeat = Atmosphere.GetSpecificHeatTerrain(HeatingDepth, Terrain[i].SoilFertility, VegetationCoverage[i]);
		Temperature[i] = LastTemperature[i] + energy / specificHeat;
	}
}

#if !UpdateTerrainJobDebug
[BurstCompile]
#endif
public struct UpdateTerrainJob : IJobParallelFor {
	public NativeArray<CellTerrain> Terrain;

	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;
	public void Execute(int i)
	{
		Terrain[i] = LastTerrain[i];
	}

}

#if !UpdateWaterSaltMassJobDebug
[BurstCompile]
#endif
public struct UpdateWaterSaltMassJob : IJobParallelFor {
	public NativeArray<WaterSaltMass> WaterSaltMass;
	[ReadOnly] public NativeArray<float> WaterLayerMass;
	[ReadOnly] public NativeArray<float> SaltLayerMass;
	public void Execute(int i)
	{
		WaterSaltMass[i] = new WaterSaltMass()
		{
			WaterMass = WaterSaltMass[i].WaterMass + WaterLayerMass[i],
			SaltMass = WaterSaltMass[i].SaltMass + SaltLayerMass[i],
		};
	}

}

#if !UpdateDependentStateJobDebug
[BurstCompile]
#endif
public struct UpdateDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceElevation;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> CloudCoverage;
	public NativeArray<float> WaterDepth;
	public NativeArray<float> SurfaceAirTemperature;
	public NativeArray<float> CloudEnergy;
	public NativeArray<float> IceEnergy;
	[ReadOnly] public NativeArray<WaterSaltMass> WaterSaltMass;
	[ReadOnly] public NativeArray<float> LowerAirTemperature;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudTemperature;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public WorldData worldData;
	[ReadOnly] public NativeArray<float> lowerAirHeight;
	public void Execute(int i)
	{
		float waterDepth = WaterSaltMass[i].WaterMass / WorldData.MassWater + WaterSaltMass[i].SaltMass / WorldData.MassSalt;
		float iceMass = IceMass[i];
		WaterDepth[i] = waterDepth;
		SurfaceElevation[i] = Terrain[i].Elevation + waterDepth + iceMass / WorldData.MassIce;
		VegetationCoverage[i] = math.saturate(Terrain[i].Vegetation * worldData.inverseFullCanopyCoverage);

		IceCoverage[i] = math.saturate(iceMass * worldData.inverseFullIceCoverage);
		IceEnergy[i] = WorldData.SpecificHeatIce * iceMass * IceTemperature[i];

		float cloudMass = CloudMass[i];
		CloudCoverage[i] = math.saturate(cloudMass * worldData.inverseCloudMassFullAbsorption);
		CloudEnergy[i] = WorldData.SpecificHeatWater * cloudMass * CloudTemperature[i];

		SurfaceAirTemperature[i] = LowerAirTemperature[i] - WorldData.TemperatureLapseRate * lowerAirHeight[i] / 2;
	}

}

#if !UpdateAirPressureJobDebug
[BurstCompile]
#endif
public struct UpdateAirPressureJob : IJobParallelFor {
	public NativeArray<float> Pressure;
	public NativeArray<float> AirMass;

	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> AirTemperature;
	[ReadOnly] public NativeArray<float2> PressureGradient;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float DewPointZero;
	[ReadOnly] public float WaterVaporMassToAirMassAtDewPoint;
	[ReadOnly] public float InverseDewPointTemperatureRange;
	[ReadOnly] public int LayerIndex;

	public void Execute(int i)
	{
	}
}

#if !UpdateDependentAirLayerJobDebug
[BurstCompile]
#endif
public struct UpdateDependentAirLayerJob : IJobParallelFor {
	public NativeArray<float> AbsoluteHumidity;
	public NativeArray<float> RelativeHumidity;
	public NativeArray<float> Pressure;
	public NativeArray<float> PotentialEnergy;
	public NativeArray<float> AirMass;
	public NativeArray<float> AirMassTotal;
	public NativeArray<float> AirMassCloud;
	public NativeArray<float> AirVaporCloud;
	public NativeArray<float> AirTemperatureCloud;
	public NativeArray<float> AirPressureCloud;
	public NativeArray<float> AirHumidityRelativeCloud;
	public NativeArray<int> AirLayerCloud;

	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> AirTemperature;
	[ReadOnly] public NativeArray<float> CloudTemperature;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float DewPointZero;
	[ReadOnly] public float WaterVaporMassToAirMassAtDewPoint;
	[ReadOnly] public float InverseDewPointTemperatureRange;
	[ReadOnly] public int LayerIndex;

	public void Execute(int i)
	{
		float layerMiddle = LayerElevation[i] + LayerHeight[i] / 2;
		float vaporMass = VaporMass[i];
		float airTemperature = AirTemperature[i];
		float airMass = Atmosphere.GetAirMass(LayerElevation[i], LayerHeight[i], airTemperature, Gravity);

		AbsoluteHumidity[i] = vaporMass / (vaporMass + airMass);
		RelativeHumidity[i] = Atmosphere.GetRelativeHumidity(airMass, vaporMass, airTemperature, DewPointZero, WaterVaporMassToAirMassAtDewPoint, InverseDewPointTemperatureRange);
		PotentialEnergy[i] = (airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * airTemperature;
		AirMass[i] = airMass;
		Pressure[i] = (AirMassTotal[i] + airMass / 2) * Gravity;

		float cloudElevation = CloudElevation[i];
		float layerElevation = LayerElevation[i];
		float layerHeight = LayerHeight[i];
		if (cloudElevation >= layerElevation && cloudElevation < layerElevation + layerHeight)
		{
			float airMassCloud = AirMass[i];
			float vaporMassCloud = VaporMass[i];
			float airTemperatureCloud = airTemperature;
			AirLayerCloud[i] = LayerIndex;
			AirMassCloud[i] = airMassCloud;
			AirTemperatureCloud[i] = airTemperatureCloud;
			AirVaporCloud[i] = VaporMass[i];
			AirPressureCloud[i] = (AirMassTotal[i] + airMass * (1.0f - (cloudElevation - layerElevation) / layerHeight)) * Gravity;
			AirHumidityRelativeCloud[i] = Atmosphere.GetRelativeHumidity(airMassCloud, vaporMassCloud, airTemperatureCloud, DewPointZero, WaterVaporMassToAirMassAtDewPoint, InverseDewPointTemperatureRange);
		}

		AirMassTotal[i] += airMass;

	}
}

#if !UpdateDependentWaterLayerJobDebug
[BurstCompile]
#endif
public struct UpdateDependentWaterLayerJob : IJobParallelFor {
	public NativeArray<float> Salinity;
	public NativeArray<float> WaterCoverage;
	public NativeArray<float> PotentialEnergy;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			float waterDepth = WaterMass[i] / WorldData.MassWater + SaltMass[i] / WorldData.MassSalt;
			WaterCoverage[i] = math.min(1, waterDepth / math.max(1, Terrain[i].Roughness));
			Salinity[i] = SaltMass[i] / (WaterMass[i] + SaltMass[i]);
			PotentialEnergy[i] = (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt) * Temperature[i];
		}
		else
		{
			WaterCoverage[i] = 0;
			Salinity[i] = 0;
			PotentialEnergy[i] = 0;
		}
	}
}

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
public struct InitDisplayJob : IJobParallelFor {
	public NativeArray<float> DisplayPressure;
	[ReadOnly] public NativeArray<float> AirLayerElevation;
	[ReadOnly] public NativeArray<float> AirLayerHeight;
	[ReadOnly] public NativeArray<float> AirTemperature;
	[ReadOnly] public NativeArray<float> AirPressure;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		DisplayPressure[i] = Atmosphere.GetPressureAtElevation(0, Gravity, AirPressure[i], AirTemperature[i], AirLayerElevation[i] + AirLayerHeight[i] / 2);
	}
}
//[BurstCompile]
//public struct TickEnergyDeltaJob : IJobParallelFor {

//	public NativeArray<CellEnergyDelta> Next;
//	public NativeArray<CellDisplay> NextDisplay;

//	[ReadOnly] public NativeArray<float> SolarRadiation;
//	[ReadOnly] public PlanetState PlanetState;
//	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;
//	[ReadOnly] public NativeArray<CellAtmosphere> LastAtmosphere;
//	[ReadOnly] public NativeArray<CellWater> LastWater;
//	[ReadOnly] public WorldData worldData;
//	[ReadOnly] public StaticState staticState;

//	public void Execute(int i)
//	{
//		var last = LastTerrain[i];
//		var lastTerrain = LastTerrain[i];
//		var next = new CellEnergyDelta();

//		float airMass = ;
//		float lastAirEnergy = ;
//		float lastWaterVaporMass = ;

//		float nextAirEnergy = ;
//		float nextWaterVaporMass = ;

//		float iceCoverage = math.min(1.0f, math.pow(last.IceMass * worldData.inverseFullIceCoverage, 0.6667f));
//		float surfaceElevation = lastTerrain.SurfaceElevation;
//		float evaporationRate = Atmosphere.GetEvaporationRate(ref worldData, iceCoverage, lastDependent.AirTemperature, lastDependent.RelativeHumidity, worldData.inverseEvapTemperatureRange);
//		float dewPoint = Atmosphere.GetDewPoint(ref worldData, lastDependent.AirTemperature, lastDependent.RelativeHumidity);

//		float waterCoverage = math.min(1.0f, math.pow(lastTerrain.WaterDepth * worldData.inverseFullWaterCoverage, 0.6667f));
//		float canopyCoverage = math.min(1.0f, math.pow(lastTerrain.Vegetation * worldData.inverseFullCanopyCoverage, 0.6667f));

//		//float groundWaterSaturation = Animals.GetGroundWaterSaturation(state.GroundWater[index], state.WaterTableDepth[index], soilFertility * world.Data.MaxSoilPorousness);


//		// SOLAR RADIATION

//		float solarRadiationAbsorbed = 0;
//		float solarRadiation = PlanetState.SolarRadiation * worldData.SecondsPerTick * sunDotSurface;

//		// get the actual atmospheric depth here based on radius of earth plus atmosphere
//		//float inverseSunAngle = PIOver2 + sunAngle;
//		//float angleFromSunToLatitudeAndAtmophereEdge = math.Asin(state.PlanetRadius * math.Sin(inverseSunAngle) / (state.PlanetRadius + world.Data.TropopauseElevation));
//		//float angleFromPlanetCenterToLatitudeAndAtmosphereEdge = math.PI - inverseSunAngle - angleFromSunToLatitudeAndAtmophereEdge;
//		//float atmosphericDepthInMeters = math.Sin(angleFromPlanetCenterToLatitudeAndAtmosphereEdge) * state.PlanetRadius / math.Sin(angleFromSunToLatitudeAndAtmophereEdge);
//		//float atmosphericDepth = math.max(1.0f, atmosphericDepthInMeters / world.Data.TropopauseElevation);

//		float atmosphericDepth = 1.0f + (1.0f - sunDotSurface);

//		// These constants obtained here, dunno if I've interpreted them correctly
//		// https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass

//		////// MAJOR TODO:
//		///// USE THIS LINK: https://www.ftexploring.com/solar-energy/sun-angle-and-insolation2.htm
//		/// With the sun 90 degrees above the horizon (SEA° = 90°), the air mass lowers the intensity of the sunlight from the 1,367 W / m2 that it is in outerspace down to about 1040 W / m2.
//		//				float consumedByAtmosphere = 1.0f - math.pow(0.7f, math.pow(atmosphericDepth, 0.678f));

//		if (solarRadiation > 0)
//		{
//			displayCell.EnergyIncoming += solarRadiation;

//			// TODO: reflect/absorb more in the atmosphere with a lower sun angle

//			// reflect some rads off atmosphere and clouds
//			// TODO: this process feels a little broken -- are we giving too much priority to reflecting/absorbing in certain layers?
//			float energyReflectedAtmosphere = solarRadiation * math.min(1, worldData.SolarReflectivityAir * (airMass + lastWaterVaporMass));
//			solarRadiation -= energyReflectedAtmosphere;
//			displayCell.EnergySolarReflectedAtmosphere += energyReflectedAtmosphere;

//			if (last.CloudMass > 0)
//			{
//				float cloudTemperatureAlbedo = WorldData.AlbedoIce + (WorldData.AlbedoWater - WorldData.AlbedoIce) * math.saturate((dewPoint - worldData.minCloudFreezingTemperature) / (worldData.maxCloudFreezingTemperature - worldData.minCloudFreezingTemperature));
//				float rainDropSizeAlbedo = math.clamp(1.0f - last.CloudDropletMass / last.CloudMass, 0, 1) * (worldData.rainDropSizeAlbedoMax - worldData.rainDropSizeAlbedoMin) + worldData.rainDropSizeAlbedoMin;
//				float cloudReflectivity = math.min(1.0f, worldData.SolarReflectivityCloud * cloudTemperatureAlbedo * last.CloudMass * rainDropSizeAlbedo / math.max(worldData.maxCloudSlopeAlbedo, 1.0f - waterSlopeAlbedo));
//				float energyReflectedClouds = solarRadiation * cloudReflectivity;
//				solarRadiation -= energyReflectedClouds;

//				float absorbedByCloudsIncoming = solarRadiation * math.min(1.0f, worldData.SolarAbsorptivityCloud * last.CloudMass);
//				solarRadiation -= absorbedByCloudsIncoming;
//				nextAirEnergy += absorbedByCloudsIncoming;
//				displayCell.EnergySolarAbsorbedCloud += absorbedByCloudsIncoming;
//				displayCell.EnergySolarAbsorbedAtmosphere += absorbedByCloudsIncoming;
//				displayCell.EnergySolarReflectedCloud += energyReflectedClouds;
//			}

//			// Absorbed by atmosphere
//			// stratosphere accounts for about a quarter of atmospheric mass
//			//	float absorbedByStratosphere = incomingRadiation * world.Data.AtmosphericHeatAbsorption * (state.StratosphereMass / massOfAtmosphericColumn);

//			float atmosphereAbsorptionRate = math.min(1, worldData.SolarAbsorptivityAir * airMass + worldData.SolarAbsorptivityWaterVapor * lastWaterVaporMass);
//			float absorbedByAtmosphereIncoming = solarRadiation * atmosphereAbsorptionRate * atmosphericDepth;

//			nextAirEnergy += absorbedByAtmosphereIncoming;
//			solarRadiation -= absorbedByAtmosphereIncoming;
//			displayCell.EnergySolarAbsorbedAtmosphere += absorbedByAtmosphereIncoming;

//			// reflection off surface
//			float energyReflected = 0;
//			{
//				if (iceCoverage > 0)
//				{
//					energyReflected += solarRadiation * iceCoverage * Atmosphere.GetAlbedo(WorldData.AlbedoIce, 0);
//				}
//				if (waterCoverage > 0)
//				{
//					energyReflected += waterCoverage * solarRadiation * Atmosphere.GetAlbedo(WorldData.AlbedoWater, waterSlopeAlbedo) * (1.0f - iceCoverage);
//				}
//				if (waterCoverage < 1 && iceCoverage < 1)
//				{
//					// reflect some incoming radiation
//					float slopeAlbedo = 0;
//					float soilReflectivity = Atmosphere.GetAlbedo(WorldData.AlbedoLand - worldData.AlbedoReductionSoilQuality * lastTerrain.SoilFertility, slopeAlbedo);
//					float heatReflectedLand = canopyCoverage * WorldData.AlbedoFoliage + math.max(0, 1.0f - canopyCoverage) * soilReflectivity;
//					energyReflected += solarRadiation * math.saturate(heatReflectedLand) * (1.0f - iceCoverage) * (1.0f - waterCoverage);
//				}
//				solarRadiation -= energyReflected;

//				// TODO: do we absorb some of this energy on the way back out of the atmosphere?
//				displayCell.EnergySolarReflectedSurface += energyReflected;
//			}

//			solarRadiationAbsorbed += solarRadiation;

//		}

//		// THERMAL RADIATION

//		float evaporation = 0;
//		float backRadiation = 0;
//		float reflected = 0;

//		// radiate heat from land
//		// TODO: deal with the fact that this also incorporates ground water
//		float shallowWaterRadiation = Atmosphere.GetRadiationRate(lastDependent.WaterTemperature, WorldData.EmissivityWater) * worldData.SecondsPerTick * waterCoverage;
//		float soilEnergy = last.GroundEnergy - last.GroundWater * worldData.maxGroundWaterTemperature * WorldData.SpecificHeatWater;
//		float radiationRate = Atmosphere.GetLandRadiationRate(ref worldData, last.GroundEnergy, last.GroundWater, lastTerrain.SoilFertility, canopyCoverage);
//		float thermalEnergyRadiatedLand = math.min(soilEnergy, radiationRate * worldData.SecondsPerTick);
//		next.GroundEnergyDelta += PlanetState.GeothermalHeat * worldData.SecondsPerTick - thermalEnergyRadiatedLand;

//		float thermalEnergyRadiatedToIce = 0;
//		float thermalEnergyRadiatedToShallowWater = 0;
//		float thermalEnergyRadiatedToAir = 0;

//		next.WaterEnergy -= shallowWaterRadiation;
//		next.GroundEnergy += shallowWaterRadiation;

//		thermalEnergyRadiatedToShallowWater += thermalEnergyRadiatedLand * waterCoverage;
//		thermalEnergyRadiatedToIce += iceCoverage * (thermalEnergyRadiatedLand - thermalEnergyRadiatedToShallowWater);
//		thermalEnergyRadiatedToAir += thermalEnergyRadiatedLand - thermalEnergyRadiatedToShallowWater - thermalEnergyRadiatedToIce;
//		next.WaterEnergy += thermalEnergyRadiatedToShallowWater;

//		// radiate from ice to air above
//		if (iceCoverage > 0)
//		{
//			float thermalEnergyRadiatedFromIce = iceCoverage * Atmosphere.GetRadiationRate(WorldData.FreezingTemperature, WorldData.EmissivityIce) * worldData.SecondsPerTick;
//			thermalEnergyRadiatedToAir += thermalEnergyRadiatedFromIce;
//			thermalEnergyRadiatedToIce -= thermalEnergyRadiatedFromIce;
//		}

//		// lose heat to air via conduction AND radiation
//		if (iceCoverage < 1)
//		{
//			// radiate heat, will be absorbed by air
//			// Net Back Radiation: The ocean transmits electromagnetic radiation into the atmosphere in proportion to the fourth power of the sea surface temperature(black-body radiation)
//			// https://eesc.columbia.edu/courses/ees/climate/lectures/o_atm.html
//			float oceanRadiation = math.min(last.WaterEnergy, shallowWaterRadiation);
//			next.WaterEnergy -= oceanRadiation;
//			thermalEnergyRadiatedToIce += oceanRadiation * iceCoverage;
//			thermalEnergyRadiatedToAir += oceanRadiation * (1.0f - iceCoverage);
//			displayCell.EnergyThermalOceanRadiation += oceanRadiation;
//		}

//		// TODO: track and emit heat from ice

//		float atmosphereEmissivity = Atmosphere.GetAtmosphericEmissivity(ref worldData, last.AirMass, last.AirMass * PlanetState.CarbonDioxide, last.AirWaterMass, last.CloudMass);
//		float surfaceEnergyReflected = 0;

//		// Thermal energy from surface to air, space, reflected off clouds
//		{
//			displayCell.EnergyThermalSurfaceRadiation += thermalEnergyRadiatedToAir;
//			float energyThroughAtmosphericWindow = thermalEnergyRadiatedToAir * worldData.EnergyLostThroughAtmosphereWindow;
//			thermalEnergyRadiatedToAir -= energyThroughAtmosphericWindow;
//			displayCell.EnergyThermalSurfaceOutAtmosphericWindow += energyThroughAtmosphericWindow;

//			float absorbed = thermalEnergyRadiatedToAir * atmosphereEmissivity;
//			thermalEnergyRadiatedToAir -= absorbed;
//			next.AirEnergy += absorbed;
//			displayCell.EnergyThermalAbsorbedAtmosphere += absorbed;

//			surfaceEnergyReflected = thermalEnergyRadiatedToAir * lastTerrain.CloudCoverage * worldData.ThermalReflectivityCloud;
//			reflected += surfaceEnergyReflected;
//			displayCell.EnergyThermalOutAtmosphere += thermalEnergyRadiatedToAir - surfaceEnergyReflected;
//		}

//		// atmosphere radiation
//		{
//			float energyEmitted = Atmosphere.GetRadiationRate(lastDependent.AirTemperature, atmosphereEmissivity) * worldData.SecondsPerTick;
//			next.AirEnergy -= 2 * energyEmitted;
//			backRadiation += energyEmitted;

//			float energyThroughAtmosphericWindow = energyEmitted * worldData.EnergyLostThroughAtmosphereWindow;
//			energyEmitted -= energyThroughAtmosphericWindow;

//			float energyReflected = energyEmitted * lastTerrain.CloudCoverage * worldData.ThermalReflectivityCloud;
//			reflected += energyReflected;
//			displayCell.EnergyThermalOutAtmosphere += energyEmitted - energyReflected + energyThroughAtmosphericWindow;
//		}

//		// reflected thermal radiation
//		{
//			float absorbed = reflected * atmosphereEmissivity;
//			next.AirEnergy += absorbed;
//			reflected -= absorbed;

//			displayCell.EnergyThermalAbsorbedAtmosphere += absorbed;

//			backRadiation += reflected;
//		}

//		displayCell.EnergyThermalBackRadiation += backRadiation;
//		displayCell.EnergySolarAbsorbedSurface += solarRadiationAbsorbed;
//		displayCell.Heat = solarRadiationAbsorbed;

//		float radiationToSurface = solarRadiationAbsorbed + backRadiation;

//		// ice
//		float remainingIceMass = last.IceMass;
//		{
//			// melt ice at surface from air temp and incoming radiation
//			if (remainingIceMass > 0)
//			{
//				float radiationAbsorbedByIce = radiationToSurface * iceCoverage;
//				radiationToSurface -= radiationAbsorbedByIce;
//				radiationAbsorbedByIce += thermalEnergyRadiatedToIce;

//				// world.Data.SpecificHeatIce * world.Data.MassIce == KJ required to raise one cubic meter by 1 degree
//				if (radiationAbsorbedByIce > 0)
//				{
//					// Remove the latent heat from the incoming energy
//					float iceMelted = math.min(remainingIceMass, radiationAbsorbedByIce / WorldData.LatentHeatWaterLiquid);
//					next.IceDelta -= iceMelted;
//					next.WaterDelta += iceMelted;
//					next.WaterEnergy += iceMelted * (WorldData.SpecificHeatWater * WorldData.FreezingTemperature);
//					remainingIceMass -= iceMelted;
//				}
//				else
//				{
//					next.WaterEnergy += radiationAbsorbedByIce * waterCoverage;
//					next.GroundEnergy += radiationAbsorbedByIce * (1.0f - waterCoverage);
//				}
//				if (lastDependent.AirTemperature > WorldData.FreezingTemperature)
//				{
//					// Remove the latent heat of water from the air
//					float temperatureDiff = lastDependent.AirTemperature - WorldData.FreezingTemperature;
//					float energyTransfer = math.min(last.AirEnergy, temperatureDiff * worldData.SecondsPerTick * iceCoverage * worldData.IceAirConductionCooling);
//					float iceMeltedFromConduction = remainingIceMass * math.saturate(energyTransfer / WorldData.LatentHeatWaterLiquid);
//					next.AirEnergy -= energyTransfer;
//					next.IceDelta -= iceMeltedFromConduction;
//					next.WaterDelta += iceMeltedFromConduction;
//					next.WaterEnergy += iceMeltedFromConduction * (WorldData.SpecificHeatWater * WorldData.FreezingTemperature);
//					remainingIceMass -= iceMeltedFromConduction;
//					displayCell.EnergySurfaceConduction -= energyTransfer;
//				}

//			}

//			// freeze the top meter based on surface temperature (plus incoming radiation)
//			if (iceCoverage < 1)
//			{
//				if (last.WaterMass > 0)
//				{
//					// world.Data.SpecificHeatIce * world.Data.MassIce == KJ required to raise one cubic meter by 1 degree
//					float specificHeatWater = Atmosphere.GetSpecificHeatOfWater(last.WaterMass, last.SaltMass);
//					float seaWaterHeatingRate = WorldData.MassWater / specificHeatWater;
//					//float surfaceTemp = lowerAirTemperature + incomingRadiation * seaWaterHeatingRate;
//					float localHeating = 0; // TODO: add in some local heating
//					float surfaceTemp = (lastDependent.AirTemperature + localHeating) * (1.0f - iceCoverage) + lastDependent.WaterTemperature * iceCoverage;
//					if (surfaceTemp < WorldData.FreezingTemperature)
//					{
//						float iceMassFrozen = math.min(last.WaterMass, math.min(math.max(0, worldData.FullIceCoverage * WorldData.MassIce - last.IceMass), (WorldData.FreezingTemperature - surfaceTemp) * seaWaterHeatingRate));
//						next.IceDelta += iceMassFrozen;
//						next.WaterDelta -= iceMassFrozen;
//						// TODO: shouldnt the latent heat be added to the air, not the water?
//						next.WaterEnergy -= iceMassFrozen * (WorldData.SpecificHeatWater * surfaceTemp - WorldData.LatentHeatWaterLiquid);
//					}

//					// TODO this should be using absolute pressure not barometric
//					float inverseLowerAirPressure = 1.0f / lastDependent.AirPressure;
//					// evaporation
//					if (evaporationRate > 0)
//					{
//						float evapotranspiration;
//						// TODO: absorb incoming radiation as latent heat (rather than from the surrounding air)
//						EvaporateWater(
//							waterCoverage,
//							evaporationRate,
//							last.GroundWaterDepth,
//							last.WaterMass,
//							last.SaltMass,
//							last.WaterEnergy,
//							lastDependent.WaterTemperature,
//							ref next.AirWaterMass,
//							ref next.AirEnergy,
//							ref next.WaterEnergy,
//							ref next.WaterMass,
//							out evaporation,
//							out evapotranspiration);
//						displayCell.EnergyEvapotranspiration += evapotranspiration;
//						displayCell.Evaporation = evaporation;
//					}
//				}
//			}
//		}


//		// absorbed by surface
//		{
//			// absorb the remainder and radiate heat
//			float absorbedByLand = (1.0f - waterCoverage) * radiationToSurface;
//			next.GroundEnergyDelta += absorbedByLand;
//			if (waterCoverage > 0)
//			{
//				// absorb remaining incoming radiation (we've already absorbed radiation in surface ice above)
//				float absorbedByWater = waterCoverage * radiationToSurface;
//				next.WaterEnergy += absorbedByWater;
//				//				displayCell.EnergySolarAbsorbedOcean += absorbedByWater;
//				//
//				// heat transfer (both ways) based on temperature differential
//				// conduction to ice from below
//				if (last.IceMass > 0 && lastDependent.WaterTemperature > WorldData.FreezingTemperature)
//				{
//					float oceanConductionRate = (lastDependent.WaterTemperature - WorldData.FreezingTemperature) * worldData.OceanIceConduction * worldData.SecondsPerTick * iceCoverage;

//					float energyToIce = math.max(0, oceanConductionRate * iceCoverage);
//					float iceMelted = math.min(remainingIceMass, energyToIce * worldData.inverseSpecificHeatIce);
//					next.IceDelta -= iceMelted;
//					next.WaterDelta += iceMelted;
//					next.WaterEnergy += iceMelted * (WorldData.SpecificHeatWater * WorldData.FreezingTemperature - WorldData.LatentHeatWaterLiquid);
//					remainingIceMass -= iceMelted;
//				}
//				// lose heat to air via conduction AND radiation
//				if (iceCoverage < 1)
//				{
//					// when ocean is warmer than air, it creates a convection current, which makes conduction more efficient)
//					float oceanConduction = (lastDependent.WaterTemperature - lastDependent.AirTemperature) * worldData.SecondsPerTick * (1.0f - iceCoverage) * math.min(1.0f, lastDependent.WaterDepth / worldData.WaterAirConductionDepth);
//					if (oceanConduction > 0)
//					{
//						oceanConduction *= worldData.OceanAirConductionWarming;
//					}
//					else
//					{
//						oceanConduction *= worldData.OceanAirConductionCooling;
//					}
//					next.AirEnergy += oceanConduction;
//					next.WaterEnergy -= oceanConduction;
//					displayCell.EnergyOceanConduction += oceanConduction;
//					displayCell.EnergySurfaceConduction += oceanConduction;
//				}

//				if (lastDependent.WaterTemperature < WorldData.FreezingTemperature)
//				{
//					float specificHeatSaltWater = (WorldData.SpecificHeatWater * last.WaterMass + WorldData.SpecificHeatWater * last.SaltMass);
//					float massFrozen = math.min(last.WaterMass,
//						specificHeatSaltWater * (lastDependent.WaterTemperature - WorldData.FreezingTemperature) /
//						(WorldData.LatentHeatWaterLiquid - WorldData.FreezingTemperature * (WorldData.SpecificHeatWater + WorldData.SpecificHeatIce)));

//					next.IceDelta += massFrozen;
//					next.WaterDelta -= massFrozen;
//					next.WaterEnergy -= massFrozen * (WorldData.SpecificHeatWater * WorldData.FreezingTemperature - WorldData.LatentHeatWaterLiquid);
//				}
//			}
//		}

//		Next[i] = next;
//	}

//}

//[BurstCompile]
//public struct TickAtmosphereJob : IJobParallelFor {

//	public NativeArray<CellAtmosphere> Cells;
//	public NativeArray<CellDisplay> DisplayCells;

//	[ReadOnly] public PlanetState PlanetState;
//	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;
//	[ReadOnly] public NativeArray<CellAtmosphere> LastAtmosphere;
//	[ReadOnly] public NativeArray<CellWater> LastWater;
//	[ReadOnly] public WorldData worldData;
//	[ReadOnly] public StaticState staticState;

//	public void Execute(int i)
//	{
//		var last = LastTerrain[i];
//		var lastAtmosphere = LastAtmosphere[i];
//		var lastTerrain = LastTerrain[i];
//		var lastWater = LastWater[i];
//		var next = new CellAtmosphere();
//		var display = new CellDisplay();

//		next.IceMass = last.IceMass;
//		next.WaterMass = last.WaterMass;
//		next.WaterEnergy = last.WaterEnergy;
//		next.SaltMass = last.SaltMass;
//		next.GroundEnergy = last.GroundEnergy;
//		next.GroundWater = last.GroundWater;
//		next.GroundWaterDepth = last.GroundWaterDepth;
//		next.CloudMass = last.CloudMass;
//		next.CloudDropletMass = last.CloudDropletMass;
//		next.AirWaterMass = last.AirWaterMass;
//		next.AirMass = last.AirMass;
//		next.AirEnergy = last.AirEnergy;

//		float iceCoverage = math.min(1.0f, math.pow(last.IceMass * worldData.inverseFullIceCoverage, 0.6667f));
//		float surfaceElevation = lastTerrain.Elevation + lastDependent.WaterAndIceDepth;
//		float evaporationRate = Atmosphere.GetEvaporationRate(ref worldData, iceCoverage, lastDependent.AirTemperature, lastDependent.RelativeHumidity, worldData.inverseEvapTemperatureRange);
//		float dewPoint = Atmosphere.GetDewPoint(ref worldData, lastDependent.AirTemperature, lastDependent.RelativeHumidity);

//		DoEnergyCycle(i, ref last, ref next, ref lastDependent, ref lastTerrain, ref display, surfaceElevation, evaporationRate, dewPoint, iceCoverage);
//		DoVerticalWaterMovement(i, ref last, ref lastDependent, ref next, ref display, surfaceElevation, dewPoint, evaporationRate);
//		DoDiffusion(i, ref last, ref lastDependent, ref lastTerrain, ref next);




//		Cells[i] = next;
//		DisplayCells[i] = display;
//	}


//	static private void EvaporateWater(
//		float waterCoverage,
//		float evapRate,
//		float waterTableDepth,
//		float shallowWaterMass,
//		float shallowSaltMass,
//		float shallowWaterEnergy,
//		float shallowWaterTemperature,
//		ref float newHumidity,
//		ref float newLowerAirEnergy,
//		ref float newShallowWaterEnergy,
//		ref float newShallowWaterMass,
//		out float evaporation,
//		out float evapotranspiration)
//	{
//		evaporation = 0;
//		evapotranspiration = 0;


//		if (waterCoverage > 0)
//		{
//			float evapMass = math.min(shallowWaterMass, waterCoverage * evapRate);
//			newHumidity += evapMass;
//			newShallowWaterMass -= evapMass;
//			evaporation += evapMass;
//			// this sucks energy out of the lower atmosphere since it uses up some energy to fill up the latent heat of water vapor
//			newShallowWaterEnergy -= evapMass * (WorldData.SpecificHeatWater * shallowWaterTemperature + WorldData.LatentHeatWaterVapor);
//			newLowerAirEnergy += evapMass * (WorldData.SpecificHeatWaterVapor * shallowWaterTemperature);
//			evapotranspiration = evapMass * (WorldData.LatentHeatWaterVapor + WorldData.SpecificHeatWaterVapor * shallowWaterTemperature);
//		}
//	}



//	private void DoVerticalWaterMovement(int i, ref CellState last, ref CellDependent lastDependent, ref CellState next, ref CellDisplay display, float surfaceElevation, float dewPoint, float evaporationRate)
//	{


//		// condensation
//		if (lastDependent.RelativeHumidity > 1)
//		{
//			float condensationMass = next.AirWaterMass * (lastDependent.RelativeHumidity - 1.0f) / lastDependent.RelativeHumidity;
//			next.AirWaterMass -= condensationMass;
//			next.AirEnergy -= condensationMass * (WorldData.SpecificHeatWaterVapor * lastDependent.AirTemperature - WorldData.LatentHeatWaterVapor);
//			if (lastDependent.AirTemperature <= WorldData.FreezingTemperature)
//			{
//				next.IceMass += condensationMass;
//			}
//			else
//			{
//				next.WaterMass += condensationMass;
//				next.WaterEnergy += condensationMass * (WorldData.SpecificHeatWater * lastDependent.AirTemperature);
//			}
//		}

//		if (lastDependent.WindVertical > 0)
//		{
//			float humidityToCloud = math.min(1.0f, lastDependent.WindVertical * worldData.SecondsPerTick / (lastDependent.CloudElevation - surfaceElevation)) * next.AirWaterMass * worldData.HumidityToCloudPercent;
//			next.CloudMass += humidityToCloud;
//			next.AirWaterMass -= humidityToCloud;

//			// TODO: figure out what to do about the 2 layers of atmosphere
//			// We're moving the latent heat of water vapor here since we want it to heat up the upper air around the cloud
//			next.AirEnergy -= humidityToCloud * WorldData.SpecificHeatWaterVapor * lastDependent.AirTemperature;
//			next.AirEnergy += humidityToCloud * (WorldData.SpecificHeatWater * lastDependent.AirTemperature + WorldData.LatentHeatWaterVapor);
//		}

//		if (last.CloudMass > 0)
//		{

//			// TODO: airDesntiy and rainDensity should probably be cleaned up (derived from other data?)
//			float rainDropVolume = math.max(worldData.rainDropMinSize, last.CloudDropletMass / (last.CloudMass * worldData.waterDensity));
//			float rainDropRadius = math.min(math.pow(rainDropVolume, 0.333f), worldData.rainDropMaxSize);
//			float rainDropVelocity = lastDependent.WindVertical - math.sqrt(8 * rainDropRadius * worldData.waterDensity * PlanetState.Gravity / (3 * worldData.airDensity * worldData.rainDropDragCoefficient));

//			next.CloudDropletMass = math.max(0, next.CloudDropletMass + last.CloudMass * (worldData.RainDropFormationSpeedTemperature / dewPoint * math.pow(math.max(0, -rainDropVelocity) * worldData.RainDropCoalescenceWind, 2)));

//			if (rainDropVelocity < 0 && last.CloudDropletMass > 0)
//			{
//				float rainfall = last.CloudMass * math.saturate(-rainDropVelocity * worldData.RainfallRate);

//				next.CloudMass -= rainfall;
//				next.CloudDropletMass = math.max(0, next.CloudDropletMass - rainfall / last.CloudMass);

//				{
//					float rainDropFallTime = -lastDependent.CloudElevation / rainDropVelocity;
//					// evap rate is based on full tile surface coverage, an occurs in the top millimeter
//					float rainDropSurfaceArea = 4 * math.PI * rainDropRadius * rainDropRadius;
//					float totalRainSurfaceArea = rainfall / (worldData.waterDensity * rainDropVolume) * rainDropSurfaceArea;
//					float rainEvapRate = evaporationRate * totalRainSurfaceArea * 1000 * staticState.InverseCellDiameter * staticState.InverseCellDiameter;
//					float rainDropMassToHumidity = math.min(rainfall, rainDropFallTime * rainEvapRate * worldData.TicksPerSecond);
//					rainfall -= rainDropMassToHumidity;
//					next.AirWaterMass += rainDropMassToHumidity;
//					// This sucks heat out of the lower atmosphere in the form of latent heat of water vapor
//					next.AirEnergy -= rainDropMassToHumidity * lastDependent.AirTemperature * WorldData.SpecificHeatWater;
//					next.AirEnergy += rainDropMassToHumidity * (lastDependent.AirTemperature * WorldData.SpecificHeatWaterVapor - WorldData.LatentHeatWaterVapor);
//				}
//				if (rainfall > 0)
//				{
//					display.Rainfall = rainfall;
//					next.WaterMass += rainfall;
//					// No real state change here
//					float energyTransfer = rainfall * lastDependent.AirTemperature * WorldData.SpecificHeatWater;
//					next.WaterEnergy += energyTransfer;
//					next.AirEnergy -= energyTransfer;
//				}
//			}

//			// dissapation
//			float dissapationSpeed = math.min(1.0f, worldData.CloudDissapationRateWind * math.max(0, -lastDependent.WindVertical) + worldData.CloudDissapationRateDryAir) * (1.0f - lastDependent.RelativeHumidity);
//			float dissapationMass = last.CloudMass * dissapationSpeed;
//			next.CloudDropletMass = math.max(0, next.CloudDropletMass - dissapationSpeed);
//			next.CloudMass -= dissapationMass;
//			next.AirWaterMass += dissapationMass;
//			next.AirEnergy -= dissapationMass * (lastDependent.AirTemperature * WorldData.SpecificHeatWater + WorldData.LatentHeatWaterVapor);
//			next.AirEnergy += dissapationMass * WorldData.SpecificHeatWaterVapor * lastDependent.AirTemperature;
//		}

//	}


//	//static private void SeepWaterIntoGround(World world, float groundWater, float groundEnergy, float shallowWaterMass, float soilFertility, float waterTableDepth, float shallowWaterTemperature, ref float newGroundWater, ref float newShallowWater, ref float newGroundEnergy, ref float newShallowEnergy)
//	//{
//	//	float maxGroundWater = soilFertility * waterTableDepth * world.Data.MaxSoilPorousness * world.Data.MassWater;
//	//	if (groundWater >= maxGroundWater && groundWater > 0)
//	//	{
//	//		float massTransfer = groundWater - maxGroundWater;
//	//		newShallowWater += massTransfer;
//	//		newGroundWater -= massTransfer;
//	//		float energyTransfer = massTransfer / groundWater * groundEnergy; // TODO: this isn't great, some of that ground energy is in the terrain, not just in the water
//	//		newShallowEnergy += energyTransfer;
//	//		newGroundEnergy -= energyTransfer;
//	//	}
//	//	else if (shallowWaterMass > 0)
//	//	{
//	//		float massTransfer = Mathf.Min(shallowWaterMass, Math.Min(soilFertility * world.Data.GroundWaterReplenishmentSpeed * world.Data.SecondsPerTick, maxGroundWater - groundWater));
//	//		newGroundWater += massTransfer;
//	//		newShallowWater -= massTransfer;
//	//		float energyTransfer = massTransfer * shallowWaterTemperature * world.Data.SpecificHeatWater;
//	//		newShallowEnergy -= energyTransfer;
//	//		newGroundEnergy += energyTransfer;
//	//	}
//	//}


//}
