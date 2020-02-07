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

[BurstCompile]
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

[BurstCompile]
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
			energyReflectedAtmosphere = incomingRadiation * (1.0f - 1.0f / math.exp10(reflectivity * beforeCloud));
			incomingRadiation -= energyReflectedAtmosphere;
			absorbed += incomingRadiation * (1.0f - 1.0f / math.exp10(absorptivity * beforeCloud));
			incomingRadiation -= absorbed;
		}


		if (isCloudLayer)
		{
			float cloudIceContent = math.saturate((CloudTemperature[i] - worldData.minCloudFreezingTemperature) / (worldData.maxCloudFreezingTemperature - worldData.minCloudFreezingTemperature));
			float cloudTemperatureAlbedo = WorldData.AlbedoIce + (WorldData.AlbedoWater - WorldData.AlbedoIce) * cloudIceContent;
			float rainDropSizeAlbedo = math.saturate(1.0f - CloudDropletMass[i] / CloudMass[i]) * (worldData.rainDropSizeAlbedoMax - worldData.rainDropSizeAlbedoMin) + worldData.rainDropSizeAlbedoMin;
			float cloudAlbedo = math.min(1.0f, worldData.SolarReflectivityCloud * cloudTemperatureAlbedo * CloudMass[i] * rainDropSizeAlbedo / math.max(worldData.maxCloudSlopeAlbedo, 1.0f - WaterSlopeAlbedo[i]));

			energyReflectedClouds = incomingRadiation * (1.0f - 1.0f / math.exp10(cloudAlbedo * cloudMass));
			incomingRadiation -= energyReflectedClouds;
			absorbedByCloudsIncoming = incomingRadiation * (1.0f - 1.0f / math.exp10(SolarAbsorptivityCloud * cloudMass));
			incomingRadiation -= absorbedByCloudsIncoming;

			SolarRadiationReflectedCloud[i] = energyReflectedClouds;
			SolarRadiationAbsorbedCloud[i] = absorbedByCloudsIncoming;
		}

		// below cloud
		if (afterCloud > 0)
		{
			float reflectedBelow = incomingRadiation * (1.0f - 1.0f / math.exp10(reflectivity * afterCloud));
			incomingRadiation -= reflectedBelow;
			energyReflectedAtmosphere += reflectedBelow;
			float absorbedBelow = incomingRadiation * (1.0f - 1.0f / math.exp10(absorptivity * afterCloud));
			absorbed += absorbedBelow;
			incomingRadiation -= absorbedBelow;
		}

		SolarRadiationAbsorbed[i] = absorbed;
		SolarRadiationIncoming[i] = incomingRadiation;
		SolarRadiationReflected[i] = energyReflectedAtmosphere;
	}
}

[BurstCompile]
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
[BurstCompile]
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
[BurstCompile]
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


[BurstCompile]
public struct EmissivityAirJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public float GreenhouseGasConcentration;
	[ReadOnly] public float AbsorptivityAir;
	[ReadOnly] public float AbsorptivityCarbonDioxide;
	[ReadOnly] public float AbsorptivityWaterVapor;
	public void Execute(int i)
	{
		Emissivity[i] = (AirMass[i] * WorldData.EmissivityAir +	VaporMass[i] * WorldData.EmissivityWater) / (AirMass[i] + VaporMass[i]);
	}
}

[BurstCompile]
public struct EmissivityWaterJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> WaterSaltMass;
	[ReadOnly] public NativeArray<float> WaterMass;
	public void Execute(int i)
	{
		Emissivity[i] = math.saturate((WorldData.EmissivityWater * WaterMass[i] + WorldData.EmissivitySalt * WaterSaltMass[i]) / (WaterMass[i] + WaterSaltMass[i]));
	}
}

[BurstCompile]
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


[BurstCompile]
public struct ThermalEnergyRadiatedJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> ThermalRadiationTransmittedDown;
	public NativeArray<float> WindowRadiationTransmittedUp;
	public NativeArray<float> WindowRadiationTransmittedDown;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		// radiate half up and half down
		float transmittedUp = math.min(Energy[i] / 2, Atmosphere.GetRadiationRate(Temperature[i], Emissivity[i]) * SecondsPerTick);
		ThermalRadiationDelta[i] = -2 * transmittedUp;

		float windowTransmittedUp = transmittedUp * PercentRadiationInAtmosphericWindow;
		transmittedUp -= windowTransmittedUp;

		WindowRadiationTransmittedUp[i] = windowTransmittedUp;
		WindowRadiationTransmittedDown[i] = windowTransmittedUp;
		ThermalRadiationTransmittedUp[i] = transmittedUp;
		ThermalRadiationTransmittedDown[i] = transmittedUp;
	}
}
[BurstCompile]
public struct ThermalEnergyRadiatedConstantEmissivityJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> ThermalRadiationTransmittedDown;
	public NativeArray<float> WindowRadiationTransmittedUp;
	public NativeArray<float> WindowRadiationTransmittedDown;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public float Emissivity;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		// radiate half up and half down
		float transmittedUp = math.min(Energy[i] / 2,  Atmosphere.GetRadiationRate(Temperature[i], Emissivity) * SecondsPerTick);
		ThermalRadiationDelta[i] = -2 * transmittedUp;

		float windowTransmittedUp = transmittedUp * PercentRadiationInAtmosphericWindow;
		transmittedUp -= windowTransmittedUp;

		WindowRadiationTransmittedUp[i] = windowTransmittedUp;
		WindowRadiationTransmittedDown[i] = windowTransmittedUp;
		ThermalRadiationTransmittedUp[i] = transmittedUp;
		ThermalRadiationTransmittedDown[i] = transmittedUp;
	}
}
[BurstCompile]
public struct ThermalEnergyRadiatedTerrainJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationDelta;
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> WindowRadiationTransmitted;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float PercentRadiationInAtmosphericWindow;
	public void Execute(int i)
	{
		// radiate half up and half down
		float emitted = math.min(Energy[i] / 2, Atmosphere.GetRadiationRate(Temperature[i], Emissivity[i]) * SecondsPerTick);
		ThermalRadiationDelta[i] = -emitted;

		float emittedOutAtmosphericWindow = emitted * PercentRadiationInAtmosphericWindow;
		emitted -= emittedOutAtmosphericWindow;
		WindowRadiationTransmitted[i] = emittedOutAtmosphericWindow;
		ThermalRadiationTransmitted[i] = emitted;
	}
}

[BurstCompile]
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

		float absorptivity = (1.0f - CarbonDioxide) * AirMass[i] * AirAbsorptivity + VaporMass[i] * VaporAbsorptivity + AirMass[i] * CarbonDioxide * CarbonAbsorptivity;

		float cloudMass = CloudMass[i];
		float cloudElevation = CloudElevation[i];
		float layerElevation = LayerElevation[LayerIndex];
		float layerHeight = LayerHeight[LayerIndex];
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
			float absorbedBeforeCloud = incoming * (1.0f - 1.0f / math.exp10(absorptivity * beforeCloud));
			incoming -= absorbedBeforeCloud;
			ThermalRadiationDelta[i] += absorbedBeforeCloud;
		}

		if (cloudMass > 0)
		{
			float absorbedByCloud = (incoming + transmitting) * (1.0f - 1.0f / math.exp10(WaterAbsorptivity * cloudMass));
			ThermalRadiationDeltaCloud[i] += absorbedByCloud;
			incoming += ThermalRadiationTransmittedCloud[i] - absorbedByCloud;
		}

		if (afterCloud > 0)
		{
			float absorbedAfterCloud = incoming * (1.0f - 1.0f / math.exp10(absorptivity * afterCloud));
			ThermalRadiationDelta[i] += absorbedAfterCloud;
			incoming -= absorbedAfterCloud;
		}
		ThermalRadiationTransmitted[i] = transmitting + incoming;
	}
}


[BurstCompile]
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

[BurstCompile]
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


[BurstCompile]
public struct ThermalEnergyAbsorbedTerrainJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationAbsorbed;
	[ReadOnly] public NativeArray<float> WindowRadiationIncoming;
	[ReadOnly] public NativeArray<float> ThermalRadiationIncoming;
	public void Execute(int i)
	{
		ThermalRadiationAbsorbed[i] += ThermalRadiationIncoming[i] + WindowRadiationIncoming[i];
	}
}

[BurstCompile]
public struct DiffusionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public float DiffusionCoefficientTemperature;
	[ReadOnly] public float DiffusionCoefficientHumidity;
	[ReadOnly] public float DiffusionCoefficientVelocity;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Humidity;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		float gradientTemperature = 0;
		float gradientWaterVapor = 0;
		float2 gradientVelocity = float2.zero;
		int neighborCount = 0;
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				gradientWaterVapor += Humidity[n];
				gradientVelocity += Velocity[n];
				gradientTemperature += Temperature[n];
				neighborCount++;
			}
		}
		gradientTemperature -= Temperature[i] * neighborCount;
		gradientWaterVapor -= Humidity[i] * neighborCount;
		gradientVelocity -= Velocity[i] * neighborCount;

		Delta[i] = new DiffusionAir()
		{
			Temperature = gradientTemperature * DiffusionCoefficientTemperature,
			Humidity = gradientWaterVapor * DiffusionCoefficientHumidity,
			Velocity = gradientVelocity * DiffusionCoefficientVelocity
		};
	}
}

[BurstCompile]
public struct DiffusionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public float DiffusionCoefficientTemperature;
	[ReadOnly] public float DiffusionCoefficientSalinity;
	[ReadOnly] public float DiffusionCoefficientVelocity;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{
		float gradientSalinity = 0;
		float gradientTemperature = 0;
		float2 gradientVelocity = float2.zero;
		int neighborCount = 0;
		for (int j = 0; j < 6; j++)
		{
			int n = Neighbors[i * 6 + j];
			if (n >= 0)
			{
				gradientSalinity += Salt[n];
				gradientTemperature += Temperature[n];
				gradientVelocity += Velocity[n];

				neighborCount++;
			}
		}
		gradientSalinity -= Salt[i] * neighborCount;
		gradientTemperature -= Temperature[i] * neighborCount;
		gradientVelocity -= Velocity[i] * neighborCount;

		Delta[i] = new DiffusionWater()
		{
			Temperature = gradientTemperature * DiffusionCoefficientTemperature,
			Salinity = gradientSalinity * DiffusionCoefficientSalinity,
			Velocity = gradientVelocity * DiffusionCoefficientVelocity
		};

	}
}

[BurstCompile]
public struct AdvectionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{

	}
}

[BurstCompile]
public struct AdvectionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{

	}
}

[BurstCompile]
public struct ConductionJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> EnergyB;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		EnergyDelta[i] = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick, -EnergyA[i], EnergyB[i]);
	}
}

[BurstCompile]
public struct ConductionPartialJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> EnergyB;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		EnergyDelta[i] = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * Coverage[i], -EnergyA[i], EnergyB[i]);
	}
}


[BurstCompile]
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
		EnergyDelta[i] = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * math.min(CoverageIce[i], CoverageWater[i]), -EnergyA[i], EnergyB[i]);
	}
}

[BurstCompile]
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
		EnergyDelta[i] = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * math.min(CoverageA[i], CoverageB[i]), -EnergyA[i], EnergyB[i]);
	}
}

[BurstCompile]
public struct ConductionIceTerrainJob : IJobParallelFor {
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
		EnergyDelta[i] = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * math.max(CoverageIce[i] - CoverageWater[i], 0), -EnergyA[i], EnergyB[i]);
	}
}

[BurstCompile]
public struct ConductionAirTerrainJob : IJobParallelFor {
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
		EnergyDelta[i] = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * (1.0f - math.max(CoverageIce[i], CoverageWater[i])), -EnergyA[i], EnergyB[i]);
	}
}

[BurstCompile]
public struct ConductionWaterTerrainJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	public NativeArray<float> EnergyDeltaWaterTotal;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> EnergyB;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public NativeArray<float> CoverageBelow;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float delta = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * Coverage[i] * (1.0f - CoverageBelow[i]), -EnergyA[i], EnergyB[i]);
		EnergyDelta[i] = delta;
		EnergyDeltaWaterTotal[i] += delta;
	}
}

[BurstCompile]
public struct ConductionWaterBottomJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	public NativeArray<float> EnergyDeltaWaterTotal;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> EnergyA;
	[ReadOnly] public NativeArray<float> EnergyB;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public float ConductionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float delta = math.clamp((TemperatureB[i] - TemperatureA[i]) * ConductionCoefficient * SecondsPerTick * Coverage[i], -EnergyA[i], EnergyB[i]);
		EnergyDelta[i] = delta;
		EnergyDeltaWaterTotal[i] += delta;
	}
}



[BurstCompile]
public struct StateChangeJob : IJobParallelFor {
	public NativeArray<float> SurfaceAirEnergy;
	public NativeArray<float> SurfaceAirVapor;
	public NativeArray<float> SurfaceWaterMass;
	public NativeArray<float> SurfaceWaterEnergy;
	public NativeArray<float> IceMass;
	public NativeArray<float> IceEnergy;
	[ReadOnly] public NativeArray<float> IceMeltedMass;
	[ReadOnly] public NativeArray<float> IceMeltedEnergy;
	[ReadOnly] public NativeArray<float> IceMeltedLatentHeatAir;
	[ReadOnly] public NativeArray<float> WaterEvaporatedMass;
	[ReadOnly] public NativeArray<float> WaterEvaporatedEnergy;
	[ReadOnly] public NativeArray<float> WaterFrozenMass;
	[ReadOnly] public NativeArray<float> WaterFrozenEnergy;
	[ReadOnly] public NativeArray<float> RainfallWaterMass;
	[ReadOnly] public NativeArray<float> RainfallWaterEnergy;
	[ReadOnly] public NativeArray<float> RainfallIceMass;
	[ReadOnly] public NativeArray<float> RainfallIceEnergy;
	public void Execute(int i)
	{
		SurfaceAirEnergy[i] += WaterEvaporatedEnergy[i] + IceMeltedLatentHeatAir[i];
		SurfaceAirVapor[i] += WaterEvaporatedMass[i];
		SurfaceWaterMass[i] += RainfallWaterMass[i] + IceMeltedMass[i];
		SurfaceWaterEnergy[i] += RainfallWaterEnergy[i] + IceMeltedEnergy[i];
		IceMass[i] += RainfallIceMass[i] + WaterFrozenMass[i];
		IceEnergy[i] += RainfallIceEnergy[i] + WaterFrozenEnergy[i];
	}
}

[BurstCompile]
public struct StateChangeAirLayerJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudEnergy;
	public NativeArray<float> CloudDropletMass;
	public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> GroundCondensationMass;
	[ReadOnly] public NativeArray<float> GroundCondensationEnergy;
	[ReadOnly] public NativeArray<float> CloudCondensationMass;
	[ReadOnly] public NativeArray<float> CloudCondensationEnergy;
	[ReadOnly] public NativeArray<float> CloudEvaporationMass;
	[ReadOnly] public NativeArray<float> CloudEvaporationEnergy;
	[ReadOnly] public int LayerIndex;
	public void Execute(int i)
	{

	}
}



[BurstCompile]
public struct EnergySurfaceAirJob : IJobParallelFor {
	public NativeArray<float> Energy;
	public NativeArray<float> Vapor;
	public NativeArray<float2> Velocity;
	public NativeArray<float> CondensationGroundMass;
	public NativeArray<float> CondensationGroundEnergy;
	public NativeArray<float> CondensationCloudMass;
	public NativeArray<float> CondensationCloudEnergy;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastEnergy;
	[ReadOnly] public NativeArray<float2> LastVelocity;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	[ReadOnly] public NativeArray<DiffusionAir> Diffusion;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyCloud;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyCloud[i]
			+ ConductionEnergyIce[i]
			+ ConductionEnergyTerrain[i]
			+ ConductionEnergyWater[i];
		float temperatureToEnergy = WorldData.SpecificHeatAtmosphere * AirMass[i] + WorldData.SpecificHeatWaterVapor * LastVapor[i];

		// TODO: account for differing Air masses
		float energy =
			LastEnergy[i] +
			SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta +
			(Diffusion[i].Temperature - Advection[i].Temperature) * temperatureToEnergy;
		float vapor = LastVapor[i] + Diffusion[i].Humidity - Advection[i].Humidity;
		float2 velocity = LastVelocity[i] + Diffusion[i].Velocity - Advection[i].Velocity;



		float condensationGroundMass = 0;
		float condensationGroundEnergy = 0;
		float condensationCloudMass = 0;
		float condensationCloudEnergy = 0;






		CondensationGroundMass[i] = condensationGroundMass;
		CondensationGroundEnergy[i] = condensationGroundEnergy;
		CondensationCloudMass[i] = condensationCloudMass;
		CondensationCloudEnergy[i] = condensationCloudEnergy;

		Energy[i] = energy;
		Vapor[i] = vapor;
		Velocity[i] = velocity;
	}
}

[BurstCompile]
public struct EnergyAirJob : IJobParallelFor {
	public NativeArray<float> Energy;
	public NativeArray<float> Vapor;
	public NativeArray<float2> Velocity;
	public NativeArray<float> CondensationGroundMass;
	public NativeArray<float> CondensationGroundEnergy;
	public NativeArray<float> CondensationCloudMass;
	public NativeArray<float> CondensationCloudEnergy;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> LastVapor;
	[ReadOnly] public NativeArray<float> LastEnergy;
	[ReadOnly] public NativeArray<float2> LastVelocity;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	[ReadOnly] public NativeArray<DiffusionAir> Diffusion;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyCloud;
	public void Execute(int i)
	{
		float conductionDelta =	-ConductionEnergyCloud[i];
		float temperatureToEnergy = WorldData.SpecificHeatAtmosphere * AirMass[i] + WorldData.SpecificHeatWaterVapor * LastVapor[i];

		// TODO: account for differing Air masses
		float energy = LastEnergy[i] +
			SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta +
			(Diffusion[i].Temperature - Advection[i].Temperature) * temperatureToEnergy;
		float vapor = LastVapor[i] + Diffusion[i].Humidity - Advection[i].Humidity;
		float2 velocity = LastVelocity[i] + Diffusion[i].Velocity - Advection[i].Velocity;

		float condensationGroundMass = 0;
		float condensationGroundEnergy = 0;
		float condensationCloudMass = 0;
		float condensationCloudEnergy = 0;




		CondensationGroundMass[i] = condensationGroundMass;
		CondensationGroundEnergy[i] = condensationGroundEnergy;
		CondensationCloudMass[i] = condensationCloudMass;
		CondensationCloudEnergy[i] = condensationCloudEnergy;

		Energy[i] = energy;
		Vapor[i] = vapor;
		Velocity[i] = velocity;
	}
}

[BurstCompile]
public struct EnergyWaterJob : IJobParallelFor {
	public NativeArray<float> Energy;
	public NativeArray<float> SaltMass;
	public NativeArray<float2> Velocity;
	public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> LastEnergy;
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
		float temperatureToEnergy = (Diffusion[i].Temperature - Advection[i].Temperature) * WorldData.SpecificHeatWater * LastMass[i] + WorldData.SpecificHeatSalt * LastSaltMass[i];
		Energy[i] = LastEnergy[i] + temperatureToEnergy + SolarRadiationIn[i] + ThermalRadiationDelta[i] + conductionDelta;
		Velocity[i] = LastVelocity[i] + Diffusion[i].Velocity - Advection[i].Velocity;
		SaltMass[i] = LastSaltMass[i];
		Mass[i] = LastMass[i];
	}
}
[BurstCompile]
public struct EnergyWaterJobSurface : IJobParallelFor {
	public NativeArray<float> Energy;
	public NativeArray<float> SaltMass;
	public NativeArray<float2> Velocity;
	public NativeArray<float> WaterMass;
	public NativeArray<float> EvaporatedWaterMass;
	public NativeArray<float> EvaporatedWaterEnergy;
	public NativeArray<float> FrozenWaterMass;
	public NativeArray<float> FrozenWaterEnergy;
	public NativeArray<float> Evapotranspiration;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastEnergy;
	[ReadOnly] public NativeArray<float> LastSaltMass;
	[ReadOnly] public NativeArray<float2> LastVelocity;
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
	public void Execute(int i)
	{
		float energyTop = -ConductionEnergyAir[i] - ConductionEnergyIce[i] + ThermalRadiationDeltaTop[i] + SolarRadiationIn[i];
		float energyBottom = ConductionEnergyTerrain[i] + ThermalRadiationDeltaBottom[i];

		float temperatureToEnergy = (Diffusion[i].Temperature - Advection[i].Temperature) * WorldData.SpecificHeatWater * LastMass[i] + WorldData.SpecificHeatSalt * LastSaltMass[i];

		float energy = LastEnergy[i] + temperatureToEnergy + energyTop + energyBottom;
		float saltMass = LastSaltMass[i] + Diffusion[i].Salinity - Advection[i].Salinity;
		float2 velocity = LastVelocity[i] + Diffusion[i].Velocity - Advection[i].Velocity;
		float waterMass = LastMass[i];

		float evapMass = 0;
		float frozenMass = 0;
		float frozenEnergy = 0;
		float evapEnergy = 0;
		float evapotranspiration = 0;
		if (waterMass > 0)
		{

			{
				float evapRate = math.saturate((1.0f - RelativeHumidity[i]) * math.max(0, WaterCoverage[i] - IceCoverage[i]) * EvaporationRate * (LastTemperature[i] - EvapTemperatureMin) / (EvapTemperatureMax - EvapTemperatureMin));
				evapMass = math.clamp(evapRate * energyTop / WorldData.LatentHeatWaterVapor, 0, waterMass);

				float freshWaterEnergy = energy * (WorldData.SpecificHeatWater * waterMass) / (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatSalt * saltMass);
				evapEnergy = freshWaterEnergy * math.saturate(evapMass / (waterMass + saltMass));
				energy -= evapEnergy;
				waterMass -= evapMass;
				evapotranspiration = WorldData.LatentHeatWaterVapor * evapMass;
				evapEnergy -= evapotranspiration;

			}


			{
				float specificHeatSaltWater = (WorldData.SpecificHeatWater * waterMass + WorldData.SpecificHeatWater * saltMass);
				float waterTemperature = energy / specificHeatSaltWater;
				frozenMass = math.clamp(specificHeatSaltWater * (waterTemperature - WorldData.FreezingTemperature) /
					(WorldData.LatentHeatWaterLiquid - WorldData.FreezingTemperature * (WorldData.SpecificHeatWater + WorldData.SpecificHeatIce)),
					0, waterMass);

				frozenEnergy = frozenMass * WorldData.SpecificHeatIce * WorldData.FreezingTemperature;
				energy -= frozenEnergy;
				waterMass -= frozenMass;
				energy += WorldData.LatentHeatWaterLiquid * frozenMass;
			}


		}
		EvaporatedWaterMass[i] = evapMass;
		EvaporatedWaterEnergy[i] = evapEnergy;
		FrozenWaterMass[i] = frozenMass;
		FrozenWaterEnergy[i] = frozenEnergy;
		Evapotranspiration[i] = evapotranspiration;

		Energy[i] = energy;
		SaltMass[i] = saltMass;
		Velocity[i] = velocity;
		WaterMass[i] = waterMass;

	}
}

[BurstCompile]
public struct EnergyCloudJob : IJobParallelFor {
	public NativeArray<float> Energy;
	public NativeArray<float> Mass;
	public NativeArray<float> Elevation;
	public NativeArray<float> DropletMass;
	public NativeArray<float2> Velocity;
	public NativeArray<float> RainfallWaterMass;
	public NativeArray<float> RainfallWaterEnergy;
	public NativeArray<float> RainfallIceMass;
	public NativeArray<float> RainfallIceEnergy;
	public NativeArray<float> EvaporationMass;
	public NativeArray<float> EvaporationEnergy;
	[ReadOnly] public NativeArray<DiffusionCloud> Advection;
	[ReadOnly] public NativeArray<DiffusionCloud> Diffusion;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> LastEnergy;
	[ReadOnly] public NativeArray<float2> LastVelocity;
	[ReadOnly] public NativeArray<float> LastElevation;
	[ReadOnly] public NativeArray<float> LastDropletMass;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> WindVertical;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float WaterDensity;
	[ReadOnly] public float AirDensity;
	[ReadOnly] public float RainDropMinSize;
	[ReadOnly] public float RainDropMaxSize;
	[ReadOnly] public float RainDropDragCoefficient;
	[ReadOnly] public float RainfallRate;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float TicksPerSecond;
	[ReadOnly] public float CloudDissapationRateDryAir;
	[ReadOnly] public float CloudDissapationRateWind;
	[ReadOnly] public float RelativeHumidityAtCloudElevation;
	public void Execute(int i)
	{

		float temperatureToEnergy = WorldData.SpecificHeatWater * LastMass[i];
		float energy = LastEnergy[i] + (Diffusion[i].Temperature - Advection[i].Temperature) / temperatureToEnergy + SolarRadiationIn[i] + ThermalRadiationDelta[i] + ConductionEnergyAir[i];
		float dropletMass = LastDropletMass[i] + Diffusion[i].DropletMass - Advection[i].DropletMass;
		float cloudMass = LastMass[i] + Diffusion[i].Mass - Advection[i].Mass;
		float cloudElevation = LastElevation[i] + Diffusion[i].Elevation - Advection[i].Elevation;
		float2 velocity = LastVelocity[i] + Diffusion[i].Velocity - Advection[i].Velocity;

		float rainfallWaterMass = 0;
		float rainfallWaterEnergy = 0;
		float rainfallIceMass = 0;
		float rainfallIceEnergy = 0;
		float evaporationMass = 0;
		float evaporationEnergy = 0;

		// TODO: airDesntiy and rainDensity should probably be cleaned up (derived from other data?)
		float rainDropVolume = math.max(RainDropMinSize, dropletMass / (cloudMass * WaterDensity));
		float rainDropRadius = math.min(math.pow(rainDropVolume, 0.333f), RainDropMaxSize);
		float rainDropVelocity = WindVertical[i] - math.sqrt(8 * rainDropRadius * WaterDensity * Gravity / (3 * AirDensity * RainDropDragCoefficient));
		if (rainDropVelocity < 0 && dropletMass > 0)
		{
			rainfallWaterMass = cloudMass * math.saturate(-rainDropVelocity * RainfallRate);
			rainfallWaterEnergy = rainfallWaterMass / cloudMass * energy;
			cloudMass -= rainfallWaterMass;
			energy -= rainfallWaterEnergy;

			float rainDropFallTime = -cloudElevation / rainDropVelocity;

			// evap rate is based on full tile surface coverage, an occurs in the top millimeter
			float evaporationRate = 0; // TODO: fill this in
			float rainDropSurfaceArea = 4 * math.PI * rainDropRadius * rainDropRadius;
			float totalRainSurfaceArea = rainfallWaterMass / (WaterDensity * rainDropVolume) * rainDropSurfaceArea;
			float rainEvapRate = evaporationRate * totalRainSurfaceArea * 1000 * InverseCellDiameter * InverseCellDiameter;
			evaporationMass = math.min(rainfallWaterMass, rainDropFallTime * rainEvapRate * TicksPerSecond);
			evaporationEnergy = evaporationMass / rainfallWaterMass * rainfallWaterEnergy;
			rainfallWaterMass -= evaporationMass;
			rainfallWaterEnergy -= evaporationEnergy;

			// This sucks heat out of the lower atmosphere in the form of latent heat of water vapor
			evaporationEnergy -= evaporationMass * WorldData.LatentHeatWaterVapor;
		}

		// dissapation
		float dissapationSpeed = math.min(1.0f, CloudDissapationRateWind * math.max(0, -WindVertical[i]) + CloudDissapationRateDryAir) * (1.0f - RelativeHumidityAtCloudElevation);
		float dissapationMass = cloudMass * dissapationSpeed;
		dropletMass = math.max(0, dropletMass - dissapationSpeed);
		float dissapationEnergy = dissapationMass / cloudMass * energy;
		energy -= dissapationEnergy;
		cloudMass -= dissapationMass;
		evaporationMass += dissapationMass;
		evaporationEnergy += dissapationEnergy - dissapationMass * WorldData.LatentHeatWaterVapor;


		Energy[i] = energy;
		DropletMass[i] = dropletMass;
		Mass[i] = cloudMass;
		Elevation[i] = cloudElevation;
		Velocity[i] = velocity;


		RainfallWaterMass[i] = rainfallWaterMass;
		RainfallWaterEnergy[i] = rainfallWaterEnergy;
		RainfallIceMass[i] = rainfallIceMass;
		RainfallIceEnergy[i] = rainfallIceEnergy;
		EvaporationMass[i] = evaporationMass;
		EvaporationEnergy[i] = evaporationEnergy;

	}
}

[BurstCompile]
public struct EnergyIceJob : IJobParallelFor {
	public NativeArray<float> Energy;
	public NativeArray<float> Mass;
	public NativeArray<float> MeltedMass;
	public NativeArray<float> MeltedEnergy;
	public NativeArray<float> LatentHeatAir;
	[ReadOnly] public NativeArray<float> LastMass;
	[ReadOnly] public NativeArray<float> LastEnergy;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaBottom;
	[ReadOnly] public NativeArray<float> ThermalRadiationDeltaTop;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyTerrain;
	public void Execute(int i)
	{
		float energyTop = -ConductionEnergyAir[i] + ThermalRadiationDeltaTop[i] + SolarRadiationIn[i];
		float energyBottom = ConductionEnergyWater[i] + ConductionEnergyTerrain[i] + ThermalRadiationDeltaBottom[i];

		float meltedMass = 0;
		float meltedEnergy = 0;
		float iceMass = LastMass[i];
		float iceEnergy = LastEnergy[i];
		float latentHeatAir = 0;
		if (iceMass > 0)
		{
			float maxHeatingDepth = 1;
			float iceDepth = iceMass / WorldData.MassIce;
			float heatingDepth = math.min(maxHeatingDepth, iceDepth);
			float heatingMass = heatingDepth * WorldData.MassIce;

			float newEnergyTop = (iceEnergy * heatingDepth / iceDepth + energyTop);
			float newTempTop = newEnergyTop * WorldData.SpecificHeatIce / heatingMass;

			meltedMass = iceMass * math.saturate(WorldData.SpecificHeatIce * (newTempTop - WorldData.FreezingTemperature) / (WorldData.SpecificHeatWater * WorldData.FreezingTemperature - WorldData.SpecificHeatIce * WorldData.FreezingTemperature));
			meltedEnergy = WorldData.SpecificHeatWater * meltedMass * WorldData.FreezingTemperature;
			energyTop -= meltedEnergy;
			iceMass -= meltedMass;
			latentHeatAir -= meltedMass * WorldData.LatentHeatWaterLiquid;

		}

		iceEnergy += energyTop;
		iceEnergy += energyBottom;

		Mass[i] = iceMass;
		Energy[i] = iceEnergy;
		MeltedMass[i] = meltedMass;
		MeltedEnergy[i] = meltedEnergy;
		LatentHeatAir[i] = latentHeatAir;
	}
}

[BurstCompile]
public struct EnergyTerrainJob : IJobParallelFor {
	public NativeArray<float> Energy;
	[ReadOnly] public NativeArray<float> LastEnergy;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> ThermalRadiationDelta;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergyAir;
	[ReadOnly] public NativeArray<float> ConductionEnergyWater;
	[ReadOnly] public NativeArray<float> ConductionEnergyIce;
	[ReadOnly] public float GeothermalEnergy;
	public void Execute(int i)
	{
		float conductionDelta =
			-ConductionEnergyAir[i]
			-ConductionEnergyWater[i]
			-ConductionEnergyIce[i];
		Energy[i] = 
			LastEnergy[i] +
			SolarRadiationIn[i]  +
			ThermalRadiationDelta[i] +
			conductionDelta +
			GeothermalEnergy
			;
	}
}

[BurstCompile]
public struct UpdateTerrainJob : IJobParallelFor {
	public NativeArray<CellTerrain> Terrain;
	public NativeArray<float> GroundWater;

	[ReadOnly] public NativeArray<CellTerrain> LastTerrain;
	[ReadOnly] public NativeArray<float> LastGroundWater;
	public void Execute(int i)
	{
		Terrain[i] = LastTerrain[i];
		GroundWater[i] = LastGroundWater[i];
	}

}

[BurstCompile]
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


[BurstCompile]
public struct UpdateDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceElevation;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> CloudCoverage;
	public NativeArray<float> WaterDepth;
	public NativeArray<float> CloudTemperature;
	public NativeArray<float> IceTemperature;
	public NativeArray<float> TerrainTemperature;
	public NativeArray<float> SurfaceAirTemperature;
	[ReadOnly] public NativeArray<WaterSaltMass> WaterSaltMass;
	[ReadOnly] public NativeArray<float> LowerAirTemperature;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> IceEnergy;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudEnergy;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> TerrainEnergy;
	[ReadOnly] public NativeArray<float> GroundWater;
	[ReadOnly] public WorldData worldData;
	[ReadOnly] public NativeArray<float> lowerAirHeight;
	public void Execute(int i)
	{
		float waterDepth = WaterSaltMass[i].WaterMass / WorldData.MassWater + WaterSaltMass[i].SaltMass / WorldData.MassSalt;
		float iceMass = IceMass[i];
		WaterDepth[i] = waterDepth;
		SurfaceElevation[i] = Terrain[i].Elevation + waterDepth + iceMass / WorldData.MassIce;
		VegetationCoverage[i] = Terrain[i].Vegetation * worldData.inverseFullCanopyCoverage;

		IceCoverage[i] = iceMass * worldData.inverseFullIceCoverage;
		if (iceMass == 0)
		{
			IceTemperature[i] = 0;
		} else
		{
			IceTemperature[i] = IceEnergy[i] / (iceMass * WorldData.SpecificHeatIce);
		}

		float cloudMass = CloudMass[i];
		CloudCoverage[i] = cloudMass * worldData.inverseCloudMassFullAbsorption;
		if (cloudMass == 0)
		{
			CloudTemperature[i] = 0;
		} else
		{
			CloudTemperature[i] = CloudEnergy[i] / (cloudMass * WorldData.SpecificHeatWater);
		}

		SurfaceAirTemperature[i] = LowerAirTemperature[i] - WorldData.TemperatureLapseRate * lowerAirHeight[i] / 2;
		TerrainTemperature[i] = Atmosphere.GetTerrainTemperature(ref worldData, TerrainEnergy[i], GroundWater[i], Terrain[i].SoilFertility, VegetationCoverage[i]);
	}

}

[BurstCompile]
public struct UpdateDependentAirLayerJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> AbsoluteHumidity;
	public NativeArray<float> RelativeHumidity;
	public NativeArray<float> Pressure;
	public NativeArray<float> LayerElevation;
	public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirVapor;
	[ReadOnly] public NativeArray<float> AirEnergy;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public WorldData WorldData;
	[ReadOnly] public float Gravity;

	[ReadOnly] public NativeArray<float> IceMass;
	public void Execute(int i)
	{
		float layerMiddle = LayerElevation[i] + LayerHeight[i] / 2;

		Temperature[i] = AirEnergy[i] / (WorldData.SpecificHeatAtmosphere * AirMass[i] + WorldData.SpecificHeatWaterVapor * VaporMass[i]);
		AbsoluteHumidity[i] = AirVapor[i] / (AirVapor[i] + AirMass[i]);
		RelativeHumidity[i] = Atmosphere.GetRelativeHumidity(ref WorldData, Temperature[i], AirVapor[i], AirMass[i], WorldData.inverseDewPointTemperatureRange);
		Pressure[i] = Atmosphere.GetAirPressure(ref WorldData, AirMass[i], LayerElevation[i], Temperature[i], Atmosphere.GetMolarMassAir(AirMass[i], VaporMass[i]), Gravity);
	}
}

[BurstCompile]
public struct UpdateDependentWaterLayerJob : IJobParallelFor {
	public NativeArray<float> Salinity;
	public NativeArray<float> Temperature;
	public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> WaterMass;
	[ReadOnly] public NativeArray<float> SaltMass;
	[ReadOnly] public NativeArray<float> Energy;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			float waterDepth = WaterMass[i] / WorldData.MassWater + SaltMass[i] / WorldData.MassSalt;
			WaterCoverage[i] = math.min(1, waterDepth / math.max(1, Terrain[i].Roughness));
			Temperature[i] = Energy[i] / (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt);
			Salinity[i] = SaltMass[i] / (WaterMass[i] + SaltMass[i]);
		}
		else
		{
			WaterCoverage[i] = 0;
			Temperature[i] = 0;
			Salinity[i] = 0;
		}
	}
}

[BurstCompile]
public struct UpdateDisplayJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbedSurface;
	public NativeArray<float> DisplayRainfall;
	public NativeArray<float> DisplayEvaporation;
	[ReadOnly] public NativeArray<float> SolarRadiationInTerrain;
	[ReadOnly] public NativeArray<float> SolarRadiationInIce;
	[ReadOnly] public NativeArray<float> SolarRadiationInWaterSurface;
	[ReadOnly] public NativeArray<float> RainfallWater;
	[ReadOnly] public NativeArray<float> RainfallIce;
	[ReadOnly] public NativeArray<float> Evaporation;
	public void Execute(int i)
	{
		SolarRadiationAbsorbedSurface[i] = SolarRadiationInTerrain[i] + SolarRadiationInIce[i] + SolarRadiationInWaterSurface[i];
		DisplayRainfall[i] = RainfallWater[i] + RainfallIce[i];
		DisplayEvaporation[i] = Evaporation[i];
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
