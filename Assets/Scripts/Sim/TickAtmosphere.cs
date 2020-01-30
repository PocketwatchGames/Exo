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
public struct DiffusionWater {
	public float Temperature;
	public float Salinity;
	public float2 Velocity;
}


[BurstCompile]
public struct SolarRadiationJob : IJobParallelFor {
	public NativeArray<float> SurfaceElevation;
	public NativeArray<float> SolarRadiation;
	public NativeArray<float> WaterSlopeAlbedo;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> WaterCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> CloudCoverage;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float3> SphericalPosition;
	[ReadOnly] public float3 SunToPlanetDir;
	[ReadOnly] public quaternion PlanetRotation;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float IncomingSolarRadiation;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		float cloudMass = CloudMass[i];
		float sunDotSurface = math.dot(SunToPlanetDir, math.rotate(PlanetRotation, -SphericalPosition[i]));
		float waterSlopeAlbedo = math.pow(1.0f - math.max(0, sunDotSurface), 9);
		SolarRadiation[i] = IncomingSolarRadiation * SecondsPerTick * sunDotSurface;
		WaterSlopeAlbedo[i] = waterSlopeAlbedo;
		VegetationCoverage[i] = Terrain[i].Vegetation * worldData.inverseFullCanopyCoverage;
		WaterCoverage[i] = Terrain[i].WaterDepth * worldData.inverseFullWaterCoverage;
		CloudCoverage[i] = CloudMass[i] * worldData.inverseCloudMassFullAbsorption;
		IceCoverage[i] = IceMass[i] * worldData.inverseFullIceCoverage;
	}
}

[BurstCompile]
public struct SolarRadiationAbsorbedAirJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationAbsorbedCloud;
	public NativeArray<float> SolarRadiationIncoming;
	[ReadOnly] public float LayerElevation;
	[ReadOnly] public float LayerHeight;
	[ReadOnly] public float SolarReflectivityAir;
	[ReadOnly] public float SolarAbsorptivityAir;
	[ReadOnly] public float SolarAbsorptivityWaterVapor;
	[ReadOnly] public float SolarAbsorptivityCloud;
	[ReadOnly] public NativeArray<float> CloudCoverage;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudTemperature;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> WaterSlopeAlbedo;
	[ReadOnly] public NativeArray<float> Humidity;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		float airMass = AirMass[i];
		float waterVaporMass = Humidity[i];
		float incomingRadiation = SolarRadiationIncoming[i];

		float aboveCloud = 1;
		float belowCloud = 0;

		float cloudMass = CloudMass[i];
		if (cloudMass > 0) {
			float cloudElevation = CloudElevation[i];
			if (cloudElevation >= LayerElevation && cloudElevation < LayerElevation + LayerHeight)
			{
				aboveCloud = (cloudElevation - LayerElevation) / LayerHeight;
				belowCloud = 1.0f - aboveCloud;
			}
		}

		float emissivity = SolarAbsorptivityAir * airMass + SolarAbsorptivityWaterVapor * waterVaporMass;
		float absorbed = 0;
		float energyReflectedAtmosphere = aboveCloud * incomingRadiation * math.min(1, SolarReflectivityAir * (airMass + waterVaporMass));
		incomingRadiation -= energyReflectedAtmosphere;
		//displayCell.EnergySolarReflectedAtmosphere += energyReflectedAtmosphere;

		absorbed += aboveCloud * incomingRadiation * emissivity;
		incomingRadiation -= absorbed;

		float absorbedByCloudsIncoming = 0;
		if (cloudMass > 0)
		{
			float cloudCoverage = CloudCoverage[i];
			float cloudIceContent = math.saturate((CloudTemperature[i] - worldData.minCloudFreezingTemperature) / (worldData.maxCloudFreezingTemperature - worldData.minCloudFreezingTemperature));
			float cloudTemperatureAlbedo = WorldData.AlbedoIce + (WorldData.AlbedoWater - WorldData.AlbedoIce) * cloudIceContent;
			float rainDropSizeAlbedo = math.saturate(1.0f - CloudDropletMass[i] / CloudMass[i]) * (worldData.rainDropSizeAlbedoMax - worldData.rainDropSizeAlbedoMin) + worldData.rainDropSizeAlbedoMin;
			float cloudAlbedo = math.min(1.0f, worldData.SolarReflectivityCloud * cloudTemperatureAlbedo * CloudMass[i] * rainDropSizeAlbedo / math.max(worldData.maxCloudSlopeAlbedo, 1.0f - WaterSlopeAlbedo[i]));

			float energyReflectedClouds = incomingRadiation * cloudAlbedo * cloudCoverage;
			incomingRadiation -= energyReflectedClouds;
			absorbedByCloudsIncoming = incomingRadiation * math.min(1.0f, SolarAbsorptivityCloud * cloudMass) * cloudCoverage;
			incomingRadiation -= absorbedByCloudsIncoming;


			//displayCell.EnergySolarAbsorbedCloud += absorbedByCloudsIncoming;
			//displayCell.EnergySolarAbsorbedAtmosphere += absorbedByCloudsIncoming;
			//displayCell.EnergySolarReflectedCloud += energyReflectedClouds;

			// below cloud
			float absorbedBelow = incomingRadiation * emissivity;
			absorbed += absorbedBelow;
			incomingRadiation -= absorbedBelow;
		}

		SolarRadiationAbsorbedCloud[i] = absorbedByCloudsIncoming;
		SolarRadiationAbsorbed[i] = absorbed;
		SolarRadiationIncoming[i] = incomingRadiation;
	}
}

[BurstCompile]
public struct SolarRadiationAbsorbedIceJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationIncoming;
	[ReadOnly] public float AlbedoIce;
	[ReadOnly] public NativeArray<float> IceCoverage;
	public void Execute(int i)
	{
		float incoming = SolarRadiationIncoming[i];
		float reflected = incoming * AlbedoIce;
		incoming -= reflected;
		float absorbed = incoming * IceCoverage[i];
		SolarRadiationAbsorbed[i] = absorbed;
		SolarRadiationIncoming[i] = incoming - absorbed;
	}
}
[BurstCompile]
public struct SolarRadiationAbsorbedWaterJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
	public NativeArray<float> SolarRadiationIncoming;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> WaterSlopeAlbedo;
	public void Execute(int i)
	{
		float incoming = SolarRadiationIncoming[i];
		float reflected = incoming * WaterSlopeAlbedo[i];
		incoming -= reflected;
		float absorbed = incoming * WaterCoverage[i];
		SolarRadiationAbsorbed[i] = absorbed;
		SolarRadiationIncoming[i] = incoming - absorbed;
	}
}
[BurstCompile]
public struct SolarRadiationAbsorbedTerrainJob : IJobParallelFor {
	public NativeArray<float> SolarRadiationAbsorbed;
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
		float reflected = incoming * vegetationCoverage * WorldData.AlbedoFoliage + math.max(0, 1.0f - vegetationCoverage) * soilReflectivity;
		incoming -= reflected;

		SolarRadiationAbsorbed[i] = incoming;
	}
}

[BurstCompile]
public struct EmissivityAirJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> Humidity;
	[ReadOnly] public float GreenhouseGasConcentration;
	[ReadOnly] public float AbsorptivityAir;
	[ReadOnly] public float AbsorptivityCarbonDioxide;
	[ReadOnly] public float AbsorptivityWaterVapor;
	public void Execute(int i)
	{
		Emissivity[i] = 
			AirMass[i] * ((1.0f - GreenhouseGasConcentration) * AbsorptivityAir + GreenhouseGasConcentration * AbsorptivityCarbonDioxide) +
			Humidity[i] * AbsorptivityWaterVapor;
	}
}

[BurstCompile]
public struct EmissivityWaterJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> WaterSaltMass;
	[ReadOnly] public NativeArray<float> WaterMass;
	public void Execute(int i)
	{
		Emissivity[i] = (WorldData.EmissivityWater * WaterMass[i] + WorldData.EmissivitySalt * WaterSaltMass[i]) / (WaterMass[i] + WaterSaltMass[i]);
	}
}

[BurstCompile]
public struct EmissivityTerrainJob : IJobParallelFor {
	public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> VegetationCoverage;
	public void Execute(int i)
	{
		Emissivity[i] = math.lerp(
			math.lerp(WorldData.EmissivitySand, WorldData.EmissivityDirt, Terrain[i].SoilFertility), 
			WorldData.EmissivityVegetation, 
			VegetationCoverage[i]);
	}
}


[BurstCompile]
public struct ThermalEnergyRadiatedJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationEmitted;
	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> ThermalRadiationTransmittedDown;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Emissivity;
	public void Execute(int i)
	{
		// radiate half up and half down
		float emitted = Atmosphere.GetRadiationRate(Temperature[i], Emissivity[i]);
		ThermalRadiationEmitted[i] = emitted;
		float transmittedUp = emitted / 2;
		ThermalRadiationTransmittedUp[i] = transmittedUp;
		ThermalRadiationTransmittedDown[i] = transmittedUp;
	}
}
[BurstCompile]
public struct ThermalEnergyRadiatedConstantEmissivityJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationEmitted;
	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> ThermalRadiationTransmittedDown;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public float Emissivity;
	public void Execute(int i)
	{
		// radiate half up and half down
		float emitted = Atmosphere.GetRadiationRate(Temperature[i], Emissivity);
		ThermalRadiationEmitted[i] = emitted;
		float transmittedUp = emitted / 2;
		ThermalRadiationTransmittedUp[i] = transmittedUp;
		ThermalRadiationTransmittedDown[i] = transmittedUp;
	}
}
[BurstCompile]
public struct ThermalEnergyRadiatedTerrainJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationEmitted;
	public NativeArray<float> ThermalRadiationTransmitted;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Emissivity;
	public void Execute(int i)
	{
		// radiate half up and half down
		float emitted = Atmosphere.GetRadiationRate(Temperature[i], Emissivity[i]) / 2;
		ThermalRadiationEmitted[i] = emitted;
		ThermalRadiationTransmitted[i] = emitted;
	}
}

[BurstCompile]
public struct ThermalEnergyAbsorbedAirJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> ThermalRadiationAbsorbed;
	public NativeArray<float> ThermalRadiationAbsorbedCloud;
	[ReadOnly] public NativeArray<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeArray<float> ThermalRadiationEmittedCloudUp;
	[ReadOnly] public NativeArray<float> ThermalRadiationEmittedCloudDown;
	[ReadOnly] public float LayerElevation;
	[ReadOnly] public float LayerHeight;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudCoverage;
	public void Execute(int i)
	{
		//float aboveCloud
		float incoming = ThermalRadiationTransmitted[i];
		float absorbed = incoming * Emissivity[i];
		ThermalRadiationEmitted[i] += incoming - absorbed;
		ThermalRadiationAbsorbed[i] += absorbed;
	}
}


[BurstCompile]
public struct ThermalEnergyAbsorbedJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> ThermalRadiationAbsorbed;
	[ReadOnly] public NativeArray<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeArray<float> Emissivity;
	public void Execute(int i)
	{
		float incoming = ThermalRadiationTransmitted[i];
		float absorbed = incoming * Emissivity[i];
		ThermalRadiationEmitted[i] += incoming - absorbed;
		ThermalRadiationAbsorbed[i] += absorbed;
	}
}

[BurstCompile]
public struct ThermalEnergyAbsorbedPartialCoverageJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> ThermalRadiationAbsorbed;
	[ReadOnly] public NativeArray<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeArray<float> Emissivity;
	[ReadOnly] public NativeArray<float> Coverage;
	public void Execute(int i)
	{
		float incoming = ThermalRadiationTransmitted[i];
		float absorbed = incoming * Emissivity[i] * Coverage[i];
		ThermalRadiationEmitted[i] += incoming - absorbed;
		ThermalRadiationAbsorbed[i] += absorbed;
	}
}

[BurstCompile]
public struct ThermalEnergyAbsorbedPartialCoverageConstantEmissivityJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationTransmitted;
	public NativeArray<float> ThermalRadiationAbsorbed;
	[ReadOnly] public NativeArray<float> ThermalRadiationEmitted;
	[ReadOnly] public NativeArray<float> Coverage;
	[ReadOnly] public float Emissivity;
	public void Execute(int i)
	{
		float incoming = ThermalRadiationTransmitted[i];
		float absorbed = incoming * Emissivity * Coverage[i];
		ThermalRadiationEmitted[i] += incoming - absorbed;
		ThermalRadiationAbsorbed[i] += absorbed;
	}
}

[BurstCompile]
public struct ThermalEnergyAbsorbedTerrainJob : IJobParallelFor {
	public NativeArray<float> ThermalRadiationAbsorbed;
	[ReadOnly] public NativeArray<float> ThermalRadiationIncoming;
	public void Execute(int i)
	{
		ThermalRadiationAbsorbed[i] += ThermalRadiationIncoming[i];
	}
}

[BurstCompile]
public struct DiffusionAirJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Humidity;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public float DiffusionCoefficientTemperature;
	[ReadOnly] public float DiffusionCoefficientHumidity;
	[ReadOnly] public float DiffusionCoefficientVelocity;
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
	[ReadOnly] public NativeArray<float> Salinity;
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
				gradientSalinity += Salinity[n];
				gradientTemperature += Temperature[n];
				gradientVelocity += Velocity[n];

				neighborCount++;
			}
		}
		gradientSalinity -= Salinity[i] * neighborCount;
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
	[ReadOnly] public NativeArray<float> Humidity;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	public void Execute(int i)
	{

	}
}

[BurstCompile]
public struct AdvectionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salinity;
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
	[ReadOnly] public NativeArray<float> MassA;
	[ReadOnly] public NativeArray<float> MassB;
	[ReadOnly] public float ConductionCoefficient;
	public void Execute(int i)
	{

	}
}

[BurstCompile]
public struct ConductionAirIceJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> MassA;
	[ReadOnly] public NativeArray<float> MassB;
	[ReadOnly] public float ConductionCoefficient;
	public void Execute(int i)
	{

	}
}


[BurstCompile]
public struct ConductionAirWaterJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> MassA;
	[ReadOnly] public NativeArray<float> MassB;
	[ReadOnly] public float ConductionCoefficientPositive;
	[ReadOnly] public float ConductionCoefficientNegative;
	public void Execute(int i)
	{

	}
}

[BurstCompile]
public struct ConductionAirTerrainJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> MassA;
	[ReadOnly] public NativeArray<float> MassB;
	[ReadOnly] public float ConductionCoefficient;
	public void Execute(int i)
	{

	}
}

[BurstCompile]
public struct ConductionIceWaterJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> MassA;
	[ReadOnly] public NativeArray<float> MassB;
	[ReadOnly] public float ConductionCoefficient;
	public void Execute(int i)
	{

	}
}

[BurstCompile]
public struct ConductionIceTerrainJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> MassA;
	[ReadOnly] public NativeArray<float> MassB;
	[ReadOnly] public float ConductionCoefficient;
	public void Execute(int i)
	{

	}
}

[BurstCompile]
public struct ConductionWaterTerrainJob : IJobParallelFor {
	public NativeArray<float> EnergyDelta;
	[ReadOnly] public NativeArray<float> TemperatureA;
	[ReadOnly] public NativeArray<float> TemperatureB;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> MassA;
	[ReadOnly] public NativeArray<float> MassB;
	[ReadOnly] public float ConductionCoefficient;
	public void Execute(int i)
	{

	}
}


[BurstCompile]
public struct EnergyAirJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Humidity;
	public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastHumidity;
	[ReadOnly] public NativeArray<float2> LastVelocity;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
	[ReadOnly] public NativeArray<DiffusionAir> Diffusion;
	[ReadOnly] public NativeArray<float> ThermalRadiationOut;
	[ReadOnly] public NativeArray<float> ThermalRadiationIn;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergy;
	public void Execute(int i)
	{
		float conductionTemperatureDelta = 0;
		Temperature[i] = LastTemperature[i] + Diffusion[i].Temperature - Advection[i].Temperature - ThermalRadiationOut[i] + ThermalRadiationIn[i] + SolarRadiationIn[i] + conductionTemperatureDelta;
		Humidity[i] = LastHumidity[i] + Diffusion[i].Humidity - Advection[i].Humidity;
		Velocity[i] = LastVelocity[i] + Diffusion[i].Velocity - Advection[i].Velocity;
	}
}

[BurstCompile]
public struct EnergyWaterJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Salinity;
	public NativeArray<float2> Velocity;
	[ReadOnly] public NativeArray<float> LastTemperature;
	[ReadOnly] public NativeArray<float> LastSalinity;
	[ReadOnly] public NativeArray<float2> LastVelocity;
	[ReadOnly] public NativeArray<DiffusionWater> Advection;
	[ReadOnly] public NativeArray<DiffusionWater> Diffusion;
	[ReadOnly] public NativeArray<float> ThermalRadiationOut;
	[ReadOnly] public NativeArray<float> ThermalRadiationIn;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergy;
	public void Execute(int i)
	{
		float conductionTemperatureDelta = 0;
		Temperature[i] = LastTemperature[i] + Diffusion[i].Temperature - Advection[i].Temperature - ThermalRadiationOut[i] + ThermalRadiationIn[i] + SolarRadiationIn[i] + conductionTemperatureDelta;
		Salinity[i] = LastSalinity[i] + Diffusion[i].Salinity - Advection[i].Salinity;
		Velocity[i] = LastVelocity[i] + Diffusion[i].Velocity - Advection[i].Velocity;
	}
}

[BurstCompile]
public struct EnergyCloudJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Elevation;
	public NativeArray<float> Evaporation;
	public NativeArray<float> Rainfall;
	public NativeArray<float> LatentHeat;
	public NativeArray<float> DropletSize;
	[ReadOnly] public NativeArray<float> ThermalRadiationOut;
	[ReadOnly] public NativeArray<float> ThermalRadiationIn;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> AirTemperature;
	[ReadOnly] public NativeArray<float> AirHumidity;
	[ReadOnly] public NativeArray<float> ConductionEnergy;
	public void Execute(int i)
	{
		// evap, ice melting/freezing, cloud sim goes here
	}
}

[BurstCompile]
public struct EnergyIceJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Melted;
	public NativeArray<float> LatentHeatBottom;
	public NativeArray<float> LatentHeatTop;
	[ReadOnly] public NativeArray<float> iceMass;
	[ReadOnly] public NativeArray<float> rainfall;
	[ReadOnly] public NativeArray<float> ThermalRadiationOut;
	[ReadOnly] public NativeArray<float> ThermalRadiationInBottom;
	[ReadOnly] public NativeArray<float> ThermalRadiationInTop;
	[ReadOnly] public NativeArray<float> SolarRadiationIn;
	[ReadOnly] public NativeArray<float> ConductionEnergy;
	public void Execute(int i)
	{
	}
}

[BurstCompile]
public struct EnergyTerrainJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> ThermalRadiationOutTerrain;
	[ReadOnly] public NativeArray<float> ThermalRadiationInTerrain;
	[ReadOnly] public NativeArray<float> SolarRadiationInTerrain;
	[ReadOnly] public NativeArray<float> ConductionEnergy;
	[ReadOnly] float GeothermalEnergy;
	public void Execute(int i)
	{
	}
}

[BurstCompile]
public struct UpdateDependentStateJob : IJobParallelFor {
	public NativeArray<float> SurfaceElevation;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> WaterCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> CloudCoverage;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public WorldData worldData;
	public void Execute(int i)
	{
		VegetationCoverage[i] = Terrain[i].Vegetation * worldData.inverseFullCanopyCoverage;
		WaterCoverage[i] = Terrain[i].WaterDepth * worldData.inverseFullWaterCoverage;
		CloudCoverage[i] = CloudMass[i] * worldData.inverseCloudMassFullAbsorption;
		IceCoverage[i] = IceMass[i] * worldData.inverseFullIceCoverage;
	}

	public void UpdateOut(ref DependentState output)
	{
		output.SurfaceElevation.CopyFrom(SurfaceElevation);
		output.IceCoverage.CopyFrom(IceCoverage);
		output.WaterCoverage.CopyFrom(WaterCoverage);
		output.VegetationCoverage.CopyFrom(VegetationCoverage);
		output.CloudCoverage.CopyFrom(CloudCoverage);
	}
	public void Dispose()
	{
		SurfaceElevation.Dispose();
		IceCoverage.Dispose();
		WaterCoverage.Dispose();
		VegetationCoverage.Dispose();
		CloudCoverage.Dispose();
		IceMass.Dispose();
		CloudMass.Dispose();
		Terrain.Dispose();
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
