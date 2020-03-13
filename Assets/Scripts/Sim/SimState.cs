using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Entities;
using Unity;
using UnityEngine;
using Unity.Mathematics;
using Unity.Collections;

public struct SimState {

	public PlanetState PlanetState;
	public NativeArray<CellTerrain> Terrain;
	public NativeArray<float> TerrainTemperature;
	public NativeArray<float> IceTemperature;
	public NativeArray<float> IceMass;
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudTemperature;
	public NativeArray<float> CloudDropletMass;
	public NativeArray<float3> CloudVelocity;
	public NativeArray<float>[] AirTemperaturePotential;
	public NativeArray<float>[] AirVapor;
	public NativeArray<float3>[] AirVelocity;
	public NativeArray<float>[] WaterTemperature;
	public NativeArray<float>[] WaterMass;
	public NativeArray<float>[] SaltMass;
	public NativeArray<float3>[] WaterVelocity;

	public void Init(int count, int airLayers, int waterLayers)
	{
		Terrain = new NativeArray<CellTerrain>(count, Allocator.Persistent);
		TerrainTemperature = new NativeArray<float>(count, Allocator.Persistent);
		IceTemperature = new NativeArray<float>(count, Allocator.Persistent);
		IceMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudTemperature = new NativeArray<float>(count, Allocator.Persistent);
		CloudDropletMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudVelocity = new NativeArray<float3>(count, Allocator.Persistent);
		AirTemperaturePotential = new NativeArray<float>[airLayers];
		AirVapor = new NativeArray<float>[airLayers];
		AirVelocity = new NativeArray<float3>[airLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirTemperaturePotential[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVapor[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVelocity[i] = new NativeArray<float3>(count, Allocator.Persistent);
		}

		WaterTemperature = new NativeArray<float>[waterLayers];
		WaterMass = new NativeArray<float>[waterLayers];
		SaltMass = new NativeArray<float>[waterLayers];
		WaterVelocity = new NativeArray<float3>[waterLayers];
		for (int i = 0; i < waterLayers; i++)
		{
			WaterTemperature[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			SaltMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterVelocity[i] = new NativeArray<float3>(count, Allocator.Persistent);
		}
	}

	public void Dispose()
	{
		Terrain.Dispose();
		TerrainTemperature.Dispose();
		IceTemperature.Dispose();
		IceMass.Dispose();
		CloudMass.Dispose();
		CloudTemperature.Dispose();
		CloudDropletMass.Dispose();
		CloudVelocity.Dispose();
		for (int i = 0; i < AirTemperaturePotential.Length; i++)
		{
			AirTemperaturePotential[i].Dispose();
			AirVapor[i].Dispose();
			AirVelocity[i].Dispose();
		}

		for (int i = 0; i < WaterTemperature.Length; i++)
		{
			WaterTemperature[i].Dispose();
			WaterMass[i].Dispose();
			SaltMass[i].Dispose();
			WaterVelocity[i].Dispose();
		}
	}
}

public struct CellTerrain {
	public float Elevation;
	public float Roughness;
	public float SoilFertility;
	public float Vegetation;
}

public struct PlanetState {
	public int Ticks;
	public float Gravity;
	public float SpinSpeed;
	public float OrbitSpeed;
	public float GeothermalHeat;
	public float SolarRadiation;
	public float CarbonDioxide;
	public float DistanceToSun;
	public float3 Rotation;
	public float3 Position;
	public float AngularSpeed;
}

public struct DependentState {

	public NativeArray<float> CloudCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> SurfaceAirTemperatureAbsolute;
	public NativeArray<float> IceEnergy;
	public NativeArray<float> WindVerticalCloud;
	public NativeArray<float>[] AirMass;
	public NativeArray<float>[] AirPressure;
	public NativeArray<float>[] AirHumidityAbsolute;
	public NativeArray<float>[] AirHumidityRelative;
	public NativeArray<float>[] WaterCoverage;
	public NativeArray<float>[] WaterDensity;
	public NativeArray<float>[] WaterPressure;
	public NativeArray<float>[] WaterLayerDepth;
	public NativeArray<float>[] WaterLayerHeight;
	public NativeArray<float>[] WaterPotentialEnergy;
	public NativeArray<float>[] AirPotentialEnergy;
	public NativeArray<float>[] LayerHeight;
	public NativeArray<float>[] LayerElevation;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> DewPoint;
	public NativeArray<float> AirDensityCloud;
	public NativeArray<float3>[] DeflectedAirVelocity;
	public NativeArray<float3>[] DeflectedWaterVelocity;
	public NativeArray<float3> DeflectedCloudVelocity;

	public void Init(int count, int airLayers, int waterLayers)
	{
		CloudCoverage = new NativeArray<float>(count, Allocator.Persistent);
		VegetationCoverage = new NativeArray<float>(count, Allocator.Persistent);
		IceCoverage = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAirTemperatureAbsolute = new NativeArray<float>(count, Allocator.Persistent);
		IceEnergy = new NativeArray<float>(count, Allocator.Persistent);
		DewPoint = new NativeArray<float>(count, Allocator.Persistent);
		WindVerticalCloud = new NativeArray<float>(count, Allocator.Persistent);
		AirDensityCloud = new NativeArray<float>(count, Allocator.Persistent);
		CloudElevation = new NativeArray<float>(count, Allocator.Persistent);
		DeflectedCloudVelocity = new NativeArray<float3>(count, Allocator.Persistent);

		AirMass = new NativeArray<float>[airLayers];
		AirPressure = new NativeArray<float>[airLayers];
		AirHumidityAbsolute = new NativeArray<float>[airLayers];
		AirHumidityRelative = new NativeArray<float>[airLayers];
		LayerHeight = new NativeArray<float>[airLayers];
		LayerElevation = new NativeArray<float>[airLayers];
		AirPotentialEnergy = new NativeArray<float>[airLayers];
		DeflectedAirVelocity = new NativeArray<float3>[waterLayers];
		WaterDensity = new NativeArray<float>[waterLayers];
		WaterPressure = new NativeArray<float>[waterLayers];
		WaterLayerDepth = new NativeArray<float>[waterLayers];
		WaterLayerHeight = new NativeArray<float>[waterLayers];
		WaterCoverage = new NativeArray<float>[waterLayers];
		WaterPotentialEnergy = new NativeArray<float>[waterLayers];
		DeflectedWaterVelocity = new NativeArray<float3>[waterLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirPressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirHumidityAbsolute[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirHumidityRelative[i] = new NativeArray<float>(count, Allocator.Persistent);
			LayerHeight[i] = new NativeArray<float>(count, Allocator.Persistent);
			LayerElevation[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirPotentialEnergy[i] = new NativeArray<float>(count, Allocator.Persistent);
			DeflectedAirVelocity[i] = new NativeArray<float3>(count, Allocator.Persistent);
		}
		for (int i = 0; i < waterLayers; i++)
		{
			WaterCoverage[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterDensity[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterPressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterLayerDepth[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterLayerHeight[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterPotentialEnergy[i] = new NativeArray<float>(count, Allocator.Persistent);
			DeflectedWaterVelocity[i] = new NativeArray<float3>(count, Allocator.Persistent);
		}
	}

	public void Dispose()
	{
		CloudCoverage.Dispose();
		VegetationCoverage.Dispose();
		IceCoverage.Dispose();
		SurfaceAirTemperatureAbsolute.Dispose();
		IceEnergy.Dispose();
		DewPoint.Dispose();
		WindVerticalCloud.Dispose();
		AirDensityCloud.Dispose();
		CloudElevation.Dispose();
		DeflectedCloudVelocity.Dispose();

		for (int i = 0; i < AirPressure.Length; i++)
		{
			AirMass[i].Dispose();
			AirPressure[i].Dispose();
			AirHumidityRelative[i].Dispose();
			AirHumidityAbsolute[i].Dispose();
			LayerHeight[i].Dispose();
			LayerElevation[i].Dispose();
			AirPotentialEnergy[i].Dispose();
			DeflectedAirVelocity[i].Dispose();
		}
		for (int i = 0; i < WaterCoverage.Length; i++)
		{
			WaterCoverage[i].Dispose();
			WaterDensity[i].Dispose();
			WaterPressure[i].Dispose();
			WaterLayerDepth[i].Dispose();
			WaterLayerHeight[i].Dispose();
			WaterPotentialEnergy[i].Dispose();
			DeflectedWaterVelocity[i].Dispose();
		}
	}
}



public struct DisplayState {
	public float SolarRadiation;
	public float GeothermalRadiation;
	public float EnergySolarReflectedAtmosphere;
	public float EnergySolarReflectedSurface;
	public float EnergySolarAbsorbedAtmosphere;
	public float EnergySolarAbsorbedSurface;
	public float EnergySolarAbsorbedOcean;
	public float EnergyThermalOceanRadiation;
	public float EnergyOceanConduction;
	public float EnergyEvapotranspiration;
	public float EnergyThermalSurfaceOutAtmosphericWindow;
	public float EnergyThermalOutAtmosphere;
	public float EnergyThermalSurfaceRadiation;
	public float EnergyThermalBackRadiation;
	public float EnergyThermalAbsorbedAtmosphere;
	public float EnergySurfaceConduction;
	public double GlobalEnthalpy;
	public double GlobalEnthalpyTerrain;
	public double GlobalEnthalpyWater;
	public double GlobalEnthalpyAir;
	public double GlobalEnthalpyDelta;
	public double GlobalEnthalpyDeltaTerrain;
	public double GlobalEnthalpyDeltaWater;
	public double GlobalEnthalpyDeltaAir;
	public double GlobalTerrainTemperature;
	public float GlobalIceMass;
	public float GlobalSurfaceTemperature;
	public double GlobalAirTemperaturePotential;
	public float GlobalOceanCoverage;
	public float GlobalOceanVolume;
	public float GlobalOceanTemperature;
	public float GlobalOceanSurfaceTemperature;
	public float GlobalSeaLevel;
	public float GlobalCloudCoverage;
	public float GlobalEvaporation;
	public float GlobalRainfall;
	public float GlobalCondensationCloud;
	public float GlobalCondensationGround;
	public double GlobalWaterVapor;
	public double GlobalAirMass;
	public float GlobalCloudMass;

	public NativeArray<float> SolarRadiationAbsorbedSurface;
	public NativeArray<float> Rainfall;
	public NativeArray<float> CondensationCloud;
	public NativeArray<float> CondensationGround;
	public NativeArray<float> Evaporation;
	public NativeArray<float> EnthalpyTerrain;
	public NativeArray<float>[] EnthalpyWater;
	public NativeArray<float>[] EnthalpyAir;
	public NativeArray<float>[] Salinity;
	public NativeArray<float>[] Pressure;
	public NativeArray<float3>[] PressureGradientForce;
	public NativeArray<float>[] ThermalDelta;
	public NativeArray<float>[] ConductionDelta;
	public NativeArray<float>[] SolarDelta;
	public NativeArray<float>[] LatentHeatDelta;

	public void Init(int count, int airLayers, int waterLayers, int totalLayers)
	{
		SolarRadiationAbsorbedSurface = new NativeArray<float>(count, Allocator.Persistent);
		Rainfall = new NativeArray<float>(count, Allocator.Persistent);
		CondensationCloud = new NativeArray<float>(count, Allocator.Persistent);
		CondensationGround = new NativeArray<float>(count, Allocator.Persistent);
		Evaporation = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyTerrain = new NativeArray<float>(count, Allocator.Persistent);

		Pressure = new NativeArray<float>[airLayers];
		PressureGradientForce = new NativeArray<float3>[airLayers];
		EnthalpyAir = new NativeArray<float>[airLayers];
		for (int i = 0; i < airLayers; i++)
		{
			Pressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			PressureGradientForce[i] = new NativeArray<float3>(count, Allocator.Persistent);
			EnthalpyAir[i] = new NativeArray<float>(count, Allocator.Persistent);
		}

		Salinity = new NativeArray<float>[waterLayers];
		EnthalpyWater = new NativeArray<float>[waterLayers];
		for (int i=0;i<waterLayers;i++)
		{
			Salinity[i] = new NativeArray<float>(count, Allocator.Persistent);
			EnthalpyWater[i] = new NativeArray<float>(count, Allocator.Persistent);
		}

		ThermalDelta = new NativeArray<float>[totalLayers];
		ConductionDelta = new NativeArray<float>[totalLayers];
		SolarDelta = new NativeArray<float>[totalLayers];
		LatentHeatDelta = new NativeArray<float>[totalLayers];
		for (int i = 0; i < totalLayers; i++)
		{
			ThermalDelta[i] = new NativeArray<float>(count, Allocator.Persistent);
			ConductionDelta[i] = new NativeArray<float>(count, Allocator.Persistent);
			SolarDelta[i] = new NativeArray<float>(count, Allocator.Persistent);
			LatentHeatDelta[i] = new NativeArray<float>(count, Allocator.Persistent);
		}

	}

	public void Dispose()
	{
		SolarRadiationAbsorbedSurface.Dispose();
		Rainfall.Dispose();
		CondensationCloud.Dispose();
		CondensationGround.Dispose();
		Evaporation.Dispose();
		EnthalpyTerrain.Dispose();
		for (int i = 0; i < Pressure.Length; i++)
		{
			Pressure[i].Dispose();
			PressureGradientForce[i].Dispose();
			EnthalpyAir[i].Dispose();
		}
		for (int i=0;i<Salinity.Length;i++)
		{
			Salinity[i].Dispose();
			EnthalpyWater[i].Dispose();
		}
		for (int i=0;i<ThermalDelta.Length;i++)
		{
			ThermalDelta[i].Dispose();
			ConductionDelta[i].Dispose();
			SolarDelta[i].Dispose();
			LatentHeatDelta[i].Dispose();
		}
	}

}


