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
	public NativeArray<float> CloudElevation;
	public NativeArray<float> CloudTemperature;
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudDropletMass;
	public NativeArray<float2> CloudVelocity;
	public NativeArray<float>[] AirTemperature;
	public NativeArray<float>[] AirHumidity;
	public NativeArray<float2>[] AirVelocity;
	public NativeArray<float>[] WaterTemperature;
	public NativeArray<float>[] WaterMass;
	public NativeArray<float>[] WaterSaltMass;
	public NativeArray<float2>[] WaterVelocity;

	public void Init(int count, int airLayers, int waterLayers)
	{
		Terrain = new NativeArray<CellTerrain>(count, Allocator.Persistent);
		TerrainTemperature = new NativeArray<float>(count, Allocator.Persistent);
		IceTemperature = new NativeArray<float>(count, Allocator.Persistent);
		IceMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudElevation = new NativeArray<float>(count, Allocator.Persistent);
		CloudTemperature = new NativeArray<float>(count, Allocator.Persistent);
		CloudMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudDropletMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudVelocity = new NativeArray<float2>(count, Allocator.Persistent);
		AirTemperature = new NativeArray<float>[waterLayers];
		AirHumidity = new NativeArray<float>[waterLayers];
		AirVelocity = new NativeArray<float2>[waterLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirTemperature[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirHumidity[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVelocity[i] = new NativeArray<float2>(count, Allocator.Persistent);
		}

		WaterTemperature = new NativeArray<float>[waterLayers];
		WaterMass = new NativeArray<float>[waterLayers];
		WaterSaltMass = new NativeArray<float>[waterLayers];
		WaterVelocity = new NativeArray<float2>[waterLayers];
		for (int i = 0; i < waterLayers; i++)
		{
			WaterTemperature[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterSaltMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterVelocity[i] = new NativeArray<float2>(count, Allocator.Persistent);
		}
	}

	public void Dispose()
	{
		Terrain.Dispose();
		TerrainTemperature.Dispose();
		IceTemperature.Dispose();
		IceMass.Dispose();
		CloudElevation.Dispose();
		CloudTemperature.Dispose();
		CloudMass.Dispose();
		CloudDropletMass.Dispose();
		CloudVelocity.Dispose();
		for (int i = 0; i < AirTemperature.Length; i++)
		{
			AirTemperature[i].Dispose();
			AirHumidity[i].Dispose();
			AirVelocity[i].Dispose();
		}

		for (int i = 0; i < WaterTemperature.Length; i++)
		{
			WaterTemperature[i].Dispose();
			WaterMass[i].Dispose();
			WaterSaltMass[i].Dispose();
			WaterVelocity[i].Dispose();
		}
	}
}

public struct CellTerrain {
	public float Elevation;
	public float Roughness;
	public float SoilFertility;
	public float Vegetation;
	public float WaterDepth;
	public float GroundWater;
	public float GroundWaterDepth;
}

public struct PlanetState {
	public int Ticks;
	public float Gravity;
	public float SpinSpeed;
	public float OrbitSpeed;
	public float GeothermalHeat;
	public float SolarRadiation;
	public float StratosphereMass;
	public float CarbonDioxide;
	public float DistanceToSun;
	public float3 Rotation;
	public float3 Position;
	public float AngularSpeed;
}

public struct DependentState {

	public NativeArray<float> WaterDepth;
	public NativeArray<float> SurfaceElevation;
	public NativeArray<float> CloudCoverage;
	public NativeArray<float> VegetationCoverage;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> WaterCoverage;
	public NativeArray<float>[] AirTemperaturePotential;
	public NativeArray<float>[] AirMass;
	public NativeArray<float>[] AirPressure;
	public NativeArray<float>[] WindVertical;
	public NativeArray<float>[] RelativeHumidity;
	public NativeArray<float> LayerHeight;
	public NativeArray<float> LayerElevation;

	// Display variables
	public PlanetDisplay DisplayPlanet;
	public NativeArray<float> Rainfall;
	public NativeArray<float> Evaporation;
	public NativeArray<CellDisplay> CellDisplays;

	public void Init(int count, int airLayers, int waterLayers)
	{
		WaterDepth = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceElevation = new NativeArray<float>(count, Allocator.Persistent);
		CloudCoverage = new NativeArray<float>(count, Allocator.Persistent);
		VegetationCoverage = new NativeArray<float>(count, Allocator.Persistent);
		IceCoverage = new NativeArray<float>(count, Allocator.Persistent);
		WaterCoverage = new NativeArray<float>(count, Allocator.Persistent);
		LayerHeight = new NativeArray<float>(airLayers + waterLayers, Allocator.Persistent);
		LayerElevation = new NativeArray<float>(airLayers + waterLayers, Allocator.Persistent);

		AirTemperaturePotential = new NativeArray<float>[airLayers];
		AirMass = new NativeArray<float>[airLayers];
		AirPressure = new NativeArray<float>[airLayers];
		RelativeHumidity = new NativeArray<float>[airLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirTemperaturePotential[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirPressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			WindVertical[i] = new NativeArray<float>(count, Allocator.Persistent);
			RelativeHumidity[i] = new NativeArray<float>(count, Allocator.Persistent);
		}
		for (int i = 0; i < waterLayers; i++)
		{
		}

		Rainfall = new NativeArray<float>(count, Allocator.Persistent);
		Evaporation = new NativeArray<float>(count, Allocator.Persistent);
		CellDisplays = new NativeArray<CellDisplay>(count, Allocator.Persistent);
	}

	public void Dispose()
	{
		WaterDepth.Dispose();
		SurfaceElevation.Dispose();
		CloudCoverage.Dispose();
		VegetationCoverage.Dispose();
		IceCoverage.Dispose();
		WaterCoverage.Dispose();
		LayerHeight.Dispose();
		LayerElevation.Dispose();

		Rainfall.Dispose();
		Evaporation.Dispose();
		CellDisplays.Dispose();
		for (int i=0;i< AirPressure.Length;i++)
		{
			AirTemperaturePotential[i].Dispose();
			AirMass[i].Dispose();
			AirPressure[i].Dispose();
			WindVertical[i].Dispose();
			RelativeHumidity[i].Dispose();
		}
	}
}




public struct PlanetDisplay {
	public float EnergyDelta;
	public float EnergyIncoming;
	public float EnergySolarReflectedCloud;
	public float EnergySolarReflectedAtmosphere;
	public float EnergySolarReflectedSurface;
	public float EnergySolarAbsorbedCloud;
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
	public float EnergyUpperAir;
	public float EnergyLowerAir;
	public float EnergyShallowWater;
	public float EnergyDeepWater;
	public float EnergyLand;
	public float IceMass;
	public float Temperature;
	public float OceanCoverage;
	public float OceanVolume;
	public float SeaLevel;
	public float CloudCoverage;
	public float Evaporation;
	public float Rainfall;
	public float WaterVapor;
	public float CloudMass;
}


public struct CellDisplay {
	public float Heat;
	public float Rainfall;
	public float Evaporation;

	public float EnergyDelta;
	public float EnergyIncoming;
	public float EnergySolarReflectedCloud;
	public float EnergySolarReflectedAtmosphere;
	public float EnergySolarReflectedSurface;
	public float EnergySolarAbsorbedCloud;
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
}

