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
	public NativeArray<float> TerrainEnergy;
	public NativeArray<float> IceEnergy;
	public NativeArray<float> IceMass;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> CloudEnergy;
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudDropletMass;
	public NativeArray<float2> CloudVelocity;
	public NativeArray<float>[] AirEnergy;
	public NativeArray<float>[] AirVapor;
	public NativeArray<float2>[] AirVelocity;
	public NativeArray<float>[] WaterEnergy;
	public NativeArray<float>[] WaterMass;
	public NativeArray<float>[] WaterSaltMass;
	public NativeArray<float2>[] WaterVelocity;

	public void Init(int count, int airLayers, int waterLayers)
	{
		Terrain = new NativeArray<CellTerrain>(count, Allocator.Persistent);
		TerrainEnergy = new NativeArray<float>(count, Allocator.Persistent);
		IceEnergy = new NativeArray<float>(count, Allocator.Persistent);
		IceMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudElevation = new NativeArray<float>(count, Allocator.Persistent);
		CloudEnergy = new NativeArray<float>(count, Allocator.Persistent);
		CloudMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudDropletMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudVelocity = new NativeArray<float2>(count, Allocator.Persistent);
		AirEnergy = new NativeArray<float>[waterLayers];
		AirVapor = new NativeArray<float>[waterLayers];
		AirVelocity = new NativeArray<float2>[waterLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirEnergy[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVapor[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVelocity[i] = new NativeArray<float2>(count, Allocator.Persistent);
		}

		WaterEnergy = new NativeArray<float>[waterLayers];
		WaterMass = new NativeArray<float>[waterLayers];
		WaterSaltMass = new NativeArray<float>[waterLayers];
		WaterVelocity = new NativeArray<float2>[waterLayers];
		for (int i = 0; i < waterLayers; i++)
		{
			WaterEnergy[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterSaltMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterVelocity[i] = new NativeArray<float2>(count, Allocator.Persistent);
		}
	}

	public void Dispose()
	{
		Terrain.Dispose();
		TerrainEnergy.Dispose();
		IceEnergy.Dispose();
		IceMass.Dispose();
		CloudElevation.Dispose();
		CloudEnergy.Dispose();
		CloudMass.Dispose();
		CloudDropletMass.Dispose();
		CloudVelocity.Dispose();
		for (int i = 0; i < AirEnergy.Length; i++)
		{
			AirEnergy[i].Dispose();
			AirVapor[i].Dispose();
			AirVelocity[i].Dispose();
		}

		for (int i = 0; i < WaterEnergy.Length; i++)
		{
			WaterEnergy[i].Dispose();
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
	public NativeArray<float>[] AirTemperature;
	public NativeArray<float>[] AirMass;
	public NativeArray<float>[] AirPressure;
	public NativeArray<float>[] AirHumidityAbsolute;
	public NativeArray<float>[] AirHumidityRelative;
	public NativeArray<float>[] WindVertical;
	public NativeArray<float>[] WaterTemperature;
	public NativeArray<float>[] WaterSalinity;
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

		AirTemperature = new NativeArray<float>[airLayers];
		AirMass = new NativeArray<float>[airLayers];
		AirPressure = new NativeArray<float>[airLayers];
		AirHumidityAbsolute = new NativeArray<float>[airLayers];
		AirHumidityRelative = new NativeArray<float>[airLayers];
		WaterTemperature = new NativeArray<float>[waterLayers];
		WaterSalinity = new NativeArray<float>[waterLayers];
		WindVertical = new NativeArray<float>[airLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirTemperature[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirPressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirHumidityAbsolute[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirHumidityRelative[i] = new NativeArray<float>(count, Allocator.Persistent);
			WindVertical[i] = new NativeArray<float>(count, Allocator.Persistent);
		}
		for (int i = 0; i < waterLayers; i++)
		{
			WaterTemperature[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterSalinity[i] = new NativeArray<float>(count, Allocator.Persistent);
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
		for (int i = 0; i < AirPressure.Length; i++)
		{
			AirTemperature[i].Dispose();
			AirMass[i].Dispose();
			AirPressure[i].Dispose();
			WindVertical[i].Dispose();
			AirHumidityRelative[i].Dispose();
			AirHumidityAbsolute[i].Dispose();
		}
		for (int i = 0; i < WaterSalinity.Length; i++)
		{
			WaterSalinity[i].Dispose();
			WaterTemperature[i].Dispose();
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

