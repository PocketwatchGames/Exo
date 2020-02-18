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
	public NativeArray<float> CloudTemperature;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudDropletMass;
	public NativeArray<float2> CloudVelocity;
	public NativeArray<float>[] AirTemperature;
	public NativeArray<float>[] AirVapor;
	public NativeArray<float2>[] Wind;
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
		AirTemperature = new NativeArray<float>[airLayers];
		AirVapor = new NativeArray<float>[airLayers];
		Wind = new NativeArray<float2>[airLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirTemperature[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVapor[i] = new NativeArray<float>(count, Allocator.Persistent);
			Wind[i] = new NativeArray<float2>(count, Allocator.Persistent);
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
			AirVapor[i].Dispose();
			Wind[i].Dispose();
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
	public NativeArray<float> SurfaceAirTemperature;
	public NativeArray<float> IceEnergy;
	public NativeArray<float> CloudEnergy;
	public NativeArray<float> WindVerticalCloud;
	public NativeArray<float>[] AirMass;
	public NativeArray<float>[] AirPressure;
	public NativeArray<float>[] AirHumidityAbsolute;
	public NativeArray<float>[] AirHumidityRelative;
	public NativeArray<float>[] WaterCoverage;
	public NativeArray<float>[] WaterSalinity;
	public NativeArray<float>[] WaterPotentialEnergy;
	public NativeArray<float>[] AirPotentialEnergy;
	public NativeArray<float>[] LayerHeight;
	public NativeArray<float>[] LayerElevation;
	public NativeArray<float2> PressureGradientAtCloudElevation;
	public NativeArray<float> AirMassCloud;
	public NativeArray<float> AirVaporCloud;
	public NativeArray<float> AirTemperatureCloud;
	public NativeArray<float> AirPressureCloud;
	public NativeArray<float> AirHumidityRelativeCloud;
	public NativeArray<int> AirLayerCloud;

	public void Init(int count, int airLayers, int waterLayers)
	{
		WaterDepth = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceElevation = new NativeArray<float>(count, Allocator.Persistent);
		CloudCoverage = new NativeArray<float>(count, Allocator.Persistent);
		VegetationCoverage = new NativeArray<float>(count, Allocator.Persistent);
		IceCoverage = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAirTemperature = new NativeArray<float>(count, Allocator.Persistent);
		IceEnergy = new NativeArray<float>(count, Allocator.Persistent);
		CloudEnergy = new NativeArray<float>(count, Allocator.Persistent);
		WindVerticalCloud = new NativeArray<float>(count, Allocator.Persistent);
		PressureGradientAtCloudElevation = new NativeArray<float2>(count, Allocator.Persistent);
		AirMassCloud = new NativeArray<float>(count, Allocator.Persistent);
		AirVaporCloud = new NativeArray<float>(count, Allocator.Persistent);
		AirTemperatureCloud = new NativeArray<float>(count, Allocator.Persistent);
		AirPressureCloud = new NativeArray<float>(count, Allocator.Persistent);
		AirHumidityRelativeCloud = new NativeArray<float>(count, Allocator.Persistent);
		AirLayerCloud = new NativeArray<int>(count, Allocator.Persistent);

		AirMass = new NativeArray<float>[airLayers];
		AirPressure = new NativeArray<float>[airLayers];
		AirHumidityAbsolute = new NativeArray<float>[airLayers];
		AirHumidityRelative = new NativeArray<float>[airLayers];
		LayerHeight = new NativeArray<float>[airLayers];
		LayerElevation = new NativeArray<float>[airLayers];
		AirPotentialEnergy = new NativeArray<float>[airLayers];
		WaterSalinity = new NativeArray<float>[waterLayers];
		WaterCoverage = new NativeArray<float>[waterLayers];
		WaterPotentialEnergy = new NativeArray<float>[waterLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirPressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirHumidityAbsolute[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirHumidityRelative[i] = new NativeArray<float>(count, Allocator.Persistent);
			LayerHeight[i] = new NativeArray<float>(count, Allocator.Persistent);
			LayerElevation[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirPotentialEnergy[i] = new NativeArray<float>(count, Allocator.Persistent);
		}
		for (int i = 0; i < waterLayers; i++)
		{
			WaterCoverage[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterSalinity[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterPotentialEnergy[i] = new NativeArray<float>(count, Allocator.Persistent);
		}
	}

	public void Dispose()
	{
		WaterDepth.Dispose();
		SurfaceElevation.Dispose();
		CloudCoverage.Dispose();
		VegetationCoverage.Dispose();
		IceCoverage.Dispose();
		SurfaceAirTemperature.Dispose();
		IceEnergy.Dispose();
		CloudEnergy.Dispose();
		WindVerticalCloud.Dispose();
		PressureGradientAtCloudElevation.Dispose();
		AirMassCloud.Dispose();
		AirVaporCloud.Dispose();
		AirTemperatureCloud.Dispose();
		AirPressureCloud.Dispose();
		AirHumidityRelativeCloud.Dispose();
		AirLayerCloud.Dispose();

		for (int i = 0; i < AirPressure.Length; i++)
		{
			AirMass[i].Dispose();
			AirPressure[i].Dispose();
			AirHumidityRelative[i].Dispose();
			AirHumidityAbsolute[i].Dispose();
			LayerHeight[i].Dispose();
			LayerElevation[i].Dispose();
			AirPotentialEnergy[i].Dispose();
		}
		for (int i = 0; i < WaterSalinity.Length; i++)
		{
			WaterCoverage[i].Dispose();
			WaterSalinity[i].Dispose();
			WaterPotentialEnergy[i].Dispose();
		}
	}
}



public struct DisplayState {
	public float SolarRadiation;
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
	public float EnergyTerrain;
	public float GlobalIceMass;
	public float GlobalTemperature;
	public float GlobalOceanCoverage;
	public float GlobalOceanVolume;
	public float GlobalSeaLevel;
	public float GlobalCloudCoverage;
	public float GlobalEvaporation;
	public float GlobalRainfall;
	public float GlobalWaterVapor;
	public float GlobalCloudMass;

	public NativeArray<float> SolarRadiationAbsorbedSurface;
	public NativeArray<float> Rainfall;
	public NativeArray<float> Evaporation;
	public NativeArray<float> Pressure;

	public void Init(int count, int airLayers, int waterLayers)
	{
		SolarRadiationAbsorbedSurface = new NativeArray<float>(count, Allocator.Persistent);
		Rainfall = new NativeArray<float>(count, Allocator.Persistent);
		Evaporation = new NativeArray<float>(count, Allocator.Persistent);
		Pressure = new NativeArray<float>(count, Allocator.Persistent);
	}

	public void Dispose()
	{
		SolarRadiationAbsorbedSurface.Dispose();
		Rainfall.Dispose();
		Evaporation.Dispose();
		Pressure.Dispose();
	}

}


