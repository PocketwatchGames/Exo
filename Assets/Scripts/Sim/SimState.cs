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
	public NativeArray<float> Elevation;
	public NativeArray<float> Roughness;
	public NativeArray<float> SoilFertility;
	public NativeArray<float> TerrainTemperature;
	public NativeArray<float> FloraMass;
	public NativeArray<float> FloraWater;
	public NativeArray<float> FloraTemperature;
	public NativeArray<float> FloraGlucose;
	public NativeArray<float> GroundWater;
	public NativeArray<float> GroundWaterTemperature;
	public NativeArray<float> IceTemperature;
	public NativeArray<float> IceMass;
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudTemperature;
	public NativeArray<float> CloudDropletMass;
	public NativeArray<float3> CloudVelocity;
	public NativeArray<float> CrustDepth;
	public NativeArray<float> MagmaMass;
	public NativeArray<float> LavaMass;
	public NativeArray<float> LavaTemperature;
	public NativeArray<float>[] AirTemperaturePotential;
	public NativeArray<float>[] AirVapor;
	public NativeArray<float>[] AirCarbonDioxide;
	public NativeArray<float>[] AirOxygen;
	public NativeArray<float3>[] AirVelocity;
	public NativeArray<float>[] Dust;
	public NativeArray<float>[] WaterTemperature;
	public NativeArray<float>[] WaterMass;
	public NativeArray<float>[] SaltMass;
	public NativeArray<float3>[] WaterVelocity;

	public void Init(int count, int airLayers, int waterLayers)
	{
		Roughness = new NativeArray<float>(count, Allocator.Persistent);
		SoilFertility = new NativeArray<float>(count, Allocator.Persistent);
		Elevation = new NativeArray<float>(count, Allocator.Persistent);
		TerrainTemperature = new NativeArray<float>(count, Allocator.Persistent);
		GroundWater = new NativeArray<float>(count, Allocator.Persistent);
		GroundWaterTemperature = new NativeArray<float>(count, Allocator.Persistent);
		FloraMass = new NativeArray<float>(count, Allocator.Persistent);
		FloraWater = new NativeArray<float>(count, Allocator.Persistent);
		FloraTemperature = new NativeArray<float>(count, Allocator.Persistent);
		FloraGlucose = new NativeArray<float>(count, Allocator.Persistent);
		IceTemperature = new NativeArray<float>(count, Allocator.Persistent);
		IceMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudTemperature = new NativeArray<float>(count, Allocator.Persistent);
		CloudDropletMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudVelocity = new NativeArray<float3>(count, Allocator.Persistent);
		LavaTemperature = new NativeArray<float>(count, Allocator.Persistent);
		LavaMass = new NativeArray<float>(count, Allocator.Persistent);
		CrustDepth = new NativeArray<float>(count, Allocator.Persistent);
		MagmaMass = new NativeArray<float>(count, Allocator.Persistent);
		AirTemperaturePotential = new NativeArray<float>[airLayers];
		AirVapor = new NativeArray<float>[airLayers];
		AirOxygen = new NativeArray<float>[airLayers];
		AirCarbonDioxide = new NativeArray<float>[airLayers];
		AirVelocity = new NativeArray<float3>[airLayers];
		Dust = new NativeArray<float>[airLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirTemperaturePotential[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVapor[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirOxygen[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirCarbonDioxide[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVelocity[i] = new NativeArray<float3>(count, Allocator.Persistent);
			Dust[i] = new NativeArray<float>(count, Allocator.Persistent);
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
		Roughness.Dispose();
		SoilFertility.Dispose();
		Elevation.Dispose();
		TerrainTemperature.Dispose();
		FloraMass.Dispose();
		FloraWater.Dispose();
		FloraGlucose.Dispose();
		FloraTemperature.Dispose();
		GroundWater.Dispose();
		GroundWaterTemperature.Dispose();
		IceTemperature.Dispose();
		IceMass.Dispose();
		CloudMass.Dispose();
		CloudTemperature.Dispose();
		CloudDropletMass.Dispose();
		CloudVelocity.Dispose();
		LavaMass.Dispose();
		LavaTemperature.Dispose();
		MagmaMass.Dispose();
		CrustDepth.Dispose();
		for (int i = 0; i < AirTemperaturePotential.Length; i++)
		{
			AirTemperaturePotential[i].Dispose();
			AirVapor[i].Dispose();
			AirOxygen[i].Dispose();
			AirCarbonDioxide[i].Dispose();
			AirVelocity[i].Dispose();
			Dust[i].Dispose();
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

public struct PlanetState {
	public int Ticks;
	public float Gravity;
	public float SpinSpeed;
	public float OrbitSpeed;
	public float GeothermalHeat;
	public float SolarRadiation;
	public float DistanceToSun;
	public float3 Rotation;
	public float3 Position;
	public float AngularSpeed;
}

public struct DependentState {

	public NativeArray<float> FloraCoverage;
	public NativeArray<float> FloraEnergy;
	public NativeArray<float> LavaEnergy;
	public NativeArray<float> LavaCoverage;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> IceEnergy;
	public NativeArray<float> SurfaceAirTemperatureAbsolute;
	public NativeArray<float> WindVerticalCloud;
	public NativeArray<float> SurfaceAreaAirIce;
	public NativeArray<float> SurfaceAreaAirWater;
	public NativeArray<float> SurfaceAreaAirFlora;
	public NativeArray<float> SurfaceAreaAirLava;
	public NativeArray<float> SurfaceAreaAirTerrain;
	public NativeArray<float> SurfaceAreaIceWater;
	public NativeArray<float> SurfaceAreaIceFlora;
	public NativeArray<float> SurfaceAreaIceLava;
	public NativeArray<float> SurfaceAreaIceTerrain;
	public NativeArray<float> SurfaceAreaWaterFlora;
	public NativeArray<float> SurfaceAreaWaterLava;
	public NativeArray<float> SurfaceAreaWaterTerrain;
	public NativeArray<float> SurfaceAreaFloraTerrain;
	public NativeArray<float> SurfaceAreaFloraLava;
	public NativeArray<float> SurfaceAreaLavaTerrain;
	public NativeArray<float>[] AirMass;
	public NativeArray<float>[] AirPressure;
	public NativeArray<float>[] AirPressureInverse;
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
	public NativeArray<float>[] LayerMiddle;
	public NativeArray<float> CloudAbsorptivity;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> DewPoint;
	public NativeArray<float> AirDensityCloud;

	public void Init(int count, int airLayers, int waterLayers)
	{
		FloraCoverage = new NativeArray<float>(count, Allocator.Persistent);
		FloraEnergy = new NativeArray<float>(count, Allocator.Persistent);
		IceCoverage = new NativeArray<float>(count, Allocator.Persistent);
		IceEnergy = new NativeArray<float>(count, Allocator.Persistent);
		LavaCoverage = new NativeArray<float>(count, Allocator.Persistent);
		LavaEnergy = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAirTemperatureAbsolute = new NativeArray<float>(count, Allocator.Persistent);
		DewPoint = new NativeArray<float>(count, Allocator.Persistent);
		WindVerticalCloud = new NativeArray<float>(count, Allocator.Persistent);
		AirDensityCloud = new NativeArray<float>(count, Allocator.Persistent);
		CloudAbsorptivity = new NativeArray<float>(count, Allocator.Persistent);
		CloudElevation = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaAirIce = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaAirWater = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaAirFlora = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaAirLava = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaAirTerrain = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaIceWater = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaIceFlora = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaIceLava = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaIceTerrain = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaWaterFlora = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaWaterLava = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaWaterTerrain = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaFloraLava = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaFloraTerrain = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaLavaTerrain = new NativeArray<float>(count, Allocator.Persistent);

		AirMass = new NativeArray<float>[airLayers];
		AirPressure = new NativeArray<float>[airLayers];
		AirPressureInverse = new NativeArray<float>[airLayers];
		AirHumidityAbsolute = new NativeArray<float>[airLayers];
		AirHumidityRelative = new NativeArray<float>[airLayers];
		LayerHeight = new NativeArray<float>[airLayers];
		LayerElevation = new NativeArray<float>[airLayers];
		LayerMiddle = new NativeArray<float>[airLayers];
		AirPotentialEnergy = new NativeArray<float>[airLayers];
		WaterDensity = new NativeArray<float>[waterLayers];
		WaterPressure = new NativeArray<float>[waterLayers];
		WaterLayerDepth = new NativeArray<float>[waterLayers];
		WaterLayerHeight = new NativeArray<float>[waterLayers];
		WaterCoverage = new NativeArray<float>[waterLayers];
		WaterPotentialEnergy = new NativeArray<float>[waterLayers];
		for (int i = 0; i < airLayers; i++)
		{
			AirMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirPressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirPressureInverse[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirHumidityAbsolute[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirHumidityRelative[i] = new NativeArray<float>(count, Allocator.Persistent);
			LayerHeight[i] = new NativeArray<float>(count, Allocator.Persistent);
			LayerElevation[i] = new NativeArray<float>(count, Allocator.Persistent);
			LayerMiddle[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirPotentialEnergy[i] = new NativeArray<float>(count, Allocator.Persistent);
		}
		for (int i = 0; i < waterLayers; i++)
		{
			WaterCoverage[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterDensity[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterPressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterLayerDepth[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterLayerHeight[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterPotentialEnergy[i] = new NativeArray<float>(count, Allocator.Persistent);
		}
	}

	public void Dispose()
	{
		FloraCoverage.Dispose();
		FloraEnergy.Dispose();
		IceCoverage.Dispose();
		IceEnergy.Dispose();
		LavaCoverage.Dispose();
		LavaEnergy.Dispose();
		SurfaceAirTemperatureAbsolute.Dispose();
		DewPoint.Dispose();
		WindVerticalCloud.Dispose();
		AirDensityCloud.Dispose();
		CloudElevation.Dispose();
		CloudAbsorptivity.Dispose();
		SurfaceAreaAirIce.Dispose();
		SurfaceAreaAirWater.Dispose();
		SurfaceAreaAirFlora.Dispose();
		SurfaceAreaAirLava.Dispose();
		SurfaceAreaAirTerrain.Dispose();
		SurfaceAreaIceWater.Dispose();
		SurfaceAreaIceFlora.Dispose();
		SurfaceAreaIceLava.Dispose();
		SurfaceAreaIceTerrain.Dispose();
		SurfaceAreaWaterFlora.Dispose();
		SurfaceAreaWaterLava.Dispose();
		SurfaceAreaWaterTerrain.Dispose();
		SurfaceAreaFloraLava.Dispose();
		SurfaceAreaFloraTerrain.Dispose();
		SurfaceAreaLavaTerrain.Dispose();

		for (int i = 0; i < AirPressure.Length; i++)
		{
			AirMass[i].Dispose();
			AirPressure[i].Dispose();
			AirPressureInverse[i].Dispose();
			AirHumidityRelative[i].Dispose();
			AirHumidityAbsolute[i].Dispose();
			LayerHeight[i].Dispose();
			LayerElevation[i].Dispose();
			LayerMiddle[i].Dispose();
			AirPotentialEnergy[i].Dispose();
		}
		for (int i = 0; i < WaterCoverage.Length; i++)
		{
			WaterCoverage[i].Dispose();
			WaterDensity[i].Dispose();
			WaterPressure[i].Dispose();
			WaterLayerDepth[i].Dispose();
			WaterLayerHeight[i].Dispose();
			WaterPotentialEnergy[i].Dispose();
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
	public double GlobalEnthalpyFlora;
	public double GlobalEnthalpyIce;
	public double GlobalEnthalpyCloud;
	public double GlobalEnthalpyWater;
	public double GlobalEnthalpyAir;
	public double GlobalEnthalpyGroundWater;
	public double GlobalEnthalpyDelta;
	public double GlobalEnthalpyDeltaTerrain;
	public double GlobalEnthalpyDeltaWater;
	public double GlobalEnthalpyDeltaAir;
	public double GlobalEnthalpyDeltaIce;
	public double GlobalEnthalpyDeltaCloud;
	public double GlobalEnthalpyDeltaFlora;
	public double GlobalEnthalpyDeltaGroundWater;
	public double GlobalTerrainTemperature;
	public float GlobalIceMass;
	public float GlobalSurfaceTemperature;
	public double GlobalAirTemperaturePotential;
	public float GlobalOceanCoverage;
	public double GlobalOceanMass;
	public float GlobalOceanTemperature;
	public float GlobalOceanSurfaceTemperature;
	public float GlobalSeaLevel;
	public float GlobalCarbonDioxide;
	public float GlobalOxygen;
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
	public NativeArray<float> EnthalpyFlora;
	public NativeArray<float> EnthalpyIce;
	public NativeArray<float> EnthalpyCloud;
	public NativeArray<float> EnthalpyGroundWater;
	public NativeArray<float> DustMass;
	public NativeArray<float>[] CarbonDioxidePercent;
	public NativeArray<float>[] OxygenPercent;
	public NativeArray<float>[] EnthalpyWater;
	public NativeArray<float>[] EnthalpyAir;
	public NativeArray<float>[] Salinity;
	public NativeArray<float>[] Pressure;
	public NativeArray<float3>[] PressureGradientForce;
	public NativeArray<float>[] ThermalDelta;
	public NativeArray<float>[] ConductionDelta;
	public NativeArray<float>[] SolarDelta;
	public NativeArray<float>[] LatentHeatDelta;
	public NativeArray<float>[] WindVertical;
	public NativeArray<SolarAbsorptivity>[] AbsorptionSolar;
	public NativeArray<ThermalAbsorptivity>[] AbsorptionThermal;

	public void Init(int count, int airLayers, int waterLayers, int totalLayers)
	{
		SolarRadiationAbsorbedSurface = new NativeArray<float>(count, Allocator.Persistent);
		Rainfall = new NativeArray<float>(count, Allocator.Persistent);
		CondensationCloud = new NativeArray<float>(count, Allocator.Persistent);
		CondensationGround = new NativeArray<float>(count, Allocator.Persistent);
		Evaporation = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyTerrain = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyFlora = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyIce = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyCloud = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyGroundWater = new NativeArray<float>(count, Allocator.Persistent);
		DustMass = new NativeArray<float>(count, Allocator.Persistent);

		Pressure = new NativeArray<float>[airLayers];
		PressureGradientForce = new NativeArray<float3>[airLayers];
		OxygenPercent = new NativeArray<float>[airLayers];
		CarbonDioxidePercent = new NativeArray<float>[airLayers];
		EnthalpyAir = new NativeArray<float>[airLayers];
		WindVertical = new NativeArray<float>[airLayers];
		AbsorptionSolar = new NativeArray<SolarAbsorptivity>[airLayers];
		AbsorptionThermal = new NativeArray<ThermalAbsorptivity>[airLayers];
		for (int i = 0; i < airLayers; i++)
		{
			Pressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			PressureGradientForce[i] = new NativeArray<float3>(count, Allocator.Persistent);
			EnthalpyAir[i] = new NativeArray<float>(count, Allocator.Persistent);
			WindVertical[i] = new NativeArray<float>(count, Allocator.Persistent);
			AbsorptionSolar[i] = new NativeArray<SolarAbsorptivity>(count, Allocator.Persistent);
			AbsorptionThermal[i] = new NativeArray<ThermalAbsorptivity>(count, Allocator.Persistent);
			OxygenPercent[i] = new NativeArray<float>(count, Allocator.Persistent);
			CarbonDioxidePercent[i] = new NativeArray<float>(count, Allocator.Persistent);
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
		EnthalpyCloud.Dispose();
		EnthalpyIce.Dispose();
		EnthalpyFlora.Dispose();
		EnthalpyGroundWater.Dispose();
		DustMass.Dispose();
		for (int i = 0; i < Pressure.Length; i++)
		{
			Pressure[i].Dispose();
			PressureGradientForce[i].Dispose();
			EnthalpyAir[i].Dispose();
			WindVertical[i].Dispose();
			AbsorptionSolar[i].Dispose();
			AbsorptionThermal[i].Dispose();
			OxygenPercent[i].Dispose();
			CarbonDioxidePercent[i].Dispose();
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


