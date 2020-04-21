
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using System.Collections.Generic;
using System;

public struct TempState {

	public NativeArray<float> FloraCoverage;
	public NativeArray<float> FloraEnergy;
	public NativeArray<float> LavaEnergy;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> IceEnergy;
	public NativeArray<float> SurfaceAirTemperatureAbsolute;
	public NativeArray<float> WindVerticalCloud;
	public NativeArray<float> SurfaceAreaAirIce;
	public NativeArray<float> SurfaceAreaAirWater;
	public NativeArray<float> SurfaceAreaAirFlora;
	public NativeArray<float> SurfaceAreaAirTerrain;
	public NativeArray<float> SurfaceAreaIceWater;
	public NativeArray<float> SurfaceAreaIceFlora;
	public NativeArray<float> SurfaceAreaIceTerrain;
	public NativeArray<float> SurfaceAreaWaterFlora;
	public NativeArray<float> SurfaceAreaWaterTerrain;
	public NativeArray<float> SurfaceAreaFloraTerrain;
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
	public NativeArray<float3> CloudVelocity;
	public NativeArray<float> DewPoint;
	public NativeArray<float> AirDensityCloud;
	public NativeArray<float> LavaDepth;


	public NativeArray<float> SolarRadiation;
	public NativeArray<float> AlbedoSlope;
	public NativeArray<float> CloudAlbedo;
	public NativeArray<SolarAbsorptivity>[] AbsorptivitySolar;
	public NativeArray<ThermalAbsorptivity>[] AbsorptivityThermal;
	public NativeArray<float>[] Emissivity;
	public NativeArray<float>[] SolarRadiationIn;
	public NativeArray<DiffusionAir>[] DiffusionAir;
	public NativeArray<DiffusionWater>[] DiffusionWater;
	public NativeArray<DiffusionAir>[] AdvectionAir;
	public NativeArray<DiffusionWater>[] AdvectionWater;
	public NativeArray<DiffusionCloud> DiffusionCloud;
	public NativeArray<DiffusionCloud> AdvectionCloud;
	public NativeArray<BarycentricValue> DestinationCloud;
	public NativeArray<BarycentricValueVertical>[] DestinationAir;
	public NativeArray<BarycentricValueVertical>[] DestinationWater;
	public NativeArray<float>[] DivergencePressureAir;
	public NativeArray<float3>[] AirAcceleration;
	public NativeArray<float3> CloudVelocityDeflected;
	public NativeArray<float> TerrainGradient;
	public NativeArray<float> WindFriction;
	public NativeArray<float3> WaterFriction;
	public NativeArray<float> ConductionAirIce;
	public NativeArray<float> ConductionAirWater;
	public NativeArray<float> ConductionAirFlora;
	public NativeArray<float> ConductionAirTerrain;
	public NativeArray<float> ConductionIceWater;
	public NativeArray<float> ConductionIceFlora;
	public NativeArray<float> ConductionIceTerrain;
	public NativeArray<float> ConductionFloraTerrain;
	public NativeArray<float>[] ConductionWaterTerrain;
	public NativeArray<float> FrozenTemperature;
	public NativeArray<float> FrozenMass;
	public NativeArray<float> SaltPlume;
	public NativeArray<float> EvaporationMassWater;
	public NativeArray<float> TemperaturePotentialFlora;
	public NativeArray<float> WaterConsumedByFlora;
	public NativeArray<float> FloraRespirationMassVapor;
	public NativeArray<float> FloraRespirationMassWater;
	public NativeArray<float> FloraMassDelta;
	public NativeArray<float> FloraWaterDelta;
	public NativeArray<float> FloraGlucoseDelta;
	public NativeArray<float> FloraDeath;
	public NativeArray<float> PlanktonMassDelta;
	public NativeArray<float> PlanktonGlucoseDelta;
	public NativeArray<float> PlanktonDeath;
	public NativeArray<float> SoilRespiration;
	public NativeArray<float> GeothermalRadiation;
	public NativeArray<float> GroundWaterFlowMass;
	public NativeArray<float> GroundWaterFlowTemperature;
	public NativeArray<float>[] DustUp;
	public NativeArray<float>[] DustDown;
	public NativeArray<float> IceMeltedMass;
	public NativeArray<float> LavaCrystalizedMass;
	public NativeArray<float> LavaEjected;
	public NativeArray<float> DustEjected;
	public NativeArray<float> CrustDelta;
	public NativeArray<float> WaterCarbonDelta;
	public NativeArray<float> AirCarbonDelta;
	public NativeArray<float> OxygenDelta;
	public NativeArray<float>[] DivergenceAir;
	public NativeArray<float> DisplaySolarRadiation;

	public NativeArray<float>[] ThermalRadiationDelta;
	public NativeArray<float>[] ThermalRadiationTransmittedUp;
	public NativeArray<float>[] ThermalRadiationTransmittedDown;
	public NativeArray<float>[] WindowRadiationTransmittedUp;
	public NativeArray<float>[] WindowRadiationTransmittedDown;
	public NativeArray<float>[] CondensationGroundMass;
	public NativeArray<float>[] CondensationCloudMass;
	public NativeArray<float>[] SolarReflected;
	public NativeArray<float>[] LatentHeat;
	public NativeArray<float> ConductionWaterTerrainTotal;
	public NativeArray<float> CloudEvaporationMass;
	public NativeArray<float> DropletDelta;
	public NativeArray<float> PrecipitationMass;
	public NativeArray<float> PrecipitationTemperature;
	public NativeArray<float> AtmosphericWindowUp;
	public NativeArray<float> AtmosphericWindowDown;

	public NativeArray<float> OutgoingFlowWater;
	public NativeArray<float> OutgoingFlowLava;
	public NativeArray<float> FlowPercentWater;
	public NativeArray<float> FlowPercentLava;
	public NativeArray<DiffusionLava> DiffusionLava;



	public void Init(int count, ref WorldData worldData)
	{
		FloraCoverage = new NativeArray<float>(count, Allocator.Persistent);
		FloraEnergy = new NativeArray<float>(count, Allocator.Persistent);
		IceCoverage = new NativeArray<float>(count, Allocator.Persistent);
		IceEnergy = new NativeArray<float>(count, Allocator.Persistent);
		LavaEnergy = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAirTemperatureAbsolute = new NativeArray<float>(count, Allocator.Persistent);
		DewPoint = new NativeArray<float>(count, Allocator.Persistent);
		WindVerticalCloud = new NativeArray<float>(count, Allocator.Persistent);
		AirDensityCloud = new NativeArray<float>(count, Allocator.Persistent);
		CloudAbsorptivity = new NativeArray<float>(count, Allocator.Persistent);
		CloudElevation = new NativeArray<float>(count, Allocator.Persistent);
		CloudVelocity = new NativeArray<float3>(count, Allocator.Persistent);
		SurfaceAreaAirIce = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaAirWater = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaAirFlora = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaAirTerrain = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaIceWater = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaIceFlora = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaIceTerrain = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaWaterFlora = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaWaterTerrain = new NativeArray<float>(count, Allocator.Persistent);
		SurfaceAreaFloraTerrain = new NativeArray<float>(count, Allocator.Persistent);
		LavaDepth = new NativeArray<float>(count, Allocator.Persistent);

		AirMass = new NativeArray<float>[worldData.AirLayers];
		AirPressure = new NativeArray<float>[worldData.AirLayers];
		AirPressureInverse = new NativeArray<float>[worldData.AirLayers];
		AirHumidityAbsolute = new NativeArray<float>[worldData.AirLayers];
		AirHumidityRelative = new NativeArray<float>[worldData.AirLayers];
		LayerHeight = new NativeArray<float>[worldData.AirLayers];
		LayerElevation = new NativeArray<float>[worldData.AirLayers];
		LayerMiddle = new NativeArray<float>[worldData.AirLayers];
		AirPotentialEnergy = new NativeArray<float>[worldData.AirLayers];
		WaterDensity = new NativeArray<float>[worldData.WaterLayers];
		WaterPressure = new NativeArray<float>[worldData.WaterLayers];
		WaterLayerDepth = new NativeArray<float>[worldData.WaterLayers];
		WaterLayerHeight = new NativeArray<float>[worldData.WaterLayers];
		WaterCoverage = new NativeArray<float>[worldData.WaterLayers];
		WaterPotentialEnergy = new NativeArray<float>[worldData.WaterLayers];
		for (int i = 0; i < worldData.AirLayers; i++)
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
		for (int i = 0; i < worldData.WaterLayers; i++)
		{
			WaterCoverage[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterDensity[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterPressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterLayerDepth[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterLayerHeight[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterPotentialEnergy[i] = new NativeArray<float>(count, Allocator.Persistent);
		}



		SolarRadiation = new NativeArray<float>(count, Allocator.Persistent);
		AlbedoSlope = new NativeArray<float>(count, Allocator.Persistent);
		CloudAlbedo = new NativeArray<float>(count, Allocator.Persistent);
		Emissivity = new NativeArray<float>[worldData.LayerCount];
		SolarRadiationIn = new NativeArray<float>[worldData.LayerCount];
		for (int i = 0; i < worldData.LayerCount; i++)
		{
			Emissivity[i] = new NativeArray<float>(count, Allocator.Persistent);
			SolarRadiationIn[i] = new NativeArray<float>(count, Allocator.Persistent);
		}
		DiffusionAir = new NativeArray<DiffusionAir>[worldData.AirLayers];
		AdvectionAir = new NativeArray<DiffusionAir>[worldData.AirLayers];
		DestinationAir = new NativeArray<BarycentricValueVertical>[worldData.AirLayers];
		DivergencePressureAir = new NativeArray<float>[worldData.AirLayers];
		AirAcceleration = new NativeArray<float3>[worldData.AirLayers];
		AbsorptivitySolar = new NativeArray<SolarAbsorptivity>[worldData.AirLayers];
		AbsorptivityThermal = new NativeArray<ThermalAbsorptivity>[worldData.AirLayers];
		DustUp = new NativeArray<float>[worldData.AirLayers];
		DustDown = new NativeArray<float>[worldData.AirLayers];
		DivergenceAir = new NativeArray<float>[worldData.AirLayers];
		for (int i = 0; i < worldData.AirLayers; i++)
		{
			DiffusionAir[i] = new NativeArray<DiffusionAir>(count, Allocator.Persistent);
			AdvectionAir[i] = new NativeArray<DiffusionAir>(count, Allocator.Persistent);
			DestinationAir[i] = new NativeArray<BarycentricValueVertical>(count, Allocator.Persistent);
			DivergencePressureAir[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirAcceleration[i] = new NativeArray<float3>(count, Allocator.Persistent);
			AbsorptivitySolar[i] = new NativeArray<SolarAbsorptivity>(count, Allocator.Persistent);
			AbsorptivityThermal[i] = new NativeArray<ThermalAbsorptivity>(count, Allocator.Persistent);
			DustUp[i] = new NativeArray<float>(count, Allocator.Persistent);
			DustDown[i] = new NativeArray<float>(count, Allocator.Persistent);
			DivergenceAir[i] = new NativeArray<float>(count, Allocator.Persistent);
		}
		DiffusionWater = new NativeArray<DiffusionWater>[worldData.WaterLayers];
		AdvectionWater = new NativeArray<DiffusionWater>[worldData.WaterLayers];
		DestinationWater = new NativeArray<BarycentricValueVertical>[worldData.WaterLayers];
		ConductionWaterTerrain = new NativeArray<float>[worldData.WaterLayers];
		for (int i = 0; i < worldData.WaterLayers; i++)
		{
			DiffusionWater[i] = new NativeArray<DiffusionWater>(count, Allocator.Persistent);
			AdvectionWater[i] = new NativeArray<DiffusionWater>(count, Allocator.Persistent);
			DestinationWater[i] = new NativeArray<BarycentricValueVertical>(count, Allocator.Persistent);
			ConductionWaterTerrain[i] = new NativeArray<float>(count, Allocator.Persistent);
		}
		FrozenTemperature = new NativeArray<float>(count, Allocator.Persistent);
		FrozenMass = new NativeArray<float>(count, Allocator.Persistent);
		SaltPlume = new NativeArray<float>(count, Allocator.Persistent);
		EvaporationMassWater = new NativeArray<float>(count, Allocator.Persistent);
		FloraRespirationMassVapor = new NativeArray<float>(count, Allocator.Persistent);
		FloraRespirationMassWater = new NativeArray<float>(count, Allocator.Persistent);
		WaterConsumedByFlora = new NativeArray<float>(count, Allocator.Persistent);
		FloraMassDelta = new NativeArray<float>(count, Allocator.Persistent);
		FloraWaterDelta = new NativeArray<float>(count, Allocator.Persistent);
		FloraGlucoseDelta = new NativeArray<float>(count, Allocator.Persistent);
		FloraDeath = new NativeArray<float>(count, Allocator.TempJob);
		PlanktonMassDelta = new NativeArray<float>(count, Allocator.Persistent);
		PlanktonGlucoseDelta = new NativeArray<float>(count, Allocator.Persistent);
		PlanktonDeath = new NativeArray<float>(count, Allocator.TempJob);
		WindFriction = new NativeArray<float>(count, Allocator.Persistent);
		WaterFriction = new NativeArray<float3>(count, Allocator.Persistent);
		DiffusionCloud = new NativeArray<DiffusionCloud>(count, Allocator.Persistent);
		AdvectionCloud = new NativeArray<DiffusionCloud>(count, Allocator.Persistent);
		DestinationCloud = new NativeArray<BarycentricValue>(count, Allocator.Persistent);
		CloudVelocityDeflected = new NativeArray<float3>(count, Allocator.Persistent);
		ConductionAirIce = new NativeArray<float>(count, Allocator.Persistent);
		ConductionAirWater = new NativeArray<float>(count, Allocator.Persistent);
		ConductionAirFlora = new NativeArray<float>(count, Allocator.Persistent);
		ConductionAirTerrain = new NativeArray<float>(count, Allocator.Persistent);
		ConductionIceWater = new NativeArray<float>(count, Allocator.Persistent);
		ConductionIceFlora = new NativeArray<float>(count, Allocator.Persistent);
		ConductionIceTerrain = new NativeArray<float>(count, Allocator.Persistent);
		ConductionFloraTerrain = new NativeArray<float>(count, Allocator.Persistent);
		DisplaySolarRadiation = new NativeArray<float>(count, Allocator.Persistent);
		GeothermalRadiation = new NativeArray<float>(count, Allocator.Persistent);
		TerrainGradient = new NativeArray<float>(count * 6, Allocator.Persistent);
		GroundWaterFlowMass = new NativeArray<float>(count, Allocator.Persistent);
		GroundWaterFlowTemperature = new NativeArray<float>(count, Allocator.Persistent);
		IceMeltedMass = new NativeArray<float>(count, Allocator.TempJob);
		LavaCrystalizedMass = new NativeArray<float>(count, Allocator.TempJob);
		LavaEjected = new NativeArray<float>(count, Allocator.TempJob);
		DustEjected = new NativeArray<float>(count, Allocator.TempJob);
		CrustDelta = new NativeArray<float>(count, Allocator.TempJob);
		AirCarbonDelta = new NativeArray<float>(count, Allocator.TempJob);
		OxygenDelta = new NativeArray<float>(count, Allocator.TempJob);
		SoilRespiration = new NativeArray<float>(count, Allocator.TempJob);
		WaterCarbonDelta = new NativeArray<float>(count, Allocator.TempJob);


		ThermalRadiationDelta = new NativeArray<float>[worldData.LayerCount];
		ThermalRadiationTransmittedUp = new NativeArray<float>[worldData.LayerCount];
		ThermalRadiationTransmittedDown = new NativeArray<float>[worldData.LayerCount];
		WindowRadiationTransmittedUp = new NativeArray<float>[worldData.LayerCount];
		WindowRadiationTransmittedDown = new NativeArray<float>[worldData.LayerCount];
		CondensationGroundMass = new NativeArray<float>[worldData.LayerCount];
		CondensationCloudMass = new NativeArray<float>[worldData.LayerCount];
		SolarReflected = new NativeArray<float>[worldData.LayerCount];
		LatentHeat = new NativeArray<float>[worldData.LayerCount];
		ConductionWaterTerrainTotal = new NativeArray<float>(count, Allocator.TempJob);
		CloudEvaporationMass = new NativeArray<float>(count, Allocator.TempJob);
		DropletDelta = new NativeArray<float>(count, Allocator.TempJob);
		PrecipitationMass = new NativeArray<float>(count, Allocator.TempJob);
		PrecipitationTemperature = new NativeArray<float>(count, Allocator.TempJob);
		AtmosphericWindowUp = new NativeArray<float>(count, Allocator.TempJob);
		AtmosphericWindowDown = new NativeArray<float>(count, Allocator.TempJob);
		for (int i = 0; i < worldData.LayerCount; i++)
		{
			LatentHeat[i] = new NativeArray<float>(count, Allocator.TempJob);
			SolarReflected[i] = new NativeArray<float>(count, Allocator.TempJob);
			ThermalRadiationDelta[i] = new NativeArray<float>(count, Allocator.TempJob);
			ThermalRadiationTransmittedUp[i] = new NativeArray<float>(count, Allocator.TempJob);
			ThermalRadiationTransmittedDown[i] = new NativeArray<float>(count, Allocator.TempJob);
			WindowRadiationTransmittedUp[i] = new NativeArray<float>(count, Allocator.TempJob);
			WindowRadiationTransmittedDown[i] = new NativeArray<float>(count, Allocator.TempJob);
			CondensationGroundMass[i] = new NativeArray<float>(count, Allocator.TempJob);
			CondensationCloudMass[i] = new NativeArray<float>(count, Allocator.TempJob);
		}


		OutgoingFlowWater = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.TempJob);
		OutgoingFlowLava = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.TempJob);
		FlowPercentWater = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.TempJob);
		FlowPercentLava = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.TempJob);
		DiffusionLava = new NativeArray<DiffusionLava>(count * StaticState.MaxNeighbors, Allocator.TempJob);

	}

	public void Dispose(ref WorldData worldData, Action completeClear)
	{
		completeClear?.Invoke();

		FloraCoverage.Dispose();
		FloraEnergy.Dispose();
		IceCoverage.Dispose();
		IceEnergy.Dispose();
		LavaEnergy.Dispose();
		SurfaceAirTemperatureAbsolute.Dispose();
		DewPoint.Dispose();
		WindVerticalCloud.Dispose();
		AirDensityCloud.Dispose();
		CloudElevation.Dispose();
		CloudVelocity.Dispose();
		CloudAbsorptivity.Dispose();
		SurfaceAreaAirIce.Dispose();
		SurfaceAreaAirWater.Dispose();
		SurfaceAreaAirFlora.Dispose();
		SurfaceAreaAirTerrain.Dispose();
		SurfaceAreaIceWater.Dispose();
		SurfaceAreaIceFlora.Dispose();
		SurfaceAreaIceTerrain.Dispose();
		SurfaceAreaWaterFlora.Dispose();
		SurfaceAreaWaterTerrain.Dispose();
		SurfaceAreaFloraTerrain.Dispose();
		LavaDepth.Dispose();

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




		SolarRadiation.Dispose();
		AlbedoSlope.Dispose();
		CloudAlbedo.Dispose();
		for (int i = 0; i < Emissivity.Length; i++)
		{
			Emissivity[i].Dispose();
			SolarRadiationIn[i].Dispose();
		}
		for (int i = 0; i < worldData.AirLayers; i++)
		{
			DiffusionAir[i].Dispose();
			AdvectionAir[i].Dispose();
			DestinationAir[i].Dispose();
			DivergencePressureAir[i].Dispose();
			AirAcceleration[i].Dispose();
			AbsorptivitySolar[i].Dispose();
			AbsorptivityThermal[i].Dispose();
			DustUp[i].Dispose();
			DustDown[i].Dispose();
			DivergenceAir[i].Dispose();
		}
		for (int i = 0; i < worldData.WaterLayers; i++)
		{
			DiffusionWater[i].Dispose();
			AdvectionWater[i].Dispose();
			DestinationWater[i].Dispose();
			ConductionWaterTerrain[i].Dispose();
		}
		FrozenTemperature.Dispose();
		FrozenMass.Dispose();
		SaltPlume.Dispose();
		EvaporationMassWater.Dispose();
		FloraRespirationMassVapor.Dispose();
		FloraRespirationMassWater.Dispose();
		WaterConsumedByFlora.Dispose();
		FloraMassDelta.Dispose();
		FloraWaterDelta.Dispose();
		FloraGlucoseDelta.Dispose();
		FloraDeath.Dispose();
		PlanktonMassDelta.Dispose();
		PlanktonGlucoseDelta.Dispose();
		PlanktonDeath.Dispose();
		WindFriction.Dispose();
		WaterFriction.Dispose();
		DiffusionCloud.Dispose();
		AdvectionCloud.Dispose();
		DestinationCloud.Dispose();
		CloudVelocityDeflected.Dispose();
		ConductionAirIce.Dispose();
		ConductionAirWater.Dispose();
		ConductionAirFlora.Dispose();
		ConductionAirTerrain.Dispose();
		ConductionIceWater.Dispose();
		ConductionIceFlora.Dispose();
		ConductionIceTerrain.Dispose();
		ConductionFloraTerrain.Dispose();
		GeothermalRadiation.Dispose();
		TerrainGradient.Dispose();
		GroundWaterFlowMass.Dispose();
		GroundWaterFlowTemperature.Dispose();
		IceMeltedMass.Dispose();
		LavaCrystalizedMass.Dispose();
		LavaEjected.Dispose();
		DustEjected.Dispose();
		CrustDelta.Dispose();
		AirCarbonDelta.Dispose();
		OxygenDelta.Dispose();
		SoilRespiration.Dispose();
		WaterCarbonDelta.Dispose();

		DisplaySolarRadiation.Dispose();

		AtmosphericWindowUp.Dispose();
		AtmosphericWindowDown.Dispose();
		DropletDelta.Dispose();
		PrecipitationMass.Dispose();
		PrecipitationTemperature.Dispose();
		CloudEvaporationMass.Dispose();
		ConductionWaterTerrainTotal.Dispose();
		for (int i = 0; i < worldData.LayerCount; i++)
		{
			ThermalRadiationDelta[i].Dispose();
			ThermalRadiationTransmittedUp[i].Dispose();
			ThermalRadiationTransmittedDown[i].Dispose();
			WindowRadiationTransmittedUp[i].Dispose();
			WindowRadiationTransmittedDown[i].Dispose();
			SolarReflected[i].Dispose();
			LatentHeat[i].Dispose();
			CondensationGroundMass[i].Dispose();
			CondensationCloudMass[i].Dispose();
		}


		OutgoingFlowWater.Dispose();
		OutgoingFlowLava.Dispose();
		FlowPercentWater.Dispose();
		FlowPercentLava.Dispose();
		DiffusionLava.Dispose();


	}

	public Action Clear(int cellCount, ref WorldData worldData)
	{
		NativeList<JobHandle> memsetHandles = new NativeList<JobHandle>(Allocator.TempJob);
		Utils.MemsetArray(memsetHandles, cellCount, ConductionWaterTerrainTotal, 0);
		Utils.MemsetArray(memsetHandles, cellCount, CloudEvaporationMass, 0);
		Utils.MemsetArray(memsetHandles, cellCount, DropletDelta, 0);
		Utils.MemsetArray(memsetHandles, cellCount, PrecipitationMass, 0);
		Utils.MemsetArray(memsetHandles, cellCount, PrecipitationTemperature, 0);
		Utils.MemsetArray(memsetHandles, cellCount, AtmosphericWindowUp, 0);
		Utils.MemsetArray(memsetHandles, cellCount, AtmosphericWindowDown, 0);
		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			Utils.MemsetArray(memsetHandles, cellCount, DivergencePressureAir[i], 0);
		}
		for (int i = 0; i < worldData.LayerCount; i++)
		{
			Utils.MemsetArray(memsetHandles, cellCount, LatentHeat[i], 0);
			Utils.MemsetArray(memsetHandles, cellCount, SolarReflected[i], 0);
			Utils.MemsetArray(memsetHandles, cellCount, ThermalRadiationDelta[i], 0);
			Utils.MemsetArray(memsetHandles, cellCount, ThermalRadiationTransmittedUp[i], 0);
			Utils.MemsetArray(memsetHandles, cellCount, ThermalRadiationTransmittedDown[i], 0);
			Utils.MemsetArray(memsetHandles, cellCount, WindowRadiationTransmittedUp[i], 0);
			Utils.MemsetArray(memsetHandles, cellCount, WindowRadiationTransmittedDown[i], 0);
			Utils.MemsetArray(memsetHandles, cellCount, CondensationGroundMass[i], 0);
			Utils.MemsetArray(memsetHandles, cellCount, CondensationCloudMass[i], 0);
		}

		return () =>
		{
			if (memsetHandles.IsCreated)
			{
				JobHandle.CompleteAll(memsetHandles);
				memsetHandles.Dispose();
			}
		};

	}

	public static JobHandle Update(JobHelper jobHelper, ref SimState state, ref TempState tempState, ref WorldData worldData, JobHandle dependencies, List<NativeArray<float>> arraysToDispose)
	{
		var waterMassTotal = new NativeArray<float>(state.IceMass, Allocator.TempJob);
		var airMassTotal = new NativeArray<float>(state.IceMass.Length, Allocator.TempJob);
		var standardLayerElevation = new NativeArray<float>(state.IceMass.Length, Allocator.TempJob);
		arraysToDispose.Add(waterMassTotal);
		arraysToDispose.Add(airMassTotal);
		arraysToDispose.Add(standardLayerElevation);
		for (int j = worldData.WaterLayers - 2; j >= 1; j--)
		{
			dependencies = jobHelper.Schedule(new UpdateWaterDepthJob()
			{
				Density = tempState.WaterDensity[j],
				LayerDepth = tempState.WaterLayerDepth[j],
				LayerHeight = tempState.WaterLayerHeight[j],
				WaterMassTotal = waterMassTotal,
				WaterCoverage = tempState.WaterCoverage[j],
				PotentialEnergy = tempState.WaterPotentialEnergy[j],
				Pressure = tempState.WaterPressure[j],

				Temperature = state.WaterTemperature[j],
				SaltMass = state.SaltMass[j],
				WaterMass = state.WaterMass[j],
				Roughness = state.Roughness,
				LayerDepthUp = tempState.WaterLayerDepth[j + 1],
				LayerHeightUp = tempState.WaterLayerHeight[j + 1],
				LayerHeightDown = tempState.WaterLayerHeight[j - 1],
				WaterDensityPerDegree = worldData.WaterDensityPerDegree,
				WaterDensityPerSalinity = worldData.WaterDensityPerSalinity,
				Gravity = state.PlanetState.Gravity,

			}, dependencies);
		}



		var surfaceElevationJob = jobHelper.Schedule(new UpdateTempStateJob()
		{
			IceEnergy = tempState.IceEnergy,
			FloraEnergy = tempState.FloraEnergy,
			LavaEnergy = tempState.LavaEnergy,
			LavaDepth = tempState.LavaDepth,
			SurfaceElevation = tempState.LayerElevation[1],

			WaterDepth = tempState.WaterLayerDepth[1],
			Elevation = state.Elevation,
			FloraMass = state.FloraMass,
			FloraWater = state.FloraWater,
			FloraTemperature = state.FloraTemperature,
			LavaMass = state.LavaMass,
			LavaTemperature = state.LavaTemperature,
			IceMass = state.IceMass,
			IceTemperature = state.IceTemperature,
			LavaToRockMassAdjustment = worldData.LavaToRockMassAdjustment,
		}, dependencies);
		dependencies = JobHandle.CombineDependencies(surfaceElevationJob, dependencies);
		dependencies = Utils.MemCopy(standardLayerElevation, tempState.LayerElevation[1], dependencies);

		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			float minumumHeight;
			float columnPercent;
			if (j == 1)
			{
				minumumHeight = worldData.BoundaryZoneElevation;
				columnPercent = 0;
			}
			else
			{
				minumumHeight = 0;
				columnPercent = 1.0f / (worldData.AirLayers - 1 - j);
			}
			dependencies = jobHelper.Schedule(new UpdateAirLayerHeightsJob()
			{
				StandardLayerElevation = standardLayerElevation,
				LayerHeight = tempState.LayerHeight[j],
				UpLayerElevation = tempState.LayerElevation[j + 1],
				AirMass = tempState.AirMass[j],
				LayerMiddle = tempState.LayerMiddle[j],

				LayerElevation = tempState.LayerElevation[j],
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				TropopauseElevation = worldData.TropopauseElevation,
				MinimumHeight = minumumHeight,
				ColumnPercent = columnPercent,
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
		}

		dependencies = jobHelper.Schedule(new UpdateStratosphereJob()
		{
			StratosphereMass = airMassTotal,

			TropopauseElevation = tempState.LayerElevation[worldData.AirLayers - 2],
			TropopauseHeight = tempState.LayerHeight[worldData.AirLayers - 2],
			Gravity = state.PlanetState.Gravity
		}, dependencies);


		for (int j = worldData.AirLayers - 2; j > 0; j--)
		{
			dependencies = jobHelper.Schedule(new UpdateAirPressureJob()
			{
				Pressure = tempState.AirPressure[j],
				PressureInverse = tempState.AirPressureInverse[j],
				AirMassTotal = airMassTotal,
				RelativeHumidity = tempState.AirHumidityRelative[j],
				AbsoluteHumidity = tempState.AirHumidityAbsolute[j],
				AirMass = tempState.AirMass[j],
				PotentialEnergy = tempState.AirPotentialEnergy[j],

				CloudMass = state.CloudMass,
				VaporMass = state.AirVapor[j],
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				LayerElevation = tempState.LayerElevation[j],
				LayerMiddle = tempState.LayerMiddle[j],
				LayerHeight = tempState.LayerHeight[j],
				SurfaceElevation = tempState.LayerElevation[1],
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
		}

		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
			dependencies = jobHelper.Schedule(new UpdateCloudJob()
			{
				CloudVelocity = tempState.CloudVelocity,
				CloudElevation = tempState.CloudElevation,
				AirDensityCloud = tempState.AirDensityCloud,
				DewPoint = tempState.DewPoint,

				AirMass = tempState.AirMass[j],
				VaporMass = state.AirVapor[j],
				SurfaceElevation = tempState.LayerElevation[1],
				RelativeHumidity = tempState.AirHumidityRelative[j],
				Pressure = tempState.AirPressure[j],
				AirTemperaturePotential = state.AirTemperaturePotential[j],
				LayerMiddle = tempState.LayerMiddle[j],
				LayerElevation = tempState.LayerElevation[j],
				LayerElevationAbove = tempState.LayerElevation[j + 1],
				LayerElevationBelow = tempState.LayerElevation[j - 1],
				LayerHeight = tempState.LayerHeight[j],
				LayerHeightAbove = tempState.LayerHeight[j + 1],
				LayerHeightBelow = tempState.LayerHeight[j - 1],
				AirVelocity = state.AirVelocity[j],
				AirVelocityBelow = state.AirVelocity[j - 1],
				AirVelocityAbove = state.AirVelocity[j + 1],
				Gravity = state.PlanetState.Gravity,
				IsTop = j == worldData.AirLayers - 2,
				IsBottom = j == 1,

			}, dependencies);
		}

		var surfaceStateJobHandle = jobHelper.Schedule(new UpdateSurfaceStateJob()
		{
			SurfaceAirTemperatureAbsolute = tempState.SurfaceAirTemperatureAbsolute,
			SurfaceAreaAirFlora = tempState.SurfaceAreaAirFlora,
			SurfaceAreaAirIce = tempState.SurfaceAreaAirIce,
			SurfaceAreaAirTerrain = tempState.SurfaceAreaAirTerrain,
			SurfaceAreaAirWater = tempState.SurfaceAreaAirWater,
			SurfaceAreaIceFlora = tempState.SurfaceAreaIceFlora,
			SurfaceAreaIceTerrain = tempState.SurfaceAreaIceTerrain,
			SurfaceAreaIceWater = tempState.SurfaceAreaIceWater,
			SurfaceAreaWaterFlora = tempState.SurfaceAreaWaterFlora,
			SurfaceAreaWaterTerrain = tempState.SurfaceAreaWaterTerrain,
			SurfaceAreaFloraTerrain = tempState.SurfaceAreaFloraTerrain,
			FloraCoverage = tempState.FloraCoverage,
			IceCoverage = tempState.IceCoverage,

			WaterCoverage = tempState.WaterCoverage[worldData.WaterLayers - 2],
			FloraMass = state.FloraMass,
			IceMass = state.IceMass,
			AirTemperaturePotential = state.AirTemperaturePotential[1],
			SurfaceLayerElevation = tempState.LayerElevation[1],
			inverseFullCoverageFloraMass = 1.0f / worldData.FullCoverageFlora,
			inverseFullCoverageIceMass = 1.0f / (worldData.FullCoverageIce * WorldData.MassIce),
			FloraAirSurfaceArea = worldData.FloraAirSurfaceArea,
			Roughness = state.Roughness
		}, surfaceElevationJob);
		dependencies = JobHandle.CombineDependencies(dependencies, surfaceStateJobHandle);


		return dependencies;
	}


	public static JobHandle UpdateWaterDepth(JobHelper jobHelper, ref SimState state, ref TempState tempState, ref WorldData worldData, JobHandle dependencies, List<NativeArray<float>> arraysToDispose)
	{
		var waterMassTotal = new NativeArray<float>(state.IceMass.Length, Allocator.TempJob);
		arraysToDispose.Add(waterMassTotal);
		for (int j = worldData.WaterLayers - 2; j >= 1; j--)
		{
			dependencies = jobHelper.Schedule(new UpdateWaterDepthJob()
			{
				Density = tempState.WaterDensity[j],
				LayerDepth = tempState.WaterLayerDepth[j],
				LayerHeight = tempState.WaterLayerHeight[j],
				WaterMassTotal = waterMassTotal,
				WaterCoverage = tempState.WaterCoverage[j],
				PotentialEnergy = tempState.WaterPotentialEnergy[j],
				Pressure = tempState.WaterPressure[j],

				Temperature = state.WaterTemperature[j],
				SaltMass = state.SaltMass[j],
				WaterMass = state.WaterMass[j],
				Roughness = state.Roughness,
				LayerDepthUp = tempState.WaterLayerDepth[j + 1],
				LayerHeightUp = tempState.WaterLayerHeight[j + 1],
				LayerHeightDown = tempState.WaterLayerHeight[j - 1],
				WaterDensityPerDegree = worldData.WaterDensityPerDegree,
				WaterDensityPerSalinity = worldData.WaterDensityPerSalinity,
				Gravity = state.PlanetState.Gravity,

			}, dependencies);
		}

		dependencies = jobHelper.Schedule(new UpdateTempStateJob()
		{
			IceEnergy = tempState.IceEnergy,
			FloraEnergy = tempState.FloraEnergy,
			LavaEnergy = tempState.LavaEnergy,
			LavaDepth = tempState.LavaDepth,
			SurfaceElevation = tempState.LayerElevation[1],

			WaterDepth = tempState.WaterLayerDepth[1],
			Elevation = state.Elevation,
			FloraMass = state.FloraMass,
			FloraWater = state.FloraWater,
			FloraTemperature = state.FloraTemperature,
			LavaMass = state.LavaMass,
			LavaTemperature = state.LavaTemperature,
			IceMass = state.IceMass,
			IceTemperature = state.IceTemperature,
			LavaToRockMassAdjustment = worldData.LavaToRockMassAdjustment,
		}, dependencies);
		dependencies.Complete();

		return dependencies;
	}


	[BurstCompile]
	public struct UpdateTempStateJob : IJobParallelFor {
		public NativeArray<float> SurfaceElevation;
		public NativeArray<float> IceEnergy;
		public NativeArray<float> FloraEnergy;
		public NativeArray<float> LavaEnergy;
		public NativeArray<float> LavaDepth;
		[ReadOnly] public NativeArray<float> WaterDepth;
		[ReadOnly] public NativeArray<float> FloraMass;
		[ReadOnly] public NativeArray<float> FloraWater;
		[ReadOnly] public NativeArray<float> FloraTemperature;
		[ReadOnly] public NativeArray<float> IceMass;
		[ReadOnly] public NativeArray<float> IceTemperature;
		[ReadOnly] public NativeArray<float> LavaMass;
		[ReadOnly] public NativeArray<float> LavaTemperature;
		[ReadOnly] public NativeArray<float> Elevation;
		[ReadOnly] public float LavaToRockMassAdjustment;
		public void Execute(int i)
		{
			float iceMass = IceMass[i];
			float lavaMass = LavaMass[i];
			SurfaceElevation[i] = Elevation[i] + WaterDepth[i] + iceMass / WorldData.MassIce + lavaMass * LavaToRockMassAdjustment / WorldData.MassLava;

			IceEnergy[i] = WorldData.SpecificHeatIce * iceMass * IceTemperature[i];
			LavaEnergy[i] = WorldData.SpecificHeatLava * lavaMass * LavaTemperature[i];
			LavaDepth[i] = lavaMass * LavaToRockMassAdjustment / WorldData.MassLava;
			FloraEnergy[i] = (WorldData.SpecificHeatFlora * FloraMass[i] + WorldData.SpecificHeatWater * FloraWater[i]) * FloraTemperature[i];
		}

	}


	[BurstCompile]
	public struct UpdateAirLayerHeightsJob : IJobParallelFor {
		public NativeArray<float> AirMass;
		public NativeArray<float> UpLayerElevation;
		public NativeArray<float> LayerHeight;
		public NativeArray<float> LayerMiddle;
		public NativeArray<float> StandardLayerElevation;

		[ReadOnly] public NativeArray<float> AirTemperaturePotential;
		[ReadOnly] public NativeArray<float> LayerElevation;
		[ReadOnly] public float Gravity;
		[ReadOnly] public float TropopauseElevation;
		[ReadOnly] public float MinimumHeight;
		[ReadOnly] public float ColumnPercent;

		public void Execute(int i)
		{
			float standardLayerHeight = math.max(MinimumHeight, (TropopauseElevation - StandardLayerElevation[i]) * ColumnPercent);
			float standardLayerElevation = StandardLayerElevation[i];
			float airMass = (Atmosphere.GetStandardPressureAtElevation(standardLayerElevation, WorldData.StandardTemperature, Gravity) - Atmosphere.GetStandardPressureAtElevation(standardLayerElevation + standardLayerHeight, WorldData.StandardTemperature, Gravity)) / Gravity;
			StandardLayerElevation[i] = standardLayerElevation + standardLayerHeight;
			AirMass[i] = airMass;

			float airTemperaturePotential = AirTemperaturePotential[i];
			float layerElevation = LayerElevation[i];
			float layerHeight = standardLayerHeight * AirTemperaturePotential[i] / WorldData.StandardTemperature;
			LayerHeight[i] = layerHeight;
			UpLayerElevation[i] = layerElevation + layerHeight;
			LayerMiddle[i] = layerElevation + layerHeight / 2;
		}
	}

	[BurstCompile]
	public struct UpdateStratosphereJob : IJobParallelFor {
		public NativeArray<float> StratosphereMass;
		[ReadOnly] public NativeArray<float> TropopauseElevation;
		[ReadOnly] public NativeArray<float> TropopauseHeight;
		[ReadOnly] public float Gravity;
		public void Execute(int i)
		{
			// TODO: rebalance stratosphere mass
			float tropopauseElevation = TropopauseElevation[i] + TropopauseHeight[i];
			float stratosphereMass = Atmosphere.GetStandardPressureAtElevation(tropopauseElevation, WorldData.StandardTemperature, Gravity) / Gravity;
			StratosphereMass[i] = stratosphereMass;
		}

	}


	[BurstCompile]
	public struct UpdateAirPressureJob : IJobParallelFor {
		public NativeArray<float> Pressure;
		public NativeArray<float> PressureInverse;
		public NativeArray<float> AirMassTotal;
		public NativeArray<float> AbsoluteHumidity;
		public NativeArray<float> RelativeHumidity;
		public NativeArray<float> PotentialEnergy;

		[ReadOnly] public NativeArray<float> AirMass;
		[ReadOnly] public NativeArray<float> AirTemperaturePotential;
		[ReadOnly] public NativeArray<float> VaporMass;
		[ReadOnly] public NativeArray<float> CloudMass;
		[ReadOnly] public NativeArray<float> LayerHeight;
		[ReadOnly] public NativeArray<float> LayerMiddle;
		[ReadOnly] public NativeArray<float> LayerElevation;
		[ReadOnly] public NativeArray<float> SurfaceElevation;
		[ReadOnly] public float Gravity;

		public void Execute(int i)
		{
			float airMass = AirMass[i];
			var pressure = (AirMassTotal[i] + airMass / 2) * Gravity;
			Pressure[i] = pressure;
			PressureInverse[i] = 1.0f / pressure;
			AirMassTotal[i] += airMass;

			float airTemperaturePotential = AirTemperaturePotential[i];
			float vaporMass = VaporMass[i];
			float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(airTemperaturePotential, LayerMiddle[i]);

			AbsoluteHumidity[i] = vaporMass / (vaporMass + airMass);
			RelativeHumidity[i] = Atmosphere.GetRelativeHumidity(airMass, vaporMass, airTemperatureAbsolute, Pressure[i]);
			PotentialEnergy[i] = (airMass * WorldData.SpecificHeatAtmosphere + vaporMass * WorldData.SpecificHeatWaterVapor) * airTemperatureAbsolute;

		}
	}


	[BurstCompile]
	public struct UpdateCloudJob : IJobParallelFor {
		public NativeArray<float> AirDensityCloud;
		public NativeArray<float> CloudElevation;
		public NativeArray<float> DewPoint;
		public NativeArray<float3> CloudVelocity;

		[ReadOnly] public NativeArray<float> RelativeHumidity;
		[ReadOnly] public NativeArray<float> SurfaceElevation;
		[ReadOnly] public NativeArray<float> Pressure;
		[ReadOnly] public NativeArray<float> AirMass;
		[ReadOnly] public NativeArray<float> VaporMass;
		[ReadOnly] public NativeArray<float> AirTemperaturePotential;
		[ReadOnly] public NativeArray<float> LayerMiddle;
		[ReadOnly] public NativeArray<float> LayerElevation;
		[ReadOnly] public NativeArray<float> LayerElevationAbove;
		[ReadOnly] public NativeArray<float> LayerElevationBelow;
		[ReadOnly] public NativeArray<float> LayerHeight;
		[ReadOnly] public NativeArray<float> LayerHeightAbove;
		[ReadOnly] public NativeArray<float> LayerHeightBelow;
		[ReadOnly] public NativeArray<float3> AirVelocity;
		[ReadOnly] public NativeArray<float3> AirVelocityAbove;
		[ReadOnly] public NativeArray<float3> AirVelocityBelow;
		[ReadOnly] public float Gravity;
		[ReadOnly] public bool IsTop;
		[ReadOnly] public bool IsBottom;

		public void Execute(int i)
		{
			float airTemperaturePotential = AirTemperaturePotential[i];
			float layerElevation = LayerElevation[i];
			float layerMiddle = LayerMiddle[i];

			float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(airTemperaturePotential, layerMiddle);
			float dewPoint = Atmosphere.GetDewPoint(RelativeHumidity[i], airTemperatureAbsolute);
			float cloudElevation = Atmosphere.GetElevationAtDewPoint(dewPoint, airTemperaturePotential);

			if ((cloudElevation >= layerElevation || IsBottom) && (cloudElevation < layerElevation + LayerHeight[i] || IsTop))
			{
				DewPoint[i] = dewPoint;

				float cloudBaseElevation = math.max(cloudElevation, SurfaceElevation[i]);
				if (cloudBaseElevation >= layerElevation)
				{
					float airMassCloud = AirMass[i];
					float vaporMassCloud = VaporMass[i];
					float airPressureCloud = Atmosphere.GetPressureAtElevation(cloudBaseElevation, Gravity, Pressure[i], airTemperaturePotential, LayerMiddle[i]);
					AirDensityCloud[i] = Atmosphere.GetAirDensity(airPressureCloud, DewPoint[i], airMassCloud, vaporMassCloud);
				}
				CloudElevation[i] = cloudElevation;

				if (cloudBaseElevation < layerMiddle)
				{
					if (IsBottom)
					{
						CloudVelocity[i] = AirVelocity[i];
					}
					else
					{
						float downLayerMidElevation = (LayerElevationBelow[i] + layerElevation) / 2;
						float t = (cloudBaseElevation - downLayerMidElevation) / (layerMiddle - downLayerMidElevation);
						CloudVelocity[i] = AirVelocity[i] * t + AirVelocityBelow[i] * (1.0f - t);
					}
				}
				else if (IsTop)
				{
					CloudVelocity[i] = AirVelocity[i];
				}
				else
				{
					float upLayerMidElevation = LayerElevationAbove[i] + LayerHeightAbove[i] / 2;
					float t = (cloudBaseElevation - layerMiddle) / (upLayerMidElevation - layerMiddle);
					CloudVelocity[i] = AirVelocityAbove[i] * t + AirVelocity[i] * (1.0f - t);
				}
			}

		}
	}




	[BurstCompile]
	public struct UpdateWaterDepthJob : IJobParallelFor {
		public NativeArray<float> LayerDepth;
		public NativeArray<float> LayerHeight;
		public NativeArray<float> Density;
		public NativeArray<float> WaterMassTotal;
		public NativeArray<float> WaterCoverage;
		public NativeArray<float> PotentialEnergy;
		public NativeArray<float> Pressure;
		[ReadOnly] public NativeArray<float> LayerDepthUp;
		[ReadOnly] public NativeArray<float> LayerHeightUp;
		[ReadOnly] public NativeArray<float> LayerHeightDown;
		[ReadOnly] public NativeArray<float> WaterMass;
		[ReadOnly] public NativeArray<float> SaltMass;
		[ReadOnly] public NativeArray<float> Temperature;
		[ReadOnly] public NativeArray<float> Roughness;
		[ReadOnly] public float WaterDensityPerDegree;
		[ReadOnly] public float WaterDensityPerSalinity;
		[ReadOnly] public float Gravity;
		public void Execute(int i)
		{
			if (WaterMass[i] > 0)
			{
				float waterMass = WaterMass[i];
				float saltMass = SaltMass[i];
				Density[i] = Atmosphere.GetWaterDensity(Atmosphere.GetWaterSalinity(waterMass, saltMass), Temperature[i]);
				float waterDepth = (waterMass + saltMass) / Density[i];
				LayerDepth[i] = LayerDepthUp[i] + waterDepth;
				LayerHeight[i] = waterDepth;
				WaterCoverage[i] = math.min(1, (LayerHeight[i] + LayerHeightDown[i]) / math.max(1, Roughness[i]));
				PotentialEnergy[i] = (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt) * Temperature[i];
				float layerMass = WaterMass[i] + SaltMass[i];
				Pressure[i] = (WaterMassTotal[i] + layerMass / 2) * Gravity;
				WaterMassTotal[i] += layerMass;
			}
			else
			{
				LayerDepth[i] = LayerDepthUp[i];
				LayerHeight[i] = 0;
				Density[i] = 0;
				Pressure[i] = 0;
				WaterCoverage[i] = 0;
				PotentialEnergy[i] = 0;
			}
		}
	}


	[BurstCompile]
	public struct UpdateSurfaceStateJob : IJobParallelFor {
		public NativeArray<float> SurfaceAirTemperatureAbsolute;
		public NativeArray<float> SurfaceAreaAirIce;
		public NativeArray<float> SurfaceAreaAirWater;
		public NativeArray<float> SurfaceAreaAirFlora;
		public NativeArray<float> SurfaceAreaAirTerrain;
		public NativeArray<float> SurfaceAreaIceWater;
		public NativeArray<float> SurfaceAreaIceFlora;
		public NativeArray<float> SurfaceAreaIceTerrain;
		public NativeArray<float> SurfaceAreaWaterFlora;
		public NativeArray<float> SurfaceAreaWaterTerrain;
		public NativeArray<float> SurfaceAreaFloraTerrain;
		public NativeArray<float> FloraCoverage;
		public NativeArray<float> IceCoverage;
		[ReadOnly] public NativeArray<float> IceMass;
		[ReadOnly] public NativeArray<float> FloraMass;
		[ReadOnly] public NativeArray<float> WaterCoverage;
		[ReadOnly] public NativeArray<float> AirTemperaturePotential;
		[ReadOnly] public NativeArray<float> SurfaceLayerElevation;
		[ReadOnly] public NativeArray<float> Roughness;
		[ReadOnly] public float inverseFullCoverageIceMass;
		[ReadOnly] public float inverseFullCoverageFloraMass;
		[ReadOnly] public float FloraAirSurfaceArea;
		public void Execute(int i)
		{
			SurfaceAirTemperatureAbsolute[i] = Atmosphere.GetAbsoluteTemperature(AirTemperaturePotential[i], SurfaceLayerElevation[i]);

			float floraCoverage = math.saturate(FloraMass[i] * inverseFullCoverageFloraMass);
			FloraCoverage[i] = floraCoverage;

			float iceCoverage = math.saturate(IceMass[i] * inverseFullCoverageIceMass);
			IceCoverage[i] = iceCoverage;

			float waterCoverage = WaterCoverage[i];

			SurfaceAreaAirIce[i] = iceCoverage;
			SurfaceAreaAirWater[i] = math.max(0, waterCoverage - iceCoverage);
			SurfaceAreaAirFlora[i] = math.min(floraCoverage, 1.0f - math.max(waterCoverage, iceCoverage)) * FloraAirSurfaceArea;
			SurfaceAreaAirTerrain[i] = math.max(0, 1.0f - (floraCoverage + math.max(waterCoverage, iceCoverage)));
			SurfaceAreaIceWater[i] = math.min(waterCoverage, iceCoverage);
			SurfaceAreaIceFlora[i] = math.max(0, (floraCoverage + math.max(0, iceCoverage - waterCoverage)) - 1.0f) * FloraAirSurfaceArea;
			SurfaceAreaIceTerrain[i] = math.max(0, iceCoverage - waterCoverage) - math.max(0, floraCoverage + iceCoverage - 1.0f);
			SurfaceAreaWaterFlora[i] = math.max(0, (floraCoverage + waterCoverage) - 1.0f) * FloraAirSurfaceArea;
			SurfaceAreaWaterTerrain[i] = waterCoverage - math.max(0, (floraCoverage + waterCoverage) - 1.0f);
			SurfaceAreaFloraTerrain[i] = floraCoverage;
		}

	}



}

