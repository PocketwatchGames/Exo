
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using System.Collections.Generic;
using System;

public struct TempState {

	public struct AirLayerHeights {
		public float MinumumHeight;
		public float ColumnPercent;
	};

	public NativeArray<float> FloraCoverage;
	public NativeArray<float> FloraEnergy;
	public NativeArray<float> LavaEnergy;
	public NativeArray<float> IceCoverage;
	public NativeArray<float> IceEnergy;
	public NativeArray<float> SurfaceAirTemperatureAbsolute;
	public NativeArray<float> WindVerticalCloud;
	public NativeArray<float> SurfaceElevation;
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

	public NativeArray<float> AirMass;
	public NativeArray<float> AirPressure;
	public NativeArray<float> AirPressureInverse;
	public NativeArray<float> AirHumidityAbsolute;
	public NativeArray<float> AirHumidityRelative;
	public NativeArray<float> AirPotentialEnergy;
	public NativeArray<float> AirLayerHeight;
	public NativeArray<float> AirLayerElevation;
	public NativeArray<float> AirLayerMiddle;

	public NativeArray<float> SolarReflectedAir;
	public NativeArray<float> SolarReflectedWater;
	public NativeArray<float> SolarReflectedTerrain;
	public NativeArray<float> SolarReflectedFlora;
	public NativeArray<float> SolarReflectedIce;
	public NativeArray<float> SolarReflectedCloud;
	public NativeArray<float> SolarRadiationInAir;
	public NativeArray<float> SolarRadiationInWater;
	public NativeArray<float> SolarRadiationInTerrain;
	public NativeArray<float> SolarRadiationInFlora;
	public NativeArray<float> SolarRadiationInIce;
	public NativeArray<float> SolarRadiationInCloud;
	public NativeArray<float> EmissivityAir;
	public NativeArray<float> EmissivityWater;
	public NativeArray<float> EmissivityTerrain;
	public NativeArray<float> EmissivityFlora;
	public NativeArray<float> EmissivityIce;
	public NativeArray<float> ThermalRadiationEmittedAir;
	public NativeArray<float> ThermalRadiationEmittedWater;
	public NativeArray<float> ThermalRadiationEmittedIce;
	public NativeArray<float> ThermalRadiationEmittedFlora;
	public NativeArray<float> ThermalRadiationEmittedTerrain;
	public NativeArray<float> ThermalRadiationDeltaAir;
	public NativeArray<float> ThermalRadiationDeltaWater;
	public NativeArray<float> ThermalRadiationDeltaIce;
	public NativeArray<float> ThermalRadiationDeltaFlora;
	public NativeArray<float> ThermalRadiationDeltaTerrain;

	public NativeArray<float>[] WaterCoverage;
	public NativeArray<float>[] WaterDensity;
	public NativeArray<float>[] WaterPressure;
	public NativeArray<float>[] WaterLayerDepth;
	public NativeArray<float>[] WaterLayerHeight;
	public NativeArray<float>[] WaterPotentialEnergy;
	public NativeArray<float> CloudAbsorptivity;
	public NativeArray<float> CloudElevation;
	public NativeArray<float3> CloudVelocity;
	public NativeArray<float> DewPoint;
	public NativeArray<float> AirDensityCloud;
	public NativeArray<float> LavaDepth;


	public NativeArray<float> SolarRadiation;
	public NativeArray<float> AlbedoSlope;
	public NativeArray<float> CloudAlbedo;
	public NativeArray<DiffusionCloud> DiffusionCloud;
	public NativeArray<DiffusionCloud> AdvectionCloud;
	public NativeArray<BarycentricValue> DestinationCloud;

	public NativeArray<SolarAbsorptivity> AbsorptivitySolar;
	public NativeArray<ThermalAbsorptivity> AbsorptivityThermal;
	public NativeArray<DiffusionAir> AdvectionAir;
	public NativeArray<DiffusionAir> DiffusionAir;
	public NativeArray<float> AirMassLeaving;
	public NativeArray<float> DestinationAir;
	public NativeArray<float> DestinationAirResolved;
	public NativeArray<float> DivergencePressureAir;
	public NativeArray<float3> AirAcceleration;
	public NativeArray<float> DivergenceAir;
	public NativeArray<float> CondensationGroundMass;
	public NativeArray<float> CondensationCloudMass;
	public NativeArray<float> LatentHeatAir;
	public NativeArray<float> LatentHeatWater;
	public NativeArray<float> LatentHeatIce;
	public NativeArray<float> LatentHeatFlora;
	public NativeArray<float> LatentHeatTerrain;
	public NativeArray<float> LatentHeatLava;

	public NativeArray<DiffusionWater>[] DiffusionWater;
	public NativeArray<DiffusionWater>[] AdvectionWater;
	public NativeArray<float>[] WaterMassLeaving;
	public NativeArray<float>[] DestinationWater;
	public NativeArray<float>[] DestinationWaterResolved;
	public NativeArray<float>[] DivergencePressureWater;
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
	public NativeArray<float> DustUp;
	public NativeArray<float> DustDown;
	public NativeArray<float> IceMeltedMass;
	public NativeArray<float> LavaCrystalizedMass;
	public NativeArray<float> LavaEjected;
	public NativeArray<float> DustEjected;
	public NativeArray<float> CrustDelta;
	public NativeArray<float> WaterCarbonDelta;
	public NativeArray<float> AirCarbonDelta;
	public NativeArray<float> OxygenDelta;
	public NativeArray<float>[] DivergenceWater;
	public NativeArray<float> DisplaySolarRadiation;

	public NativeArray<float> ConductionWaterTerrainTotal;
	public NativeArray<float> CloudEvaporationMass;
	public NativeArray<float> DropletDelta;
	public NativeArray<float> PrecipitationMass;
	public NativeArray<float> PrecipitationTemperature;
	public NativeArray<float> AtmosphericWindowUp;
	public NativeArray<float> AtmosphericWindowDown;
	public NativeArray<float> AirMassTotal;
	public NativeArray<float> WaterMassTotal;
	public NativeArray<float> StandardLayerElevation;

	public NativeArray<float> OutgoingFlowWater;
	public NativeArray<float> OutgoingFlowLava;
	public NativeArray<float> FlowPercentWater;
	public NativeArray<float> FlowPercentLava;
	public NativeArray<DiffusionLava> DiffusionLava;

	public NativeArray<float> ThermalRadiationTransmittedUp;
	public NativeArray<float> ThermalRadiationTransmittedDown;
	public NativeArray<float> WindowRadiationTransmittedUp;
	public NativeArray<float> WindowRadiationTransmittedDown;


	private NativeArray<AirLayerHeights> _airLayerHeights;

	private JobHelper _jobHelper, _jobHelperAir;


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
		SurfaceElevation = new NativeArray<float>(count, Allocator.Persistent);
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

		AirMass = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirPressure = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirPressureInverse = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirHumidityAbsolute = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirHumidityRelative = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirLayerHeight = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirLayerElevation = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirLayerMiddle = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirPotentialEnergy = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);

		SolarRadiationInAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		SolarRadiationInWater = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		SolarRadiationInIce = new NativeArray<float>(count, Allocator.Persistent);
		SolarRadiationInFlora = new NativeArray<float>(count, Allocator.Persistent);
		SolarRadiationInTerrain = new NativeArray<float>(count, Allocator.Persistent);
		SolarRadiationInCloud = new NativeArray<float>(count, Allocator.Persistent);
		SolarReflectedAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		SolarReflectedWater = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		SolarReflectedIce = new NativeArray<float>(count, Allocator.Persistent);
		SolarReflectedFlora = new NativeArray<float>(count, Allocator.Persistent);
		SolarReflectedTerrain = new NativeArray<float>(count, Allocator.Persistent);
		SolarReflectedCloud = new NativeArray<float>(count, Allocator.Persistent);
		EmissivityAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		EmissivityWater = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		EmissivityIce = new NativeArray<float>(count, Allocator.Persistent);
		EmissivityFlora = new NativeArray<float>(count, Allocator.Persistent);
		EmissivityTerrain = new NativeArray<float>(count, Allocator.Persistent);
		ThermalRadiationEmittedAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		ThermalRadiationEmittedWater = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		ThermalRadiationEmittedIce = new NativeArray<float>(count, Allocator.Persistent);
		ThermalRadiationEmittedFlora = new NativeArray<float>(count, Allocator.Persistent);
		ThermalRadiationEmittedTerrain = new NativeArray<float>(count, Allocator.Persistent);
		ThermalRadiationDeltaAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		ThermalRadiationDeltaWater = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		ThermalRadiationDeltaIce = new NativeArray<float>(count, Allocator.Persistent);
		ThermalRadiationDeltaFlora = new NativeArray<float>(count, Allocator.Persistent);
		ThermalRadiationDeltaTerrain = new NativeArray<float>(count, Allocator.Persistent);

		WaterDensity = new NativeArray<float>[worldData.WaterLayers];
		WaterPressure = new NativeArray<float>[worldData.WaterLayers];
		WaterLayerDepth = new NativeArray<float>[worldData.WaterLayers];
		WaterLayerHeight = new NativeArray<float>[worldData.WaterLayers];
		WaterCoverage = new NativeArray<float>[worldData.WaterLayers];
		WaterPotentialEnergy = new NativeArray<float>[worldData.WaterLayers];
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

		AirMassLeaving = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		DiffusionAir = new NativeArray<DiffusionAir>(count * worldData.AirLayers, Allocator.Persistent);
		AdvectionAir = new NativeArray<DiffusionAir>(count * worldData.AirLayers, Allocator.Persistent);
		DestinationAir = new NativeArray<float>(count * worldData.AirLayers * StaticState.MaxNeighborsVert, Allocator.Persistent);
		DestinationAirResolved = new NativeArray<float>(count * worldData.AirLayers * StaticState.MaxNeighborsVert, Allocator.Persistent);
		DivergencePressureAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirAcceleration = new NativeArray<float3>(count * worldData.AirLayers, Allocator.Persistent);
		AbsorptivitySolar = new NativeArray<SolarAbsorptivity>(count * worldData.AirLayers, Allocator.Persistent);
		AbsorptivityThermal = new NativeArray<ThermalAbsorptivity>(count * worldData.AirLayers, Allocator.Persistent);
		DustUp = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		DustDown = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		DivergenceAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		CondensationGroundMass = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		CondensationCloudMass = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		StandardLayerElevation = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		LatentHeatAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		LatentHeatWater = new NativeArray<float>(count, Allocator.Persistent);
		LatentHeatIce = new NativeArray<float>(count, Allocator.Persistent);
		LatentHeatFlora = new NativeArray<float>(count, Allocator.Persistent);
		LatentHeatLava = new NativeArray<float>(count, Allocator.Persistent);
		LatentHeatTerrain = new NativeArray<float>(count, Allocator.Persistent);

		DiffusionWater = new NativeArray<DiffusionWater>[worldData.WaterLayers];
		AdvectionWater = new NativeArray<DiffusionWater>[worldData.WaterLayers];
		WaterMassLeaving = new NativeArray<float>[worldData.WaterLayers];
		DestinationWater = new NativeArray<float>[worldData.WaterLayers];
		DestinationWaterResolved = new NativeArray<float>[worldData.WaterLayers];
		ConductionWaterTerrain = new NativeArray<float>[worldData.WaterLayers];
		DivergencePressureWater = new NativeArray<float>[worldData.WaterLayers];
		DivergenceWater = new NativeArray<float>[worldData.WaterLayers];
		for (int i = 0; i < worldData.WaterLayers; i++)
		{
			WaterMassLeaving[i] = new NativeArray<float>(count, Allocator.Persistent);
			DiffusionWater[i] = new NativeArray<DiffusionWater>(count, Allocator.Persistent);
			AdvectionWater[i] = new NativeArray<DiffusionWater>(count, Allocator.Persistent);
			DestinationWater[i] = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.Persistent);
			DestinationWaterResolved[i] = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.Persistent);
			ConductionWaterTerrain[i] = new NativeArray<float>(count, Allocator.Persistent);
			DivergencePressureWater[i] = new NativeArray<float>(count, Allocator.Persistent);
			DivergenceWater[i] = new NativeArray<float>(count, Allocator.Persistent);
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
		FloraDeath = new NativeArray<float>(count, Allocator.Persistent);
		PlanktonMassDelta = new NativeArray<float>(count, Allocator.Persistent);
		PlanktonGlucoseDelta = new NativeArray<float>(count, Allocator.Persistent);
		PlanktonDeath = new NativeArray<float>(count, Allocator.Persistent);
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
		IceMeltedMass = new NativeArray<float>(count, Allocator.Persistent);
		LavaCrystalizedMass = new NativeArray<float>(count, Allocator.Persistent);
		LavaEjected = new NativeArray<float>(count, Allocator.Persistent);
		DustEjected = new NativeArray<float>(count, Allocator.Persistent);
		CrustDelta = new NativeArray<float>(count, Allocator.Persistent);
		AirCarbonDelta = new NativeArray<float>(count, Allocator.Persistent);
		OxygenDelta = new NativeArray<float>(count, Allocator.Persistent);
		SoilRespiration = new NativeArray<float>(count, Allocator.Persistent);
		WaterCarbonDelta = new NativeArray<float>(count, Allocator.Persistent);


		ConductionWaterTerrainTotal = new NativeArray<float>(count, Allocator.Persistent);
		CloudEvaporationMass = new NativeArray<float>(count, Allocator.Persistent);
		DropletDelta = new NativeArray<float>(count, Allocator.Persistent);
		PrecipitationMass = new NativeArray<float>(count, Allocator.Persistent);
		PrecipitationTemperature = new NativeArray<float>(count, Allocator.Persistent);
		AtmosphericWindowUp = new NativeArray<float>(count, Allocator.Persistent);
		AtmosphericWindowDown = new NativeArray<float>(count, Allocator.Persistent);
		AirMassTotal = new NativeArray<float>(count, Allocator.Persistent);
		WaterMassTotal = new NativeArray<float>(count, Allocator.Persistent);


		OutgoingFlowWater = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.Persistent);
		OutgoingFlowLava = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.Persistent);
		FlowPercentWater = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.Persistent);
		FlowPercentLava = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.Persistent);
		DiffusionLava = new NativeArray<DiffusionLava>(count * StaticState.MaxNeighbors, Allocator.Persistent);

		ThermalRadiationTransmittedUp = new NativeArray<float>(count, Allocator.Persistent);
		ThermalRadiationTransmittedDown = new NativeArray<float>(count, Allocator.Persistent);
		WindowRadiationTransmittedUp = new NativeArray<float>(count, Allocator.Persistent);
		WindowRadiationTransmittedDown = new NativeArray<float>(count, Allocator.Persistent);



		_airLayerHeights = new NativeArray<AirLayerHeights>(worldData.AirLayers, Allocator.Persistent);
		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			if (i == 1)
			{
				_airLayerHeights[i] = new AirLayerHeights
				{
					MinumumHeight = worldData.BoundaryZoneElevation,
					ColumnPercent = 0,
				};
			}
			else
			{
				_airLayerHeights[i] = new AirLayerHeights
				{
					MinumumHeight = 0,
					ColumnPercent = 1.0f / (worldData.AirLayers - 1 - i),
				};
			}
		}
		_jobHelper = new JobHelper(count);
		_jobHelperAir = new JobHelper(count * (worldData.AirLayers - 2));
	}

	public void Dispose(ref WorldData worldData)
	{

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
		SurfaceElevation.Dispose();
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

		AirMass.Dispose();
		AirPressure.Dispose();
		AirPressureInverse.Dispose();
		AirHumidityRelative.Dispose();
		AirHumidityAbsolute.Dispose();
		AirLayerHeight.Dispose();
		AirLayerElevation.Dispose();
		AirLayerMiddle.Dispose();
		AirPotentialEnergy.Dispose();
		AirMassLeaving.Dispose();
		DiffusionAir.Dispose();
		AdvectionAir.Dispose();
		DestinationAir.Dispose();
		DestinationAirResolved.Dispose();
		AirAcceleration.Dispose();
		AbsorptivitySolar.Dispose();
		AbsorptivityThermal.Dispose();
		DustUp.Dispose();
		DustDown.Dispose();
		DivergencePressureAir.Dispose();
		DivergenceAir.Dispose();
		CondensationGroundMass.Dispose();
		CondensationCloudMass.Dispose();
		LatentHeatAir.Dispose();
		LatentHeatWater.Dispose();
		LatentHeatIce.Dispose();
		LatentHeatTerrain.Dispose();
		LatentHeatFlora.Dispose();
		LatentHeatLava.Dispose();

		for (int i = 0; i < WaterCoverage.Length; i++)
		{
			WaterCoverage[i].Dispose();
			WaterDensity[i].Dispose();
			WaterPressure[i].Dispose();
			WaterLayerDepth[i].Dispose();
			WaterLayerHeight[i].Dispose();
			WaterPotentialEnergy[i].Dispose();
			WaterMassLeaving[i].Dispose();
			DiffusionWater[i].Dispose();
			AdvectionWater[i].Dispose();
			DestinationWater[i].Dispose();
			DestinationWaterResolved[i].Dispose();
			ConductionWaterTerrain[i].Dispose();
			DivergencePressureWater[i].Dispose();
			DivergenceWater[i].Dispose();
		}

		SolarRadiationInAir.Dispose();
		SolarRadiationInWater.Dispose();
		SolarRadiationInIce.Dispose();
		SolarRadiationInFlora.Dispose();
		SolarRadiationInTerrain.Dispose();
		SolarRadiationInCloud.Dispose();
		SolarReflectedAir.Dispose();
		SolarReflectedWater.Dispose();
		SolarReflectedIce.Dispose();
		SolarReflectedFlora.Dispose();
		SolarReflectedTerrain.Dispose();
		SolarReflectedCloud.Dispose();
		EmissivityAir.Dispose();
		EmissivityWater.Dispose();
		EmissivityIce.Dispose();
		EmissivityFlora.Dispose();
		EmissivityTerrain.Dispose();
		ThermalRadiationEmittedAir.Dispose();
		ThermalRadiationEmittedWater.Dispose();
		ThermalRadiationEmittedIce.Dispose();
		ThermalRadiationEmittedFlora.Dispose();
		ThermalRadiationEmittedTerrain.Dispose();
		ThermalRadiationDeltaAir.Dispose();
		ThermalRadiationDeltaWater.Dispose();
		ThermalRadiationDeltaIce.Dispose();
		ThermalRadiationDeltaFlora.Dispose();
		ThermalRadiationDeltaTerrain.Dispose();



		ThermalRadiationTransmittedUp.Dispose();
		ThermalRadiationTransmittedDown.Dispose();
		WindowRadiationTransmittedUp.Dispose();
		WindowRadiationTransmittedDown.Dispose();


		SolarRadiation.Dispose();
		AlbedoSlope.Dispose();
		CloudAlbedo.Dispose();

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
		AirMassTotal.Dispose();
		WaterMassTotal.Dispose();
		StandardLayerElevation.Dispose();
		DropletDelta.Dispose();
		PrecipitationMass.Dispose();
		PrecipitationTemperature.Dispose();
		CloudEvaporationMass.Dispose();
		ConductionWaterTerrainTotal.Dispose();

		OutgoingFlowWater.Dispose();
		OutgoingFlowLava.Dispose();
		FlowPercentWater.Dispose();
		FlowPercentLava.Dispose();
		DiffusionLava.Dispose();

		_airLayerHeights.Dispose();
	}

	public JobHandle Clear(int cellCount, ref WorldData worldData, JobHandle dependency)
	{
		var h = default(JobHandle);
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, ConductionWaterTerrainTotal, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, CloudEvaporationMass, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, DropletDelta, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, PrecipitationMass, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, PrecipitationTemperature, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, AtmosphericWindowUp, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, AtmosphericWindowDown, 0));

		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount * worldData.AirLayers, dependency, DivergencePressureAir, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount * worldData.AirLayers, dependency, CondensationGroundMass, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount * worldData.AirLayers, dependency, CondensationCloudMass, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount * worldData.AirLayers, dependency, LatentHeatAir, 0));

		for (int i = 1; i < worldData.WaterLayers - 1; i++)
		{
			h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, DivergencePressureWater[i], 0));
		}
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, LatentHeatWater, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, LatentHeatIce, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, ThermalRadiationTransmittedUp, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, ThermalRadiationTransmittedDown, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, WindowRadiationTransmittedUp, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, WindowRadiationTransmittedDown, 0));

		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount * worldData.AirLayers, dependency, ThermalRadiationDeltaAir, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount * worldData.WaterLayers, dependency, ThermalRadiationDeltaWater, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, ThermalRadiationDeltaIce, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, ThermalRadiationDeltaFlora, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, ThermalRadiationDeltaTerrain, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount * worldData.AirLayers, dependency, SolarReflectedAir, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount * worldData.WaterLayers, dependency, SolarReflectedWater, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, SolarReflectedIce, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, SolarReflectedFlora, 0));
		h = JobHandle.CombineDependencies(h, Utils.MemsetArray(cellCount, dependency, SolarReflectedTerrain, 0));

		return h;
	}

	public JobHandle Update(ref SimState state, ref WorldData worldData, ref StaticState staticState, JobHandle dependencies)
	{
		dependencies = Utils.MemCopy(WaterMassTotal, state.IceMass, dependencies);
		dependencies = Utils.MemsetArray(staticState.Count, dependencies, AirMassTotal, 0);
		for (int j = worldData.WaterLayers - 2; j >= 1; j--)
		{
			dependencies = _jobHelper.Schedule(new UpdateWaterDepthJob()
			{
				Density = WaterDensity[j],
				LayerDepth = WaterLayerDepth[j],
				LayerHeight = WaterLayerHeight[j],
				WaterMassTotal = WaterMassTotal,
				WaterCoverage = WaterCoverage[j],
				PotentialEnergy = WaterPotentialEnergy[j],
				Pressure = WaterPressure[j],

				Temperature = state.WaterTemperature[j],
				SaltMass = state.SaltMass[j],
				WaterMass = state.WaterMass[j],
				Roughness = state.Roughness,
				LayerDepthUp = WaterLayerDepth[j + 1],
				LayerHeightUp = WaterLayerHeight[j + 1],
				LayerHeightDown = WaterLayerHeight[j - 1],
				WaterDensityPerDegree = worldData.WaterDensityPerDegree,
				WaterDensityPerSalinity = worldData.WaterDensityPerSalinity,
				Gravity = state.PlanetState.Gravity,

			}, dependencies);
		}



		dependencies = _jobHelper.Schedule(new UpdateTempStateJob()
		{
			IceEnergy = IceEnergy,
			FloraEnergy = FloraEnergy,
			LavaEnergy = LavaEnergy,
			LavaDepth = LavaDepth,
			SurfaceElevation = SurfaceElevation,

			WaterDepth = WaterLayerDepth[1],
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
		dependencies = Utils.MemCopy(staticState.GetSliceLayer(AirLayerElevation, worldData.SurfaceAirLayer), SurfaceElevation, dependencies);
		dependencies = Utils.MemCopy(staticState.GetSliceLayer(StandardLayerElevation, worldData.SurfaceAirLayer), SurfaceElevation, dependencies);

		dependencies = _jobHelperAir.Run(new UpdateAirLayerHeightsJob()
		{
			StandardLayerElevation = staticState.GetSliceAir(StandardLayerElevation),
			LayerHeight = staticState.GetSliceAir(AirLayerHeight),
			LayerElevation = staticState.GetSliceAir(AirLayerElevation),
			AirMass = staticState.GetSliceAir(AirMass),
			LayerMiddle = staticState.GetSliceAir(AirLayerMiddle),

			AirTemperaturePotential = staticState.GetSliceAir(state.AirTemperaturePotential),
			AirLayerHeights = new NativeSlice<AirLayerHeights>(_airLayerHeights, 1, worldData.AirLayers-2),
			SurfaceElevation = SurfaceElevation,
			TropopauseElevation = worldData.TropopauseElevation,
			Gravity = state.PlanetState.Gravity,
			Count = staticState.Count
		}, dependencies);

		dependencies = _jobHelper.Schedule(new UpdateStratosphereJob()
		{
			StratosphereMass = AirMassTotal,

			TropopauseElevation = staticState.GetSliceLayer(AirLayerElevation,worldData.AirLayers - 2),
			TropopauseHeight = staticState.GetSliceLayer(AirLayerHeight,worldData.AirLayers - 2),
			Gravity = state.PlanetState.Gravity
		}, dependencies);


		for (int j = worldData.AirLayers - 2; j > 0; j--)
		{
			dependencies = _jobHelper.Schedule(new UpdateAirPressureJob()
			{
				Pressure = staticState.GetSliceLayer(AirPressure, j),
				PressureInverse = staticState.GetSliceLayer(AirPressureInverse, j),
				AirMassTotal = AirMassTotal,
				RelativeHumidity = staticState.GetSliceLayer(AirHumidityRelative, j),
				AbsoluteHumidity = staticState.GetSliceLayer(AirHumidityAbsolute, j),
				AirMass = staticState.GetSliceLayer(AirMass, j),
				PotentialEnergy = staticState.GetSliceLayer(AirPotentialEnergy, j),

				CloudMass = state.CloudMass,
				VaporMass = staticState.GetSliceLayer(state.AirVapor, j),
				AirTemperaturePotential = staticState.GetSliceLayer(state.AirTemperaturePotential, j),
				LayerElevation = staticState.GetSliceLayer(AirLayerElevation, j),
				LayerMiddle = staticState.GetSliceLayer(AirLayerMiddle, j),
				LayerHeight = staticState.GetSliceLayer(AirLayerHeight, j),
				SurfaceElevation = SurfaceElevation,
				Gravity = state.PlanetState.Gravity,
			}, dependencies);
		}

		dependencies = _jobHelper.Schedule(new UpdateCloudJob()
		{
			CloudVelocity = CloudVelocity,
			CloudElevation = CloudElevation,
			AirDensityCloud = AirDensityCloud,
			DewPoint = DewPoint,

			AirMass = AirMass,
			VaporMass = state.AirVapor,
			RelativeHumidity = AirHumidityRelative,
			Pressure = AirPressure,
			AirTemperaturePotential = state.AirTemperaturePotential,
			LayerMiddle = AirLayerMiddle,
			LayerElevation = AirLayerElevation,
			LayerHeight = AirLayerHeight,
			AirVelocity = state.AirVelocity,
			Gravity = state.PlanetState.Gravity,
			AirLayerCount = worldData.AirLayers,
			Count = staticState.Count

		}, dependencies);

		var surfaceStateJobHandle = _jobHelper.Schedule(new UpdateSurfaceStateJob()
		{
			SurfaceAirTemperatureAbsolute = SurfaceAirTemperatureAbsolute,
			SurfaceAreaAirFlora = SurfaceAreaAirFlora,
			SurfaceAreaAirIce = SurfaceAreaAirIce,
			SurfaceAreaAirTerrain = SurfaceAreaAirTerrain,
			SurfaceAreaAirWater = SurfaceAreaAirWater,
			SurfaceAreaIceFlora = SurfaceAreaIceFlora,
			SurfaceAreaIceTerrain = SurfaceAreaIceTerrain,
			SurfaceAreaIceWater = SurfaceAreaIceWater,
			SurfaceAreaWaterFlora = SurfaceAreaWaterFlora,
			SurfaceAreaWaterTerrain = SurfaceAreaWaterTerrain,
			SurfaceAreaFloraTerrain = SurfaceAreaFloraTerrain,
			FloraCoverage = FloraCoverage,
			IceCoverage = IceCoverage,

			WaterCoverage = WaterCoverage[worldData.SurfaceWaterLayer],
			FloraMass = state.FloraMass,
			IceMass = state.IceMass,
			AirTemperaturePotential = staticState.GetSliceLayer(state.AirTemperaturePotential, worldData.SurfaceAirLayer),
			SurfaceLayerElevation = staticState.GetSliceLayer(AirLayerElevation,worldData.SurfaceAirLayer),
			inverseFullCoverageFloraMass = 1.0f / worldData.FullCoverageFlora,
			inverseFullCoverageIceMass = 1.0f / (worldData.FullCoverageIce * WorldData.MassIce),
			FloraAirSurfaceArea = worldData.FloraAirSurfaceArea,
			Roughness = state.Roughness
		}, dependencies);
		dependencies = JobHandle.CombineDependencies(dependencies, surfaceStateJobHandle);


		return dependencies;
	}


	[BurstCompile]
	public struct UpdateTempStateJob : IJobParallelFor {
		public NativeSlice<float> SurfaceElevation;
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


	//[BurstCompile]
	public struct UpdateAirLayerHeightsJob : IJobParallelFor {
		public NativeSlice<float> AirMass;
		public NativeSlice<float> LayerElevation;
		public NativeSlice<float> LayerHeight;
		public NativeSlice<float> LayerMiddle;
		public NativeSlice<float> StandardLayerElevation;

		[ReadOnly] public NativeSlice<float> AirTemperaturePotential;
		[ReadOnly] public NativeSlice<AirLayerHeights> AirLayerHeights;
		[ReadOnly] public NativeArray<float> SurfaceElevation;
		[ReadOnly] public float Gravity;
		[ReadOnly] public float TropopauseElevation;
		[ReadOnly] public int Count;

		public void Execute(int i)
		{
			int layer = i / Count;
			int index = i % Count;
			float layerElevation = SurfaceElevation[index];
			float standardLayerElevation = layerElevation;
			for (int l = 0; l <= layer; l++)
			{
				float standardLayerHeight = math.max(AirLayerHeights[l].MinumumHeight, (TropopauseElevation - standardLayerElevation) * AirLayerHeights[l].ColumnPercent);
				float airTemperaturePotential = AirTemperaturePotential[l * Count + index];
				float layerHeight = standardLayerHeight * airTemperaturePotential / WorldData.StandardTemperature;

				if (l == layer)
				{
					float airMass = (Atmosphere.GetStandardPressureAtElevation(standardLayerElevation, WorldData.StandardTemperature, Gravity) - Atmosphere.GetStandardPressureAtElevation(standardLayerElevation + standardLayerHeight, WorldData.StandardTemperature, Gravity)) / Gravity;
					StandardLayerElevation[i] = standardLayerElevation + standardLayerHeight;
					AirMass[i] = airMass;
					LayerHeight[i] = layerHeight;
					LayerElevation[i] = layerElevation;
					LayerMiddle[i] = layerElevation + layerHeight / 2;
					return;
				}

				layerElevation += layerHeight;
				standardLayerElevation += standardLayerHeight;
			}
		}
	}

	[BurstCompile]
	public struct UpdateStratosphereJob : IJobParallelFor {
		public NativeArray<float> StratosphereMass;
		[ReadOnly] public NativeSlice<float> TropopauseElevation;
		[ReadOnly] public NativeSlice<float> TropopauseHeight;
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
		public NativeSlice<float> Pressure;
		public NativeSlice<float> PressureInverse;
		public NativeSlice<float> AirMassTotal;
		public NativeSlice<float> AbsoluteHumidity;
		public NativeSlice<float> RelativeHumidity;
		public NativeSlice<float> PotentialEnergy;

		[ReadOnly] public NativeSlice<float> AirMass;
		[ReadOnly] public NativeSlice<float> AirTemperaturePotential;
		[ReadOnly] public NativeSlice<float> VaporMass;
		[ReadOnly] public NativeSlice<float> CloudMass;
		[ReadOnly] public NativeSlice<float> LayerHeight;
		[ReadOnly] public NativeSlice<float> LayerMiddle;
		[ReadOnly] public NativeSlice<float> LayerElevation;
		[ReadOnly] public NativeSlice<float> SurfaceElevation;
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
		[ReadOnly] public NativeArray<float> Pressure;
		[ReadOnly] public NativeArray<float> AirMass;
		[ReadOnly] public NativeArray<float> VaporMass;
		[ReadOnly] public NativeArray<float> AirTemperaturePotential;
		[ReadOnly] public NativeArray<float> LayerMiddle;
		[ReadOnly] public NativeArray<float> LayerElevation;
		[ReadOnly] public NativeArray<float> LayerHeight;
		[ReadOnly] public NativeArray<float3> AirVelocity;
		[ReadOnly] public float Gravity;
		[ReadOnly] public int AirLayerCount;
		[ReadOnly] public int Count;

		public void Execute(int i)
		{
			for (int layer = 1; layer < AirLayerCount - 1; layer++)
			{
				int airIndex = layer * Count + i;
				float airTemperaturePotential = AirTemperaturePotential[airIndex];
				float layerElevation = LayerElevation[airIndex];
				float layerMiddle = LayerMiddle[airIndex];

				float airTemperatureAbsolute = Atmosphere.GetAbsoluteTemperature(airTemperaturePotential, layerMiddle);
				float dewPoint = Atmosphere.GetDewPoint(RelativeHumidity[airIndex], airTemperatureAbsolute);
				float cloudElevation = Atmosphere.GetElevationAtDewPoint(dewPoint, airTemperaturePotential);

				bool isTop = layer == AirLayerCount - 2;
				bool isBottom = layer == 1;

				if ((cloudElevation >= layerElevation || isBottom) && (cloudElevation < layerElevation + LayerHeight[airIndex] || isTop))
				{
					DewPoint[i] = dewPoint;

					float cloudBaseElevation = math.max(cloudElevation, LayerElevation[Count + i]);
					if (cloudBaseElevation >= layerElevation)
					{
						float airMassCloud = AirMass[airIndex];
						float vaporMassCloud = VaporMass[airIndex];
						float airPressureCloud = Atmosphere.GetPressureAtElevation(cloudBaseElevation, Gravity, Pressure[airIndex], airTemperaturePotential, LayerMiddle[airIndex]);
						AirDensityCloud[i] = Atmosphere.GetAirDensity(airPressureCloud, DewPoint[i], airMassCloud, vaporMassCloud);
					}
					CloudElevation[i] = cloudElevation;

					if (cloudBaseElevation < layerMiddle)
					{
						if (isBottom)
						{
							CloudVelocity[i] = AirVelocity[airIndex];
						}
						else
						{
							int airIndexDown = (layer - 1) * Count + i;
							float downLayerMidElevation = (LayerElevation[airIndexDown] + layerElevation) / 2;
							float t = (cloudBaseElevation - downLayerMidElevation) / (layerMiddle - downLayerMidElevation);
							CloudVelocity[i] = AirVelocity[airIndex] * t + AirVelocity[airIndexDown] * (1.0f - t);
						}
					}
					else if (isTop)
					{
						CloudVelocity[i] = AirVelocity[airIndex];
					}
					else
					{
						int airIndexUp = (layer + 1) * Count + i;
						float upLayerMidElevation = LayerElevation[airIndexUp] + LayerHeight[airIndexUp] / 2;
						float t = (cloudBaseElevation - layerMiddle) / (upLayerMidElevation - layerMiddle);
						CloudVelocity[i] = AirVelocity[airIndexUp] * t + AirVelocity[airIndex] * (1.0f - t);
					}
					return;
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
		[ReadOnly] public NativeSlice<float> AirTemperaturePotential;
		[ReadOnly] public NativeSlice<float> SurfaceLayerElevation;
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

