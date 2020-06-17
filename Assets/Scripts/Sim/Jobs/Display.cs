#define LayerRefactor

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

public struct DisplayState {
	public static JobHelper DisplayJob;
	public static JobHelper DisplayJobAir;
	public static JobHelper DisplayJobWater;

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
	public float GlobalAirCarbon;
	public float GlobalWaterCarbon;
	public float GlobalCloudCoverage;
	public float GlobalEvaporation;
	public float GlobalRainfall;
	public float GlobalSoilFertility;
	public float GlobalCondensationCloud;
	public float GlobalCondensationGround;
	public double GlobalWaterVapor;
	public double GlobalAirMass;
	public float GlobalCloudMass;

	public NativeArray<float> Plate;
	public NativeArray<float> SolarRadiationAbsorbedSurface;
	public NativeArray<float> Rainfall;
	public NativeArray<float> CondensationCloud;
	public NativeArray<float> CondensationGround;
	public NativeArray<float> Evaporation;
	public NativeArray<float> EnthalpyTerrain;
	public NativeArray<float> EnthalpyIce;
	public NativeArray<float> EnthalpyCloud;
	public NativeArray<float> EnthalpyGroundWater;
	public NativeArray<float> DustMass;
	public NativeArray<float> CarbonDioxidePercent;
	public NativeArray<float> EnthalpyAir;
	public NativeArray<float> Pressure;
	public NativeArray<float> DivergenceAir;
	public NativeArray<float3> PressureGradientForce;
	public NativeArray<float> WindVertical;
	public NativeArray<SolarAbsorptivity> AbsorptionSolar;
	public NativeArray<ThermalAbsorptivity> AbsorptionThermal;
	public NativeArray<float> LatentHeatDelta;
	public NativeArray<float> EnthalpyWater;
	public NativeArray<float> Salinity;
	public NativeArray<float> WaterCarbonDioxidePercent;
	public NativeArray<float> DivergenceWater;
	public NativeArray<float>[] ConductionDelta;

	private bool _initialized;

	public void Init(int count, ref WorldData worldData)
	{
		_initialized = true;
		Plate = new NativeArray<float>(count, Allocator.Persistent);
		SolarRadiationAbsorbedSurface = new NativeArray<float>(count, Allocator.Persistent);
		Rainfall = new NativeArray<float>(count, Allocator.Persistent);
		Evaporation = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyTerrain = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyIce = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyCloud = new NativeArray<float>(count, Allocator.Persistent);
		EnthalpyGroundWater = new NativeArray<float>(count, Allocator.Persistent);

		DivergenceAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		Pressure = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		PressureGradientForce = new NativeArray<float3>(count * worldData.AirLayers, Allocator.Persistent);
		EnthalpyAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		WindVertical = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AbsorptionSolar = new NativeArray<SolarAbsorptivity>(count * worldData.AirLayers, Allocator.Persistent);
		AbsorptionThermal = new NativeArray<ThermalAbsorptivity>(count * worldData.AirLayers, Allocator.Persistent);
		CarbonDioxidePercent = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		CondensationCloud = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		CondensationGround = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		DustMass = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);

		WaterCarbonDioxidePercent = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		Salinity = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		EnthalpyWater = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		DivergenceWater = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);

		ConductionDelta = new NativeArray<float>[worldData.LayerCount];
		for (int i = 0; i < worldData.LayerCount; i++)
		{
			ConductionDelta[i] = new NativeArray<float>(count, Allocator.Persistent);
		}

		LatentHeatDelta = new NativeArray<float>(count* worldData.LayerCount, Allocator.Persistent);

	}

	public void Dispose()
	{
		if (!_initialized)
		{
			return;
		}
		Plate.Dispose();
		SolarRadiationAbsorbedSurface.Dispose();
		Rainfall.Dispose();
		CondensationCloud.Dispose();
		CondensationGround.Dispose();
		Evaporation.Dispose();
		EnthalpyTerrain.Dispose();
		EnthalpyCloud.Dispose();
		EnthalpyIce.Dispose();
		EnthalpyGroundWater.Dispose();
		DustMass.Dispose();

		DivergenceAir.Dispose();
		Pressure.Dispose();
		PressureGradientForce.Dispose();
		EnthalpyAir.Dispose();
		WindVertical.Dispose();
		AbsorptionSolar.Dispose();
		AbsorptionThermal.Dispose();
		CarbonDioxidePercent.Dispose();

		LatentHeatDelta.Dispose();

		WaterCarbonDioxidePercent.Dispose();
		Salinity.Dispose();
		EnthalpyWater.Dispose();
		DivergenceWater.Dispose();

		for (int i = 0; i < ConductionDelta.Length; i++)
		{
			ConductionDelta[i].Dispose();
		}
	}

	public static JobHandle Update(
		ref DisplayState display,
		ref DisplayState lastDisplay,
		ref WorldData worldData, 
		ref TempState tempState, 
		ref SimState nextState, 
		ref StaticState staticState,
		ref SimSettings settings
		)
	{
		if (DisplayJob == null)
		{
			DisplayJob = new JobHelper(staticState.Count);
			DisplayJobAir = new JobHelper(staticState.Count * (worldData.AirLayers - 2));
			DisplayJobWater = new JobHelper(staticState.Count * (worldData.WaterLayers - 2));
		}

		JobHandle initDisplayAirHandle = default(JobHandle);
		JobHandle initDisplayWaterHandle = default(JobHandle);

		initDisplayAirHandle = JobHandle.CombineDependencies(initDisplayAirHandle, DisplayJobAir.Schedule(
			JobType.Schedule, 64,
			new GetDivergenceJob()
			{
				Divergence = staticState.GetSliceAir(display.DivergenceAir),
				Destination = staticState.GetSliceAirNeighbors(tempState.DestinationAirResolved),
			}));

		initDisplayWaterHandle = JobHandle.CombineDependencies(initDisplayWaterHandle, DisplayJobWater.Schedule(
			JobType.Schedule, 64,
			new GetDivergenceJob()
			{
				Divergence = staticState.GetSliceWater(display.DivergenceWater),
				Destination = staticState.GetSliceAirNeighbors(tempState.DestinationWater),
			}));


		tempState.AbsorptivitySolar.CopyTo(display.AbsorptionSolar);
		tempState.AbsorptivityThermal.CopyTo(display.AbsorptionThermal);
		initDisplayAirHandle = JobHandle.CombineDependencies(initDisplayAirHandle, DisplayJobAir.Schedule(
			JobType.Schedule, 64,
			new InitDisplayAirLayerJob()
			{
				DisplayPressure = staticState.GetSliceAir(display.Pressure),
				DisplayPressureGradientForce = staticState.GetSliceAir(display.PressureGradientForce),
				DisplayCondensationGround = staticState.GetSliceAir(display.CondensationGround),
				DisplayCondensationCloud = staticState.GetSliceAir(display.CondensationCloud),
				Enthalpy = staticState.GetSliceAir(display.EnthalpyAir),
				DustCoverage = staticState.GetSliceAir(display.DustMass),
				CarbonDioxidePercent = staticState.GetSliceAir(display.CarbonDioxidePercent),
				Divergence = staticState.GetSliceAir(display.DivergenceAir),

				CarbonDioxide = staticState.GetSliceAir(nextState.AirCarbonDioxide),
				AirTemperaturePotential = staticState.GetSliceAir(nextState.AirTemperaturePotential),
				AirPressure = staticState.GetSliceAir(tempState.AirPressure),
				LayerMiddle = staticState.GetSliceAir(tempState.AirLayerMiddle),
				PressureGradientForce = staticState.GetSliceAir(tempState.AirAcceleration),
				CondensationCloud = staticState.GetSliceAir(tempState.CondensationCloudMass),
				CondensationGround = staticState.GetSliceAir(tempState.CondensationGroundMass),
				AirMass = staticState.GetSliceAir(tempState.AirMass),
				VaporMass = staticState.GetSliceAir(nextState.AirVapor),
				DustMass = staticState.GetSliceAir(nextState.AirDust),
				CloudElevation = tempState.CloudElevation,
				CloudMass = nextState.CloudMass,
				LayerElevation = staticState.GetSliceAir(tempState.AirLayerElevation),
				LayerHeight = staticState.GetSliceAir(tempState.AirLayerHeight),
				Gravity = nextState.PlanetState.Gravity,
				Count = staticState.Count,
			}, initDisplayAirHandle));

		initDisplayWaterHandle = JobHandle.CombineDependencies(initDisplayWaterHandle, DisplayJobWater.Schedule(
			JobType.Schedule, 64,
			new InitDisplayWaterLayerJob()
			{
				Enthalpy = staticState.GetSliceWater(display.EnthalpyWater),
				Salinity = staticState.GetSliceWater(display.Salinity),
				CarbonPercent = staticState.GetSliceWater(display.WaterCarbonDioxidePercent),

				WaterTemperature = nextState.WaterTemperature,
				SaltMass = nextState.WaterSaltMass,
				WaterMass = nextState.WaterMass,
				WaterCarbon = nextState.WaterCarbonDioxide,
				CountPerLayer = staticState.Count,
			}, initDisplayWaterHandle));

		var updateDisplayJobHandle = DisplayJob.Schedule(
			JobType.Schedule, 64,
			new UpdateDisplayJob()
			{
				DisplayPlate = display.Plate,
				SolarRadiationAbsorbedSurface = display.SolarRadiationAbsorbedSurface,
				DisplayEvaporation = display.Evaporation,
				DisplayPrecipitation = display.Rainfall,
				EnthalpyTerrain = display.EnthalpyTerrain,
				EnthalpyCloud = display.EnthalpyCloud,
				EnthalpyIce = display.EnthalpyIce,
				EnthalpyGroundWater = display.EnthalpyGroundWater,

				Plate = nextState.Plate,
				SolarRadiationInTerrain = tempState.SolarRadiationInTerrain,
				SolarRadiationInIce = tempState.SolarRadiationInIce,
				SolarRadiationInWaterSurface = tempState.SolarRadiationInWater,
				EvaporationWater = tempState.EvaporationMassWater,
				EvaporationFlora = tempState.FloraRespirationMassVapor,
				Precipitation = tempState.PrecipitationMass,
				SoilFertility = tempState.SoilFertility,
				TerrainTemperature = nextState.GroundTemperature,
				HeatingDepth = worldData.SoilHeatDepth,
				CloudMass = nextState.CloudMass,
				IceMass = nextState.IceMass,
				IceTemperature = nextState.IceTemperature,
				GroundWaterMass = nextState.GroundWater,
				GroundWaterTemperature = nextState.GroundWaterTemperature
			});

		updateDisplayJobHandle = JobHandle.CombineDependencies(initDisplayAirHandle, initDisplayWaterHandle, updateDisplayJobHandle);
		UpdateGlobals(ref display, ref worldData, ref tempState, ref nextState, ref staticState, default(JobHandle));
		if (settings.CollectGlobalsDebug)
		{
			updateDisplayJobHandle.Complete();

			UpdateGlobalsDebug(ref display, ref worldData, ref tempState, ref nextState, ref staticState);

			display.GlobalEnthalpyDelta = display.GlobalEnthalpy - lastDisplay.GlobalEnthalpy;
			display.GlobalEnthalpyDeltaTerrain = display.GlobalEnthalpyTerrain - lastDisplay.GlobalEnthalpyTerrain;
			display.GlobalEnthalpyDeltaAir = display.GlobalEnthalpyAir - lastDisplay.GlobalEnthalpyAir;
			display.GlobalEnthalpyDeltaWater = display.GlobalEnthalpyWater - lastDisplay.GlobalEnthalpyWater;
			display.GlobalEnthalpyDeltaCloud = display.GlobalEnthalpyCloud - lastDisplay.GlobalEnthalpyCloud;
			display.GlobalEnthalpyDeltaIce = display.GlobalEnthalpyIce - lastDisplay.GlobalEnthalpyIce;
			display.GlobalEnthalpyDeltaGroundWater = display.GlobalEnthalpyGroundWater - lastDisplay.GlobalEnthalpyGroundWater;
		}


		return updateDisplayJobHandle;
	}

	public static void UpdateGlobals(
		ref DisplayState display,
		ref WorldData worldData,
		ref TempState tempState,
		ref SimState nextState,
		ref StaticState staticState,
		JobHandle dependency
		)
	{
		dependency.Complete();
		display.GlobalSurfaceTemperature = 0;
		for (int i = 0; i < staticState.Count; i++)
		{
			display.GlobalSurfaceTemperature += tempState.SurfaceAirTemperatureAbsolute[i];
		}
		display.GlobalSurfaceTemperature /= staticState.Count;
	}
	public static void UpdateGlobalsDebug(
		ref DisplayState display,
		ref WorldData worldData,
		ref TempState tempState,
		ref SimState nextState,
		ref StaticState staticState
		)
	{
		float globalWaterMass = 0;
		float globalWaterSurfaceMass = 0;
		for (int i = 0; i < staticState.Count; i++)
		{
			int surfaceWaterIndex = staticState.GetWaterIndex(worldData.SurfaceWaterLayer, i);
			float waterMassSurface = nextState.WaterMass[surfaceWaterIndex];
			globalWaterSurfaceMass += waterMassSurface;
			display.GlobalSoilFertility += nextState.GroundCarbonDioxide[i];
			display.GlobalOceanSurfaceTemperature += nextState.WaterTemperature[surfaceWaterIndex] * waterMassSurface;
			display.SolarRadiation += tempState.DisplaySolarRadiation[i];
			display.GeothermalRadiation += tempState.GeothermalRadiation[i];
			display.GlobalCloudMass += nextState.CloudMass[i];
			display.GlobalIceMass += nextState.IceMass[i];
			display.GlobalOceanCoverage += tempState.WaterCoverage[staticState.GetWaterIndex(worldData.SurfaceWaterLayer, i)];
			display.GlobalSeaLevel += tempState.AirLayerElevation[staticState.GetLayerIndexAir(worldData.SurfaceAirLayer, i)];
			display.GlobalEvaporation += display.Evaporation[i];
			display.GlobalRainfall += display.Rainfall[i];
			display.GlobalEnthalpyTerrain += display.EnthalpyTerrain[i];
			display.GlobalEnthalpyIce += display.EnthalpyIce[i];
			display.GlobalEnthalpyCloud += display.EnthalpyCloud[i];
			display.GlobalEnthalpyGroundWater += display.EnthalpyGroundWater[i];
			display.GlobalTerrainTemperature += nextState.GroundTemperature[i];
			for (int j = 1; j < worldData.AirLayers - 1; j++)
			{
				int index = staticState.GetLayerIndexAir(j, i);
				display.GlobalAirTemperaturePotential += nextState.AirTemperaturePotential[index] * tempState.AirMass[index];
				display.GlobalAirMass += tempState.AirMass[index];
				display.GlobalWaterVapor += nextState.AirVapor[index];
				display.EnergySolarReflectedAtmosphere += tempState.SolarReflectedAir[j * staticState.Count + i];
				display.EnergySolarAbsorbedAtmosphere += tempState.SolarRadiationInAir[j * staticState.Count + i];
				display.GlobalEnthalpyAir += display.EnthalpyAir[index];
				display.GlobalCloudCoverage += math.min(1, tempState.AbsorptivitySolar[index].AbsorptivityCloud * 100);
				display.GlobalAirCarbon += nextState.AirCarbonDioxide[index];
				display.GlobalCondensationCloud += display.CondensationCloud[index];
				display.GlobalCondensationGround += display.CondensationGround[index];
			}
			display.EnergySolarAbsorbedOcean += tempState.SolarRadiationInWater[i];
			display.EnergySolarAbsorbedSurface += tempState.SolarRadiationInWater[i] + tempState.SolarRadiationInTerrain[i] + tempState.SolarRadiationInIce[i];
			display.EnergySolarReflectedSurface += tempState.SolarReflectedWater[i] + tempState.SolarReflectedTerrain[i] + tempState.SolarReflectedIce[i];
			for (int j = 1; j < worldData.WaterLayers - 1; j++)
			{
				int index = staticState.GetLayerIndexWater(j, i);
				float waterMass = nextState.WaterMass[index];
				globalWaterMass += waterMass;
				display.GlobalOceanTemperature += nextState.WaterTemperature[index] * waterMass;
				display.GlobalEnthalpyWater += display.EnthalpyWater[index];
				display.GlobalOceanMass += nextState.WaterMass[index];
				display.GlobalWaterCarbon += nextState.WaterCarbonDioxide[index];
			}
			display.EnergySurfaceConduction += tempState.ConductionAirIce[i] + tempState.ConductionAirTerrain[i] + tempState.ConductionAirWater[i];
			display.EnergyOceanConduction += tempState.ConductionAirWater[i];
			display.EnergyEvapotranspiration += (tempState.EvaporationMassWater[i] + tempState.FloraRespirationMassVapor[i]) * WorldData.LatentHeatWaterVapor;
#if !LayerRefactor
			display.EnergyThermalBackRadiation += tempState.WindowRadiationTransmittedDown[(worldData.AirLayer0 + 1) * staticState.Count + i] + tempState.ThermalRadiationTransmittedDown[(worldData.AirLayer0 + 1) * staticState.Count + i];
			display.EnergyThermalOceanRadiation += (tempState.WindowRadiationTransmittedUp[worldData.SurfaceWaterLayerGlobal * staticState.Count + i] + tempState.ThermalRadiationTransmittedUp[worldData.SurfaceWaterLayerGlobal * staticState.Count + i]) * tempState.WaterCoverage[worldData.SurfaceWaterLayer][i];

			float surfaceRadiation = tempState.WindowRadiationTransmittedUp[worldData.IceLayer * staticState.Count + i] + tempState.ThermalRadiationTransmittedUp[worldData.IceLayer * staticState.Count + i];
			float surfaceRadiationOutWindow = tempState.WindowRadiationTransmittedUp[worldData.IceLayer * staticState.Count + i];
			float radiationToSpace = tempState.ThermalRadiationTransmittedUp[(worldData.AirLayer0 + worldData.AirLayers - 2) * staticState.Count + i] + tempState.WindowRadiationTransmittedUp[(worldData.AirLayer0 + worldData.AirLayers - 2) * staticState.Count + i] - surfaceRadiationOutWindow;
			display.EnergyThermalAbsorbedAtmosphere += surfaceRadiation - surfaceRadiationOutWindow;
			display.EnergyThermalOutAtmosphere += radiationToSpace;
			display.EnergyThermalSurfaceOutAtmosphericWindow += surfaceRadiationOutWindow;
			display.EnergyThermalSurfaceRadiation += surfaceRadiation;
#endif
		}
		display.GlobalAirTemperaturePotential /= display.GlobalAirMass;
		display.GlobalOceanSurfaceTemperature /= globalWaterSurfaceMass;
		display.GlobalOceanTemperature /= globalWaterMass;
		display.GlobalTerrainTemperature /= staticState.Count;
		display.GlobalEnthalpy = display.GlobalEnthalpyTerrain + display.GlobalEnthalpyAir + display.GlobalEnthalpyWater + display.GlobalEnthalpyCloud + display.GlobalEnthalpyDeltaTerrain + display.GlobalEnthalpyDeltaIce;
	}


	[BurstCompile]
	private struct UpdateDisplayJob : IJobParallelFor {
		public NativeArray<float> DisplayPlate;
		public NativeArray<float> SolarRadiationAbsorbedSurface;
		public NativeArray<float> DisplayPrecipitation;
		public NativeArray<float> DisplayEvaporation;
		public NativeArray<float> EnthalpyTerrain;
		public NativeArray<float> EnthalpyCloud;
		public NativeArray<float> EnthalpyIce;
		public NativeArray<float> EnthalpyGroundWater;
		[ReadOnly] public NativeArray<short> Plate;
		[ReadOnly] public NativeArray<float> SolarRadiationInTerrain;
		[ReadOnly] public NativeArray<float> SolarRadiationInIce;
		[ReadOnly] public NativeArray<float> SolarRadiationInWaterSurface;
		[ReadOnly] public NativeArray<float> Precipitation;
		[ReadOnly] public NativeArray<float> EvaporationWater;
		[ReadOnly] public NativeArray<float> EvaporationFlora;
		[ReadOnly] public NativeArray<float> TerrainTemperature;
		[ReadOnly] public NativeArray<float> CloudMass;
		[ReadOnly] public NativeArray<float> IceMass;
		[ReadOnly] public NativeArray<float> IceTemperature;
		[ReadOnly] public NativeArray<float> GroundWaterMass;
		[ReadOnly] public NativeArray<float> GroundWaterTemperature;
		[ReadOnly] public NativeArray<float> SoilFertility;
		[ReadOnly] public float HeatingDepth;
		public void Execute(int i)
		{
			DisplayPlate[i] = Plate[i];
			SolarRadiationAbsorbedSurface[i] = SolarRadiationInTerrain[i] + SolarRadiationInIce[i] + SolarRadiationInWaterSurface[i];
			DisplayPrecipitation[i] = Precipitation[i];
			DisplayEvaporation[i] = EvaporationWater[i] + EvaporationFlora[i];
			float vegetation = 0;
			float vegetationWater = 0;
			EnthalpyTerrain[i] = TerrainTemperature[i] * Atmosphere.GetSpecificHeatTerrain(HeatingDepth, SoilFertility[i], vegetation, vegetationWater) + vegetationWater * WorldData.LatentHeatWaterLiquid;
			EnthalpyCloud[i] = CloudMass[i] * WorldData.LatentHeatWaterLiquid;
			EnthalpyIce[i] = IceMass[i] * IceTemperature[i] * WorldData.SpecificHeatIce;
			EnthalpyGroundWater[i] = GroundWaterMass[i] * (WorldData.LatentHeatWaterLiquid + WorldData.SpecificHeatWater * GroundWaterTemperature[i]);
		}
	}

	[BurstCompile]
	private struct InitDisplayAirLayerJob : IJobParallelFor {
		public NativeSlice<float> DisplayPressure;
		public NativeSlice<float3> DisplayPressureGradientForce;
		public NativeSlice<float> DisplayCondensationGround;
		public NativeSlice<float> DisplayCondensationCloud;
		public NativeSlice<float> Enthalpy;
		public NativeSlice<float> DustCoverage;
		public NativeSlice<float> CarbonDioxidePercent;
		public NativeSlice<float> Divergence;
		[ReadOnly] public NativeSlice<float> AirTemperaturePotential;
		[ReadOnly] public NativeSlice<float> AirPressure;
		[ReadOnly] public NativeSlice<float> LayerMiddle;
		[ReadOnly] public NativeSlice<float> CondensationCloud;
		[ReadOnly] public NativeSlice<float> CondensationGround;
		[ReadOnly] public NativeSlice<float3> PressureGradientForce;
		[ReadOnly] public NativeSlice<float> AirMass;
		[ReadOnly] public NativeSlice<float> VaporMass;
		[ReadOnly] public NativeSlice<float> DustMass;
		[ReadOnly] public NativeSlice<float> CarbonDioxide;
		[ReadOnly] public NativeArray<float> CloudMass;
		[ReadOnly] public NativeArray<float> CloudElevation;
		[ReadOnly] public NativeSlice<float> LayerElevation;
		[ReadOnly] public NativeSlice<float> LayerHeight;
		[ReadOnly] public float Gravity;
		[ReadOnly] public int Count;
		public void Execute(int i)
		{
			DisplayPressure[i] = Atmosphere.GetPressureAtElevation(0, Gravity, AirPressure[i], AirTemperaturePotential[i], LayerMiddle[i]);
			DisplayPressureGradientForce[i] = PressureGradientForce[i];
			DisplayCondensationGround[i] = CondensationCloud[i];
			DisplayCondensationCloud[i] = CondensationGround[i];
			DustCoverage[i] += DustMass[i];
			CarbonDioxidePercent[i] = CarbonDioxide[i] / AirMass[i];
			Divergence[i] /= AirMass[i];
			if (AirMass[i] > 0)
			{
				int columnIndex = i % Count;
				float cloudMass = Atmosphere.GetCloudMassInLayer(CloudMass[columnIndex], CloudElevation[columnIndex], LayerElevation[i], LayerHeight[i]);
				Enthalpy[i] = AirTemperaturePotential[i] * Atmosphere.GetSpecificHeatAir(AirMass[i], VaporMass[i], cloudMass) 
					+ VaporMass[i] * (WorldData.LatentHeatWaterLiquid + WorldData.LatentHeatWaterVapor);
			}
		}
	}
	[BurstCompile]
	private struct InitDisplayWaterLayerJob : IJobParallelFor {
		public NativeSlice<float> Enthalpy;
		public NativeSlice<float> Salinity;
		public NativeSlice<float> CarbonPercent;
		[ReadOnly] public NativeArray<float> WaterTemperature;
		[ReadOnly] public NativeArray<float> WaterMass;
		[ReadOnly] public NativeArray<float> SaltMass;
		[ReadOnly] public NativeArray<float> WaterCarbon;
		[ReadOnly] public int CountPerLayer;
		public void Execute(int i)
		{
			Debug.Assert(CountPerLayer > 0);
			int index = i + CountPerLayer;
			Salinity[i] = Atmosphere.GetWaterSalinity(WaterMass[index], SaltMass[index]);
			if (WaterMass[index] > 0)
			{
				Enthalpy[i] = WaterTemperature[index] * (WorldData.SpecificHeatWater * WaterMass[index] + WorldData.SpecificHeatSalt * SaltMass[index]) + WaterMass[index] * WorldData.LatentHeatWaterLiquid;
				CarbonPercent[i] = WaterCarbon[index] / (WaterMass[index] + SaltMass[index] + WaterCarbon[index]);
			}
		}
	}


}

