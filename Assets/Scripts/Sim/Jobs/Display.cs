#define LayerRefactor

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


public struct DisplayState {
	public static JobHelper DisplayJob;
	public static JobHelper DisplayJobAir;

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
	public float GlobalAirCarbon;
	public float GlobalWaterCarbon;
	public float GlobalCloudCoverage;
	public float GlobalEvaporation;
	public float GlobalRainfall;
	public float GlobalFloraMass;
	public float GlobalPlanktonMass;
	public float GlobalSoilFertility;
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
	public NativeArray<float> CarbonDioxidePercent;
	public NativeArray<float> EnthalpyAir;
	public NativeArray<float> Pressure;
	public NativeArray<float> DivergenceAir;
	public NativeArray<float3> PressureGradientForce;
	public NativeArray<float> WindVertical;
	public NativeArray<SolarAbsorptivity> AbsorptionSolar;
	public NativeArray<ThermalAbsorptivity> AbsorptionThermal;
	public NativeArray<float>[] EnthalpyWater;
	public NativeArray<float>[] Salinity;
	public NativeArray<float>[] WaterCarbonDioxidePercent;
	public NativeArray<float>[] DivergenceWater;
	public NativeArray<float>[] ThermalDelta;
	public NativeArray<float>[] ConductionDelta;
	public NativeArray<float>[] SolarDelta;
	public NativeArray<float>[] LatentHeatDelta;

	private bool _initialized;

	public void Init(int count, ref WorldData worldData)
	{
		_initialized = true;
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

		DivergenceAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		Pressure = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		PressureGradientForce = new NativeArray<float3>(count * worldData.AirLayers, Allocator.Persistent);
		EnthalpyAir = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		WindVertical = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AbsorptionSolar = new NativeArray<SolarAbsorptivity>(count * worldData.AirLayers, Allocator.Persistent);
		AbsorptionThermal = new NativeArray<ThermalAbsorptivity>(count * worldData.AirLayers, Allocator.Persistent);
		CarbonDioxidePercent = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);

		WaterCarbonDioxidePercent = new NativeArray<float>[worldData.WaterLayers];
		Salinity = new NativeArray<float>[worldData.WaterLayers];
		EnthalpyWater = new NativeArray<float>[worldData.WaterLayers];
		DivergenceWater = new NativeArray<float>[worldData.WaterLayers];
		for (int i = 0; i < worldData.WaterLayers; i++)
		{
			WaterCarbonDioxidePercent[i] = new NativeArray<float>(count, Allocator.Persistent);
			Salinity[i] = new NativeArray<float>(count, Allocator.Persistent);
			EnthalpyWater[i] = new NativeArray<float>(count, Allocator.Persistent);
			DivergenceWater[i] = new NativeArray<float>(count, Allocator.Persistent);
		}

		ThermalDelta = new NativeArray<float>[worldData.LayerCount];
		ConductionDelta = new NativeArray<float>[worldData.LayerCount];
		SolarDelta = new NativeArray<float>[worldData.LayerCount];
		LatentHeatDelta = new NativeArray<float>[worldData.LayerCount];
		for (int i = 0; i < worldData.LayerCount; i++)
		{
			ThermalDelta[i] = new NativeArray<float>(count, Allocator.Persistent);
			ConductionDelta[i] = new NativeArray<float>(count, Allocator.Persistent);
			SolarDelta[i] = new NativeArray<float>(count, Allocator.Persistent);
			LatentHeatDelta[i] = new NativeArray<float>(count, Allocator.Persistent);
		}

	}

	public void Dispose()
	{
		if (!_initialized)
		{
			return;
		}
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

		DivergenceAir.Dispose();
		Pressure.Dispose();
		PressureGradientForce.Dispose();
		EnthalpyAir.Dispose();
		WindVertical.Dispose();
		AbsorptionSolar.Dispose();
		AbsorptionThermal.Dispose();
		CarbonDioxidePercent.Dispose();

		for (int i = 0; i < Salinity.Length; i++)
		{
			WaterCarbonDioxidePercent[i].Dispose();
			Salinity[i].Dispose();
			EnthalpyWater[i].Dispose();
			DivergenceWater[i].Dispose();
		}
		for (int i = 0; i < ThermalDelta.Length; i++)
		{
			ThermalDelta[i].Dispose();
			ConductionDelta[i].Dispose();
			SolarDelta[i].Dispose();
			LatentHeatDelta[i].Dispose();
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
			DisplayJobAir = new JobHelper(staticState.Count * worldData.AirLayers);
		}

		JobHandle initDisplayAirHandle = default(JobHandle);
		JobHandle initDisplayWaterHandle = default(JobHandle);

#if !LayerRefactor
		initDisplayAirHandle = JobHandle.CombineDependencies(initDisplayAirHandle, DisplayJobAir.Schedule(new GetDivergenceJob()
		{
			Divergence = display.DivergenceAir,
			Destination = tempState.DestinationAirResolved,
			Neighbors = staticState.NeighborsAir,
			Mass = tempState.AirMass,
		}));

		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			int layer = worldData.WaterLayer0 + j;
			initDisplayWaterHandle = JobHandle.CombineDependencies(initDisplayWaterHandle, DisplayJob.Schedule(new GetDivergenceJob()
			{
				Divergence = display.DivergenceWater[j],
				Destination = tempState.DestinationWater[j],
				Neighbors = staticState.Neighbors,
				Mass = nextState.WaterMass[j],
			}));
		}
#endif



		tempState.AbsorptivitySolar.CopyTo(display.AbsorptionSolar);
		tempState.AbsorptivityThermal.CopyTo(display.AbsorptionThermal);
		initDisplayAirHandle = JobHandle.CombineDependencies(initDisplayAirHandle, DisplayJobAir.Schedule(new InitDisplayAirLayerJob()
		{
			DisplayPressure = staticState.GetSliceAir(display.Pressure),
			DisplayPressureGradientForce = staticState.GetSliceAir(display.PressureGradientForce),
			DisplayCondensationGround = staticState.GetSliceAir(display.CondensationGround),
			DisplayCondensationCloud = staticState.GetSliceAir(display.CondensationCloud),
			Enthalpy = staticState.GetSliceAir(display.EnthalpyAir),
			DustCoverage = staticState.GetSliceAir(display.DustMass),
			CarbonDioxidePercent = staticState.GetSliceAir(display.CarbonDioxidePercent),
			Divergence = staticState.GetSliceAir(display.DivergenceAir),

			CarbonDioxide = staticState.GetSliceAir(nextState.AirCarbon),
			AirTemperaturePotential = staticState.GetSliceAir(nextState.AirTemperaturePotential),
			AirPressure = staticState.GetSliceAir(tempState.AirPressure),
			LayerMiddle = staticState.GetSliceAir(tempState.AirLayerMiddle),
			PressureGradientForce = staticState.GetSliceAir(tempState.AirAcceleration),
			CondensationCloud = staticState.GetSliceAir(tempState.CondensationCloudMass),
			CondensationGround = staticState.GetSliceAir(tempState.CondensationGroundMass),
			AirMass = staticState.GetSliceAir(tempState.AirMass),
			VaporMass = staticState.GetSliceAir(nextState.AirVapor),
			DustMass = staticState.GetSliceAir(nextState.Dust),
			Gravity = nextState.PlanetState.Gravity,
		}, initDisplayAirHandle));

		for (int i = 1; i < worldData.WaterLayers - 1; i++)
		{
			initDisplayWaterHandle = JobHandle.CombineDependencies(initDisplayWaterHandle, DisplayJob.Schedule(new InitDisplayWaterLayerJob()
			{
				Enthalpy = display.EnthalpyWater[i],
				Salinity = display.Salinity[i],
				CarbonPercent = display.WaterCarbonDioxidePercent[i],

				WaterTemperature = nextState.WaterTemperature[i],
				SaltMass = nextState.SaltMass[i],
				WaterMass = nextState.WaterMass[i],
				WaterCarbon = nextState.WaterCarbon[i]
			}, initDisplayWaterHandle));
		}
		var updateDisplayJobHandle = DisplayJob.Schedule(new UpdateDisplayJob()
		{
			SolarRadiationAbsorbedSurface = display.SolarRadiationAbsorbedSurface,
			DisplayEvaporation = display.Evaporation,
			DisplayPrecipitation = display.Rainfall,
			EnthalpyTerrain = display.EnthalpyTerrain,
			EnthalpyCloud = display.EnthalpyCloud,
			EnthalpyFlora = display.EnthalpyFlora,
			EnthalpyIce = display.EnthalpyIce,
			EnthalpyGroundWater = display.EnthalpyGroundWater,

			SolarRadiationInTerrain = tempState.SolarRadiationIn[worldData.TerrainLayer],
			SolarRadiationInIce = tempState.SolarRadiationIn[worldData.IceLayer],
			SolarRadiationInWaterSurface = tempState.SolarRadiationIn[worldData.SurfaceWaterLayerGlobal],
			EvaporationWater = tempState.EvaporationMassWater,
			EvaporationFlora = tempState.FloraRespirationMassVapor,
			Precipitation = tempState.PrecipitationMass,
			SoilFertility = nextState.GroundCarbon,
			TerrainTemperature = nextState.GroundTemperature,
			Flora = nextState.FloraMass,
			FloraWater = nextState.FloraWater,
			FloraTemperature = nextState.FloraTemperature,
			HeatingDepth = worldData.SoilHeatDepth,
			CloudMass = nextState.CloudMass,
			IceMass = nextState.IceMass,
			IceTemperature = nextState.IceTemperature,
			GroundWaterMass = nextState.GroundWater,
			GroundWaterTemperature = nextState.GroundWaterTemperature
		});

		for (int i = 0; i < worldData.LayerCount; i++)
		{
			tempState.SolarRadiationIn[i].CopyTo(display.SolarDelta[i]);
			tempState.ThermalRadiationDelta[i].CopyTo(display.ThermalDelta[i]);
		}

		updateDisplayJobHandle = JobHandle.CombineDependencies(initDisplayAirHandle, initDisplayWaterHandle, updateDisplayJobHandle);
		if (settings.CollectGlobals)
		{
			updateDisplayJobHandle.Complete();

			UpdateGlobals(ref display, ref worldData, ref tempState, ref nextState, ref staticState);

			display.GlobalEnthalpyDelta = display.GlobalEnthalpy - lastDisplay.GlobalEnthalpy;
			display.GlobalEnthalpyDeltaTerrain = display.GlobalEnthalpyTerrain - lastDisplay.GlobalEnthalpyTerrain;
			display.GlobalEnthalpyDeltaAir = display.GlobalEnthalpyAir - lastDisplay.GlobalEnthalpyAir;
			display.GlobalEnthalpyDeltaWater = display.GlobalEnthalpyWater - lastDisplay.GlobalEnthalpyWater;
			display.GlobalEnthalpyDeltaFlora = display.GlobalEnthalpyFlora - lastDisplay.GlobalEnthalpyFlora;
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
		ref StaticState staticState
		)
	{
		float globalWaterMass = 0;
		float globalWaterSurfaceMass = 0;
		for (int i = 0; i < staticState.Count; i++)
		{
			float waterMassSurface = nextState.WaterMass[worldData.SurfaceWaterLayer][i];
			globalWaterSurfaceMass += waterMassSurface;
			display.GlobalFloraMass += nextState.FloraMass[i];
			display.GlobalPlanktonMass += nextState.PlanktonMass[worldData.SurfaceWaterLayer][i];
			display.GlobalSoilFertility += nextState.GroundCarbon[i];
			display.GlobalOceanSurfaceTemperature += nextState.WaterTemperature[worldData.SurfaceWaterLayer][i] * waterMassSurface;
			display.SolarRadiation += tempState.DisplaySolarRadiation[i];
			display.GeothermalRadiation += tempState.GeothermalRadiation[i];
			display.GlobalCloudMass += nextState.CloudMass[i];
			display.GlobalIceMass += nextState.IceMass[i];
			display.GlobalOceanCoverage += tempState.WaterCoverage[worldData.SurfaceWaterLayer][i];
			display.GlobalSurfaceTemperature += tempState.SurfaceAirTemperatureAbsolute[i];
			display.GlobalSeaLevel += tempState.AirLayerElevation[staticState.GetLayerIndexAir(worldData.SurfaceAirLayer, i)];
			display.GlobalEvaporation += display.Evaporation[i];
			display.GlobalRainfall += display.Rainfall[i];
			display.GlobalCondensationCloud += display.CondensationCloud[i];
			display.GlobalCondensationGround += display.CondensationGround[i];
			display.GlobalEnthalpyTerrain += display.EnthalpyTerrain[i];
			display.GlobalEnthalpyFlora += display.EnthalpyFlora[i];
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
				display.EnergySolarReflectedAtmosphere += tempState.SolarReflected[j + worldData.AirLayer0][i];
				display.EnergySolarAbsorbedAtmosphere += tempState.SolarRadiationIn[j + worldData.AirLayer0][i];
				display.GlobalEnthalpyAir += display.EnthalpyAir[index];
				display.GlobalCloudCoverage += math.min(1, tempState.AbsorptivitySolar[index].AbsorptivityCloud * 100);
				display.GlobalAirCarbon += nextState.AirCarbon[index];
			}
			display.EnergySolarAbsorbedSurface += tempState.SolarRadiationIn[worldData.TerrainLayer][i] + tempState.SolarRadiationIn[worldData.IceLayer][i];
			display.EnergySolarReflectedSurface += tempState.SolarReflected[worldData.TerrainLayer][i] + tempState.SolarReflected[worldData.IceLayer][i];
			for (int j = 1; j < worldData.WaterLayers - 1; j++)
			{
				float waterMass = nextState.WaterMass[j][i];
				globalWaterMass += waterMass;
				display.GlobalOceanTemperature += nextState.WaterTemperature[j][i] * waterMass;

				float absorbed = tempState.SolarRadiationIn[j + worldData.WaterLayer0][i];
				display.EnergySolarAbsorbedOcean += absorbed;
				display.EnergySolarAbsorbedSurface += absorbed;
				display.EnergySolarReflectedSurface += tempState.SolarReflected[j + worldData.WaterLayer0][i];
				display.GlobalEnthalpyWater += display.EnthalpyWater[j][i];
				display.GlobalOceanMass += nextState.WaterMass[j][i];
				display.GlobalWaterCarbon += nextState.WaterCarbon[j][i];
			}
			display.EnergySurfaceConduction += tempState.ConductionAirIce[i] + tempState.ConductionAirTerrain[i] + tempState.ConductionAirWater[i];
			display.EnergyOceanConduction += tempState.ConductionAirWater[i];
			display.EnergyEvapotranspiration += (tempState.EvaporationMassWater[i] + tempState.FloraRespirationMassVapor[i]) * WorldData.LatentHeatWaterVapor;
			display.EnergyThermalBackRadiation += tempState.WindowRadiationTransmittedDown[worldData.AirLayer0 + 1][i] + tempState.ThermalRadiationTransmittedDown[worldData.AirLayer0 + 1][i];
			display.EnergyThermalOceanRadiation += (tempState.WindowRadiationTransmittedUp[worldData.SurfaceWaterLayerGlobal][i] + tempState.ThermalRadiationTransmittedUp[worldData.SurfaceWaterLayerGlobal][i]) * tempState.WaterCoverage[worldData.SurfaceWaterLayer][i];

			float surfaceRadiation = tempState.WindowRadiationTransmittedUp[worldData.IceLayer][i] + tempState.ThermalRadiationTransmittedUp[worldData.IceLayer][i];
			float surfaceRadiationOutWindow = tempState.WindowRadiationTransmittedUp[worldData.IceLayer][i];
			float radiationToSpace = tempState.ThermalRadiationTransmittedUp[worldData.AirLayer0 + worldData.AirLayers - 2][i] + tempState.WindowRadiationTransmittedUp[worldData.AirLayer0 + worldData.AirLayers - 2][i] - surfaceRadiationOutWindow;
			display.EnergyThermalAbsorbedAtmosphere += surfaceRadiation - surfaceRadiationOutWindow;
			display.EnergyThermalOutAtmosphere += radiationToSpace;
			display.EnergyThermalSurfaceOutAtmosphericWindow += surfaceRadiationOutWindow;
			display.EnergyThermalSurfaceRadiation += surfaceRadiation;

		}
		display.GlobalAirTemperaturePotential /= display.GlobalAirMass;
		display.GlobalOceanSurfaceTemperature /= globalWaterSurfaceMass;
		display.GlobalOceanTemperature /= globalWaterMass;
		display.GlobalTerrainTemperature /= staticState.Count;
		display.GlobalEnthalpy = display.GlobalEnthalpyTerrain + display.GlobalEnthalpyAir + display.GlobalEnthalpyWater + display.GlobalEnthalpyCloud + display.GlobalEnthalpyDeltaFlora + display.GlobalEnthalpyDeltaTerrain + display.GlobalEnthalpyDeltaIce;
	}

	[BurstCompile]
	private struct UpdateDisplayJob : IJobParallelFor {
		public NativeArray<float> SolarRadiationAbsorbedSurface;
		public NativeArray<float> DisplayPrecipitation;
		public NativeArray<float> DisplayEvaporation;
		public NativeArray<float> EnthalpyTerrain;
		public NativeArray<float> EnthalpyFlora;
		public NativeArray<float> EnthalpyCloud;
		public NativeArray<float> EnthalpyIce;
		public NativeArray<float> EnthalpyGroundWater;
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
		[ReadOnly] public NativeArray<float> Flora;
		[ReadOnly] public NativeArray<float> FloraWater;
		[ReadOnly] public NativeArray<float> FloraTemperature;
		[ReadOnly] public NativeArray<float> GroundWaterMass;
		[ReadOnly] public NativeArray<float> GroundWaterTemperature;
		[ReadOnly] public NativeArray<float> SoilFertility;
		[ReadOnly] public float HeatingDepth;
		public void Execute(int i)
		{
			SolarRadiationAbsorbedSurface[i] = SolarRadiationInTerrain[i] + SolarRadiationInIce[i] + SolarRadiationInWaterSurface[i];
			DisplayPrecipitation[i] = Precipitation[i];
			DisplayEvaporation[i] = EvaporationWater[i] + EvaporationFlora[i];
			EnthalpyTerrain[i] = TerrainTemperature[i] * Atmosphere.GetSpecificHeatTerrain(HeatingDepth, SoilFertility[i]);
			EnthalpyFlora[i] = FloraTemperature[i] * (Flora[i] * WorldData.SpecificHeatFlora + FloraWater[i] * (WorldData.SpecificHeatWater + WorldData.LatentHeatWaterLiquid));
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
		[ReadOnly] public float Gravity;
		public void Execute(int i)
		{
			DisplayPressure[i] = Atmosphere.GetPressureAtElevation(0, Gravity, AirPressure[i], AirTemperaturePotential[i], LayerMiddle[i]);
			DisplayPressureGradientForce[i] = PressureGradientForce[i];
			DisplayCondensationGround[i] += CondensationCloud[i];
			DisplayCondensationCloud[i] += CondensationGround[i];
			DustCoverage[i] += DustMass[i];
			CarbonDioxidePercent[i] += CarbonDioxide[i] / AirMass[i];
			Divergence[i] /= AirMass[i];
			if (AirMass[i] > 0)
			{
				Enthalpy[i] = AirTemperaturePotential[i] * (WorldData.SpecificHeatAtmosphere * AirMass[i] + WorldData.SpecificHeatWaterVapor * VaporMass[i]) + VaporMass[i] * (WorldData.LatentHeatWaterLiquid + WorldData.LatentHeatWaterVapor);
			}
		}
	}
	[BurstCompile]
	private struct InitDisplayWaterLayerJob : IJobParallelFor {
		public NativeArray<float> Enthalpy;
		public NativeArray<float> Salinity;
		public NativeArray<float> CarbonPercent;
		[ReadOnly] public NativeArray<float> WaterTemperature;
		[ReadOnly] public NativeArray<float> WaterMass;
		[ReadOnly] public NativeArray<float> SaltMass;
		[ReadOnly] public NativeArray<float> WaterCarbon;
		public void Execute(int i)
		{
			Salinity[i] = Atmosphere.GetWaterSalinity(WaterMass[i], SaltMass[i]);
			if (WaterMass[i] > 0)
			{
				Enthalpy[i] = WaterTemperature[i] * (WorldData.SpecificHeatWater * WaterMass[i] + WorldData.SpecificHeatSalt * SaltMass[i]) + WaterMass[i] * WorldData.LatentHeatWaterLiquid;
				CarbonPercent[i] = WaterCarbon[i] / (WaterMass[i] + SaltMass[i] + WaterCarbon[i]);
			}
		}
	}


}

