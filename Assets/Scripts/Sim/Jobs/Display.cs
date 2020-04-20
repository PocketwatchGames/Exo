using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


public struct DisplayState {
	public static JobHelper DisplayJob;

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
	public NativeArray<float>[] WaterCarbonDioxidePercent;
	public NativeArray<float>[] CarbonDioxidePercent;
	public NativeArray<float>[] EnthalpyWater;
	public NativeArray<float>[] EnthalpyAir;
	public NativeArray<float>[] Salinity;
	public NativeArray<float>[] Pressure;
	public NativeArray<float>[] Divergence;
	public NativeArray<float3>[] PressureGradientForce;
	public NativeArray<float>[] ThermalDelta;
	public NativeArray<float>[] ConductionDelta;
	public NativeArray<float>[] SolarDelta;
	public NativeArray<float>[] LatentHeatDelta;
	public NativeArray<float>[] WindVertical;
	public NativeArray<SolarAbsorptivity>[] AbsorptionSolar;
	public NativeArray<ThermalAbsorptivity>[] AbsorptionThermal;

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

		Divergence = new NativeArray<float>[worldData.AirLayers];
		Pressure = new NativeArray<float>[worldData.AirLayers];
		PressureGradientForce = new NativeArray<float3>[worldData.AirLayers];
		CarbonDioxidePercent = new NativeArray<float>[worldData.AirLayers];
		EnthalpyAir = new NativeArray<float>[worldData.AirLayers];
		WindVertical = new NativeArray<float>[worldData.AirLayers];
		AbsorptionSolar = new NativeArray<SolarAbsorptivity>[worldData.AirLayers];
		AbsorptionThermal = new NativeArray<ThermalAbsorptivity>[worldData.AirLayers];
		for (int i = 0; i < worldData.AirLayers; i++)
		{
			Divergence[i] = new NativeArray<float>(count, Allocator.Persistent);
			Pressure[i] = new NativeArray<float>(count, Allocator.Persistent);
			PressureGradientForce[i] = new NativeArray<float3>(count, Allocator.Persistent);
			EnthalpyAir[i] = new NativeArray<float>(count, Allocator.Persistent);
			WindVertical[i] = new NativeArray<float>(count, Allocator.Persistent);
			AbsorptionSolar[i] = new NativeArray<SolarAbsorptivity>(count, Allocator.Persistent);
			AbsorptionThermal[i] = new NativeArray<ThermalAbsorptivity>(count, Allocator.Persistent);
			CarbonDioxidePercent[i] = new NativeArray<float>(count, Allocator.Persistent);
		}

		WaterCarbonDioxidePercent = new NativeArray<float>[worldData.WaterLayers];
		Salinity = new NativeArray<float>[worldData.WaterLayers];
		EnthalpyWater = new NativeArray<float>[worldData.WaterLayers];
		for (int i = 0; i < worldData.WaterLayers; i++)
		{
			WaterCarbonDioxidePercent[i] = new NativeArray<float>(count, Allocator.Persistent);
			Salinity[i] = new NativeArray<float>(count, Allocator.Persistent);
			EnthalpyWater[i] = new NativeArray<float>(count, Allocator.Persistent);
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
		for (int i = 0; i < Pressure.Length; i++)
		{
			Divergence[i].Dispose();
			Pressure[i].Dispose();
			PressureGradientForce[i].Dispose();
			EnthalpyAir[i].Dispose();
			WindVertical[i].Dispose();
			AbsorptionSolar[i].Dispose();
			AbsorptionThermal[i].Dispose();
			CarbonDioxidePercent[i].Dispose();
		}
		for (int i = 0; i < Salinity.Length; i++)
		{
			WaterCarbonDioxidePercent[i].Dispose();
			Salinity[i].Dispose();
			EnthalpyWater[i].Dispose();
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
		ref WorldData worldData, 
		ref TempState tempState, 
		ref SimState nextState, 
		ref StaticState staticState
		)
	{
		if (DisplayJob == null)
		{
			DisplayJob = new JobHelper(staticState.Count);
		}

		JobHandle initDisplayAirHandle = default(JobHandle);
		JobHandle initDisplayWaterHandle = default(JobHandle);
		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			tempState.AbsorptivitySolar[i].CopyTo(display.AbsorptionSolar[i]);
			tempState.AbsorptivityThermal[i].CopyTo(display.AbsorptionThermal[i]);

			initDisplayAirHandle = JobHandle.CombineDependencies(initDisplayAirHandle, DisplayJob.Schedule(new InitDisplayAirLayerJob()
			{
				DisplayPressure = display.Pressure[i],
				DisplayPressureGradientForce = display.PressureGradientForce[i],
				DisplayCondensationGround = display.CondensationGround,
				DisplayCondensationCloud = display.CondensationCloud,
				Enthalpy = display.EnthalpyAir[i],
				DustCoverage = display.DustMass,
				CarbonDioxidePercent = display.CarbonDioxidePercent[i],

				CarbonDioxide = nextState.AirCarbon[i],
				Gravity = nextState.PlanetState.Gravity,
				AirTemperaturePotential = nextState.AirTemperaturePotential[i],
				AirPressure = tempState.AirPressure[i],
				LayerMiddle = tempState.LayerMiddle[i],
				PressureGradientForce = tempState.AirAcceleration[i],
				CondensationCloud = tempState.CondensationCloudMass[i],
				CondensationGround = tempState.CondensationGroundMass[i],
				AirMass = tempState.AirMass[i],
				VaporMass = nextState.AirVapor[i],
				DustMass = nextState.Dust[i],
			}, initDisplayAirHandle));
		}

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
			}));
		}
		for (int i = 0; i < worldData.LayerCount; i++)
		{
			tempState.SolarRadiationIn[i].CopyTo(display.SolarDelta[i]);
			tempState.ThermalRadiationDelta[i].CopyTo(display.ThermalDelta[i]);
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
			SolarRadiationInWaterSurface = tempState.SolarRadiationIn[worldData.WaterLayer0 + worldData.SurfaceWaterLayer],
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
		return JobHandle.CombineDependencies(initDisplayAirHandle, initDisplayWaterHandle, updateDisplayJobHandle);
	}

	public static void UpdateGlobals(
		ref DisplayState display, 
		ref WorldData worldData, 
		ref TempState tempState, 
		ref SimState nextState, 
		ref StaticState staticState,
		NativeArray<float> displaySolarRadiation,
		NativeArray<float> geothermalRadiation,
		NativeArray<float>[] solarReflected,
		NativeArray<float>[] solarRadiationIn,
		NativeArray<SolarAbsorptivity>[] absorptivitySolar,
		NativeArray<float> conductionAirIce,
		NativeArray<float> conductionAirTerrain,
		NativeArray<float> conductionAirWater,
		NativeArray<float> evaporationMassWater,
		NativeArray<float> floraRespirationMassVapor,
		NativeArray<float>[] windowRadiationTransmittedUp,
		NativeArray<float>[] windowRadiationTransmittedDown,
		NativeArray<float>[] thermalRadiationTransmittedUp,
		NativeArray<float>[] thermalRadiationTransmittedDown
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
			display.SolarRadiation += displaySolarRadiation[i];
			display.GeothermalRadiation += geothermalRadiation[i];
			display.GlobalCloudMass += nextState.CloudMass[i];
			display.GlobalIceMass += nextState.IceMass[i];
			display.GlobalOceanCoverage += tempState.WaterCoverage[worldData.SurfaceWaterLayer][i];
			display.GlobalSurfaceTemperature += tempState.SurfaceAirTemperatureAbsolute[i];
			display.GlobalSeaLevel += tempState.LayerElevation[1][i];
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
				display.GlobalAirTemperaturePotential += nextState.AirTemperaturePotential[j][i] * tempState.AirMass[j][i];
				display.GlobalAirMass += tempState.AirMass[j][i];
				display.GlobalWaterVapor += nextState.AirVapor[j][i];
				display.EnergySolarReflectedAtmosphere += solarReflected[j + worldData.AirLayer0][i];
				display.EnergySolarAbsorbedAtmosphere += solarRadiationIn[j + worldData.AirLayer0][i];
				display.GlobalEnthalpyAir += display.EnthalpyAir[j][i];
				display.GlobalCloudCoverage += math.min(1, absorptivitySolar[j][i].AbsorptivityCloud * 100);
				display.GlobalAirCarbon += nextState.AirCarbon[j][i];
			}
			display.EnergySolarAbsorbedSurface += solarRadiationIn[worldData.TerrainLayer][i] + solarRadiationIn[worldData.IceLayer][i];
			display.EnergySolarReflectedSurface += solarReflected[worldData.TerrainLayer][i] + solarReflected[worldData.IceLayer][i];
			for (int j = 1; j < worldData.WaterLayers - 1; j++)
			{
				float waterMass = nextState.WaterMass[j][i];
				globalWaterMass += waterMass;
				display.GlobalOceanTemperature += nextState.WaterTemperature[j][i] * waterMass;

				float absorbed = solarRadiationIn[j + worldData.WaterLayer0][i];
				display.EnergySolarAbsorbedOcean += absorbed;
				display.EnergySolarAbsorbedSurface += absorbed;
				display.EnergySolarReflectedSurface += solarReflected[j + worldData.WaterLayer0][i];
				display.GlobalEnthalpyWater += display.EnthalpyWater[j][i];
				display.GlobalOceanMass += nextState.WaterMass[j][i];
				display.GlobalWaterCarbon += nextState.WaterCarbon[j][i];
			}
			display.EnergySurfaceConduction += conductionAirIce[i] + conductionAirTerrain[i] + conductionAirWater[i];
			display.EnergyOceanConduction += conductionAirWater[i];
			display.EnergyEvapotranspiration += (evaporationMassWater[i] + floraRespirationMassVapor[i]) * WorldData.LatentHeatWaterVapor;
			display.EnergyThermalBackRadiation += windowRadiationTransmittedDown[worldData.AirLayer0 + 1][i] + thermalRadiationTransmittedDown[worldData.AirLayer0 + 1][i];
			display.EnergyThermalOceanRadiation += (windowRadiationTransmittedUp[worldData.WaterLayer0 + worldData.SurfaceWaterLayer][i] + thermalRadiationTransmittedUp[worldData.WaterLayer0 + worldData.SurfaceWaterLayer][i]) * tempState.WaterCoverage[worldData.SurfaceWaterLayer][i];

			float surfaceRadiation = windowRadiationTransmittedUp[worldData.IceLayer][i] + thermalRadiationTransmittedUp[worldData.IceLayer][i];
			float surfaceRadiationOutWindow = windowRadiationTransmittedUp[worldData.IceLayer][i];
			float radiationToSpace = thermalRadiationTransmittedUp[worldData.AirLayer0 + worldData.AirLayers - 2][i] + windowRadiationTransmittedUp[worldData.AirLayer0 + worldData.AirLayers - 2][i] - surfaceRadiationOutWindow;
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
		public NativeArray<float> DisplayPressure;
		public NativeArray<float3> DisplayPressureGradientForce;
		public NativeArray<float> DisplayCondensationGround;
		public NativeArray<float> DisplayCondensationCloud;
		public NativeArray<float> Enthalpy;
		public NativeArray<float> DustCoverage;
		public NativeArray<float> CarbonDioxidePercent;
		[ReadOnly] public NativeArray<float> AirTemperaturePotential;
		[ReadOnly] public NativeArray<float> AirPressure;
		[ReadOnly] public NativeArray<float> LayerMiddle;
		[ReadOnly] public NativeArray<float> CondensationCloud;
		[ReadOnly] public NativeArray<float> CondensationGround;
		[ReadOnly] public NativeArray<float3> PressureGradientForce;
		[ReadOnly] public NativeArray<float> AirMass;
		[ReadOnly] public NativeArray<float> VaporMass;
		[ReadOnly] public NativeArray<float> DustMass;
		[ReadOnly] public NativeArray<float> CarbonDioxide;
		[ReadOnly] public float Gravity;
		public void Execute(int i)
		{
			DisplayPressure[i] = Atmosphere.GetPressureAtElevation(0, Gravity, AirPressure[i], AirTemperaturePotential[i], LayerMiddle[i]);
			DisplayPressureGradientForce[i] = PressureGradientForce[i];
			DisplayCondensationGround[i] += CondensationCloud[i];
			DisplayCondensationCloud[i] += CondensationGround[i];
			DustCoverage[i] += DustMass[i];
			CarbonDioxidePercent[i] += CarbonDioxide[i] / AirMass[i];
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

