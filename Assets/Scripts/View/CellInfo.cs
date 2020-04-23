using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity;
using UnityEngine;
using System.Globalization;
using Unity.Mathematics;

public static class CellInfo {

	public enum TemperatureUnits {
		Celsius,
		Farenheit,
		Kelvin,
	}

	public enum CellInfoType {
		Global,
		Enthalpy,
		Energy,
		Cell,
		Atmosphere,
		Water,
		Ground
	}


	public static string GetCellInfoGlobal(TemperatureUnits ActiveTemperatureUnits, float inverseCellCount, ref WorldData worldData, ref SimState state, ref TempState dependent, ref DisplayState display)
	{
		StringBuilder s = new StringBuilder();
		NumberFormatInfo nfi1 = new NumberFormatInfo() { NumberDecimalDigits = 1 };
		NumberFormatInfo nfi2 = new NumberFormatInfo() { NumberDecimalDigits = 2 };
		s.AppendFormat("CO2: {0:N0} ppm", display.GlobalAirCarbon * 1000000 / display.GlobalAirMass);
		s.AppendFormat("\nCO2 Mass Air: {0:N1} kg", display.GlobalAirCarbon);
		s.AppendFormat("\nCO2 Mass Water: {0:N1} kg", display.GlobalWaterCarbon);
		s.AppendFormat("\nFlora: {0:N1} kg", display.GlobalFloraMass);
		s.AppendFormat("\nPlankton: {0:N1} kg", display.GlobalPlanktonMass);
		s.AppendFormat("\nSoil Carbon: {0:N1} kg", display.GlobalSoilFertility);
		s.AppendFormat("\nCloud Coverage: {0:N1}%", display.GlobalCloudCoverage * 100 * inverseCellCount);
		s.AppendFormat("\nSurface Temp Air: {0}", GetTemperatureString(display.GlobalSurfaceTemperature * inverseCellCount, ActiveTemperatureUnits, 2));
		s.AppendFormat("\nSurface Temp Ocean: {0:N0}", GetTemperatureString(display.GlobalOceanSurfaceTemperature, ActiveTemperatureUnits, 2));
		s.AppendFormat("\nTemperature Air: {0}", GetTemperatureString((float)(display.GlobalAirTemperaturePotential), ActiveTemperatureUnits, 2));
		s.AppendFormat("\nTemperature Ocean: {0:N0}", GetTemperatureString(display.GlobalOceanTemperature, ActiveTemperatureUnits, 4));
		s.AppendFormat("\nTemperature Terrain: {0:N0}", GetTemperatureString((float)display.GlobalTerrainTemperature, ActiveTemperatureUnits, 4));
		s.AppendFormat("\nWater Vapor: {0:N0}", display.GlobalWaterVapor);
		s.AppendFormat("\nRainfall: {0:N3}", display.GlobalRainfall * worldData.TicksPerYear * inverseCellCount / WorldData.MassWater);
		s.AppendFormat("\nCondensationCloud: {0:N3}", display.GlobalCondensationCloud * worldData.TicksPerYear * inverseCellCount / WorldData.MassWater);
		s.AppendFormat("\nCondensationGround: {0:N3}", display.GlobalCondensationGround * worldData.TicksPerYear * inverseCellCount / WorldData.MassWater);
		s.AppendFormat("\nEvaporation: {0:N3}", display.GlobalEvaporation * worldData.TicksPerYear * inverseCellCount / WorldData.MassWater);
		s.AppendFormat("\nCloud Mass: {0:N2}", display.GlobalCloudMass);
		s.AppendFormat("\nGlobal Sea Level: {0:N2}", display.GlobalSeaLevel * inverseCellCount);
		s.AppendFormat("\nOcean Coverage: {0:N1}%", display.GlobalOceanCoverage * 100 * inverseCellCount);
		s.AppendFormat("\nOcean Mass: {0:N} M", display.GlobalOceanMass / 1000000);

		return s.ToString();
	}
	public static string GetCellInfoEnthalpy(TemperatureUnits ActiveTemperatureUnits, float inverseCellCount, ref WorldData worldData, ref SimState state, ref TempState dependent, ref DisplayState display)
	{
		StringBuilder s = new StringBuilder();
		s.AppendFormat("Enthalpy Delta: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDelta * inverseCellCount), worldData.SecondsPerTick));
		s.AppendFormat("\nEnthalpy Delta Terrain: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaTerrain * inverseCellCount), worldData.SecondsPerTick));
		s.AppendFormat("\nEnthalpy Delta Water: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaWater * inverseCellCount), worldData.SecondsPerTick));
		s.AppendFormat("\nEnthalpy Delta Air: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaAir * inverseCellCount), worldData.SecondsPerTick));
		s.AppendFormat("\nEnthalpy Delta Ice: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaIce * inverseCellCount), worldData.SecondsPerTick));
		s.AppendFormat("\nEnthalpy Delta Cloud: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaCloud * inverseCellCount), worldData.SecondsPerTick));
		s.AppendFormat("\nEnthalpy Delta Terrain: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaTerrain * inverseCellCount), worldData.SecondsPerTick));
		s.AppendFormat("\nEnthalpy Delta Flora: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaFlora * inverseCellCount), worldData.SecondsPerTick));
		s.AppendFormat("\nEnthalpy Delta Ground Water: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaGroundWater * inverseCellCount), worldData.SecondsPerTick));

		return s.ToString();
	}
	public static string GetCellInfoEnergy(TemperatureUnits ActiveTemperatureUnits, float inverseCellCount, ref WorldData worldData, ref SimState state, ref TempState dependent, ref DisplayState display)
	{
		StringBuilder s = new StringBuilder();
		NumberFormatInfo nfi1 = new NumberFormatInfo() { NumberDecimalDigits = 1 };
		NumberFormatInfo nfi2 = new NumberFormatInfo() { NumberDecimalDigits = 2 };

		var totalReflected = display.EnergySolarReflectedAtmosphere + display.EnergySolarReflectedSurface;
		var totalOutgoing = display.EnergyThermalSurfaceOutAtmosphericWindow + display.EnergyThermalOutAtmosphere;
		var emittedByAtmosphere = display.EnergyThermalOutAtmosphere - display.EnergyThermalSurfaceOutAtmosphericWindow;
		s.AppendFormat("Delta: {0:N1}", ConvertTileEnergyToWatts((display.SolarRadiation + display.GeothermalRadiation - totalReflected - totalOutgoing) * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nS Incoming: {0:N1}", ConvertTileEnergyToWatts(display.SolarRadiation * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nS Reflected: {0:N1}", ConvertTileEnergyToWatts((totalReflected) * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nS Reflected Atmos: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarReflectedAtmosphere * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nS Reflected Surf: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarReflectedSurface * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nS Abs Atm: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarAbsorbedAtmosphere * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nS Abs Surface Total: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarAbsorbedSurface * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nS Abs Ocean: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarAbsorbedOcean * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nT Outgoing: {0:N1}", ConvertTileEnergyToWatts(totalOutgoing * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nT Out Atm Window: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalSurfaceOutAtmosphericWindow * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nT Out Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalOutAtmosphere * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nT Surface Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalSurfaceRadiation * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nT Atm Absorbed: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalAbsorbedAtmosphere * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nT Atm Emitted: {0:N1}", ConvertTileEnergyToWatts(emittedByAtmosphere * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nT Back Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalBackRadiation * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nGeothermal Incoming: {0:N1}", ConvertTileEnergyToWatts(display.GeothermalRadiation * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nEvapotranspiration: {0:N1}", ConvertTileEnergyToWatts(display.EnergyEvapotranspiration * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nSurface Conduction: {0:N1}", ConvertTileEnergyToWatts(display.EnergySurfaceConduction * inverseCellCount, worldData.SecondsPerTick));
		s.AppendFormat("\nOcean Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalOceanRadiation * inverseCellCount / display.GlobalOceanCoverage, worldData.SecondsPerTick));
		s.AppendFormat("\nOcean Conduction: {0:N1}", ConvertTileEnergyToWatts(display.EnergyOceanConduction * inverseCellCount / display.GlobalOceanCoverage, worldData.SecondsPerTick));

		return s.ToString();
	}
	public static string GetCellInfoCell(TemperatureUnits ActiveTemperatureUnits, int ActiveCellIndex, ref StaticState staticState, ref SimState state, ref TempState dependent, ref DisplayState display)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		var coord = staticState.Coordinate[ActiveCellIndex];
		var pos = staticState.SphericalPosition[ActiveCellIndex];
		s.AppendFormat("INDEX: {0}", ActiveCellIndex);
		s.AppendFormat("\nCOORD: ({0:N1}, {1:N1})", math.degrees(coord.x), math.degrees(coord.y));
		s.AppendFormat("\nPOS: ({0:N2}, {1:N2}, {2:N2})", pos.x, pos.y, pos.z);
		//s.AppendFormat("\nSolar Terrain:   {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[0][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Terrain: {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[0][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Water0:    {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[2][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Water0:  {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[2][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Water1:    {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[3][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Water1:  {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[3][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Water2:    {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[4][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Water2:  {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[4][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Ice:       {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[6][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Ice:     {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[6][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Air0:       {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[8][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Air0:     {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[8][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Air1:       {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[1][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Air1:     {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[1][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Air2:       {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[2][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Air2:     {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[2][ActiveCellIndex]));

		return s.ToString();
	}
	public static string GetCellInfoAtmosphere(TemperatureUnits ActiveTemperatureUnits, int ActiveCellIndex, ref WorldData worldData, ref SimState state, ref TempState dependent, ref StaticState staticState, ref DisplayState display)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		float cloudMass = state.CloudMass[ActiveCellIndex];
		int upperAtmosphereLayerIndex = worldData.AirLayers - 2;
		s.AppendFormat("RAIN: {0:N3} kg", display.Rainfall[ActiveCellIndex]);
		s.AppendFormat("\nEVAP: {0:N3} kg", display.Evaporation[ActiveCellIndex]);
		s.AppendFormat("\nSURFACE TEMP: {0:N3}", GetTemperatureString(dependent.SurfaceAirTemperatureAbsolute[ActiveCellIndex], ActiveTemperatureUnits, 1));
		s.AppendFormat("\nTROPOPAUSE ELE: {0:N0}m", dependent.LayerElevation[worldData.AirLayers - 1][ActiveCellIndex]);
		s.AppendFormat("\nDUST: {0:N3} kg", display.DustMass[ActiveCellIndex]);

		if (cloudMass > 0)
		{
			float dropletSize = 1000 * Atmosphere.GetDropletRadius(state.CloudDropletMass[ActiveCellIndex], Atmosphere.GetWaterDensityAtElevation(dependent.DewPoint[ActiveCellIndex], dependent.CloudElevation[ActiveCellIndex]));
			float3 vel = Utils.GetPolarCoordinates(staticState.SphericalPosition[ActiveCellIndex], dependent.CloudVelocity[ActiveCellIndex]);
			s.AppendFormat("\nCLOUD: {0:N3} kg ELE: {1:N0} m R: {2:N3} mm VEL: ({3:N1}, {4:N1}, {5:N1})",
				state.CloudMass[ActiveCellIndex],
				dependent.CloudElevation[ActiveCellIndex],
				state.CloudDropletMass[ActiveCellIndex],
				vel.x, vel.y, vel.z);
		}
		else
		{
			s.AppendFormat("\nCLOUD: 0 kg");
		}
		s.AppendLine();

		for (int i = 1; i < worldData.AirLayers - 1; i++)
		{
			var wind = Utils.GetPolarCoordinates(staticState.SphericalPosition[ActiveCellIndex], state.AirVelocity[i][ActiveCellIndex]);
			s.AppendFormat("\nLAYER {0} | TEMP: {1} RH: {2:P1}",
				i,
				GetTemperatureString(Atmosphere.GetAbsoluteTemperature(state.AirTemperaturePotential[i][ActiveCellIndex], dependent.LayerMiddle[i][ActiveCellIndex]), ActiveTemperatureUnits, 1),
				(dependent.AirHumidityRelative[i][ActiveCellIndex]));
			s.AppendFormat("\nELE: {0:N0} m P: {1:N0} Pa WIND: ({2:N1}, {3:N1}, {4:N1})",
				dependent.LayerElevation[i][ActiveCellIndex],
				display.Pressure[i][ActiveCellIndex],
				wind.x, wind.y, wind.z);
			s.AppendFormat("\nMASS: {0:N0} kg VAPOR: {1:N0} kg CO2: {2:N0} ppm", dependent.AirMass[i][ActiveCellIndex], state.AirVapor[i][ActiveCellIndex], display.CarbonDioxidePercent[i][ActiveCellIndex] * 1000000);
			s.AppendFormat("\nSOLAR ABSORB: {0:P0} REFLECT: {1:P0}", display.AbsorptionSolar[i][ActiveCellIndex].AbsorptivityAirAbove + (1.0f - display.AbsorptionSolar[i][ActiveCellIndex].AbsorptivityAirAbove) * display.AbsorptionSolar[i][ActiveCellIndex].AbsorptivityAirBelow, display.AbsorptionSolar[i][ActiveCellIndex].ReflectivityAirBelow + (1.0f - display.AbsorptionSolar[i][ActiveCellIndex].ReflectivityAirBelow) * display.AbsorptionSolar[i][ActiveCellIndex].ReflectivityAirBelow);
			s.AppendFormat("\nCLOUD ABSORB: {0:P0} REFLECT: {1:P0}", display.AbsorptionSolar[i][ActiveCellIndex].AbsorptivityCloud, display.AbsorptionSolar[i][ActiveCellIndex].ReflectivityCloud);
			s.AppendFormat("\nTHERMAL ABSORB: {0:P0} CLOUD: {1:P0}", display.AbsorptionThermal[i][ActiveCellIndex].AbsorptivityAirAbove + (1.0f - display.AbsorptionThermal[i][ActiveCellIndex].AbsorptivityAirAbove) * display.AbsorptionThermal[i][ActiveCellIndex].AbsorptivityAirBelow, display.AbsorptionThermal[i][ActiveCellIndex].AbsorptivityCloud);
			s.AppendLine();
		}
		return s.ToString();
	}
	public static string GetCellInfoGround(TemperatureUnits ActiveTemperatureUnits, int ActiveCellIndex, ref SimState state, ref TempState dependent)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		s.AppendFormat("ELE: {0:N0} m", state.Elevation[ActiveCellIndex]);
		s.AppendFormat("\nROUGH: {0:N0} m", state.Roughness[ActiveCellIndex]);
		s.AppendFormat("\nTEMP: {0}", GetTemperatureString(state.GroundTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		s.AppendFormat("\nSOIL: {0:N2}", state.GroundCarbon[ActiveCellIndex]);
		s.AppendFormat("\nFLORA: {0:N2} kg TEMP: {1:N2}",
			state.FloraMass[ActiveCellIndex],
			GetTemperatureString(state.FloraTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		s.AppendFormat("\nFWATER: {0:N2} kg GLUCOSE: {1:N2} kg",
			state.FloraWater[ActiveCellIndex],
			state.FloraGlucose[ActiveCellIndex]);
		s.AppendFormat("\nGWATER: {0:N2} TEMP: {1:N2}",
			state.GroundWater[ActiveCellIndex],
			GetTemperatureString(state.GroundWaterTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		s.AppendFormat("\nLAVA: {0:N2} TEMP: {1:N0}",
			state.LavaMass[ActiveCellIndex],
			GetTemperatureString(state.LavaTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		return s.ToString();
	}
	public static string GetCellInfoWater(TemperatureUnits ActiveTemperatureUnits, int ActiveCellIndex, ref WorldData worldData, ref SimState state, ref TempState dependent, ref StaticState staticState, ref DisplayState display)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		if (state.IceMass[ActiveCellIndex] > 0)
		{
			s.AppendFormat("ICE: {0:N3} m TEMP: {1}",
				(state.IceMass[ActiveCellIndex] / WorldData.MassIce),
				GetTemperatureString(state.IceTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		}
		else
		{
			s.AppendFormat("ICE: 0 m");
		}
		float depth = dependent.WaterLayerDepth[1][ActiveCellIndex];
		var nfi = new NumberFormatInfo() { NumberDecimalDigits = (depth >= 1) ? 0 : 3 };
		s.AppendFormat(nfi, "\nDEPTH: {0:N} m", depth);

		if (state.WaterMass[worldData.WaterLayers - 2][ActiveCellIndex] > 0)
		{
			s.AppendFormat(nfi, "\nPLANKTON: {0:N2} kg", state.PlanktonMass[worldData.WaterLayers - 2][ActiveCellIndex]);
		}
		s.AppendLine();


		for (int i = worldData.WaterLayers - 2; i >= 1; i--)
		{
			int layerIndex = (worldData.WaterLayers - 2 - i);
			if (state.WaterMass[i][ActiveCellIndex] > 0)
			{
				var current = Utils.GetPolarCoordinates(staticState.SphericalPosition[ActiveCellIndex], state.WaterVelocity[i][ActiveCellIndex]);
				s.AppendFormat("\nLAYER {0} | TEMP: {1} SALT: {2:P4}",
					layerIndex,
					GetTemperatureString(state.WaterTemperature[i][ActiveCellIndex], ActiveTemperatureUnits, 2),
					display.Salinity[i][ActiveCellIndex]);
				s.AppendFormat("\nVEL: ({0:N3}, {1:N3}, {2:N3})",
					current.x, current.y, current.z);
				s.AppendFormat("\nCO2: {0:N3}",
					state.WaterCarbon[i][ActiveCellIndex]);
				s.AppendFormat("\nP: {0} D: {1}",
					dependent.WaterPressure[i][ActiveCellIndex],
					Atmosphere.GetWaterDensity(display.Salinity[i][ActiveCellIndex], state.WaterTemperature[i][ActiveCellIndex]));
				s.AppendFormat(nfi, "\nDEPTH: {0:N} m HEIGHT: {1:N} m",
					dependent.WaterLayerDepth[i][ActiveCellIndex],
					dependent.WaterLayerHeight[i][ActiveCellIndex]
					);
				s.AppendLine();
			}
		}
		return s.ToString();
	}


	public static float ConvertTemperature(float kelvin, TemperatureUnits units)
	{
		switch (units)
		{
			case TemperatureUnits.Celsius:
				return kelvin - WorldData.FreezingTemperature;
			case TemperatureUnits.Farenheit:
				return (kelvin - WorldData.FreezingTemperature) * 9 / 5 + 32;
			case TemperatureUnits.Kelvin:
			default:
				return kelvin;
		}
	}

	public static string GetTemperatureString(float kelvin, TemperatureUnits units, int decimals)
	{
		string tFormat = "0";
		if (decimals > 0)
		{
			tFormat += ".";
		}
		for (int i = 0; i < decimals; i++)
		{
			tFormat += "0";
		}
		string t = ConvertTemperature(kelvin, units).ToString(tFormat);
		switch (units)
		{
			case TemperatureUnits.Celsius:
				return t + " C";
			case TemperatureUnits.Farenheit:
				return t + " F";
			case TemperatureUnits.Kelvin:
			default:
				return t + " K";
		}
	}

	public static float ConvertTileEnergyToWatts(float energy, float secondsPerTick)
	{
		return energy * 1000 / secondsPerTick;
	}

	public static void PrintState(string title, int i, ref StaticState staticState, ref SimState state, ref WorldData worldData, List<string> degenVarNames)
	{
		StringBuilder s = new StringBuilder();
		s.AppendFormat("{0} Index: {1} Time: {2}", title, i, state.PlanetState.Ticks);
		foreach (var n in degenVarNames)
		{
			s.AppendFormat(" | {0}", n);
		}
		s.AppendLine("");
		s.AppendFormat("X: {0} Y: {1}\n", staticState.Coordinate[i].x, staticState.Coordinate[i].y);
		s.AppendFormat("Elevation: {0}\n", state.Elevation[i]);
		s.AppendFormat("Roughness: {0}\n", state.Roughness[i]);
		s.AppendFormat("SoilFertility: {0}\n", state.GroundCarbon[i]);
		s.AppendFormat("Ground Water: {0} kg\n", state.GroundWater[i]);
		s.AppendFormat("TerrainTemperature: {0}\n", state.GroundTemperature[i]);
		s.AppendFormat("IceMass: {0}\n", state.IceMass[i]);
		s.AppendFormat("IceTemperature: {0}\n", state.IceTemperature[i]);
		for (int j = 0; j < StaticState.GetMaxNeighbors(i, staticState.Neighbors); j++)
		{
			s.AppendFormat("Flow Velocity {0}: {1}\n", j, state.FlowWater[i * StaticState.MaxNeighbors + j]);
		}

		s.AppendFormat("\nLAVA\n");
		s.AppendFormat("Mass: {0}\n", state.LavaMass[i]);
		s.AppendFormat("Temperature: {0}\n", state.LavaTemperature[i]);
		s.AppendFormat("MagmaMass: {0}\n", state.MagmaMass[i]);
		s.AppendFormat("CrustDepth: {0}\n", state.CrustDepth[i]);

		s.AppendFormat("\nFLORA\n");
		s.AppendFormat("Mass: {0}\n", state.FloraMass[i]);
		s.AppendFormat("Glucose: {0}\n", state.FloraGlucose[i]);
		s.AppendFormat("Water: {0}\n", state.FloraWater[i]);
		s.AppendFormat("Temperature: {0}\n", state.FloraTemperature[i]);

		s.AppendFormat("\nPLANKTON\n");
		s.AppendFormat("Mass: {0}\n", state.PlanktonMass[worldData.SurfaceWaterLayer][i]);
		s.AppendFormat("Glucose: {0}\n", state.PlanktonGlucose[worldData.SurfaceWaterLayer][i]);

		s.AppendFormat("\nCLOUD\n");
		s.AppendFormat("CloudMass: {0}\n", state.CloudMass[i]);
		s.AppendFormat("CloudDropletMass: {0}\n", state.CloudDropletMass[i]);

		for (int j = worldData.AirLayers - 2; j >= 1; j--)
		{
			s.AppendFormat("\nAIR LAYER {0}\n", j);
			s.AppendFormat("Temperature: {0}\n", state.AirTemperaturePotential[j][i]);
			s.AppendFormat("Vapor: {0}\n", state.AirVapor[j][i]);
			s.AppendFormat("CarbonDioxide: {0}\n", state.AirCarbon[j][i]);
			s.AppendFormat("Velocity: {0}\n", state.AirVelocity[j][i]);
		}

		for (int j = worldData.WaterLayers - 2; j >= 1; j--)
		{
			s.AppendFormat("\nWATER LAYER {0}\n", j);
			s.AppendFormat("WaterMass: {0}\n", state.WaterMass[j][i]);
			s.AppendFormat("SaltMass: {0}\n", state.SaltMass[j][i]);
			s.AppendFormat("Temperature: {0}\n", state.WaterTemperature[j][i]);
			s.AppendFormat("Carbon: {0}\n", state.WaterCarbon[j][i]);
			s.AppendFormat("Velocity: {0}\n", state.WaterVelocity[j][i]);
		}
		Debug.Log(s);
	}
	public static void PrintDependentState(string title, int i, ref TempState dependent, ref WorldData worldData)
	{
		StringBuilder s = new StringBuilder();
		s.AppendFormat("{0} Index: {1}\n", title, i);
		s.AppendFormat("Surface Elevation: {0}\n", dependent.LayerElevation[worldData.SurfaceAirLayer][i]);
		s.AppendFormat("Water Depth: {0}\n", dependent.WaterLayerDepth[1][i]);
		s.AppendFormat("Ice Coverage: {0}\n", dependent.IceCoverage[i]);
		s.AppendFormat("Flora Coverage: {0}\n", dependent.FloraCoverage[i]);
		s.AppendFormat("Ice Terrain SA: {0}\n", dependent.SurfaceAreaIceTerrain[i]);
		s.AppendFormat("Ice Flora SA: {0}\n", dependent.SurfaceAreaIceFlora[i]);
		s.AppendFormat("Ice Water SA: {0}\n", dependent.SurfaceAreaIceWater[i]);
		s.AppendFormat("Air Ice SA: {0}\n", dependent.SurfaceAreaAirIce[i]);
		s.AppendFormat("Cloud Elevation: {0}\n", dependent.CloudElevation[i]);
		for (int j = 1; j < worldData.WaterLayers - 1; j++)
		{
			s.AppendFormat("Water Coverage: {0}\n", dependent.WaterCoverage[j][i]);
		}
		for (int j = 1; j < worldData.AirLayers - 1; j++)
		{
		}
		Debug.Log(s);
	}


}
