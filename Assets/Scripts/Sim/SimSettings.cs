using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

[Serializable]
public struct SimSettings {
	public bool MakeAirIncompressible;
	public bool MakeWaterIncompressible;
	public bool WaterSurfaceFlowEnabled;
	public bool Condensation;
	public bool Evaporation;
	public bool Freezing;
	public bool Plankton;
	public bool Precipitation;
	public bool Flora;
	public bool IceMelting;
	public bool SoilRespiration;
	public bool GroundWater;
	public bool ConductionAirIce;
	public bool ConductionAirWater;
	public bool ConductionAirFlora;
	public bool ConductionAirTerrain;
	public bool ConductionIceWater;
	public bool ConductionIceFlora;
	public bool ConductionIceTerrain;
	public bool ConductionFloraTerrain;
	public bool ConductionWaterTerrain;
	public int IncompressibilityIterations;
	public bool CheckForDegeneracy;
	public bool CollectGlobals;
	public bool CollectOverlay;
	public bool LogState;
	public int LogStateIndex;
}
