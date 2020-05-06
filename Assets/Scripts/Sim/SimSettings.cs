using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

[Serializable]
public struct SimSettings {
	[Header("Debug")]
	public bool CheckForDegeneracy;
	public bool LogState;
	public int LogStateIndex;
	public SynchronousOverride SynchronousOverrides;

	[Header("Sim")]
	public bool MakeAirIncompressible;
	public bool MakeWaterIncompressible;
	public bool WaterSurfaceFlowEnabled;
	public bool RebalanceWaterLayers;
	public bool AdvectionAir;
	public bool DiffusionAir;
	public bool AdvectionWater;
	public bool DiffusionWater;
	public bool AdvectionCloud;
	public bool DiffusionCloud;
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

	[HideInInspector] public bool CollectGlobals;
	[HideInInspector] public bool CollectOverlay;

	[Serializable]
	public struct SynchronousOverride {
		public bool FluxDust;
		public bool FluxFreeze;
		public bool FluxIceMelt;
		public bool FluxCondensation;
		public bool FluxEvaporation;
		public bool FluxPlankton;
		public bool FluxCloud;
		public bool FluxFlora;
		public bool FluxLava;
		public bool FluxTerrain;
	}

}


