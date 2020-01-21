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

	public SimPlanetState PlanetState;
	public SimWind[] Wind;
	public SimCell[] Cells;
	public NativeArray<DisplayCell> DisplayCells;

	public void Init(int count)
	{
		Cells = new SimCell[count];
		Wind = new SimWind[count];
		DisplayCells = new NativeArray<DisplayCell>(count, Allocator.Persistent);
	}

	public void Dispose()
	{
		DisplayCells.Dispose();
	}
}

public struct SimPlanetState {
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

public struct SimCell {
	public float Elevation;
	public float Roughness;
	public float SoilFertility;
	public float Vegetation;
	public float IceMass;
	public float WaterMass;
	public float WaterEnergy;
	public float SaltMass;
	public float GroundEnergy;
	public float GroundWater;
	public float GroundWaterDepth;
	public float AirMass;
	public float AirEnergy;
	public float AirWaterMass;
	public float CloudMass;
	public float CloudDropletMass;
	public float CloudCoverage;
	public float CloudElevation;
	public float WaterDepth;
	public float WaterAndIceDepth;
	public float WaterDensity;
	public float RelativeHumidity;
	public float AirTemperature;
	public float AirPressure;
	public float WaterTemperature;
}

public struct SimWind {
	public float WindVertical;
	public float2 WindSurface;
	public float2 WindTropopause;
	public float2 CurrentSurface;
	public float2 CurrentDeep;
	public float CurrentVertical;
}

public struct DisplayCell {
	public float Heat;
	public float Rainfall;
	public float Evaporation;
}


public struct RenderState {

	public float Ticks;
	public Vector3 Position;
	public Vector3 Rotation;
	public Color32[] TerrainColor;
	public Color32[] WaterColor;
	public Color32[] CloudColor;
	public Vector3[] TerrainPosition;
	public Vector3[] WaterPosition;
	public Vector3[] CloudPosition;
	public Vector3[] TerrainNormal;
	public Vector3[] WaterNormal;
	public Vector3[] CloudNormal;

	public void Init(int count)
	{
		TerrainColor = new Color32[count];
		WaterColor = new Color32[count];
		CloudColor = new Color32[count];
		TerrainPosition = new Vector3[count];
		WaterPosition = new Vector3[count];
		CloudPosition = new Vector3[count];
		TerrainNormal = new Vector3[count];
		WaterNormal = new Vector3[count];
		CloudNormal = new Vector3[count];
	}
}
