using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Entities;
using Unity;
using UnityEngine;
using Unity.Mathematics;

public struct SimState {

	public int Count;
	public int Ticks;
	public float Gravity;
	public float SpinAngle;
	public float SpinSpeed;
	public float OrbitSpeed;
	public float TiltAngle;
	public float GeothermalHeat;
	public float SolarRadiation;
	public float StratosphereMass;
	public float CarbonDioxide;

	public SimStateCell[] Cells;

	public void Init(int count)
	{
		Count = count;
		Cells = new SimStateCell[Count];
	}
}

public struct SimStateCell {
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
	public float RelativeHumidity;
	public float AirTemperature;
	public float AirPressure;
	public float WaterTemperature;
	public float WindVertical;
	public float2 WindSurface;
	public float2 WindTropopause;
}


public struct RenderState {

	public float Ticks;
	public float SpinAngle;
	public float SpinSpeed;
	public float OrbitSpeed;
	public float TiltAngle;
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
