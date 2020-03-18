using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity;
using Unity.Jobs;
using Unity.Burst;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;


public struct RenderStateCell {
	public Color32 TerrainColor;
	public Color32 WaterColor;
	public Color32 CloudColor;
	public float3 SurfacePosition;
	public float3 TerrainPosition;
	public float3 WaterPosition;
	public float3 CloudPosition;
	public float3 TerrainNormal;
	public float3 WaterNormal;
	public float3 CloudNormal;
	public float2 VelocityArrow;
}
public struct RenderState {

	public float Ticks;
	public float3 Position;
	public float3 Rotation;

	public NativeArray<Color32> TerrainColor;
	public NativeArray<Color32> WaterColor;
	public NativeArray<Color32> CloudColor;
	public NativeArray<Vector3> TerrainPosition;
	public NativeArray<Vector3> WaterPosition;
	public NativeArray<Vector3> CloudPosition;
	public NativeArray<Vector3> TerrainNormal;
	public NativeArray<Vector3> WaterNormal;
	public NativeArray<Vector3> CloudNormal;
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> VelocityArrow;

	public void Init(int count)
	{
		TerrainColor = new NativeArray<Color32>(count, Allocator.Persistent);
		TerrainNormal = new NativeArray<Vector3>(count, Allocator.Persistent);
		TerrainPosition = new NativeArray<Vector3>(count, Allocator.Persistent);
		WaterColor = new NativeArray<Color32>(count, Allocator.Persistent);
		WaterNormal = new NativeArray<Vector3>(count, Allocator.Persistent);
		WaterPosition = new NativeArray<Vector3>(count, Allocator.Persistent);
		CloudColor = new NativeArray<Color32>(count, Allocator.Persistent);
		CloudNormal = new NativeArray<Vector3>(count, Allocator.Persistent);
		CloudPosition = new NativeArray<Vector3>(count, Allocator.Persistent);
		VelocityArrow = new NativeArray<float3>(count, Allocator.Persistent);
		SurfacePosition = new NativeArray<float3>(count, Allocator.Persistent);
	}

	public void Dispose()
	{
		TerrainColor.Dispose();
		TerrainNormal.Dispose();
		TerrainPosition.Dispose();
		WaterColor.Dispose();
		WaterNormal.Dispose();
		WaterPosition.Dispose();
		CloudColor.Dispose();
		CloudNormal.Dispose();
		CloudPosition.Dispose();
		VelocityArrow.Dispose();
		SurfacePosition.Dispose();
	}
}

[BurstCompile]
public struct BuildRenderStateJob : IJobParallelFor {
	public NativeArray<Color32> TerrainColor;
	public NativeArray<Color32> WaterColor;
	public NativeArray<Color32> CloudColor;
	public NativeArray<Vector3> TerrainPosition;
	public NativeArray<Vector3> WaterPosition;
	public NativeArray<Vector3> CloudPosition;
	public NativeArray<Vector3> TerrainNormal;
	public NativeArray<Vector3> WaterNormal;
	public NativeArray<Vector3> CloudNormal;
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> VelocityArrow;

	[ReadOnly] public NativeArray<float> Roughness;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public NativeArray<float> Elevation;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> GroundWater;
	[ReadOnly] public NativeArray<float> FloraCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> WaterDepth;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float3> Icosphere;
	[ReadOnly] public NativeArray<float> MeshOverlayData;
	[ReadOnly] public NativeArray<CVP> MeshOverlayColors;
	[ReadOnly] public NativeArray<float3> WindOverlayData;
	[ReadOnly] public bool MeshOverlayActive;
	[ReadOnly] public float MeshOverlayMin;
	[ReadOnly] public float MeshOverlayInverseRange;
	[ReadOnly] public bool WindOverlayActive;
	[ReadOnly] public bool WindMaskedByLand;
	[ReadOnly] public float WindVelocityMax;
	[ReadOnly] public float TerrainScale;
	[ReadOnly] public float AtmosphereScale;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float CloudDropletSizeMin;
	[ReadOnly] public float InverseCloudDropletSizeRange;
	[ReadOnly] public float InverseCloudMass;
	[ReadOnly] public float GroundWaterMax;

	public void Execute(int i)
	{
		Color32 terrainColor;
		Color32 waterColor;
		Color32 cloudColor;
		float3 terrainNormal;
		float3 waterNormal;
		float3 cloudNormal;
		float3 terrainPosition;
		float3 waterPosition;
		float3 cloudPosition;
		float3 velocityArrow;
		float3 surfacePosition;


		float waterDepth = WaterDepth[i];
		float cloudCoverage = CloudMass[i]* InverseCloudMass;
		float waterCoverage = WaterCoverage[i];
		float iceCoverage = IceCoverage[i];
		float floraCoverage = FloraCoverage[i];
		var icosphere = Icosphere[i];
		var elevation = Elevation[i];
		float roughness = Roughness[i];
		float surfaceElevation = elevation + math.max(roughness, waterDepth);

		if (MeshOverlayActive)
		{
			var overlayColor = CVP.Lerp(MeshOverlayColors, (MeshOverlayData[i]- MeshOverlayMin) * MeshOverlayInverseRange);
			terrainColor = overlayColor;
			waterColor = overlayColor;
		}
		else
		{
			terrainColor = GetTerrainColor(roughness, SoilFertility[i], waterDepth, iceCoverage, floraCoverage, GroundWater[i] / GroundWaterMax);
			waterColor = GetWaterColor(iceCoverage);
		}
		cloudColor = GetCloudColor(math.saturate((CloudDropletMass[i] - CloudDropletSizeMin) * InverseCloudDropletSizeRange), cloudCoverage);

		terrainNormal = icosphere;
		waterNormal = icosphere;
		cloudNormal = icosphere;
		terrainPosition = icosphere * ((elevation + roughness) * TerrainScale + PlanetRadius) / PlanetRadius;
		waterPosition = icosphere * ((elevation + waterDepth) * TerrainScale + PlanetRadius) / PlanetRadius * math.saturate(waterDepth / roughness);
		surfacePosition = icosphere * (surfaceElevation * TerrainScale + PlanetRadius) / PlanetRadius;
		cloudPosition = icosphere * (((math.max(0, CloudElevation[i] - surfaceElevation)) * AtmosphereScale) + ((surfaceElevation + 100) * TerrainScale + PlanetRadius)) / PlanetRadius;

		if (WindOverlayActive)
		{
			velocityArrow = WindOverlayData[i] / WindVelocityMax;
			if (math.length(velocityArrow) < 0.01f)
			{
				velocityArrow = float3.zero;
			}
			else if (WindMaskedByLand && waterDepth < roughness)
			{
				velocityArrow = float3.zero;
			}
		} else
		{
			velocityArrow = float3.zero;
		}

		TerrainColor[i] = terrainColor;
		TerrainNormal[i] = terrainNormal;
		TerrainPosition[i] = terrainPosition;
		WaterColor[i] = waterColor;
		WaterNormal[i] = waterNormal;
		WaterPosition[i] = waterPosition;
		CloudColor[i] = cloudColor;
		CloudNormal[i] = cloudNormal;
		CloudPosition[i] = cloudPosition;
		VelocityArrow[i] = velocityArrow;
		SurfacePosition[i] = surfacePosition;
	}

	private Color32 GetTerrainColor(float roughness, float soilFertility, float waterDepth, float iceCoverage, float floraCoverage, float groundWaterSaturation)
	{
		if (iceCoverage > 0.5f)
		{
			return new Color32(255, 255, 255, 255);
		//} else if (waterDepth / roughness > 0.5f)
		//{
		//	return new Color32(0, 0, 255, 255);
		}
		else if (groundWaterSaturation >= 1)
		{
			return new Color32(0, 0, 255, 255);
		}
		else if (floraCoverage > 0.3f)
		{
			return new Color32(0, 60, 10, 255);
		}
		else if (floraCoverage > 0.15f)
		{
			return new Color32(40, 100, 30, 255);
		}
		else if (soilFertility > 0.5f)
		{
			if (groundWaterSaturation > 0.5f)
			{
				return new Color32(40, 20, 20, 255);
			}
			return new Color32(50, 30, 20, 255);
		}
		else if (soilFertility > 0.25f)
		{
			if (groundWaterSaturation > 0.5f)
			{
				return new Color32(70, 60, 30, 255);
			}
			return new Color32(80, 70, 30, 255);
		}
		else
		{
			if (groundWaterSaturation > 0.5f)
			{
				return new Color32(120, 120, 60, 255);
			}
			return new Color32(130, 130, 60, 255);
		}
		//var groundColor = Color32.Lerp(new Color32(50, 50, 80, 255), new Color32(100, 60, 20, 255), soilFertility);
		//var waterColor = Color32.Lerp(groundColor, new Color32(0, 0, 255, 255), math.saturate(math.pow(waterDepth / roughness, 2)));
		//var iceColor = Color32.Lerp(waterColor, new Color32(255, 255, 255, 255), iceCoverage);
		//var floraColor = Color32.Lerp(groundColor, new Color32(0, 220, 30, 255), floraCoverage);
		//return floraColor;
	}

	private Color32 GetWaterColor(float iceCoverage)
	{
		var waterColor = new Color32(0, 0, 255, 255);
		var iceColor = Color32.Lerp(waterColor, new Color32(255, 255, 255, 255), iceCoverage);
		return iceColor;
	}

	private Color32 GetCloudColor(float dropletSize, float cloudCoverage)
	{
		var c = Color32.Lerp(new Color32(255, 255, 255, 255), new Color32(0, 0, 0, 255), dropletSize);
		float opacity = -math.cos(cloudCoverage * math.PI) / 2 + 0.5f;
		c.a = (byte)(255 * opacity);
		return c;
	}


}

[BurstCompile]
struct LerpJobVector3 : IJobParallelFor {
	public NativeArray<Vector3> Out;
	[ReadOnly] public NativeArray<Vector3> Start;
	[ReadOnly] public NativeArray<Vector3> End;
	[ReadOnly] public float Progress;
	public void Execute(int i)
	{
		Out[i] = Vector3.Lerp(Start[i], End[i], Progress);
	}
}

[BurstCompile]
struct LerpJobColor32 : IJobParallelFor {
	public NativeArray<Color32> Out;
	[ReadOnly] public NativeArray<Color32> Start;
	[ReadOnly] public NativeArray<Color32> End;
	[ReadOnly] public float Progress;
	public void Execute(int i)
	{
		Out[i] = Color32.Lerp(Start[i], End[i], Progress);
	}
}


[BurstCompile]
struct LerpJobfloat3 : IJobParallelFor {
	public NativeArray<float3> Out;
	[ReadOnly] public NativeArray<float3> Start;
	[ReadOnly] public NativeArray<float3> End;
	[ReadOnly] public float Progress;
	public void Execute(int i)
	{
		Out[i] = math.lerp(Start[i], End[i], Progress);
	}
}

[BurstCompile]
struct LerpJobfloat2 : IJobParallelFor {
	public NativeArray<float2> Out;
	[ReadOnly] public NativeArray<float2> Start;
	[ReadOnly] public NativeArray<float2> End;
	[ReadOnly] public float Progress;
	public void Execute(int i)
	{
		Out[i] = math.lerp(Start[i], End[i], Progress);
	}
}

