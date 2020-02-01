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
	public NativeArray<float2> VelocityArrow;

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
		VelocityArrow = new NativeArray<float2>(count, Allocator.Persistent);
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
	public NativeArray<float2> VelocityArrow;

	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> CloudCoverage;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> VegetationCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> WaterDepth;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float3> Icosphere;
	[ReadOnly] public NativeArray<float> MeshOverlayData;
	[ReadOnly] public NativeArray<CVP> MeshOverlayColors;
	[ReadOnly] public NativeArray<float2> WindOverlayData;
	[ReadOnly] public bool MeshOverlayActive;
	[ReadOnly] public float MeshOverlayMin;
	[ReadOnly] public float MeshOverlayInverseRange;
	[ReadOnly] public bool WindOverlayActive;
	[ReadOnly] public bool WindMaskedByLand;
	[ReadOnly] public float WindVelocityMax;
	[ReadOnly] public float TerrainScale;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float CloudDropletSizeMin;
	[ReadOnly] public float InverseCloudDropletSizeRange;

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
		float2 velocityArrow;
		float3 surfacePosition;


		float waterDepth = WaterDepth[i];
		float cloudCoverage = CloudCoverage[i];
		float waterCoverage = WaterCoverage[i];
		float iceCoverage = IceCoverage[i];
		float vegetationCoverage = VegetationCoverage[i];
		var icosphere = Icosphere[i];
		var elevation = Terrain[i].Elevation;
		float roughness = Terrain[i].Roughness;

		if (MeshOverlayActive)
		{
			var overlayColor = CVP.Lerp(MeshOverlayColors, (MeshOverlayData[i]- MeshOverlayMin) * MeshOverlayInverseRange);
			terrainColor = overlayColor;
			waterColor = overlayColor;
		}
		else
		{
			terrainColor = GetTerrainColor(roughness, Terrain[i].SoilFertility, waterDepth, iceCoverage, vegetationCoverage);
			waterColor = GetWaterColor(iceCoverage);
		}
		cloudColor = GetCloudColor(math.saturate((CloudDropletMass[i] - CloudDropletSizeMin) * InverseCloudDropletSizeRange), cloudCoverage);

		terrainNormal = icosphere;
		waterNormal = icosphere;
		cloudNormal = icosphere;
		terrainPosition = icosphere * ((elevation + roughness) * TerrainScale + PlanetRadius) / PlanetRadius;
		waterPosition = icosphere * ((elevation + waterDepth) * TerrainScale + PlanetRadius) / PlanetRadius * math.saturate(waterDepth / roughness);
		surfacePosition = icosphere * ((elevation + math.max(roughness, waterDepth)) * TerrainScale + PlanetRadius) / PlanetRadius;
		cloudPosition = icosphere * (CloudElevation[i] * TerrainScale + PlanetRadius) / PlanetRadius;

		if (WindOverlayActive)
		{
			velocityArrow = WindOverlayData[i] / WindVelocityMax;
			if (math.length(velocityArrow) < 0.01f)
			{
				velocityArrow = float2.zero;
			}
			else if (WindMaskedByLand && waterDepth < roughness)
			{
				velocityArrow = float2.zero;
			}
		} else
		{
			velocityArrow = float2.zero;
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

	private Color32 GetTerrainColor(float roughness, float soilFertility, float waterDepth, float iceCoverage, float vegetationCoverage)
	{
		var groundColor = Color32.Lerp(new Color32(50, 50, 80, 255), new Color32(100, 60, 20, 255), soilFertility);
		var waterColor = Color32.Lerp(groundColor, new Color32(0, 0, 255, 255), math.saturate(math.pow(waterDepth / roughness, 2)));
		var iceColor = Color32.Lerp(waterColor, new Color32(255, 255, 255, 255), iceCoverage);
		var vegetationColor = Color32.Lerp(groundColor, new Color32(0, 220, 30, 255), vegetationCoverage);
		return vegetationColor;
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

[BurstCompile]
struct LerpJobfloat : IJobParallelFor {
	public NativeArray<float> Out;
	[ReadOnly] public NativeArray<float> Start;
	[ReadOnly] public NativeArray<float> End;
	[ReadOnly] public float Progress;
	public void Execute(int i)
	{
		Out[i] = math.lerp(Start[i], End[i], Progress);
	}
	public LerpJobfloat(NativeArray<float> a, NativeArray<float> b, float p)
	{
		Out = new NativeArray<float>(a.Length, Allocator.TempJob);
		Start = new NativeArray<float>(a, Allocator.TempJob);
		End = new NativeArray<float>(b, Allocator.TempJob);
		Progress = p;
	}
	public void CopyAndDispose(ref NativeArray<float> c)
	{
		Out.CopyTo(c);
		Out.Dispose();
		Start.Dispose();
		End.Dispose();
	}
}



