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

public struct RenderState {

	public float Ticks;
	public float3 Position;
	public float3 Rotation;
	public NativeArray<Color32> TerrainColor;
	public NativeArray<Color32> WaterColor;
	public NativeArray<Color32> CloudColor;
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> TerrainPosition;
	public NativeArray<float3> WaterPosition;
	public NativeArray<float3> CloudPosition;
	public NativeArray<float3> TerrainNormal;
	public NativeArray<float3> WaterNormal;
	public NativeArray<float3> CloudNormal;
	public NativeArray<float2> VelocityArrow;

	public void Init(int count)
	{
		TerrainColor = new NativeArray<Color32>(count, Allocator.Persistent);
		WaterColor = new NativeArray<Color32>(count, Allocator.Persistent);
		CloudColor = new NativeArray<Color32>(count, Allocator.Persistent);
		SurfacePosition = new NativeArray<float3>(count, Allocator.Persistent);
		TerrainPosition = new NativeArray<float3>(count, Allocator.Persistent);
		WaterPosition = new NativeArray<float3>(count, Allocator.Persistent);
		CloudPosition = new NativeArray<float3>(count, Allocator.Persistent);
		TerrainNormal = new NativeArray<float3>(count, Allocator.Persistent);
		WaterNormal = new NativeArray<float3>(count, Allocator.Persistent);
		CloudNormal = new NativeArray<float3>(count, Allocator.Persistent);
		VelocityArrow = new NativeArray<float2>(count, Allocator.Persistent);
	}

	public void Dispose()
	{
		TerrainColor.Dispose();
		WaterColor.Dispose();
		CloudColor.Dispose();
		SurfacePosition.Dispose();
		TerrainPosition.Dispose();
		WaterPosition.Dispose();
		CloudPosition.Dispose();
		TerrainNormal.Dispose();
		WaterNormal.Dispose();
		CloudNormal.Dispose();
		VelocityArrow.Dispose();
	}
}

[BurstCompile]
public struct BuildRenderStateJob : IJobParallelFor {
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> TerrainPosition;
	public NativeArray<float3> TerrainNormal;
	public NativeArray<Color32> TerrainColor;
	public NativeArray<float3> WaterPosition;
	public NativeArray<float3> WaterNormal;
	public NativeArray<Color32> WaterColor;
	public NativeArray<float3> CloudPosition;
	public NativeArray<float3> CloudNormal;
	public NativeArray<Color32> CloudColor;
	public NativeArray<float2> VelocityArrow;

	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float3> Icosphere;
	[ReadOnly] public NativeArray<float2> Velocity;
	[ReadOnly] public float MaxVelocity;
	[ReadOnly] public float TerrainScale;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public bool UseOverlay;
	[ReadOnly] public bool VelocityMaskedByLand;
	[ReadOnly] public WorldView.Overlay Overlay;

	public void Execute(int i)
	{
		float cloudCoverage =;
		float waterDepth =;
		float iceCoverage =;
		float vegetationCoverage =;
		float waterAndIceDepth =;
		float relativeHumidity =;
		var icosphere = Icosphere[i];
		var elevation = Terrain[i].Elevation;
		float roughness = Terrain[i].Roughness;

		if (UseOverlay)
		{
			var overlayColor = CVP.Lerp(Overlay.ColorValuePairs, (Overlay.Values[i]-Overlay.Min)*Overlay.InverseRange);
			TerrainColor[i] = overlayColor;
			WaterColor[i] = overlayColor;
		}
		else
		{
			TerrainColor[i] = WorldView.GetTerrainColor(roughness, Terrain[i].SoilFertility, waterDepth, iceCoverage, vegetationCoverage);
			WaterColor[i] = WorldView.GetWaterColor(iceCoverage);
		}
		CloudColor[i] = WorldView.GetCloudColor(relativeHumidity, cloudCoverage);

		TerrainNormal[i] = icosphere;
		WaterNormal[i] = icosphere;
		CloudNormal[i] = icosphere;
		TerrainPosition[i] = icosphere * ((elevation + roughness) * TerrainScale + PlanetRadius) / PlanetRadius;
		WaterPosition[i] = icosphere * ((elevation + waterDepth) * TerrainScale + PlanetRadius) / PlanetRadius * math.saturate(waterDepth / roughness);
		SurfacePosition[i] = icosphere * ((elevation + math.max(roughness, waterAndIceDepth)) * TerrainScale + PlanetRadius) / PlanetRadius;
		CloudPosition[i] = icosphere * (CloudElevation[i] * TerrainScale + PlanetRadius) / PlanetRadius;

		float2 velocity = Velocity[i] / MaxVelocity;
		if (math.length(velocity) < 0.01f)
		{
			velocity = float2.zero;
		}
		else if (VelocityMaskedByLand && waterDepth < roughness)
		{
			velocity = float2.zero;
		}
		VelocityArrow[i] = velocity;		
	}

}

[BurstCompile]
public struct LerpColor32Job : IJobParallelFor {
	public NativeArray<Color32> Out;
	[ReadOnly] public NativeArray<Color32> Start;
	[ReadOnly] public NativeArray<Color32> End;
	[ReadOnly] public float T;
	public void Execute(int i)
	{
		Out[i] = Color32.Lerp(Start[i], End[i], T);
	}
}
[BurstCompile]
public struct Lerpfloat3Job : IJobParallelFor {
	public NativeArray<float3> Out;
	[ReadOnly] public NativeArray<float3> Start;
	[ReadOnly] public NativeArray<float3> End;
	[ReadOnly] public float T;
	public void Execute(int i)
	{
		Out[i] = math.lerp(Start[i], End[i], T);
	}
}
[BurstCompile]
public struct Lerpfloat2Job : IJobParallelFor {
	public NativeArray<float2> Out;
	[ReadOnly] public NativeArray<float2> Start;
	[ReadOnly] public NativeArray<float2> End;
	[ReadOnly] public float T;
	public void Execute(int i)
	{
		Out[i] = math.lerp(Start[i], End[i], T);
	}
}
[BurstCompile]
public struct LerpfloatJob : IJobParallelFor {
	public NativeArray<float> Out;
	[ReadOnly] public NativeArray<float> Start;
	[ReadOnly] public NativeArray<float> End;
	[ReadOnly] public float T;
	public void Execute(int i)
	{
		Out[i] = math.lerp(Start[i], End[i], T);
	}
}


