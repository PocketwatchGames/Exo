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
	public NativeArray<float> TerrainElevation;
	public NativeArray<Color32> WaterColor;
	public NativeArray<float> WaterElevation;
	public NativeArray<Color32> CloudColor;
	public NativeArray<float> CloudElevation;
	public NativeArray<Color32> LavaColor;
	public NativeArray<float> LavaElevation;
	public NativeArray<Color32> DustColor;
	public NativeArray<float> DustElevation;
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> VelocityArrow;

	public void Init(int count)
	{
		TerrainColor = new NativeArray<Color32>(count, Allocator.Persistent);
		TerrainElevation = new NativeArray<float>(count, Allocator.Persistent);
		WaterColor = new NativeArray<Color32>(count, Allocator.Persistent);
		WaterElevation = new NativeArray<float>(count, Allocator.Persistent);
		CloudColor = new NativeArray<Color32>(count, Allocator.Persistent);
		CloudElevation = new NativeArray<float>(count, Allocator.Persistent);
		LavaColor = new NativeArray<Color32>(count, Allocator.Persistent);
		LavaElevation = new NativeArray<float>(count, Allocator.Persistent);
		DustColor = new NativeArray<Color32>(count, Allocator.Persistent);
		DustElevation = new NativeArray<float>(count, Allocator.Persistent);
		VelocityArrow = new NativeArray<float3>(count, Allocator.Persistent);
		SurfacePosition = new NativeArray<float3>(count, Allocator.Persistent);
	}

	public void Dispose()
	{
		TerrainColor.Dispose();
		TerrainElevation.Dispose();
		WaterColor.Dispose();
		WaterElevation.Dispose();
		CloudColor.Dispose();
		CloudElevation.Dispose();
		LavaColor.Dispose();
		LavaElevation.Dispose();
		DustColor.Dispose();
		DustElevation.Dispose();
		VelocityArrow.Dispose();
		SurfacePosition.Dispose();
	}
}

[BurstCompile]
public struct BuildRenderStateCellJob : IJobParallelFor {
	public NativeArray<Color32> TerrainColor;
	public NativeArray<float> TerrainElevation;
	public NativeArray<Color32> WaterColor;
	public NativeArray<float> WaterElevation;
	public NativeArray<float> CloudElevation;
	public NativeArray<Color32> CloudColor;
	public NativeArray<Color32> LavaColor;
	public NativeArray<float> LavaElevation;
	public NativeArray<Color32> DustColor;
	public NativeArray<float> DustElevation;
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> VelocityArrow;

	[ReadOnly] public NativeArray<float> Roughness;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public NativeArray<float> Elevation;
	[ReadOnly] public NativeArray<float> CloudElevationSim;
	[ReadOnly] public NativeArray<float> CloudDropletMass;
	[ReadOnly] public NativeArray<float> CloudAbsorption;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> GroundWater;
	[ReadOnly] public NativeArray<float> FloraCoverage;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> WaterDepth;
	[ReadOnly] public NativeArray<float> WaterTemperature;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float3> Icosphere;
	[ReadOnly] public NativeArray<float> MeshOverlayData;
	[ReadOnly] public NativeArray<CVP> MeshOverlayColors;
	[ReadOnly] public NativeArray<float3> WindOverlayData;
	[ReadOnly] public NativeArray<float> LavaMass;
	[ReadOnly] public NativeArray<float> LavaTemperature;
	[ReadOnly] public NativeArray<float> DustCoverage;
	[ReadOnly] public NativeArray<float> PlanktonMass;
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
	[ReadOnly] public float GroundWaterMax;
	[ReadOnly] public float SoilFertilityMax;
	[ReadOnly] public float DustHeight;
	[ReadOnly] public float LavaCrystalizationTemperature;
	[ReadOnly] public float LavaTemperatureRangeInverse;
	[ReadOnly] public float DustMaxInverse;
	[ReadOnly] public float LavaDensityAdjustment;
	[ReadOnly] public float DisplaySoilWeight;
	[ReadOnly] public float DisplaySandWeight;
	[ReadOnly] public float DisplayFloraWeight;

	public void Execute(int i)
	{
		Color32 terrainColor;
		float terrainElevation;
		Color32 waterColor;
		float waterElevation;
		Color32 cloudColor;
		float cloudElevation;
		Color32 lavaColor;
		float lavaElevation;
		Color32 dustColor;
		float dustElevation;
		float3 velocityArrow;
		float3 surfacePosition;


		float waterDepth = WaterDepth[i];
		float cloudCoverage = Utils.Sqr(CloudAbsorption[i]) * 0.25f;
		float waterCoverage = WaterCoverage[i];
		float iceCoverage = IceCoverage[i];
		float floraCoverage = FloraCoverage[i];
		var icosphere = Icosphere[i];
		var elevation = Elevation[i];
		float roughness = math.max(1, Roughness[i]);
		float surfaceElevation = elevation + math.max(roughness, waterDepth);

		if (MeshOverlayActive)
		{
			var overlayColor = CVP.Lerp(MeshOverlayColors, (MeshOverlayData[i] - MeshOverlayMin) * MeshOverlayInverseRange);
			terrainColor = overlayColor;
			waterColor = overlayColor;
			lavaColor = overlayColor;
		}
		else
		{
			// Terrain color
			float fertility = math.saturate(SoilFertility[i] / SoilFertilityMax);
			var t = math.normalize(
				new float4(
				0,
				(1 - fertility) * DisplaySandWeight,
				fertility * DisplaySoilWeight,
				floraCoverage * DisplayFloraWeight));
			t *= (1 - (waterDepth < 1 ? iceCoverage : 0));
			t *= 255;
			terrainColor = new Color32((byte)t.x, (byte)t.y, (byte)t.z, (byte)t.w);

			waterColor = GetWaterColor(iceCoverage, WaterTemperature[i], waterDepth, PlanktonMass[i]);
			lavaColor = GetLavaColor(LavaTemperature[i], LavaCrystalizationTemperature, LavaTemperatureRangeInverse);
		}
		cloudColor = GetCloudColor(cloudCoverage);
		dustColor = GetDustColor(DustCoverage[i], DustMaxInverse);

		terrainElevation = ((elevation + roughness) * TerrainScale + PlanetRadius) / PlanetRadius;
		waterElevation = ((waterDepth == 0) ? 0.99f : ((elevation + waterDepth) * TerrainScale + PlanetRadius) / PlanetRadius);
		surfacePosition = icosphere * (surfaceElevation * TerrainScale + PlanetRadius) / PlanetRadius;
		cloudElevation = (((math.max(0, CloudElevationSim[i] - surfaceElevation)) * AtmosphereScale) + ((surfaceElevation + 100) * TerrainScale + PlanetRadius)) / PlanetRadius;


		dustElevation = (((math.max(0, DustHeight - surfaceElevation)) * AtmosphereScale) + ((surfaceElevation + 100) * TerrainScale + PlanetRadius)) / PlanetRadius;

		float lavaDepth = LavaMass[i] / (WorldData.MassLava * LavaDensityAdjustment);
		lavaElevation = ((lavaDepth == 0) ? 0.75f : ((elevation + roughness + lavaDepth) * TerrainScale + PlanetRadius) / PlanetRadius);

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
		}
		else
		{
			velocityArrow = float3.zero;
		}

		TerrainColor[i] = terrainColor;
		TerrainElevation[i] = terrainElevation;
		WaterColor[i] = waterColor;
		WaterElevation[i] = waterElevation;
		CloudColor[i] = cloudColor;
		CloudElevation[i] = cloudElevation;
		LavaColor[i] = lavaColor;
		LavaElevation[i] = lavaElevation;
		DustColor[i] = dustColor;
		DustElevation[i] = dustElevation;
		VelocityArrow[i] = velocityArrow;
		SurfacePosition[i] = surfacePosition;
	}


	private Color32 GetWaterColor(float iceCoverage, float waterTemperature, float depth, float plankton)
	{
		int blue = 150;
		int red = 0;
		int green = 0;
		int alpha = 200;
		if (plankton > 1)
		{
			green += 80;
		} else if (plankton > 0.01f)
		{
			green += 40;
		} else
		{
			green += 0;
		}
		if (iceCoverage > 0.95f)
		{
			red += 220;
			green += 220;
			blue += 220;
			alpha = 240;
		}
		else if (iceCoverage > 0.1f)
		{
			red += 40;
			green += 40;
			blue += 40;
			alpha = 220;
		}
		return new Color32((byte)math.min(255, red), (byte)math.min(255, green), (byte)math.min(255, blue), (byte)alpha);
	}					   
	private Color32 GetCloudColor(float cloudCoverage)
	{
		return new Color32(255, 255, 255, (byte)(255 * cloudCoverage));
	}

	private Color32 GetDustColor(float dustCoverage, float dustMaxInverse)
	{
		return Color32.Lerp(new Color32(0, 0, 0, 0), new Color32(0, 0, 0, 255), dustCoverage*dustMaxInverse);
	}

	private Color32 GetLavaColor(float temperature, float minTemp, float inverseTempRange)
	{
		return new Color32((byte)(255 * (temperature - minTemp) * inverseTempRange), 0, 0, 255);
	}

}



[BurstCompile]
public struct BuildHexVertsJob : IJobParallelFor {
	public NativeArray<Vector3> VTerrainPosition;
	public NativeArray<Color32> VTerrainColor;
	public NativeArray<Vector3> VWaterPosition;
	public NativeArray<Vector3> VWaterNormal;
	public NativeArray<Color32> VWaterColor;
	public NativeArray<Vector3> VCloudPosition;
	public NativeArray<Color32> VCloudColor;
	public NativeArray<Vector3> VLavaPosition;
	public NativeArray<Color32> VLavaColor;
	public NativeArray<Vector3> VDustPosition;
	public NativeArray<Color32> VDustColor;

	[ReadOnly] public NativeArray<float> TerrainElevation;
	[ReadOnly] public NativeArray<float> WaterElevation;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> LavaElevation;
	[ReadOnly] public NativeArray<float> DustElevation;
	[ReadOnly] public NativeArray<Color32> TerrainColor;
	[ReadOnly] public NativeArray<Color32> WaterColor;
	[ReadOnly] public NativeArray<Color32> CloudColor;
	[ReadOnly] public NativeArray<Color32> LavaColor;
	[ReadOnly] public NativeArray<Color32> DustColor;
	[ReadOnly] public NativeArray<float3> HexVerts;
	[ReadOnly] public NativeArray<float3> IcosphereVerts;

	public void Execute(int i)
	{
		float3 v = HexVerts[i];
		int j = (int)(i / WorldView.VertsPerCell);

		VTerrainPosition[i] = v * TerrainElevation[j];
		VWaterPosition[i] = v * WaterElevation[j];
		VCloudPosition[i] = v * CloudElevation[j];
		VLavaPosition[i] = v * LavaElevation[j];
		VDustPosition[i] = v * DustElevation[j];
		VTerrainColor[i] = TerrainColor[j];
		VWaterColor[i] = WaterColor[j];
		VCloudColor[i] = CloudColor[j];
		VLavaColor[i] = LavaColor[j];
		VDustColor[i] = DustColor[j];

		VWaterNormal[i] = v;
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
}

