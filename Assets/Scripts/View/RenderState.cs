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

	public NativeArray<Vector4> TerrainColor1;
	public NativeArray<Vector4> TerrainColor2;
	public NativeArray<float> TerrainElevation;
	public NativeArray<Color32> WaterColor;
	public NativeArray<float> WaterElevation;
	public NativeArray<Color32> CloudColor;
	public NativeArray<Color32> OverlayColor;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> CloudHeight;
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> VelocityArrow;

	public void Init(int count)
	{
		TerrainColor1 = new NativeArray<Vector4>(count, Allocator.Persistent);
		TerrainColor2 = new NativeArray<Vector4>(count, Allocator.Persistent);
		TerrainElevation = new NativeArray<float>(count, Allocator.Persistent);
		WaterColor = new NativeArray<Color32>(count, Allocator.Persistent);
		WaterElevation = new NativeArray<float>(count, Allocator.Persistent);
		OverlayColor = new NativeArray<Color32>(count, Allocator.Persistent);
		CloudColor = new NativeArray<Color32>(count, Allocator.Persistent);
		CloudElevation = new NativeArray<float>(count, Allocator.Persistent);
		CloudHeight = new NativeArray<float>(count, Allocator.Persistent);
		VelocityArrow = new NativeArray<float3>(count, Allocator.Persistent);
		SurfacePosition = new NativeArray<float3>(count, Allocator.Persistent);
	}

	public void Dispose()
	{
		TerrainColor1.Dispose();
		TerrainColor2.Dispose();
		TerrainElevation.Dispose();
		WaterColor.Dispose();
		WaterElevation.Dispose();
		OverlayColor.Dispose();
		CloudColor.Dispose();
		CloudHeight.Dispose();
		CloudElevation.Dispose();
		VelocityArrow.Dispose();
		SurfacePosition.Dispose();
	}
}

[BurstCompile]
public struct BuildRenderStateCellJob : IJobParallelFor {
	public NativeArray<Vector4> TerrainColor1;
	public NativeArray<Vector4> TerrainColor2;
	public NativeArray<float> TerrainElevation;
	public NativeArray<Color32> WaterColor;
	public NativeArray<Color32> OverlayColor;
	public NativeArray<float> WaterElevation;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> CloudHeight;
	public NativeArray<Color32> CloudColor;
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> VelocityArrow;

	[ReadOnly] public NativeArray<float> Roughness;
	[ReadOnly] public NativeArray<float> SoilFertility;
	[ReadOnly] public NativeArray<float> Elevation;
	[ReadOnly] public NativeArray<float> CloudElevationSim;
	[ReadOnly] public NativeArray<float> CloudMass;
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
	[ReadOnly] public float LavaCrystalizationTemperature;
	[ReadOnly] public float LavaTemperatureRangeInverse;
	[ReadOnly] public float DustMaxInverse;
	[ReadOnly] public float LavaToRockMassAdjustment;
	[ReadOnly] public float DisplaySoilWeight;
	[ReadOnly] public float DisplaySandWeight;
	[ReadOnly] public float DisplayFloraWeight;
	[ReadOnly] public float DisplayCloudHeight;
	[ReadOnly] public float DisplayCloudMin;
	[ReadOnly] public float DisplayCloudRangeInverse;

	public void Execute(int i)
	{
		Vector4 terrainColor1;
		Vector4 terrainColor2;
		float terrainElevation;
		Color32 waterColor;
		Color32 overlayColor;
		float waterElevation;
		Color32 cloudColor;
		float cloudElevation;
		float cloudHeight;
		float3 velocityArrow;
		float3 surfacePosition;


		float waterDepth = WaterDepth[i];
		var icosphere = Icosphere[i];
		var elevation = Elevation[i];
		float roughness = math.max(1, Roughness[i]);
		// TODO: draw roughness correctly
		roughness = 0.1f;

		float surfaceElevation = elevation + math.max(roughness, waterDepth);

		if (MeshOverlayActive)
		{
			overlayColor = CVP.Lerp(MeshOverlayColors, (MeshOverlayData[i] - MeshOverlayMin) * MeshOverlayInverseRange);
		} else
		{
			overlayColor = new Color32();
		}
		// Terrain color
		float fertility = math.saturate(SoilFertility[i] / SoilFertilityMax);
		float lavaCoverage = math.min(1, LavaMass[i] / WorldData.MassLava) * 0.1f;
		float iceCoverage = IceCoverage[i];
		float waterCoverage = WaterCoverage[i];
		float floraCoverage = FloraCoverage[i] * DisplayFloraWeight;
		float dirtCoverage = 1.0f;

		terrainColor1 = new Vector4(
			0.5f,
			dirtCoverage * (1 - fertility),
			dirtCoverage * fertility,
			floraCoverage
			);

		terrainColor2 = new Vector4(
			iceCoverage,
			waterCoverage,
			lavaCoverage,
			0
			);

		waterColor = GetWaterColor(iceCoverage, WaterTemperature[i], waterDepth, PlanktonMass[i]);

		float cloudVolume = CloudMass[i] * 2;
		float cloudCoverage = math.saturate((math.pow(cloudVolume, 0.6667f) - DisplayCloudMin) * DisplayCloudRangeInverse);
		float dustCoverage = math.saturate(math.sqrt(DustCoverage[i] * DustMaxInverse));
		cloudColor = Color32.Lerp(new Color32(255, 255, 255, 255), new Color32(0,0,0,255), dustCoverage);
		cloudColor.a = (byte)(255 * math.max(dustCoverage, cloudCoverage * 0.75f));

		terrainElevation = ((elevation + roughness) * TerrainScale + PlanetRadius) / PlanetRadius;
		waterElevation = ((waterDepth == 0) ? 0.99f : ((elevation + waterDepth) * TerrainScale + PlanetRadius) / PlanetRadius);
		surfacePosition = icosphere * ((surfaceElevation + 50) * TerrainScale + PlanetRadius) / PlanetRadius;
		cloudElevation = (5000 * TerrainScale + PlanetRadius) / PlanetRadius;
		cloudHeight = cloudCoverage == 0 ? 0 : (cloudVolume / cloudCoverage * DisplayCloudHeight);


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

		TerrainColor1[i] = terrainColor1;
		TerrainColor2[i] = terrainColor2;
		TerrainElevation[i] = terrainElevation;
		WaterColor[i] = waterColor;
		WaterElevation[i] = waterElevation;
		CloudColor[i] = cloudColor;
		CloudElevation[i] = cloudElevation;
		CloudHeight[i] = cloudHeight;
		VelocityArrow[i] = velocityArrow;
		SurfacePosition[i] = surfacePosition;
		OverlayColor[i] = overlayColor;
	}


	private Color32 GetWaterColor(float iceCoverage, float waterTemperature, float depth, float plankton)
	{
		int blue = 150;
		int red = 0;
		int green = 20;
		int alpha = 200;
		if (plankton > 1)
		{
			green += 20;
		}
		else if (plankton > 0.25f)
		{
			green += 15;
		}
		else if (plankton > 0.05f)
		{
			green += 10;
		}
		else if (plankton > 0.01f)
		{
			green += 5;
		}
		else
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
		if (depth == 0)
		{
			alpha = 0;
		}
		return new Color32((byte)math.min(255, red), (byte)math.min(255, green), (byte)math.min(255, blue), (byte)alpha);
	}					   


}



[BurstCompile]
public struct BuildTerrainVertsJob : IJobParallelFor {
	public NativeArray<Vector3> VTerrainPosition;
	public NativeArray<Vector4> VTerrainColor1;
	public NativeArray<Vector4> VTerrainColor2;
	public NativeArray<Vector3> VWaterPosition;
	public NativeArray<Vector3> VWaterNormal;
	public NativeArray<Color32> VWaterColor;
	public NativeArray<Color32> VOverlayColor;

	[ReadOnly] public NativeArray<float> TerrainElevation;
	[ReadOnly] public NativeArray<float> WaterElevation;
	[ReadOnly] public NativeArray<Vector4> TerrainColor1;
	[ReadOnly] public NativeArray<Vector4> TerrainColor2;
	[ReadOnly] public NativeArray<Color32> WaterColor;
	[ReadOnly] public NativeArray<float3> StandardVerts;
	[ReadOnly] public NativeArray<Color32> OverlayColor;

	public void Execute(int i)
	{
		float3 v = StandardVerts[i];
		int j = (int)(i / WorldView.VertsPerCell);

		VTerrainPosition[i] = v * TerrainElevation[j];
		VWaterPosition[i] = v * WaterElevation[j];
		VTerrainColor1[i] = TerrainColor1[j];
		VTerrainColor2[i] = TerrainColor2[j];
		VWaterColor[i] = WaterColor[j];
		VOverlayColor[i] = OverlayColor[j];

		VWaterNormal[i] = v;
	}

}


[BurstCompile]
public struct BuildCloudVertsJob : IJobParallelFor {
	public NativeArray<Vector3> VCloudPosition;
	public NativeArray<Vector3> VCloudNormal;
	public NativeArray<Color32> VCloudColor;

	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudHeight;
	[ReadOnly] public NativeArray<Color32> CloudColor;
	[ReadOnly] public NativeArray<float3> StandardVerts;

	public void Execute(int i)
	{
		float3 v = StandardVerts[i];
		int j = (int)(i / WorldView.VertsPerCell);

		float cloudElevation = CloudElevation[j];
		//if (i % WorldView.VertsPerCell < 7)
		//{
		//	cloudElevation += CloudHeight[j];
		//}
		VCloudPosition[i] = v * cloudElevation;
		VCloudNormal[i] = v;
		VCloudColor[i] = CloudColor[j];
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
struct LerpJobVector4 : IJobParallelFor {
	public NativeArray<Vector4> Out;
	[ReadOnly] public NativeArray<Vector4> Start;
	[ReadOnly] public NativeArray<Vector4> End;
	[ReadOnly] public float Progress;
	public void Execute(int i)
	{
		Out[i] = Vector4.Lerp(Start[i], End[i], Progress);
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

