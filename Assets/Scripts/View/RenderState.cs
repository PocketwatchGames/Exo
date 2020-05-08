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

	public NativeArray<float4> TerrainColor1;
	public NativeArray<float4> TerrainColor2;
	public NativeArray<float> TerrainElevation;
	public NativeArray<float4> WaterColor;
	public NativeArray<float> WaterElevation;
	public NativeArray<float3> WaterCurrent;
	public NativeArray<Color32> CloudColor;
	public NativeArray<Color32> OverlayColor;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> CloudHeight;
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> VelocityHorizontal;
	public NativeArray<float> VelocityVertical;

	public void Init(int count)
	{
		TerrainColor1 = new NativeArray<float4>(count, Allocator.Persistent);
		TerrainColor2 = new NativeArray<float4>(count, Allocator.Persistent);
		TerrainElevation = new NativeArray<float>(count, Allocator.Persistent);
		WaterColor = new NativeArray<float4>(count, Allocator.Persistent);
		WaterElevation = new NativeArray<float>(count, Allocator.Persistent);
		WaterCurrent = new NativeArray<float3>(count, Allocator.Persistent);
		OverlayColor = new NativeArray<Color32>(count, Allocator.Persistent);
		CloudColor = new NativeArray<Color32>(count, Allocator.Persistent);
		CloudElevation = new NativeArray<float>(count, Allocator.Persistent);
		CloudHeight = new NativeArray<float>(count, Allocator.Persistent);
		VelocityHorizontal = new NativeArray<float3>(count, Allocator.Persistent);
		VelocityVertical = new NativeArray<float>(count, Allocator.Persistent);
		SurfacePosition = new NativeArray<float3>(count, Allocator.Persistent);
	}

	public void Dispose()
	{
		TerrainColor1.Dispose();
		TerrainColor2.Dispose();
		TerrainElevation.Dispose();
		WaterColor.Dispose();
		WaterElevation.Dispose();
		WaterCurrent.Dispose();
		OverlayColor.Dispose();
		CloudColor.Dispose();
		CloudHeight.Dispose();
		CloudElevation.Dispose();
		VelocityHorizontal.Dispose();
		VelocityVertical.Dispose();
		SurfacePosition.Dispose();
	}
}

[BurstCompile]
public struct BuildRenderStateCellJob : IJobParallelFor {
	public NativeArray<float4> TerrainColor1;
	public NativeArray<float4> TerrainColor2;
	public NativeArray<float> TerrainElevation;
	public NativeArray<float4> WaterColor;
	public NativeArray<Color32> OverlayColor;
	public NativeArray<float> WaterElevation;
	public NativeArray<float3> RWaterCurrent;
	public NativeArray<float> CloudElevation;
	public NativeArray<float> CloudHeight;
	public NativeArray<Color32> CloudColor;
	public NativeArray<float3> SurfacePosition;
	public NativeArray<float3> VelocityHorizontal;
	public NativeArray<float> VelocityVertical;

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
	[ReadOnly] public NativeSlice<float> WaterCoverage;
	[ReadOnly] public NativeSlice<float> WaterDepth;
	[ReadOnly] public NativeSlice<float> WaterTemperature;
	[ReadOnly] public NativeSlice<float3> WaterCurrent;
	[ReadOnly] public NativeSlice<float> PlanktonMass;
	[ReadOnly] public NativeArray<float> SurfaceElevation;
	[ReadOnly] public NativeArray<float3> Icosphere;
	[ReadOnly] public NativeSlice<float> MeshOverlayData;
	[ReadOnly] public NativeArray<CVP> MeshOverlayColors;
	[ReadOnly] public NativeArray<CVP> WindColors;
	[ReadOnly] public NativeSlice<float3> WindOverlayData;
	[ReadOnly] public NativeArray<float> LavaMass;
	[ReadOnly] public NativeArray<float> LavaTemperature;
	[ReadOnly] public NativeSlice<float> DustCoverage;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public bool MeshOverlayActive;
	[ReadOnly] public float MeshOverlayMin;
	[ReadOnly] public float MeshOverlayInverseRange;
	[ReadOnly] public bool WindOverlayActive;
	[ReadOnly] public bool WindMaskedByLand;
	[ReadOnly] public float WindVelocityMax;
	[ReadOnly] public float WindVerticalMax;
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
	[ReadOnly] public float PlanktonPower;
	[ReadOnly] public float PlanktonMax;
	[ReadOnly] public int PlanktonLevels;
	[ReadOnly] public int IceLevels;
	[ReadOnly] public float WaterTemperatureMax;
	[ReadOnly] public int WaterTemperatureLevels;
	[ReadOnly] public int AirLayers;
	[ReadOnly] public int Count;

	public void Execute(int i)
	{
		float4 terrainColor1;
		float4 terrainColor2;
		float terrainElevation;
		float4 waterColor;
		Color32 overlayColor;
		float waterElevation;
		Color32 cloudColor;
		float cloudElevation;
		float cloudHeight;
		float3 velocityArrow;
		float3 surfacePosition;
		float velocityVertical = 0;

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
		}
		else
		{
			overlayColor = new Color32();
		}
		// Terrain color
		float fertility = math.saturate(SoilFertility[i] / SoilFertilityMax);
		float lavaCoverage = math.min(1, LavaMass[i] / WorldData.MassLava) * 0.1f;
		float iceCoverage = math.floor(IceLevels * math.pow(IceCoverage[i], 0.4f)) / IceLevels;
		float waterCoverage = WaterCoverage[i];
		float floraCoverage = FloraCoverage[i] * DisplayFloraWeight;
		float dirtCoverage = 1.0f;
		float plankton = math.floor(PlanktonLevels * math.saturate(math.pow(PlanktonMass[i] / PlanktonMax, PlanktonPower))) / PlanktonLevels;

		terrainColor1 = new float4(
			0.5f,
			dirtCoverage * (1 - fertility),
			dirtCoverage * fertility,
			floraCoverage
			);

		terrainColor2 = new float4(
			iceCoverage * (1.0f - waterCoverage),
			waterCoverage,
			lavaCoverage,
			0
			);

		float waterTemperatureColor = math.floor(WaterTemperatureLevels * math.saturate((WaterTemperature[i] - WorldData.FreezingTemperature) / WaterTemperatureMax)) / WaterTemperatureLevels;
		waterColor = new float4(plankton, iceCoverage, waterTemperatureColor, waterDepth);

		float cloudVolume = CloudMass[i] * 2;
		float cloudCoverage = math.saturate((math.pow(cloudVolume, 0.6667f) - DisplayCloudMin) * DisplayCloudRangeInverse);

		float dust = 0;
		for (int j = 0; j < AirLayers; j++)
		{
			dust += DustCoverage[i + j * Count];
		}
		float dustCoverage = math.saturate(math.sqrt(dust * DustMaxInverse));
		cloudColor = Color32.Lerp(new Color32(255, 255, 255, 255), new Color32(0, 0, 0, 255), dustCoverage);
		cloudColor.a = (byte)(255 * math.max(dustCoverage, cloudCoverage * 0.75f));

		terrainElevation = ((elevation + roughness) * TerrainScale + PlanetRadius) / PlanetRadius;
		waterElevation = ((waterDepth == 0) ? 0.99f : ((elevation + waterDepth) * TerrainScale + PlanetRadius) / PlanetRadius);
		surfacePosition = icosphere * ((surfaceElevation + 50) * TerrainScale + PlanetRadius) / PlanetRadius;
		cloudElevation = (5000 * TerrainScale + PlanetRadius) / PlanetRadius;
		cloudHeight = cloudCoverage == 0 ? 0 : (cloudVolume / cloudCoverage * DisplayCloudHeight);


		if (WindOverlayActive)
		{
			velocityArrow = Utils.GetHorizontalComponent(WindOverlayData[i], Positions[i]) / WindVelocityMax;
			if (math.length(velocityArrow) < 0.01f)
			{
				velocityArrow = float3.zero;
			}
			else if (WindMaskedByLand && waterDepth < roughness)
			{
				velocityArrow = float3.zero;
			}
			velocityVertical = math.dot(WindOverlayData[i], Positions[i]) / WindVerticalMax;
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
		VelocityHorizontal[i] = velocityArrow;
		VelocityVertical[i] = velocityVertical;
		SurfacePosition[i] = surfacePosition;
		OverlayColor[i] = overlayColor;
		RWaterCurrent[i] = Utils.GetHorizontalComponent(WaterCurrent[i], Positions[i]);
	}


}


[BurstCompile]
public struct BuildTerrainVertsJob : IJobParallelFor {
	public NativeArray<Vector3> VTerrainPosition;
	public NativeArray<Vector4> VTerrainColor1;
	public NativeArray<Vector4> VTerrainColor2;
	public NativeArray<Vector3> VWaterPosition;
	public NativeArray<Vector3> VWaterNormal;
	public NativeArray<Vector4> VWaterColor;
	public NativeArray<Vector4> VWaterCurrent;
	public NativeArray<Color32> VOverlayColor;

	[ReadOnly] public NativeArray<float> TerrainElevation;
	[ReadOnly] public NativeArray<float> WaterElevation;
	[ReadOnly] public NativeArray<float4> TerrainColor1;
	[ReadOnly] public NativeArray<float4> TerrainColor2;
	[ReadOnly] public NativeArray<float4> WaterColor;
	[ReadOnly] public NativeArray<float3> WaterCurrent;
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
		VOverlayColor[i] = OverlayColor[j];
		var c = WaterCurrent[j];
		int cellVert = (i % WorldView.VertsPerCell);
		VWaterCurrent[i] = new float4(c.x, c.y, c.z, (cellVert == 0) ? 1 : 0);
		var w = WaterColor[j];
		VWaterColor[i] = new Vector4(
			w.x,
			w.y,
			w.z,
			(w.w > 0) ? ((cellVert == 0) ? 2 : 1) : 0
			);

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
struct LerpJobfloat4 : IJobParallelFor {
	public NativeArray<float4> Out;
	[ReadOnly] public NativeArray<float4> Start;
	[ReadOnly] public NativeArray<float4> End;
	[ReadOnly] public float Progress;
	public void Execute(int i)
	{
		Out[i] = math.lerp(Start[i], End[i], Progress);
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

