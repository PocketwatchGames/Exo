//#define PressureGradientForceAirJobDebug
//#define WaterFrictionJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !WindFrictionJobDebug
[BurstCompile]
#endif
public struct WindFrictionJob : IJobParallelFor {
	public NativeArray<float> Friction;
	[ReadOnly] public NativeArray<CellTerrain> Terrain;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> VegetationCoverage;
	[ReadOnly] public float IceFriction;
	[ReadOnly] public float WaterFriction;
	[ReadOnly] public float TerrainFrictionMin;
	[ReadOnly] public float TerrainFrictionMax;
	[ReadOnly] public float VegetationFriction;
	[ReadOnly] public float MaxTerrainRoughness;
	public void Execute(int i)
	{
		float exposedIce = IceCoverage[i];
		float exposedWater = math.max(0, WaterCoverage[i] - exposedIce);
		float exposedVegetation = math.max(0, VegetationCoverage[i] - exposedIce - exposedWater);
		float exposedTerrain = 1.0f - exposedWater - exposedVegetation - exposedIce;
		Friction[i] = exposedIce * IceFriction +
			exposedWater * WaterFriction +
			exposedVegetation * VegetationFriction +
			exposedTerrain * (TerrainFrictionMin + (TerrainFrictionMax - TerrainFrictionMin) * math.saturate(Terrain[i].Roughness / MaxTerrainRoughness));
		//Debug.Log("I: " + i + " F: " + Friction[i]);
	}
}

#if !PressureGradientForceAirJobDebug
[BurstCompile]
#endif
public struct PressureGradientForceAirJob : IJobParallelFor {
	public NativeArray<float3> Delta;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<float> UpTemperature;
	[ReadOnly] public NativeArray<float> UpHumidity;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> DownTemperature;
	[ReadOnly] public NativeArray<float> DownHumidity;
	[ReadOnly] public NativeArray<float> DownAirMass;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float Gravity;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	public void Execute(int i)
	{
		float3 gradientPressure = 0;
		float3 position = Positions[i];
		float elevation = LayerElevation[i] + LayerHeight[i] / 2;
		float pressure = Pressure[i];
		float3 force = 0;

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float3 diff = math.normalize(position - Positions[n]);

				float neighborMidElevation = LayerElevation[n] + LayerHeight[n] / 2;
				float neighborElevationAtPressure = Atmosphere.GetElevationAtPressure(pressure, Temperature[n], Pressure[n], neighborMidElevation, Gravity);
				gradientPressure += diff * (neighborElevationAtPressure - elevation);
			}
		}
		float inverseDensity = Atmosphere.GetInverseAirDensity(pressure, Temperature[i], AirMass[i], VaporMass[i]);
		force = gradientPressure * Gravity * InverseCellDiameter * inverseDensity;


		float buoyancy = 0;
		if (!IsTop)
		{
			float heightDiff = (UpLayerElevation[i] + UpLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
			float potentialTemperatureUp = UpTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
			buoyancy += Temperature[i] / potentialTemperatureUp - 1;
		}
		if (!IsBottom)
		{
			float heightDiff = (DownLayerElevation[i] + DownLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
			float potentialTemperatureDown = DownTemperature[i] - WorldData.TemperatureLapseRate * heightDiff;
			buoyancy -= potentialTemperatureDown / Temperature[i] - 1;
		}

		Delta[i] = force + buoyancy * Gravity * Positions[i];
	}
}

#if !PressureGradientForceAirJobDebug
[BurstCompile]
#endif
public struct PressureGradientForceCloudJob : IJobParallelFor {
	public NativeArray<float3> Force;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float3> LayerForce;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float3> UpForce;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float3> DownForce;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	public void Execute(int i)
	{
		float cloudElevation = CloudElevation[i];
		float layerElevation = LayerElevation[i];
		float layerMidHeight = layerElevation + LayerHeight[i] / 2;
		if (cloudElevation >= layerElevation || IsBottom)
		{
			if (cloudElevation < layerElevation + LayerHeight[i] || IsTop)
			{
				if (cloudElevation < layerMidHeight)
				{
					if (IsBottom)
					{
						Force[i] = LayerForce[i];
					}else
					{
						float downLayerMidElevation = (DownLayerElevation[i] + layerElevation) / 2;
						float t = (cloudElevation - downLayerMidElevation) / (layerMidHeight - downLayerMidElevation);
						Force[i] = LayerForce[i] * t + DownForce[i] * (1.0f - t);
					}
				} else if (IsTop)
				{
					Force[i] = LayerForce[i];
				} else
				{
					float upLayerMidElevation = UpLayerElevation[i] + UpLayerHeight[i] / 2;
					float t = (cloudElevation - layerMidHeight) / (upLayerMidElevation - layerMidHeight);
					Force[i] = UpForce[i] * t + LayerForce[i] * (1.0f - t);
				}
			}
		}
	}
}

#if !WaterFrictionJobDebug
[BurstCompile]
#endif
public struct WaterFrictionForceJob : IJobParallelFor {
	public NativeArray<float3> Force;
	[ReadOnly] public NativeArray<float3> Current;
	[ReadOnly] public NativeArray<float3> WindUp;
	[ReadOnly] public NativeArray<float3> WindDown;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public float FrictionCoefficientUp;
	[ReadOnly] public float FrictionCoefficientDown;
	public void Execute(int i)
	{
		var horizontalWindUp = math.cross(math.cross(Positions[i], WindUp[i]), Positions[i]);
		var horizontalWindDown = math.cross(math.cross(Positions[i], WindDown[i]), Positions[i]);
		Force[i] = (horizontalWindUp - Current[i]) * FrictionCoefficientUp + (horizontalWindDown - Current[i]) * FrictionCoefficientDown;
	}
}

#if !WaterDensityGradientForceJobDebug
[BurstCompile]
#endif
public struct WaterDensityGradientForceJob : IJobParallelFor {
	public NativeArray<float3> Force;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<float> WaterDensity;
	[ReadOnly] public NativeArray<float> UpWaterDensity;
	[ReadOnly] public NativeArray<float> DownWaterDensity;
	[ReadOnly] public float InverseCellDiameter;
	[ReadOnly] public float Gravity;
	public void Execute(int i)
	{
		float3 densityGradient = 0;
		var pos = Positions[i];
		var density = WaterDensity[i];
		for (int j=0;j<6;j++)
		{
			var n = Neighbors[i * 6 + j];
			if (n >= 0)
			{
				densityGradient += (WaterDensity[n] - density) * math.normalize(pos - Positions[n]);
			}
		}

		float3 buoyancy = pos * (DownWaterDensity[i] - UpWaterDensity[i]) * Gravity;

		Force[i] = densityGradient * InverseCellDiameter * Gravity;
	}
}

