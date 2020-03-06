//#define AccelerationAirJobDebug
//#define WaterFrictionJobDebug
//#define WaterDensityGradientForceJobDebug

using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !WindFrictionJobDebug
[BurstCompile]
#endif
public struct AirTerrainFrictionJob : IJobParallelFor {
	public NativeArray<float> Force;
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
		Force[i] = exposedIce * IceFriction +
			exposedWater * WaterFriction +
			exposedVegetation * VegetationFriction +
			exposedTerrain * (TerrainFrictionMin + (TerrainFrictionMax - TerrainFrictionMin) * math.saturate(Terrain[i].Roughness / MaxTerrainRoughness));
	}
}

#if !AccelerationAirJobDebug
[BurstCompile]
#endif
public struct AccelerationAirJob : IJobParallelFor {
	public NativeArray<float3> Velocity;
	public NativeArray<float3> Force;
	[ReadOnly] public NativeArray<float> Friction;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> TemperaturePotential;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<float> UpTemperaturePotential;
	[ReadOnly] public NativeArray<float> UpHumidity;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> DownTemperaturePotential;
	[ReadOnly] public NativeArray<float> DownHumidity;
	[ReadOnly] public NativeArray<float> DownAirMass;
	[ReadOnly] public NativeArray<float> DownLayerElevation;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float Gravity;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float FrictionCoefficient;
	public void Execute(int i)
	{
		float3 gradientPressure = 0;
		float3 position = Positions[i];
		float elevation = LayerElevation[i] + LayerHeight[i] / 2;
		float pressure = Pressure[i];
		float3 force = 0;
		int neighborCount = 0;
		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float3 diff = (position - Positions[n]) * PlanetRadius;
				// TODO: this should only be using horizontal component, which we should cache

				float neighborElevationAtPressure = Atmosphere.GetElevationAtPressure(pressure, TemperaturePotential[n], Pressure[n], LayerElevation[n] + LayerHeight[n] / 2, Gravity);
				gradientPressure += diff / math.lengthsq(diff) * (neighborElevationAtPressure - elevation);
				neighborCount++;
			}
		}
		float inverseDensity = Atmosphere.GetInverseAirDensity(pressure, Atmosphere.GetAbsoluteTemperature(TemperaturePotential[i], elevation), AirMass[i], VaporMass[i]);
		force = gradientPressure * Gravity * inverseDensity / neighborCount;


		float buoyancy = 0;
		//if (!IsTop)
		//{
		//	float heightDiff = (UpLayerElevation[i] + UpLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
		//	buoyancy += Temperature[i] / potentialTemperatureUp - 1;
		//}
		//if (!IsBottom)
		//{
		//	float heightDiff = (DownLayerElevation[i] + DownLayerHeight[i] / 2) - (LayerElevation[i] + LayerHeight[i] / 2);
		//	buoyancy -= potentialTemperatureDown / Temperature[i] - 1;
		//}

		force += buoyancy * Gravity * Positions[i];
		Force[i] = force;
		Velocity[i] += force * SecondsPerTick - Velocity[i] * Friction[i] * FrictionCoefficient;
	}
}

#if !WaterSurfaceFrictionJobDebug
[BurstCompile]
#endif
public struct WaterSurfaceFrictionJob : IJobParallelFor {
	public NativeArray<float3> Force;
	[ReadOnly] public NativeArray<float3> Current;
	[ReadOnly] public NativeArray<float3> AirVelocityUp;
	[ReadOnly] public NativeArray<float3> AirVelocityDown;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float FrictionCoefficientUp;
	[ReadOnly] public float FrictionCoefficientDown;
	[ReadOnly] public float WaterSurfaceFrictionDepth;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		var horizontalWindUp = math.cross(math.cross(Position[i], AirVelocityUp[i]), Position[i]);
		var horizontalWindDown = math.cross(math.cross(Position[i], AirVelocityDown[i]), Position[i]);

		// TODO: the ekman depth is calculable!
		// https://en.wikipedia.org/wiki/Ekman_transport

		// Surface current is generally at about 45 degrees
		// this averages out to about a 90 degree turn if it has the full depth of the ekman spiral (200 meters)
		// http://oceanmotion.org/html/background/ocean-in-motion.htm

		var ekmanCurrent = math.cross(Position[i], horizontalWindUp) * CoriolisMultiplier[i] + horizontalWindUp * (1.0f - CoriolisMultiplier[i]);

		var force = (ekmanCurrent - Current[i]) * FrictionCoefficientUp + (horizontalWindDown - Current[i]) * FrictionCoefficientDown;

		Force[i] = force;
	}
}

#if !WaterDensityGradientForceJobDebug
[BurstCompile]
#endif
public struct AccelerationWaterJob : IJobParallelFor {
	public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> Friction;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<float> WaterDensity;
	[ReadOnly] public NativeArray<float> WaterPressure;
	[ReadOnly] public NativeArray<float> LayerDepth;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> UpLayerDepth;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float> UpWaterDensity;
	[ReadOnly] public NativeArray<float> UpWaterPressure;
	[ReadOnly] public NativeArray<float> DownLayerDepth;
	[ReadOnly] public NativeArray<float> DownLayerHeight;
	[ReadOnly] public NativeArray<float> DownWaterDensity;
	[ReadOnly] public NativeArray<float> DownWaterPressure;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float FrictionCoefficient;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		var density = WaterDensity[i];
		if (density > 0)
		{
			float3 pressureGradient = 0;
			var pos = Positions[i];
			float midDepth = LayerDepth[i] + LayerHeight[i] / 2;
			float pressure = WaterPressure[i];
			for (int j = 0; j < 6; j++)
			{
				var n = Neighbors[i * 6 + j];
				if (n >= 0 && WaterDensity[n] > 0)
				{
					// TODO: we should cache this value
					float3 diff = (pos - Positions[n]) * PlanetRadius;

					float neighborDepthAtPressure = Atmosphere.GetDepthAtPressure(pressure, WaterPressure[n], LayerDepth[n] + LayerHeight[n] / 2, WaterDensity[n], Gravity);
					pressureGradient += diff / math.lengthsq(diff) * (neighborDepthAtPressure - midDepth);
				}
			}

			var force = pressureGradient * Gravity / density;

			//	if (DownWaterDensity[i] > 0)
			//	{
			//		Force[i] += pos * Gravity * (DownWaterDensity[i] / density - 1)/* / (LayerHeight[i] + DownLayerHeight[i]) * 0.5f*/;
			//	}
			//	if (UpWaterDensity[i] > 0)
			//	{
			//		Force[i] += pos * Gravity * (1 - UpWaterDensity[i] / density)/* / (LayerHeight[i] + UpLayerHeight[i]) * 0.5f*/;
			//	}
			Velocity[i] += force * SecondsPerTick - Velocity[i] * Friction[i] * FrictionCoefficient;
		}

	}
}

#if !UpdateVelocityCloudJobDebug
[BurstCompile]
#endif
public struct UpdateVelocityCloudJob : IJobParallelFor {
	public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float3> LayerVelocity;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float3> UpVelocity;
	[ReadOnly] public NativeArray<float> UpLayerElevation;
	[ReadOnly] public NativeArray<float> UpLayerHeight;
	[ReadOnly] public NativeArray<float3> DownVelocity;
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
						Velocity[i] = LayerVelocity[i];
					}
					else
					{
						float downLayerMidElevation = (DownLayerElevation[i] + layerElevation) / 2;
						float t = (cloudElevation - downLayerMidElevation) / (layerMidHeight - downLayerMidElevation);
						Velocity[i] = LayerVelocity[i] * t + DownVelocity[i] * (1.0f - t);
					}
				}
				else if (IsTop)
				{
					Velocity[i] = LayerVelocity[i];
				}
				else
				{
					float upLayerMidElevation = UpLayerElevation[i] + UpLayerHeight[i] / 2;
					float t = (cloudElevation - layerMidHeight) / (upLayerMidElevation - layerMidHeight);
					Velocity[i] = UpVelocity[i] * t + LayerVelocity[i] * (1.0f - t);
				}
			}
		}
	}
}

