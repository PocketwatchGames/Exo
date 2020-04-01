
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


#if !WindFrictionJobDebug
[BurstCompile]
#endif
public struct AirTerrainFrictionJob : IJobParallelFor {
	public NativeArray<float> Force;
	[ReadOnly] public NativeArray<float> Roughness;
	[ReadOnly] public NativeArray<float> WaterCoverage;
	[ReadOnly] public NativeArray<float> IceCoverage;
	[ReadOnly] public NativeArray<float> FloraCoverage;
	[ReadOnly] public float IceFriction;
	[ReadOnly] public float WaterFriction;
	[ReadOnly] public float TerrainFrictionMin;
	[ReadOnly] public float TerrainFrictionMax;
	[ReadOnly] public float FloraFriction;
	[ReadOnly] public float MaxTerrainRoughness;
	public void Execute(int i)
	{
		float exposedIce = IceCoverage[i];
		float exposedWater = math.max(0, WaterCoverage[i] - exposedIce);
		float exposedFlora = math.max(0, FloraCoverage[i] - exposedIce - exposedWater);
		float exposedTerrain = 1.0f - exposedWater - exposedFlora - exposedIce;
		Force[i] = exposedIce * IceFriction +
			exposedWater * WaterFriction +
			exposedFlora * FloraFriction +
			exposedTerrain * (TerrainFrictionMin + (TerrainFrictionMax - TerrainFrictionMin) * math.saturate(Roughness[i] / MaxTerrainRoughness));
	}
}

[BurstCompile]
public struct AccelerationAirJob : IJobParallelFor {
	public NativeArray<float3> Velocity;
	public NativeArray<float3> Force;
	[ReadOnly] public NativeArray<float> Friction;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> TemperaturePotential;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> LayerMiddle;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> NeighborDiffInverse;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<float> UpTemperaturePotential;
	[ReadOnly] public NativeArray<float> UpHumidity;
	[ReadOnly] public NativeArray<float> UpAirMass;
	[ReadOnly] public NativeArray<float> UpLayerMiddle;
	[ReadOnly] public NativeArray<float> DownTemperaturePotential;
	[ReadOnly] public NativeArray<float> DownHumidity;
	[ReadOnly] public NativeArray<float> DownAirMass;
	[ReadOnly] public NativeArray<float> DownLayerMiddle;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float GravityInverse;
	[ReadOnly] public float Gravity;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float FrictionCoefficient;
	public void Execute(int i)
	{

		float3 gradientPressure = 0;
		float3 position = Positions[i];
		float pressure = Pressure[i];
		float3 force = 0;
		int neighborCount = 0;

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float neighborElevationAtPressure = Atmosphere.GetElevationAtPressure(pressure, TemperaturePotential[n], Pressure[n], LayerMiddle[n], GravityInverse);
				gradientPressure += NeighborDiffInverse[neighborIndex] * (neighborElevationAtPressure - LayerMiddle[i]);
				neighborCount++;
			}
		}
		float inverseDensity = Atmosphere.GetInverseAirDensity(pressure, Atmosphere.GetAbsoluteTemperature(TemperaturePotential[i], LayerMiddle[i]), AirMass[i], VaporMass[i]);
		force = gradientPressure * Gravity * inverseDensity / neighborCount;

		Force[i] = force;
		var vel = Velocity[i] + force * SecondsPerTick - Velocity[i] * Friction[i] * FrictionCoefficient;
		vel -= Positions[i] * math.dot(Positions[i], vel);

		float buoyancy = 0;
		// Cold air moves into warm air via gravity
		// Note that if warm air is above cold air, it is stratified but stable, so there's no force
		if (!IsBottom)
		{
			buoyancy = math.min(0, SecondsPerTick * Gravity * (DownTemperaturePotential[i] / TemperaturePotential[i] - 1));

			// TODO: this is temp, what's a reasonable way to apply a buoyant force over 3600 seconds?
			buoyancy = buoyancy * 0.001f;

			// TODO: it might be a good idea to just separate buoyancy velocity since we aren't preserving momentum, then we could get rid of this
			vel += Positions[i] * (buoyancy - math.dot(Positions[i], vel));
		}


		Velocity[i] = vel;
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
	[ReadOnly] public NativeArray<float3> NeighborDiffInverse;
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
	[ReadOnly] public NativeArray<float> SurfaceElevation;
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
			float midDepthElevation = SurfaceElevation[i] - (LayerDepth[i] - LayerHeight[i] / 2);
			float pressure = WaterPressure[i];
			for (int j = 0; j < 6; j++)
			{
				int neighborIndex = i * 6 + j;
				var n = Neighbors[neighborIndex];
				if (n >= 0 && WaterDensity[n] > 0)
				{
					float neighborElevationAtPressure = SurfaceElevation[n] - Atmosphere.GetDepthAtPressure(pressure, WaterPressure[n], LayerDepth[n] - LayerHeight[n] / 2, WaterDensity[n], Gravity);
					pressureGradient += NeighborDiffInverse[neighborIndex] * (neighborElevationAtPressure - midDepthElevation);
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

