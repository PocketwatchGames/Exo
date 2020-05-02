
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


[BurstCompile]
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
	public NativeSlice<float3> Velocity;
	public NativeSlice<float3> Force;
	[ReadOnly] public NativeArray<float> Friction;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeSlice<float> VaporMass;
	[ReadOnly] public NativeSlice<float> TemperaturePotential;
	[ReadOnly] public NativeSlice<float> NewTemperaturePotential;
	[ReadOnly] public NativeSlice<float> Pressure;
	[ReadOnly] public NativeSlice<float> LayerMiddle;
	[ReadOnly] public NativeSlice<float3> LastVelocity;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float3> NeighborDiffInverse;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public float PlanetRadius;
	[ReadOnly] public float GravityInverse;
	[ReadOnly] public float Gravity;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public int LayerCount;
	[ReadOnly] public int Count;
	public void Execute(int i)
	{

		float3 gradientPressure = 0;
		float pressure = Pressure[i];
		float layerMiddle = LayerMiddle[i];
		float3 force = 0;
		int neighborCount = 0;

		int layer = i / Count;
		bool isTop = layer == LayerCount - 1;
		bool isBottom = layer == 0;
		int cellIndex = i - layer * Count;

		for (int j = 0; j < StaticState.MaxNeighbors; j++)
		{
			int neighborIndex = cellIndex * StaticState.MaxNeighbors + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				n += layer * Count;
				float neighborElevationAtPressure = Atmosphere.GetElevationAtPressure(pressure, TemperaturePotential[n], Pressure[n], LayerMiddle[n], GravityInverse);
				gradientPressure += NeighborDiffInverse[neighborIndex] * (neighborElevationAtPressure - layerMiddle);
				neighborCount++;
			}
		}
		float inverseDensity = Atmosphere.GetInverseAirDensity(pressure, Atmosphere.GetAbsoluteTemperature(TemperaturePotential[i], layerMiddle), AirMass[i], VaporMass[i]);
		force = gradientPressure * Gravity * inverseDensity / neighborCount;

		Force[i] = force;
		var vel = LastVelocity[i] + force * SecondsPerTick;

		// we are leaving out the SecondsPerTick term here
		// TODO: is there a better way to integrate this over the full time step?
		if (isBottom)
		{
			vel -= LastVelocity[i] * Friction[cellIndex];
		}

		float3 pos = Positions[cellIndex];

		// Remove vertical component
		// TODO: we shouldn't need to do this if we are correctly turning the velocity
		vel -= Utils.GetVerticalComponent(vel, pos);

		float buoyancy = 0;
		if (!isTop)
		{
			int upIndex = i + Count;
			// When air is heated, it pushes surrounding air out of the way
			buoyancy += math.max(0, 1 - NewTemperaturePotential[upIndex] / NewTemperaturePotential[i]);
		}
		if (!isBottom)
		{
			int downIndex = i - Count;
			// Cold air sitting on top of warm air produces a downdraft
			// Note that if warm air is above cold air, it is stratified but stable, so there's no force
			buoyancy += math.min(0, 1 - TemperaturePotential[downIndex] / TemperaturePotential[i]);
		}

		// TODO: I've removed the secondsPerTick term here because it produces massive vertical velocities, but it's not physically accurate
		// We ARE accelerating over the entire time of the tick, so figure out how the pro ACGMs deal with buoyancy
		vel += pos * Gravity * buoyancy;
		//vel += Positions[i] * SecondsPerTick * Gravity * buoyancy;

		// TODO: is there a better way to deal with boundary conditions than this?
		if (isTop)
		{
			float dotUp = math.dot(pos, vel);
			if (dotUp > 0)
			{
				vel -= pos * dotUp;
			}
		}
		if (isBottom)
		{
			float dotUp = math.dot(pos, vel);
			if (dotUp < 0)
			{
				vel -= pos * dotUp;
			}
		}

		Velocity[i] = vel;
	}
}

[BurstCompile]
public struct WaterSurfaceFrictionJob : IJobParallelFor {
	public NativeArray<float3> Force;
	[ReadOnly] public NativeArray<float3> Current;
	[ReadOnly] public NativeSlice<float3> AirVelocityUp;
	[ReadOnly] public NativeSlice<float3> AirVelocityDown;
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
		var horizontalWindUp = Utils.GetHorizontalComponent(AirVelocityUp[i], Position[i]);
		var horizontalWindDown = Utils.GetHorizontalComponent(AirVelocityDown[i], Position[i]);

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

[BurstCompile]
public struct AccelerationWaterJob : IJobParallelFor {
	public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> LastVelocity;
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
			Velocity[i] = LastVelocity[i] + (force + Friction[i] * FrictionCoefficient) * SecondsPerTick;
		}

	}
}

