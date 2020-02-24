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
	public NativeArray<float2> Delta;
	public NativeArray<float> Buoyancy;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> VaporMass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Pressure;
	[ReadOnly] public NativeArray<float> LayerHeight;
	[ReadOnly] public NativeArray<float> LayerElevation;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float2> Coords;
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
	[ReadOnly] public float InverseCoordDiff;
	[ReadOnly] public float Gravity;
	[ReadOnly] public bool IsTop;
	[ReadOnly] public bool IsBottom;
	public void Execute(int i)
	{
		float2 gradientPressure = 0;
		float2 coord = Coords[i];
		float elevation = LayerElevation[i] + LayerHeight[i] / 2;
		float pressure = Pressure[i];

		for (int j = 0; j < 6; j++)
		{
			int neighborIndex = i * 6 + j;
			int n = Neighbors[neighborIndex];
			if (n >= 0)
			{
				float2 coordDiff = Coords[n] - coord;
				float2 diff = math.normalize(math.float2(Utils.WrapAngle(coordDiff.x), Utils.WrapAngle(coordDiff.y))); // TODO: cache the normalized vectors in staticstate

				float neighborMidElevation = LayerElevation[n] + LayerHeight[n] / 2;
				float neighborElevationAtPressure = Atmosphere.GetElevationAtPressure(pressure, Temperature[n], Pressure[n], neighborMidElevation, Gravity);
				gradientPressure += diff * (neighborElevationAtPressure - elevation);
				//	Debug.Log("I: " + i + " D: " + diff + " P: " + pressure + " E: " + (neighborElevationAtPressure - elevation) + " NT: " + Temperature[n] + " NP: " + Pressure[n] + " NE: " + neighborMidElevation + " NEAP: " + neighborElevationAtPressure);
			}
		}
		float inverseDensity = Atmosphere.GetInverseAirDensity(pressure, Temperature[i], AirMass[i], VaporMass[i]);
		Delta[i] = gradientPressure * Gravity * InverseCellDiameter * inverseDensity;


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
		//		float moveToNeutralBuoyancy = (UpTemperature[i] - Temperature[i]) / WorldData.TemperatureLapseRate - heightDiff;
		//		float vertMovement = math.min(MaxVerticalMovement, math.clamp(moveToNeutralBuoyancy + DiffusionCoefficient, 0, 1));

		//Debug.Log("I: " + i + " G: " + gradientWaterVapor + " HA: " + Humidity[i] + " HB: " + UpHumidity[i] + " AA: " + atmosphereMass + " AB: " + atmosphereMassUp);

		Buoyancy[i] = Gravity * buoyancy;

	}
}

