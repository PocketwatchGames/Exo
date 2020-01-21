using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


[BurstCompile]
public struct WindJob : IJobParallelFor {

	public NativeArray<SimWind> Wind;

	[ReadOnly] public SimPlanetState PlanetState;
	[ReadOnly] public NativeArray<SimCell> Last;
	[ReadOnly] public WorldData worldData;
	[ReadOnly] public StaticState staticState;

	public void Execute(int i)
	{
		SimWind next = new SimWind();
		next.WindTropopause = float2.zero;
		next.WindSurface = float2.zero;
		next.CurrentDeep = float2.zero;
		next.CurrentSurface = float2.zero;

		var last = Last[i];
		float inverseFullIceCoverage = 1.0f / worldData.FullIceCoverage;
		var windInfo = staticState.WindInfo[i];

		float latitude = staticState.Coordinate[i].y;
		float lowerPressure = last.AirPressure;
		float lowerTemperature = last.AirTemperature;
		float elevation = last.Elevation;
		float waterDepth = last.WaterDepth;
		float waterAndIceDepth = last.WaterAndIceDepth;
		float elevationOrSeaLevel = elevation + waterAndIceDepth;
		float iceCoverage = last.IceMass * inverseFullIceCoverage;
		float friction;
		if (waterDepth > 0)
		{
			friction = worldData.WindOceanFriction;
		}
		else
		{
			friction = worldData.WindLandFriction;
		}
		friction = math.saturate(math.lerp(friction, worldData.WindIceFriction, iceCoverage));
		float lowerTemperatureAtSeaLevel = lowerTemperature - WorldData.TemperatureLapseRate * elevationOrSeaLevel;
		float molarMassLowerAir = Atmosphere.GetMolarMassAir(last.AirMass, last.AirWaterMass);
		float surfaceAirDensity = Atmosphere.GetAirDensity(lowerPressure, lowerTemperatureAtSeaLevel, molarMassLowerAir);
		float boundaryElevation = elevationOrSeaLevel + worldData.BoundaryZoneElevation;

		var lowerWindH = GetHorizontalWind(
			i, 
			latitude,
			PlanetState.AngularSpeed, 
			windInfo.coriolisParam, 
			windInfo.inverseCoriolisParam, 
			worldData.GlobalCoriolisInfluenceWindLower, 
			last.AirPressure, 
			elevationOrSeaLevel, 
			elevationOrSeaLevel,
			friction,
			surfaceAirDensity, 
			molarMassLowerAir);

		//// within 1 km of the ground, frictional forces slow wind down
		float neighborPressureDifferential = 0;
		float neighborElevationDifferential = 0;
		//if (lowerWindH.x < 0)
		//{
		//	int neighborIndex = world.GetNeighborIndex(x, y, 0);
		//	neighborPressureDifferential += -lowerWindH.x * (lowerPressure - state.LowerAirPressure[neighborIndex]);
		//	neighborElevationDifferential += -lowerWindH.x * (state.Elevation[neighborIndex] + state.WaterAndIceDepth[neighborIndex] - elevationOrSeaLevel);
		//}
		//else
		//{
		//	var neighborIndex = world.GetNeighborIndex(x, y, 1);
		//	neighborPressureDifferential += lowerWindH.x * (lowerPressure - state.LowerAirPressure[neighborIndex]);
		//	neighborElevationDifferential += lowerWindH.x * (state.Elevation[neighborIndex] + state.WaterAndIceDepth[neighborIndex] - elevationOrSeaLevel);
		//}
		//if (lowerWindH.y < 0)
		//{
		//	var neighborIndex = world.GetNeighborIndex(x, y, 3);
		//	neighborPressureDifferential += -lowerWindH.y * (lowerPressure - state.LowerAirPressure[neighborIndex]);
		//	neighborElevationDifferential += -lowerWindH.y * (state.Elevation[neighborIndex] + state.WaterAndIceDepth[neighborIndex] - elevationOrSeaLevel);
		//}
		//else
		//{
		//	var neighborIndex = world.GetNeighborIndex(x, y, 2);
		//	neighborPressureDifferential += lowerWindH.y * (lowerPressure - state.LowerAirPressure[neighborIndex]);
		//	neighborElevationDifferential += lowerWindH.y * (state.Elevation[neighborIndex] + state.WaterAndIceDepth[neighborIndex] - elevationOrSeaLevel);
		//}
		//var verticalTemperatureDifferential = lowerTemperatureAtSeaLevel - upperTemperatureAtSeaLevel;

		//float lowerWindSpeedH = lowerWindH.magnitude;
		float lowerWindV = neighborElevationDifferential * worldData.MountainUpdraftWindSpeed;
		lowerWindV += neighborPressureDifferential * worldData.DestinationPressureDifferentialToVerticalWindSpeed; // thermal
																												   //lowerWindV += (lowerPressure - upperPressure) * worldData.PressureToVerticalWindSpeed; // thermal

		next.WindTropopause = float2.zero;
		next.WindSurface = lowerWindH;
		next.WindVertical = lowerWindV;

		if (last.WaterDepth > 0)
		{
			float2 densityDifferential = float2.zero;
			float2 shallowCurrentH;
			float shallowCurrentV = 0;

			if (iceCoverage < 1)
			{
				shallowCurrentH = GetCurrentHorizontal(latitude, PlanetState.AngularSpeed, windInfo.coriolisParam, windInfo.inverseCoriolisParam, lowerWindH);
				if (iceCoverage > 0)
				{
					shallowCurrentH *= 1.0f - iceCoverage;
				}
			}
			else
			{
				shallowCurrentH = float2.zero;
			}

			//if (state.DeepWaterMass[index] > 0)
			//{
			//	float density = state.DeepWaterDensity[index];
			//	for (int i = 0; i < 4; i++)
			//	{
			//		var neighbor = world.GetNeighbor(x, y, i);
			//		int nIndex = world.GetIndex(neighbor.x, neighbor.y);
			//		if (last.DeepWaterMass[nIndex] > 0)
			//		{
			//			//var neighborWind = state.Wind[nIndex];
			//			//nWind += neighborWind;

			//			switch (i)
			//			{
			//				case 0:
			//					densityDifferential.x += last.DeepWaterDensity[nIndex] - density;
			//					break;
			//				case 1:
			//					densityDifferential.x -= last.DeepWaterDensity[nIndex] - density;
			//					break;
			//				case 2:
			//					densityDifferential.y -= last.DeepWaterDensity[nIndex] - density;
			//					break;
			//				case 3:
			//					densityDifferential.y += last.DeepWaterDensity[nIndex] - density;
			//					break;
			//			}
			//		}
			//		else
			//		{
			//			//switch (i)
			//			//{
			//			//	case 0:
			//			//		shallowCurrentV += shallowCurrentH.x;
			//			//		if (shallowCurrentH.x < 0)
			//			//		{
			//			//			shallowCurrentH.x = 0;
			//			//		}
			//			//		break;
			//			//	case 1:
			//			//		shallowCurrentV -= shallowCurrentH.x;
			//			//		if (shallowCurrentH.x > 0)
			//			//		{
			//			//			shallowCurrentH.x = 0;
			//			//		}
			//			//		break;
			//			//	case 2:
			//			//		shallowCurrentV -= shallowCurrentH.y;
			//			//		if (shallowCurrentH.y > 0)
			//			//		{
			//			//			shallowCurrentH.y = 0;
			//			//		}
			//			//		break;
			//			//	case 3:
			//			//		shallowCurrentV += shallowCurrentH.y;
			//			//		if (shallowCurrentH.y < 0)
			//			//		{
			//			//			shallowCurrentH.y = 0;
			//			//		}
			//			//		break;
			//			//}
			//		}
			//	}
			//}
			densityDifferential *= worldData.OceanDensityCurrentSpeed;
			next.CurrentDeep = densityDifferential;
			next.CurrentSurface = shallowCurrentH;
			next.CurrentVertical = shallowCurrentV;
		}
		else
		{
			next.CurrentDeep = float2.zero;
			next.CurrentSurface = float2.zero;
			next.CurrentVertical = 0;
		}

		Wind[i] = next;
	}


	public float2 GetHorizontalWind(int index, float latitude, float planetRotationSpeed, float coriolisParam, float inverseCoriolisParam, float coriolisInfluenceAtElevation, float thisPressure, float landElevation, float windElevation, float friction, float density, float molarMass)
	{
		float inverseDensity = 1.0f / density;
		float altitude = math.max(0, windElevation - landElevation);
		float complementFrictionAtElevation = 1.0f - friction * math.max(0, (worldData.BoundaryZoneElevation - altitude) / worldData.BoundaryZoneElevation);
		var pressureGradientForce = GetPressureGradient(index, thisPressure, windElevation, molarMass);
		pressureGradientForce.x *= PlanetState.Gravity * staticState.InverseCellDiameter;
		pressureGradientForce.y *= PlanetState.Gravity * staticState.InverseCellDiameter;
		float2 wind = float2.zero;

		//for (int i = 0; i < 6; i++)
		//{
		//	var nIndex = staticState.Neighbors[index * 6 + i];
		//	var neighborWind = state.UpperWind[nIndex];
		//	float nWindSpeed = Mathf.Sqrt(neighborWind.x * neighborWind.x + neighborWind.y * neighborWind.y);
		//	switch (i)
		//	{
		//		case 0:
		//			if (neighborWind.x > 0)
		//			{
		//				wind += neighborWind.x / nWindSpeed * new Vector3(neighborWind.x, neighborWind.y, 0);
		//			}
		//			break;
		//		case 1:
		//			if (neighborWind.x < 0)
		//			{
		//				wind += -neighborWind.x / nWindSpeed * new Vector3(neighborWind.x, neighborWind.y, 0);
		//			}
		//			break;
		//		case 2:
		//			if (neighborWind.y < 0)
		//			{
		//				wind += -neighborWind.y / nWindSpeed * new Vector3(neighborWind.x, neighborWind.y, 0);
		//			}
		//			break;
		//		case 3:
		//			if (neighborWind.y > 0)
		//			{
		//				wind += neighborWind.y / nWindSpeed * new Vector3(neighborWind.x, neighborWind.y, 0);
		//			}
		//			break;
		//	}
		//}

		var pressureWind = pressureGradientForce * worldData.PressureGradientWindMultiplier;
		if (coriolisParam != 0)
		{
			float geostrophicInfluence = math.saturate(math.pow(math.abs(coriolisParam) * 2, 2)) * complementFrictionAtElevation * coriolisInfluenceAtElevation;
			var geostrophicWind = new float2(-pressureGradientForce.y, pressureGradientForce.x) * inverseCoriolisParam / planetRotationSpeed;
			wind = (geostrophicWind * geostrophicInfluence + pressureWind * (1.0f - geostrophicInfluence));
		}
		else
		{
			wind = pressureWind;
		}
		wind *= complementFrictionAtElevation * inverseDensity;

		//wind += inertialWind * worldData.windInertia;
		return wind;
	}
	public float2 GetCurrentHorizontal(float latitude, float planetRotationSpeed, float coriolisParam, float inverseCoriolisParam, float2 pressureGradientForce)
	{
		float coriolisPower = math.abs(coriolisParam);
		float2 w = pressureGradientForce * worldData.WindToOceanCurrentFactor;
		if (coriolisPower > 0)
		{
			float geostrophicInfluence = math.saturate(math.pow(coriolisPower * 2, 2)) * worldData.GlobalCoriolisInfluenceOcean;
			var geostrophicCurrent = new float2(-w.y, w.x) * PlanetState.Gravity * staticState.InverseCellDiameter * inverseCoriolisParam / planetRotationSpeed;
			w = geostrophicCurrent * geostrophicInfluence + (1.0f - geostrophicInfluence) * w;
		}
		return w;
	}


	private float2 GetPressureGradient(int index, float pressureAtWindElevation, float windElevation, float molarMass)
	{
		float2 pressureDifferential = float2.zero;
		float2 nWind = float2.zero;
		var coord = staticState.Coordinate[index];
		int nCount = 0;
		for (int i = 0; i < 6; i++)
		{
			var nIndex = staticState.Neighbors[index*6+i];
			if (nIndex < 0)
			{
				continue;
			}

			var nCell = Last[nIndex];

			// see bottom of: https://en.wikipedia.org/wiki/Vertical_pressure_variation
			float neighborElevation = nCell.Elevation + nCell.WaterAndIceDepth;
			float neighborTemperatureElevation;
			neighborTemperatureElevation = neighborElevation;
			float neighborTemperatureAtSeaLevel = nCell.AirTemperature - neighborTemperatureElevation * WorldData.TemperatureLapseRate;
			float neighborElevationAtPressure = neighborTemperatureAtSeaLevel / WorldData.TemperatureLapseRate * (math.pow(pressureAtWindElevation / nCell.AirPressure, -1.0f / (worldData.PressureExponent * molarMass)) - 1);

			var coordDiff = staticState.Coordinate[nIndex] - coord;
			if (coordDiff.x > math.PI / 2)
			{
				coordDiff.x -= math.PI * 2;
			} else if (coordDiff.x < -math.PI / 2)
			{
				coordDiff.x += math.PI * 2;
			}
			pressureDifferential += (neighborElevationAtPressure - windElevation) * coordDiff;
			nCount++;
		}
		return pressureDifferential / nCount;

	}

}
