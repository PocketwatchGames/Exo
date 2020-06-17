using System.Collections.Generic;
using Unity.Collections;
using Unity.Mathematics;

public static class SimTest {
	public static bool CheckForDegeneracy(int cellCount, ref SimState state, ref SimState lastState, ref StaticState staticState, ref TempState dependent, ref WorldData worldData, ref SortedSet<int> degenIndices, ref List<string> degenVarNames)
	{
		bool degen = false;
		degen |= CheckDegenMinMaxValues(cellCount, degenIndices, "TerrainTemperature", state.GroundTemperature, 0, 1200, degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "GroundWater", state.GroundWater, degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "CloudMass", state.CloudMass, degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "CloudDropletMass", state.CloudDropletMass, degenVarNames);
		degen |= CheckDegenMinMaxValues(cellCount, degenIndices, "CloudElevation", dependent.CloudElevation, -100000, 100000, degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "IceMass", state.IceMass, degenVarNames);
		degen |= CheckDegenMinMaxValues(cellCount, degenIndices, "IceTemperature", state.IceTemperature, 0, 1200, degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "CrustDepth", state.LavaMass, degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "MagmaMass", state.LavaMass, degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "LavaMass", state.LavaMass, degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "LavaTemperature", state.LavaTemperature, degenVarNames);

		degen |= CheckDegenMinMaxValues(cellCount, degenIndices, "AirTemperature", staticState.GetSliceAir(state.AirTemperaturePotential), 0, 1200, degenVarNames);
		degen |= CheckDegenMinMaxValues(cellCount, degenIndices, "AirVapor", staticState.GetSliceAir(state.AirVapor), 0, 10000, degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "CarbonDioxide", staticState.GetSliceAir(state.AirCarbonDioxide), degenVarNames);
		degen |= CheckDegen(cellCount, degenIndices, "AirVelocity", staticState.GetSliceAir(state.AirVelocity), degenVarNames);

		degen |= CheckDegenPosValues(cellCount, degenIndices, "WaterMass", staticState.GetSliceWater(state.WaterMass), degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "SaltMass", staticState.GetSliceWater(state.WaterSaltMass), degenVarNames);
		degen |= CheckDegenPosValues(cellCount, degenIndices, "CarbonMass", staticState.GetSliceWater(state.WaterCarbonDioxide), degenVarNames);
		degen |= CheckDegenMinMaxValues(cellCount, degenIndices, "WaterTemperature", staticState.GetSliceWater(state.WaterTemperature), 0, 1200, degenVarNames);
		degen |= CheckDegen(cellCount, degenIndices, "Current", staticState.GetSliceWater(state.WaterVelocity), degenVarNames);

		return degen;
	}

	private static bool CheckDegen(int count, SortedSet<int> degenIndices, string name, NativeSlice<float> values, List<string> degenVarNames)
	{
		for (int i = 0; i < count; i++)
		{
			float v = values[i];
			if (!math.isfinite(v))
			{
				degenVarNames.Add(name);
				degenIndices.Add(i);
				return true;
			}
		}
		return false;
	}
	private static bool CheckDegen(int count, SortedSet<int> degenIndices, string name, NativeSlice<float3> values, List<string> degenVarNames)
	{
		for (int i = 0; i < count; i++)
		{
			if (!math.isfinite(values[i].x) || !math.isfinite(values[i].y) || !math.isfinite(values[i].z))
			{
				degenVarNames.Add(name);
				degenIndices.Add(i);
				return true;
			}
		}
		return false;
	}

	private static bool CheckDegenPosValues(int count, SortedSet<int> degenIndices, string name, NativeSlice<float> values, List<string> degenVarNames)
	{
		for (int i = 0; i < count; i++)
		{
			if (!math.isfinite(values[i]) || values[i] < 0)
			{
				degenVarNames.Add(name);
				degenIndices.Add(i);
				return true;
			}
		}
		return false;
	}

	private static bool CheckDegenMinMaxValues(int count, SortedSet<int> degenIndices, string name, NativeSlice<float> values, float min, float max, List<string> degenVarNames)
	{
		for (int i = 0; i < count; i++)
		{
			if (!math.isfinite(values[i]) || values[i] < min || values[i] > max)
			{
				degenVarNames.Add(name);
				degenIndices.Add(i);
				return true;
			}
		}
		return false;
	}

}
