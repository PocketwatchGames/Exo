using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

public static class WorldGen {
	static System.Random _random;
	public static void Generate(SimState state, StaticState staticState, int seed, WorldGenData worldGenData)
	{
		_random = new Random(seed);
		staticState.Radius = worldGenData.Radius;
		for (int i = 0; i < state.Count; i++)
		{
			state.Elevation[i] = (float)(_random.NextDouble() * (worldGenData.MaxElevation - worldGenData.MinElevation) + worldGenData.MinElevation);
			state.Vegetation[i] = (float)_random.NextDouble();
		}
	}
}
