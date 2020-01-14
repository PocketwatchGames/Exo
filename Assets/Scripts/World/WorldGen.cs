﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;

public static class WorldGen {
	static System.Random _random;
	static FastNoise _noise;

	private static float GetPerlinMinMax(float x, float y, float frequency, float hash, float min, float max)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency) + 1.0f) * (max - min) / 2 + min;
	}
	private static float GetPerlinNormalized(float x, float y, float frequency, float hash)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency) + 1.0f) / 2;
	}
	private static float GetPerlin(float x, float y, float frequency, float hash)
	{
		return _noise.GetPerlin(x * frequency + hash, y * frequency);
	}
	public static void Generate(SimState state, StaticState staticState, int seed, WorldGenData worldGenData)
	{
		float inversePI = 1.0f / math.PI;
		_noise = new FastNoise(seed);
		_noise.SetFrequency(10);
		_random = new System.Random(seed);
		staticState.Radius = worldGenData.Radius;
		for (int i = 0; i < state.Count; i++)
		{
			var coord = staticState.Coordinate[i] * 2 * inversePI;
			var coordNormalized = (coord + 1) / 2;
			float elevation =
				0.4f * GetPerlinMinMax(coordNormalized.x, coordNormalized.y, 0.1f, 0, worldGenData.MinElevation, worldGenData.MaxElevation) +
				0.3f * GetPerlinMinMax(coordNormalized.x, coordNormalized.y, 0.25f, 0, worldGenData.MinElevation, worldGenData.MaxElevation) +
				0.2f * GetPerlinMinMax(coordNormalized.x, coordNormalized.y, 1f, 0, worldGenData.MinElevation, worldGenData.MaxElevation) +
				0.1f * GetPerlinMinMax(coordNormalized.x, coordNormalized.y, 2.5f, 0, worldGenData.MinElevation, worldGenData.MaxElevation);
			state.Elevation[i] = elevation;
			state.WaterElevation[i] = 0;
			state.Vegetation[i] = elevation < 0 ? 0 : (float)_random.NextDouble();
			state.Ice[i] = GetPerlinNormalized(coordNormalized.x, coordNormalized.y, 10, 100) * coord.y * coord.y;

			state.CloudElevation[i] =
				0.5f * GetPerlinMinMax(coordNormalized.x, coordNormalized.y, 0.1f, 20, 0, worldGenData.MaxElevation) +
				0.5f * GetPerlinMinMax(coordNormalized.x, coordNormalized.y, 1.0f, 20, 0, worldGenData.MaxElevation);

			state.CloudCoverage[i] =
				0.5f * GetPerlinNormalized(coordNormalized.x, coordNormalized.y, 0.1f, 30) +
				0.5f * GetPerlinNormalized(coordNormalized.x, coordNormalized.y, 1.0f, 30);

			state.RelativeHumidity[i] =
				0.5f * GetPerlinNormalized(coordNormalized.x, coordNormalized.y, 0.1f, 40) +
				0.5f * GetPerlinNormalized(coordNormalized.x, coordNormalized.y, 1.0f, 40);

		}
	}
}
