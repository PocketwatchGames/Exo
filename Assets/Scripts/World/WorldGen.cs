using System;
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
	private static float GetPerlinMinMax(float x, float y, float z, float frequency, float hash, float min, float max)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency, z * frequency) + 1.0f) * (max - min) / 2 + min;
	}
	private static float GetPerlinNormalized(float x, float y, float frequency, float hash)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency) + 1.0f) / 2;
	}
	private static float GetPerlinNormalized(float x, float y, float z, float frequency, float hash)
	{
		return (_noise.GetPerlin(x * frequency + hash, y * frequency, z * frequency) + 1.0f) / 2;
	}
	private static float GetPerlin(float x, float y, float frequency, float hash)
	{
		return _noise.GetPerlin(x * frequency + hash, y * frequency);
	}
	private static float GetPerlin(float x, float y, float z, float frequency, float hash)
	{
		return _noise.GetPerlin(x * frequency + hash, y * frequency, z * frequency);
	}
	public static void Generate(StaticState staticState, int seed, WorldGenData worldGenData, WorldData worldData, ref SimState state)
	{
		float inversePI = 1.0f / math.PI;
		_noise = new FastNoise(seed);
		_noise.SetFrequency(10);
		_random = new System.Random(seed);
		staticState.Radius = worldGenData.Radius;
		state.TiltAngle = worldGenData.TiltAngle;
		state.SpinSpeed = 1.0f / worldGenData.SpinTime;
		state.OrbitSpeed = 1.0f / worldGenData.OrbitTime;
		for (int i = 0; i < state.Count; i++)
		{
			SimStateCell cell;
			var coord = staticState.Coordinate[i] * 2 * inversePI;
			var pos = staticState.SphericalPosition[i];
			var coordNormalized = (coord + 1) / 2;
			cell.Roughness = math.max(1, math.pow(GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 6643230), 3) * worldGenData.MaxRoughness);
			float elevation =
				0.6f * GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 0, worldGenData.MinElevation, worldGenData.MaxElevation) +
				0.3f * GetPerlinMinMax(pos.x, pos.y, pos.z, 0.5f, 0, worldGenData.MinElevation, worldGenData.MaxElevation) +
				0.1f * GetPerlinMinMax(pos.x, pos.y, pos.z, 2.0f, 0, worldGenData.MinElevation, worldGenData.MaxElevation);
			cell.Elevation = elevation;
			cell.WaterDepth = math.max(0, -elevation) + math.max(0, GetPerlinMinMax(pos.x, pos.y, pos.z, 1.0f, 25430, -100, 100));
			//cell.WaterDepth = GetPerlinNormalized(coordNormalized.x, coordNormalized.y, 2.5f, 110) * 2000;

			float soil =
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.2f, 6630) +
				0.3f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 43560) +
				0.2f * GetPerlinNormalized(pos.x, pos.y, pos.z, 1f, 8740);
			cell.SoilFertility = soil;

			cell.Vegetation = elevation < 0 ? 0 : math.pow(math.sqrt(soil) * math.sqrt(GetPerlinNormalized(pos.x, pos.y, pos.z, 1f, 6630423)), 3);
			cell.Ice = 3 * math.pow(GetPerlinNormalized(pos.x, pos.y, pos.z, 10, 100), 0.25f) * math.abs(math.pow(coord.y, 3));

			cell.CloudElevation =
				0.5f * GetPerlinMinMax(pos.x, pos.y, pos.z, 0.1f, 20, 0, worldGenData.MaxElevation) +
				0.5f * GetPerlinMinMax(pos.x, pos.y, pos.z, 1.0f, 20, 0, worldGenData.MaxElevation);

			cell.CloudCoverage =
				math.pow(
					0.75f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 30) +
					0.25f * GetPerlinNormalized(pos.x, pos.y, pos.z, 1.0f, 30),
					2);

			cell.RelativeHumidity =
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.1f, 40) +
				0.5f * GetPerlinNormalized(pos.x, pos.y, pos.z, 0.5f, 40);

			state.Cells[i] = cell;
		}
	}
}
