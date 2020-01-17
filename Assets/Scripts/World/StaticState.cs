using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;
using Unity.Collections;

public struct StaticState {

	public int Count;
	public float PlanetRadius;
	public float CellSurfaceArea;
	public float CellDiameter;
	public float InverseCellDiameter;
	public NativeArray<float2> Coordinate;
	public NativeArray<float3> SphericalPosition;
	public NativeArray<int> Neighbors;


	public void Init(float radius, Icosphere icosphere, ref WorldData worldData)
	{
		PlanetRadius = radius;
		Count = icosphere.Vertices.Count;
		Coordinate = new NativeArray<float2>(Count, Allocator.Persistent);
		SphericalPosition = new NativeArray<float3>(Count, Allocator.Persistent);;
		Neighbors = new NativeArray<int>(Count * 6, Allocator.Persistent);
		float surfaceArea = 4 * math.PI * PlanetRadius * PlanetRadius;
		CellSurfaceArea = surfaceArea / Count;
		CellDiameter = 2 * math.sqrt(CellSurfaceArea / (4 * math.PI));
		InverseCellDiameter = 1.0f / CellDiameter;

		var neighborList = new List<int>[Count];
		for (int i = 0; i < Count; i++)
		{
			neighborList[i] = new List<int>();
		}

		for (int i = 0; i < Count; i++)
		{
			var v = icosphere.Vertices[i];
			Coordinate[i] = new float2(math.atan2(v.x, v.z), math.asin(v.y));
			SphericalPosition[i] = new float3(v.x, v.y, v.z);
		}

		for (int i = 0; i < icosphere.Polygons.Count; i++)
		{
			var p = icosphere.Polygons[i];
			for (int j = 0; j < 3; j++)
			{
				neighborList[p.m_Vertices[j]].Add(p.m_Vertices[(j + 1) % 3]);
			}
		}
		for (int i = 0; i < Count; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				if (j < neighborList[i].Count)
				{
					Neighbors[i * 6 + j] = neighborList[i][j];
				}
				else
				{
					Neighbors[i * 6 + j] = -1;
				}
			}
		}




		//for (int y = 0; y < size; y++)
		//{
		//	float latitude = ((float)y / size) * 2 - 1.0f;
		//	float yaw = (float)(latitude * Math.PI * 1.5f);
		//	float pitch = (float)(latitude * Math.PI * 3f);
		//	float absSinPitch = (float)(Math.Abs(Math.Sin(pitch)));
		//	float cosYaw = (float)Math.Cos(yaw);
		//	float cosPitch = (float)Math.Cos(pitch);

		//	float tropopauseElevation = (1.0f - Math.Abs(latitude)) * (MaxTropopauseElevation - MinTropopauseElevation) + MinTropopauseElevation + TropopauseElevationSeason * latitude;
		//	windInfo[y] = new WindInfo()
		//	{
		//		latitude = latitude,
		//		yaw = yaw,
		//		tropopauseElevationMax = tropopauseElevation,
		//		coriolisParam = Mathf.Sin(latitude * Mathf.PI / 2),
		//		inverseCoriolisParam = 1.0f / Mathf.Sin(latitude * Mathf.PI / 2)
		//	};
		//}

	}

	public void Dispose()
	{
		Neighbors.Dispose();
		Coordinate.Dispose();
		SphericalPosition.Dispose();
	}
}
