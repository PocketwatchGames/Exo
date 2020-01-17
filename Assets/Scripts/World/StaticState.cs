using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;
using Unity.Collections;

public class StaticState {

	public int Count;
	public float Radius;
	public float2[] Coordinate;
	public float3[] SphericalPosition;
	public List<int>[] NeighborList;
	public NativeArray<int> Neighbors;

	public float TicksPerHour;
	public float TicksPerSecond;
	public float SecondsPerTick;
	public float InverseMetersPerTile;


	public void Init(int count, Icosphere icosphere)
	{
		Count = count;
		Coordinate = new float2[count];
		SphericalPosition = new float3[count];
		NeighborList = new List<int>[count];
		Neighbors = new NativeArray<int>(count * 6, Allocator.Persistent);
		for (int i = 0; i < count; i++)
		{
			NeighborList[i] = new List<int>();
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
				NeighborList[p.m_Vertices[j]].Add(p.m_Vertices[(j + 1) % 3]);
			}
		}
		for (int i = 0; i < Count; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				if (j < NeighborList[i].Count)
				{
					Neighbors[i * 6 + j] = NeighborList[i][j];
				}
				else
				{
					Neighbors[i * 6 + j] = -1;
				}
			}
		}


		//TicksPerHour = (float)TicksPerYear / (365 * 24);
		//int secondsPerYear = 365 * 24 * 60 * 60;
		//TicksPerSecond = (float)TicksPerYear / secondsPerYear;
		//SecondsPerTick = 1.0f / TicksPerSecond;

		//InverseMetersPerTile = 1.0f / MetersPerTile;

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
	}
}
