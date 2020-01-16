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

	}
}
