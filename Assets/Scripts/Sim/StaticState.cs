using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;
using Unity.Collections;

public struct StaticState {

	public const int MaxNeighbors = 6;

	public int Count;
	public float PlanetRadius;
	public float CellSurfaceArea;
	public float CellRadius;
	public float CellCircumference;
	public NativeArray<float2> Coordinate;
	public NativeArray<float3> SphericalPosition;
	public NativeArray<int> Neighbors;
	public NativeArray<int> ReverseNeighbors;
	public NativeArray<float> NeighborDistInverse;
	public NativeArray<float> NeighborDist;
	public NativeArray<float3> NeighborDir;
	public NativeArray<float3> NeighborTangent;
	public NativeArray<float3> NeighborDiffInverse;
	public NativeArray<float> CoriolisMultiplier;


	public void Init(float radius, Icosphere icosphere, ref WorldData worldData)
	{
		PlanetRadius = radius;
		Count = icosphere.Vertices.Length;
		Coordinate = new NativeArray<float2>(Count, Allocator.Persistent);
		SphericalPosition = new NativeArray<float3>(Count, Allocator.Persistent);
		CoriolisMultiplier = new NativeArray<float>(Count, Allocator.Persistent);
		Neighbors = new NativeArray<int>(Count * MaxNeighbors, Allocator.Persistent);
		ReverseNeighbors = new NativeArray<int>(Count * MaxNeighbors, Allocator.Persistent);
		NeighborDir = new NativeArray<float3>(Count * MaxNeighbors, Allocator.Persistent);
		NeighborDistInverse = new NativeArray<float>(Count * MaxNeighbors, Allocator.Persistent);
		NeighborTangent = new NativeArray<float3>(Count * MaxNeighbors, Allocator.Persistent);
		NeighborDiffInverse = new NativeArray<float3>(Count * MaxNeighbors, Allocator.Persistent);
		NeighborDist = new NativeArray<float>(Count * MaxNeighbors, Allocator.Persistent);
		float surfaceArea = 4 * math.PI * PlanetRadius * PlanetRadius;
		CellSurfaceArea = surfaceArea / Count;
		CellRadius = math.sqrt(CellSurfaceArea / math.PI);
		CellCircumference = math.PI * 2 * CellRadius;

		var neighborList = new List<Tuple<int, float3>>[Count];
		for (int i = 0; i < Count; i++)
		{
			neighborList[i] = new List<Tuple<int, float3>>();
		}

		for (int i = 0; i < Count; i++)
		{
			var v = icosphere.Vertices[i];
			Coordinate[i] = new float2(-math.atan2(v.x, v.z), math.asin(v.y));
			SphericalPosition[i] = new float3(v.x, v.y, v.z);
		}

		for (int i = 0; i < icosphere.Polygons.Count; i++)
		{
			var p = icosphere.Polygons[i];
			for (int j = 0; j < 3; j++)
			{
				int vertIndex = p.m_Vertices[(j + 1) % 3];
				neighborList[p.m_Vertices[j]].Add(new Tuple<int, float3>(vertIndex, SphericalPosition[vertIndex]));
			}
		}
		for (int i = 0; i < Count; i++)
		{
			var pos = SphericalPosition[i];
			var forward = neighborList[i][0].Item2 - pos;
			float forwardLength = math.length(forward);

			neighborList[i].Sort(delegate (Tuple<int, float3> a, Tuple<int, float3> b)
			{
				float3 diffA = a.Item2 - pos;
				float3 diffB = b.Item2 - pos;
				float dotA = math.dot(diffA, forward);
				float dotB = math.dot(diffB, forward);
				float angleA = diffA.Equals(forward) ? 0 : math.acos(dotA / (math.length(diffA) * forwardLength));
				float angleB = diffB.Equals(forward) ? 0 : math.acos(dotB / (math.length(diffB) * forwardLength));
				angleA *= math.sign(math.dot(pos, math.cross(forward, diffA)));
				angleB *= math.sign(math.dot(pos, math.cross(forward, diffB)));
				return (int)math.sign(angleB - angleA);
			});
			for (int j = 0; j < MaxNeighbors; j++)
			{
				int index = i * MaxNeighbors + j;
				if (j < neighborList[i].Count)
				{
					int n = neighborList[i][j].Item1;
					Neighbors[index] = n;
					var diff = pos - SphericalPosition[n];
					float dist = math.length(diff * PlanetRadius);
					NeighborDist[index] = dist;
					NeighborDistInverse[index] = 1.0f / dist;
					NeighborDir[index] = math.normalize(math.cross(math.cross(pos, diff), pos));
					NeighborDiffInverse[index] = NeighborDir[index] * NeighborDistInverse[index];
					NeighborTangent[index] = NeighborDir[index] * dist;

				}
				else
				{
					Neighbors[index] = -1;
				}
			}
		}

		for (int i = 0; i < Count * MaxNeighbors; i++)
		{
			ReverseNeighbors[i] = -1;
			int cellIndex = i / MaxNeighbors;
			int nIndex = Neighbors[i];
			if (nIndex >= 0)
			{
				for (int j = 0; j < MaxNeighbors; j++)
				{
					if (Neighbors[nIndex * MaxNeighbors + j] == cellIndex)
					{
						ReverseNeighbors[i] = nIndex * MaxNeighbors + j;
					}
				}
			}
		}

		SortedDictionary<float, SortedDictionary<float, int>> vertsByCoord = new SortedDictionary<float, SortedDictionary<float, int>>();
		for (int i = 0; i < Coordinate.Length; i++)
		{
			SortedDictionary<float, int> vertsAtLatitude;
			float latitude = Coordinate[i].y;
			if (!vertsByCoord.TryGetValue(latitude, out vertsAtLatitude))
			{
				vertsAtLatitude = new SortedDictionary<float, int>();
				vertsByCoord.Add(latitude, vertsAtLatitude);
			}
			vertsAtLatitude.Add(Coordinate[i].x, i);

			CoriolisMultiplier[i] = math.sin(latitude);

		}


	}

	public void Dispose()
	{
		Neighbors.Dispose();
		ReverseNeighbors.Dispose();
		NeighborDir.Dispose();
		NeighborDist.Dispose();
		NeighborDistInverse.Dispose();
		NeighborDiffInverse.Dispose();
		NeighborTangent.Dispose();
		Coordinate.Dispose();
		SphericalPosition.Dispose();
		CoriolisMultiplier.Dispose();
	}

	public int GetWaterIndex(int layer, int i)
	{
		return Count * layer + i;
	}

	public static int GetMaxNeighbors(int cell, NativeArray<int> neighbors)
	{
		return (neighbors[(cell + 1) * MaxNeighbors - 1] >= 0) ? MaxNeighbors : (MaxNeighbors - 1);
	}
}
