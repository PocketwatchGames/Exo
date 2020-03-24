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
	public float CellRadius;
	public float CellCircumference;
	public NativeArray<float2> Coordinate;
	public NativeArray<float3> SphericalPosition;
	public NativeArray<int> Neighbors;
	public NativeArray<float> NeighborDistInverse;
	public NativeArray<float> NeighborDist;
	public NativeArray<float3> NeighborDir;
	public NativeArray<float3> NeighborDiffInverse;
	public NativeArray<float> CoriolisMultiplier;


	public void Init(float radius, Icosphere icosphere, ref WorldData worldData)
	{
		PlanetRadius = radius;
		Count = icosphere.Vertices.Length;
		Coordinate = new NativeArray<float2>(Count, Allocator.Persistent);
		SphericalPosition = new NativeArray<float3>(Count, Allocator.Persistent);
		CoriolisMultiplier = new NativeArray<float>(Count, Allocator.Persistent);
		Neighbors = new NativeArray<int>(Count * 6, Allocator.Persistent);
		NeighborDir = new NativeArray<float3>(Count * 6, Allocator.Persistent);
		NeighborDistInverse = new NativeArray<float>(Count * 6, Allocator.Persistent);
		NeighborDiffInverse = new NativeArray<float3>(Count * 6, Allocator.Persistent);
		NeighborDist = new NativeArray<float>(Count * 6, Allocator.Persistent);
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
			for (int j = 0; j < 6; j++)
			{
				int index = i * 6 + j;
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

				}
				else
				{
					Neighbors[index] = -1;
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
		NeighborDir.Dispose();
		NeighborDist.Dispose();
		NeighborDistInverse.Dispose();
		NeighborDiffInverse.Dispose();
		Coordinate.Dispose();
		SphericalPosition.Dispose();
		CoriolisMultiplier.Dispose();
	}

	public int GetWaterIndex(int layer, int i)
	{
		return Count * layer + i;
	}
}
