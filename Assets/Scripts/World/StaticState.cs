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
	public NativeArray<StaticWindInfo> WindInfo;


	private NativeArray<PolarCoordinateLookup> _polarLookup;
	private NativeArray<PolarCoordinateLookup> _polarLookupLatitudeIndex;

	private struct PolarCoordinateLookup {
		public float angle;
		public int index;
	}
	public struct StaticWindInfo {
		public float coriolisParam;
		public float inverseCoriolisParam;
	}


	public void Init(float radius, Icosphere icosphere, ref WorldData worldData)
	{
		PlanetRadius = radius;
		Count = icosphere.Vertices.Count;
		Coordinate = new NativeArray<float2>(Count, Allocator.Persistent);
		SphericalPosition = new NativeArray<float3>(Count, Allocator.Persistent); ;
		WindInfo = new NativeArray<StaticWindInfo>(Count, Allocator.Persistent); ;
		Neighbors = new NativeArray<int>(Count * 6, Allocator.Persistent);
		_polarLookup = new NativeArray<PolarCoordinateLookup>(Count, Allocator.Persistent);
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
			Coordinate[i] = new float2(-math.atan2(v.x, v.z), math.asin(v.y));
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

			WindInfo[i] = new StaticWindInfo()
			{
				coriolisParam = math.sin(latitude),
				inverseCoriolisParam = 1.0f / math.sin(latitude)
			};

		}

		_polarLookupLatitudeIndex = new NativeArray<PolarCoordinateLookup>(vertsByCoord.Count, Allocator.Persistent);
		int curIndex = 0;
		for (int i = 0; i< vertsByCoord.Count;i++)
		{
			var vertsByLat = vertsByCoord.ElementAt(i);
			_polarLookupLatitudeIndex[i] = new PolarCoordinateLookup { angle = vertsByLat.Key, index = curIndex };
			for (int j = 0; j < vertsByLat.Value.Count; j++)
			{
				var vertByLong = vertsByLat.Value.ElementAt(j);
				_polarLookup[curIndex] = new PolarCoordinateLookup { angle = vertByLong.Key, index = vertByLong.Value };
				curIndex++;
			}
		}


	}

	public void Dispose()
	{
		Neighbors.Dispose();
		Coordinate.Dispose();
		SphericalPosition.Dispose();
		WindInfo.Dispose();
		_polarLookup.Dispose();
		_polarLookupLatitudeIndex.Dispose();
	}

	public int GetClosestVertByPolarCoord(float2 pos)
	{
		int latitudeStartIndex = 0;
		int latitudeEndIndex = 0;
		for (int i=0;i<_polarLookupLatitudeIndex.Length;i++)
		{
			if (_polarLookupLatitudeIndex[i].angle > pos.y)
			{
				int latitudeSpan = math.max(0, i - 1);
				latitudeStartIndex = _polarLookupLatitudeIndex[latitudeSpan].index;
				latitudeEndIndex = ((latitudeSpan < _polarLookupLatitudeIndex.Length - 1) ? _polarLookupLatitudeIndex[i].index : _polarLookup.Length);
				for (int j = latitudeStartIndex; j< latitudeEndIndex;j++)
				{
					if (_polarLookup[j].angle > pos.x)
					{
						return _polarLookup[j].index;
					}
				}
			}
		}

		return 0;

	}
}
