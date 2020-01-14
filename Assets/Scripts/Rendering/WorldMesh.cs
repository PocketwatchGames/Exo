using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;

public struct MeshBuilder {


	List<Polygon> m_Polygons;
	List<Vector3> m_Vertices;
	Vector3[] terrainVertices;
	Vector3[] terrainNormals;
	Color32[] terrainColors;
	Vector3[] waterVertices;
	Vector3[] waterNormals;
	Color32[] waterColors;
	Vector3[] cloudVertices;
	Vector3[] cloudNormals;
	Color32[] cloudColors;

	static Color32 green = new Color32(20, 255, 30, 255);
	static Color32 brown = new Color32(220, 150, 70, 255);
	static Color32 blue = new Color32(0, 0, 255, 255);
	static Color32 white = new Color32(255, 255, 255, 255);
	static Color32 black = new Color32(0, 0, 0, 255);


	public void Init(int subdivisions)
	{
		InitAsIcosohedron();

		// TODO: this subdivides recursively, not by a fraction of the edge length, which allows for more precision -- see wikipedia page
		Subdivide(subdivisions);
	}

	public List<Vector3> GetVertices()
	{
		return m_Vertices;
	}


	public void InitMesh(Mesh terrainMesh, Mesh waterMesh, Mesh cloudMesh)
	{

		int vertexCount = m_Polygons.Count * 3;

		int[] indices = new int[vertexCount];
		terrainVertices = new Vector3[m_Vertices.Count];
		terrainNormals = new Vector3[m_Vertices.Count];
		terrainColors = new Color32[m_Vertices.Count];
		waterVertices = new Vector3[m_Vertices.Count];
		waterNormals = new Vector3[m_Vertices.Count];
		waterColors = new Color32[m_Vertices.Count];
		cloudVertices = new Vector3[m_Vertices.Count];
		cloudNormals = new Vector3[m_Vertices.Count];
		cloudColors = new Color32[m_Vertices.Count];

		for (int i = 0; i < m_Polygons.Count; i++)
		{
			var poly = m_Polygons[i];

			indices[i * 3 + 0] = poly.m_Vertices[0];
			indices[i * 3 + 1] = poly.m_Vertices[1];
			indices[i * 3 + 2] = poly.m_Vertices[2];
		}

		terrainMesh.vertices = terrainVertices;
		terrainMesh.normals = terrainNormals;
		terrainMesh.colors32 = terrainColors;
		terrainMesh.SetTriangles(indices, 0);

		waterMesh.vertices = waterVertices;
		waterMesh.normals = waterNormals;
		waterMesh.colors32 = waterColors;
		waterMesh.SetTriangles(indices, 0);

		cloudMesh.vertices = cloudVertices;
		cloudMesh.normals = cloudNormals;
		cloudMesh.colors32 = cloudColors;
		cloudMesh.SetTriangles(indices, 0);

	}

	public void UpdateMesh(Mesh terrainMesh, Mesh waterMesh, Mesh cloudMesh, StaticState staticState, SimState lastState, SimState nextState, float t, float scale)
	{
		for (int i = 0; i < m_Vertices.Count; i++)
		{
			var lastStateCell = lastState.Cells[i];
			var nextStateCell = nextState.Cells[i];
			float lastElevation = lastStateCell.Elevation;
			float nextElevation = nextStateCell.Elevation;
			terrainVertices[i] = m_Vertices[i] * ((nextElevation - lastElevation) * t + lastElevation + staticState.Radius) * scale;
			terrainColors[i] = Color32.Lerp(GetTerrainColor(staticState, lastState, i), GetTerrainColor(staticState, nextState, i), t);
			terrainNormals[i] = m_Vertices[i];

			float lastWaterElevation = lastStateCell.WaterElevation;
			float nextWaterElevation = nextStateCell.WaterElevation;
			waterVertices[i] = m_Vertices[i] * ((nextWaterElevation - lastWaterElevation) * t + lastWaterElevation + staticState.Radius) * scale;
			waterColors[i] = Color32.Lerp(GetWaterColor(staticState, lastState, i), GetWaterColor(staticState, nextState, i), t);
			waterNormals[i] = m_Vertices[i];

			float lastCloudElevation = lastStateCell.CloudElevation;
			float nextCloudElevation = nextStateCell.CloudElevation;
			cloudVertices[i] = m_Vertices[i] * ((nextCloudElevation - lastCloudElevation) * t + lastCloudElevation + staticState.Radius) * scale;
			cloudColors[i] = Color32.Lerp(GetCloudColor(staticState, lastState, i), GetCloudColor(staticState, nextState, i), t);
			cloudNormals[i] = m_Vertices[i];

		}

		terrainMesh.vertices = terrainVertices;
		terrainMesh.normals = terrainNormals;
		terrainMesh.colors32 = terrainColors;
		terrainMesh.RecalculateBounds();

		waterMesh.vertices = waterVertices;
		waterMesh.normals = waterNormals;
		waterMesh.colors32 = waterColors;
		waterMesh.RecalculateBounds();

		cloudMesh.vertices = cloudVertices;
		cloudMesh.normals = cloudNormals;
		cloudMesh.colors32 = cloudColors;
		cloudMesh.RecalculateBounds();

	}

	#region private functions

	private Color32 GetTerrainColor(StaticState staticState, SimState state, int index)
	{
		return Color32.Lerp(brown, green, state.Cells[index].Vegetation);
	}

	private Color32 GetWaterColor(StaticState staticState, SimState state, int index)
	{
		return Color32.Lerp(blue, white, state.Cells[index].Ice);
	}

	private Color32 GetCloudColor(StaticState staticState, SimState state, int index)
	{
		var humidity = state.Cells[index].RelativeHumidity;
		var c = Color32.Lerp(white, black, humidity);
		float opacity = humidity > 0.5f ? math.pow(humidity, 0.25f) * 0.75f : math.pow(humidity, 4);
		c.a = (byte)(255 * opacity);
		return c;
	}

	private void InitAsIcosohedron()
	{
		m_Polygons = new List<Polygon>();
		m_Vertices = new List<Vector3>();

		// An icosahedron has 12 vertices, and
		// since they're completely symmetrical the
		// formula for calculating them is kind of
		// symmetrical too:

		float t = (1.0f + Mathf.Sqrt(5.0f)) / 2.0f;

		m_Vertices.Add(new Vector3(-1, t, 0).normalized);
		m_Vertices.Add(new Vector3(1, t, 0).normalized);
		m_Vertices.Add(new Vector3(-1, -t, 0).normalized);
		m_Vertices.Add(new Vector3(1, -t, 0).normalized);
		m_Vertices.Add(new Vector3(0, -1, t).normalized);
		m_Vertices.Add(new Vector3(0, 1, t).normalized);
		m_Vertices.Add(new Vector3(0, -1, -t).normalized);
		m_Vertices.Add(new Vector3(0, 1, -t).normalized);
		m_Vertices.Add(new Vector3(t, 0, -1).normalized);
		m_Vertices.Add(new Vector3(t, 0, 1).normalized);
		m_Vertices.Add(new Vector3(-t, 0, -1).normalized);
		m_Vertices.Add(new Vector3(-t, 0, 1).normalized);

		// And here's the formula for the 20 sides,
		// referencing the 12 vertices we just created.

		m_Polygons.Add(new Polygon(0, 11, 5));
		m_Polygons.Add(new Polygon(0, 5, 1));
		m_Polygons.Add(new Polygon(0, 1, 7));
		m_Polygons.Add(new Polygon(0, 7, 10));
		m_Polygons.Add(new Polygon(0, 10, 11));
		m_Polygons.Add(new Polygon(1, 5, 9));
		m_Polygons.Add(new Polygon(5, 11, 4));
		m_Polygons.Add(new Polygon(11, 10, 2));
		m_Polygons.Add(new Polygon(10, 7, 6));
		m_Polygons.Add(new Polygon(7, 1, 8));
		m_Polygons.Add(new Polygon(3, 9, 4));
		m_Polygons.Add(new Polygon(3, 4, 2));
		m_Polygons.Add(new Polygon(3, 2, 6));
		m_Polygons.Add(new Polygon(3, 6, 8));
		m_Polygons.Add(new Polygon(3, 8, 9));
		m_Polygons.Add(new Polygon(4, 9, 5));
		m_Polygons.Add(new Polygon(2, 4, 11));
		m_Polygons.Add(new Polygon(6, 2, 10));
		m_Polygons.Add(new Polygon(8, 6, 7));
		m_Polygons.Add(new Polygon(9, 8, 1));
	}

	private void Subdivide(int recursions)
	{
		var midPointCache = new Dictionary<int, int>();

		for (int i = 0; i < recursions; i++)
		{
			var newPolys = new List<Polygon>();
			foreach (var poly in m_Polygons)
			{
				int a = poly.m_Vertices[0];
				int b = poly.m_Vertices[1];
				int c = poly.m_Vertices[2];

				// Use GetMidPointIndex to either create a
				// new vertex between two old vertices, or
				// find the one that was already created.

				int ab = GetMidPointIndex(midPointCache, a, b);
				int bc = GetMidPointIndex(midPointCache, b, c);
				int ca = GetMidPointIndex(midPointCache, c, a);

				// Create the four new polygons using our original
				// three vertices, and the three new midpoints.
				newPolys.Add(new Polygon(a, ab, ca));
				newPolys.Add(new Polygon(b, bc, ab));
				newPolys.Add(new Polygon(c, ca, bc));
				newPolys.Add(new Polygon(ab, bc, ca));
			}
			// Replace all our old polygons with the new set of
			// subdivided ones.
			m_Polygons = newPolys;
		}
	}
	private int GetMidPointIndex(Dictionary<int, int> cache, int indexA, int indexB)
	{
		// We create a key out of the two original indices
		// by storing the smaller index in the upper two bytes
		// of an integer, and the larger index in the lower two
		// bytes. By sorting them according to whichever is smaller
		// we ensure that this function returns the same result
		// whether you call
		// GetMidPointIndex(cache, 5, 9)
		// or...
		// GetMidPointIndex(cache, 9, 5)

		int smallerIndex = Mathf.Min(indexA, indexB);
		int greaterIndex = Mathf.Max(indexA, indexB);
		int key = (smallerIndex << 16) + greaterIndex;

		// If a midpoint is already defined, just return it.

		int ret;
		if (cache.TryGetValue(key, out ret))
			return ret;

		// If we're here, it's because a midpoint for these two
		// vertices hasn't been created yet. Let's do that now!

		Vector3 p1 = m_Vertices[indexA];
		Vector3 p2 = m_Vertices[indexB];
		Vector3 middle = Vector3.Lerp(p1, p2, 0.5f).normalized;

		ret = m_Vertices.Count;
		m_Vertices.Add(middle);

		// Add our new midpoint to the cache so we don't have
		// to do this again. =)

		cache.Add(key, ret);
		return ret;
	}
	#endregion
}
