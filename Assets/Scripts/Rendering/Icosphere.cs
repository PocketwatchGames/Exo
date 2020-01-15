﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

public class Icosphere {
	public List<Polygon> Polygons { get; private set; }
	public List<Vector3> Vertices { get; private set; }

	public Icosphere(int subdivisions)
	{
		Polygons = new List<Polygon>();
		Vertices = new List<Vector3>();

		// An icosahedron has 12 vertices, and
		// since they're completely symmetrical the
		// formula for calculating them is kind of
		// symmetrical too:

		float t = (1.0f + Mathf.Sqrt(5.0f)) / 2.0f;

		Vertices.Add(new Vector3(-1, t, 0).normalized);
		Vertices.Add(new Vector3(1, t, 0).normalized);
		Vertices.Add(new Vector3(-1, -t, 0).normalized);
		Vertices.Add(new Vector3(1, -t, 0).normalized);
		Vertices.Add(new Vector3(0, -1, t).normalized);
		Vertices.Add(new Vector3(0, 1, t).normalized);
		Vertices.Add(new Vector3(0, -1, -t).normalized);
		Vertices.Add(new Vector3(0, 1, -t).normalized);
		Vertices.Add(new Vector3(t, 0, -1).normalized);
		Vertices.Add(new Vector3(t, 0, 1).normalized);
		Vertices.Add(new Vector3(-t, 0, -1).normalized);
		Vertices.Add(new Vector3(-t, 0, 1).normalized);

		// And here's the formula for the 20 sides,
		// referencing the 12 vertices we just created.

		Polygons.Add(new Polygon(0, 11, 5));
		Polygons.Add(new Polygon(0, 5, 1));
		Polygons.Add(new Polygon(0, 1, 7));
		Polygons.Add(new Polygon(0, 7, 10));
		Polygons.Add(new Polygon(0, 10, 11));
		Polygons.Add(new Polygon(1, 5, 9));
		Polygons.Add(new Polygon(5, 11, 4));
		Polygons.Add(new Polygon(11, 10, 2));
		Polygons.Add(new Polygon(10, 7, 6));
		Polygons.Add(new Polygon(7, 1, 8));
		Polygons.Add(new Polygon(3, 9, 4));
		Polygons.Add(new Polygon(3, 4, 2));
		Polygons.Add(new Polygon(3, 2, 6));
		Polygons.Add(new Polygon(3, 6, 8));
		Polygons.Add(new Polygon(3, 8, 9));
		Polygons.Add(new Polygon(4, 9, 5));
		Polygons.Add(new Polygon(2, 4, 11));
		Polygons.Add(new Polygon(6, 2, 10));
		Polygons.Add(new Polygon(8, 6, 7));
		Polygons.Add(new Polygon(9, 8, 1));

		Subdivide(subdivisions);
	}

	private void Subdivide(int recursions)
	{
		var midPointCache = new Dictionary<int, int>();

		for (int i = 0; i < recursions; i++)
		{
			var newPolys = new List<Polygon>();
			foreach (var poly in Polygons)
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
			Polygons = newPolys;
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

		Vector3 p1 = Vertices[indexA];
		Vector3 p2 = Vertices[indexB];
		Vector3 middle = Vector3.Lerp(p1, p2, 0.5f).normalized;

		ret = Vertices.Count;
		Vertices.Add(middle);

		// Add our new midpoint to the cache so we don't have
		// to do this again. =)

		cache.Add(key, ret);
		return ret;
	}

}
