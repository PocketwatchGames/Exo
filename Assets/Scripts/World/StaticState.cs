using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;

public class StaticState {

	public int Count;
	public float Radius;
	public float2[] Coordinate;
	public float3[] SphericalPosition;
	public int[] Neighbors;

	public void Init(int count)
	{
		Count = count;
		Neighbors = new int[count * 6];
		Coordinate = new float2[count];
		SphericalPosition = new float3[count];
	}

	public void ExtractCoordinates(List<UnityEngine.Vector3> vertices)
	{
		for (int i=0;i<Count;i++)
		{
			var v = vertices[i];
			Coordinate[i] = new float2(math.asin(v.x), math.asin(v.y));
			SphericalPosition[i] = new float3(v.x, v.y, v.z);
		}
	}
}
