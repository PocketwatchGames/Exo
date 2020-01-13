using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

public class SimState {

	public int Count;
	public float[] Elevation;
	public float[] Vegetation;

	public void Init(int count)
	{
		Count = count;
		Elevation = new float[Count];
		Vegetation = new float[Count];
	}
}
