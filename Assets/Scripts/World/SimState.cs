using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

public class SimState {

	public int Count;
	public float[] Elevation;
	public float[] WaterElevation;
	public float[] Vegetation;
	public float[] Ice;
	public float[] CloudCoverage;
	public float[] CloudElevation;
	public float[] RelativeHumidity;

	public void Init(int count)
	{
		Count = count;
		Elevation = new float[Count];
		WaterElevation = new float[Count];
		Vegetation = new float[Count];
		Ice = new float[Count];
		CloudCoverage = new float[Count];
		CloudElevation = new float[Count];
		RelativeHumidity = new float[Count];
	}
}
