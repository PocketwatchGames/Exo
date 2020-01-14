using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Entities;

public struct SimState {

	public int Count;
	public int Ticks;
	public float Gravity;
	public float SpinAngle;
	public float SpinSpeed;
	public float OrbitSpeed;
	public float TiltAngle;
	public SimStateCell[] Cells;

	public void Init(int count)
	{
		Count = count;
		Cells = new SimStateCell[Count];
	}
}

public struct SimStateCell {
	public float Elevation;
	public float WaterElevation;
	public float Vegetation;
	public float Ice;
	public float CloudCoverage;
	public float CloudElevation;
	public float RelativeHumidity;
}
