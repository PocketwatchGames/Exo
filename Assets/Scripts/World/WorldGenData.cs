using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

[Serializable]
public class WorldGenData {
	public float Radius;
	public float MinElevation;
	public float MaxElevation;
	public float StratosphereMass;
	public float TroposphereMass;
	public float MinTemperature;
	public float MaxTemperature;
	public float Gravity;
	public float CarbonDioxide;
	public float SolarRadiation; // extraterrestrial solar radiation // https://en.wikipedia.org/wiki/Sunlight (1367 w/m^2)
	public float TiltAngle;
	public float SpinAngle;
	public float SpinSpeed;
	public float GeothermalHeat;
	public float MaxSalinity;
	public float MinSalinity;
	public int NumPlates;

}
