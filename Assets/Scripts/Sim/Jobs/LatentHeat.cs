using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;


[BurstCompile]
public struct ApplyLatentHeatIceJob : IJobParallelFor {
	public NativeArray<float> IceTemperature;
	[ReadOnly] public NativeArray<float> IceMass;
	[ReadOnly] public NativeSlice<float> LatentHeat;
	public void Execute(int i)
	{
		if (IceMass[i] > 0)
		{
			IceTemperature[i] = math.max(0, IceTemperature[i] + LatentHeat[i] / (WorldData.SpecificHeatIce * IceMass[i]));
		}
	}
}
[BurstCompile]
public struct ApplyLatentHeatTerrainJob : IJobParallelFor {
	public NativeArray<float> TerrainTemperature;
	[ReadOnly] public NativeSlice<float> LatentHeat;
	[ReadOnly] public NativeArray<float> SpecificHeatTerrain;
	public void Execute(int i)
	{
		TerrainTemperature[i] = math.max(0, TerrainTemperature[i] + LatentHeat[i] / SpecificHeatTerrain[i]);
	}
}
[BurstCompile]
public struct ApplyLatentHeatFloraJob : IJobParallelFor {
	public NativeArray<float> FloraTemperature;
	[ReadOnly] public NativeArray<float> LatentHeat;
	[ReadOnly] public NativeArray<float> FloraWater;
	[ReadOnly] public NativeArray<float> FloraMass;
	public void Execute(int i)
	{
		float specificHeat = (FloraMass[i] * WorldData.SpecificHeatFlora + FloraWater[i] * WorldData.SpecificHeatWater);
		FloraTemperature[i] = math.max(0, FloraTemperature[i] + LatentHeat[i] / specificHeat);
	}
}

[BurstCompile]
public struct ApplyLatentHeatAirJob : IJobParallelFor {
	public NativeSlice<float> AirTemperaturePotential;
	[ReadOnly] public NativeSlice<float> LatentHeatAir;
	[ReadOnly] public NativeSlice<float> LatentHeatCloud;
	[ReadOnly] public NativeArray<float> CloudElevation;
	[ReadOnly] public NativeArray<float> CloudMass;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeSlice<float> VaporMass;
	[ReadOnly] public NativeSlice<float> LayerElevation;
	[ReadOnly] public NativeSlice<float> LayerHeight;
	[ReadOnly] public int Count;
	[ReadOnly] public int LayerCount;
	public void Execute(int i)
	{
		int columnIndex = i % Count;
		float cloudMass = Atmosphere.GetCloudMassInLayer(CloudMass[columnIndex], CloudElevation[columnIndex], LayerElevation[i], LayerHeight[i]);
		float latentHeatCloud = 0;
		if (CloudElevation[columnIndex] >= LayerElevation[i] && CloudElevation[columnIndex] < LayerElevation[i] + LayerHeight[i])
		{
			for (int j = 0; j < LayerCount; j++)
			{
				latentHeatCloud += LatentHeatCloud[columnIndex + j * Count];
			}
		}
		float sh = Atmosphere.GetSpecificHeatAir(
			AirMass[i],
			VaporMass[i],
			cloudMass);
		AirTemperaturePotential[i] += (LatentHeatAir[i] + latentHeatCloud) / sh;
	}
}

[BurstCompile]
public struct ApplyLatentHeatWaterJob : IJobParallelFor {
	public NativeSlice<float> WaterTemperature;
	[ReadOnly] public NativeSlice<float> WaterMass;
	[ReadOnly] public NativeSlice<float> SaltMass;
	[ReadOnly] public NativeSlice<float> LatentHeat;
	public void Execute(int i)
	{
		if (WaterMass[i] > 0)
		{
			float newTemp = math.max(0, WaterTemperature[i] + LatentHeat[i] / (WaterMass[i] * WorldData.SpecificHeatWater + SaltMass[i] * WorldData.SpecificHeatSalt));
			WaterTemperature[i] = newTemp;
		}
	}
}

[BurstCompile]
public struct ApplyLatentHeatLavaJob : IJobParallelFor {
	public NativeArray<float> LavaTemperature;
	[ReadOnly] public NativeSlice<float> LatentHeat;
	[ReadOnly] public NativeArray<float> LavaMass;
	public void Execute(int i)
	{
		if (LavaMass[i] > 0)
		{
			LavaTemperature[i] = math.max(0, LavaTemperature[i] + LatentHeat[i] / (LavaMass[i] * WorldData.SpecificHeatLava));
		}
	}
}
