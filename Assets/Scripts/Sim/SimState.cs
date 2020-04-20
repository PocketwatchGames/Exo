using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Entities;
using Unity;
using UnityEngine;
using Unity.Mathematics;
using Unity.Collections;

public struct SimState {

	public PlanetState PlanetState;
	public NativeArray<float> Elevation;
	public NativeArray<float> Roughness;
	public NativeArray<float> GroundCarbon;
	public NativeArray<float> GroundTemperature;
	public NativeArray<float> GroundWater;
	public NativeArray<float> GroundWaterTemperature;
	public NativeArray<float> FloraMass;
	public NativeArray<float> FloraWater;
	public NativeArray<float> FloraTemperature;
	public NativeArray<float> FloraGlucose;
	public NativeArray<float> IceTemperature;
	public NativeArray<float> IceMass;
	public NativeArray<float> CloudMass;
	public NativeArray<float> CloudTemperature;
	public NativeArray<float> CloudDropletMass;
	public NativeArray<float> CrustDepth;
	public NativeArray<float> MagmaMass;
	public NativeArray<float> LavaMass;
	public NativeArray<float> LavaTemperature;
	public NativeArray<float> FlowWater;
	public NativeArray<float> FlowLava;
	public NativeArray<float>[] AirTemperaturePotential;
	public NativeArray<float>[] AirVapor;
	public NativeArray<float>[] AirCarbon;
	public NativeArray<float3>[] AirVelocity;
	public NativeArray<float>[] Dust;
	public NativeArray<float>[] WaterTemperature;
	public NativeArray<float>[] WaterMass;
	public NativeArray<float>[] WaterCarbon;
	public NativeArray<float3>[] WaterVelocity;
	public NativeArray<float>[] SaltMass;
	public NativeArray<float>[] PlanktonMass;
	public NativeArray<float>[] PlanktonGlucose;

	public void Init(int count, ref WorldData worldData)
	{
		Roughness = new NativeArray<float>(count, Allocator.Persistent);
		GroundCarbon = new NativeArray<float>(count, Allocator.Persistent);
		Elevation = new NativeArray<float>(count, Allocator.Persistent);
		GroundTemperature = new NativeArray<float>(count, Allocator.Persistent);
		GroundWater = new NativeArray<float>(count, Allocator.Persistent);
		GroundWaterTemperature = new NativeArray<float>(count, Allocator.Persistent);
		FloraMass = new NativeArray<float>(count, Allocator.Persistent);
		FloraWater = new NativeArray<float>(count, Allocator.Persistent);
		FloraTemperature = new NativeArray<float>(count, Allocator.Persistent);
		FloraGlucose = new NativeArray<float>(count, Allocator.Persistent);
		IceTemperature = new NativeArray<float>(count, Allocator.Persistent);
		IceMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudMass = new NativeArray<float>(count, Allocator.Persistent);
		CloudTemperature = new NativeArray<float>(count, Allocator.Persistent);
		CloudDropletMass = new NativeArray<float>(count, Allocator.Persistent);
		LavaTemperature = new NativeArray<float>(count, Allocator.Persistent);
		LavaMass = new NativeArray<float>(count, Allocator.Persistent);
		CrustDepth = new NativeArray<float>(count, Allocator.Persistent);
		MagmaMass = new NativeArray<float>(count, Allocator.Persistent);

		FlowWater = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.Persistent);
		FlowLava = new NativeArray<float>(count * StaticState.MaxNeighbors, Allocator.Persistent);

		AirTemperaturePotential = new NativeArray<float>[worldData.AirLayers];
		AirVapor = new NativeArray<float>[worldData.AirLayers];
		AirCarbon = new NativeArray<float>[worldData.AirLayers];
		AirVelocity = new NativeArray<float3>[worldData.AirLayers];
		Dust = new NativeArray<float>[worldData.AirLayers];
		for (int i = 0; i < worldData.AirLayers; i++)
		{
			AirTemperaturePotential[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVapor[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirCarbon[i] = new NativeArray<float>(count, Allocator.Persistent);
			AirVelocity[i] = new NativeArray<float3>(count, Allocator.Persistent);
			Dust[i] = new NativeArray<float>(count, Allocator.Persistent);
		}

		WaterTemperature = new NativeArray<float>[worldData.WaterLayers];
		WaterMass = new NativeArray<float>[worldData.WaterLayers];
		SaltMass = new NativeArray<float>[worldData.WaterLayers];
		WaterCarbon = new NativeArray<float>[worldData.WaterLayers];
		PlanktonMass = new NativeArray<float>[worldData.WaterLayers];
		PlanktonGlucose = new NativeArray<float>[worldData.WaterLayers];
		WaterVelocity = new NativeArray<float3>[worldData.WaterLayers];
		for (int i = 0; i < worldData.WaterLayers; i++)
		{
			WaterTemperature[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterCarbon[i] = new NativeArray<float>(count, Allocator.Persistent);
			PlanktonMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			PlanktonGlucose[i] = new NativeArray<float>(count, Allocator.Persistent);
			SaltMass[i] = new NativeArray<float>(count, Allocator.Persistent);
			WaterVelocity[i] = new NativeArray<float3>(count, Allocator.Persistent);
		}
	}

	public void CopyFrom(ref SimState from)
	{
		PlanetState = from.PlanetState;
		Roughness.CopyFrom(from.Roughness);
		GroundCarbon.CopyFrom(from.GroundCarbon);
		Elevation.CopyFrom(from.Elevation);
		GroundTemperature.CopyFrom(from.GroundTemperature);
		GroundWater.CopyFrom(from.GroundWater);
		GroundWaterTemperature.CopyFrom(from.GroundWaterTemperature);
		FloraMass.CopyFrom(from.FloraMass);
		FloraWater.CopyFrom(from.FloraWater);
		FloraTemperature.CopyFrom(from.FloraTemperature);
		FloraGlucose.CopyFrom(from.FloraGlucose);
		IceTemperature.CopyFrom(from.IceTemperature);
		IceMass.CopyFrom(from.IceMass);
		CloudMass.CopyFrom(from.CloudMass);
		CloudTemperature.CopyFrom(from.CloudTemperature);
		CloudDropletMass.CopyFrom(from.CloudDropletMass);
		LavaTemperature.CopyFrom(from.LavaTemperature);
		LavaMass.CopyFrom(from.LavaMass);
		CrustDepth.CopyFrom(from.CrustDepth);
		MagmaMass.CopyFrom(from.MagmaMass);
		FlowWater.CopyFrom(from.FlowWater);
		FlowLava.CopyFrom(from.FlowLava);
		for (int i = 0; i < from.AirTemperaturePotential.Length; i++)
		{
			AirTemperaturePotential[i].CopyFrom(from.AirTemperaturePotential[i]);
			AirVapor[i].CopyFrom(from.AirVapor[i]);
			AirCarbon[i].CopyFrom(from.AirCarbon[i]);
			AirVelocity[i].CopyFrom(from.AirVelocity[i]);
			Dust[i].CopyFrom(from.Dust[i]);
		}

		for (int i = 0; i < from.WaterTemperature.Length; i++)
		{
			WaterTemperature[i].CopyFrom(from.WaterTemperature[i]);
			WaterMass[i].CopyFrom(from.WaterMass[i]);
			WaterCarbon[i].CopyFrom(from.WaterCarbon[i]);
			PlanktonMass[i].CopyFrom(from.PlanktonMass[i]);
			PlanktonGlucose[i].CopyFrom(from.PlanktonGlucose[i]);
			SaltMass[i].CopyFrom(from.SaltMass[i]);
			WaterVelocity[i].CopyFrom(from.WaterVelocity[i]);
		}
	}

	public void Dispose()
	{
		Roughness.Dispose();
		GroundCarbon.Dispose();
		Elevation.Dispose();
		GroundTemperature.Dispose();
		FloraMass.Dispose();
		FloraWater.Dispose();
		FloraGlucose.Dispose();
		FloraTemperature.Dispose();
		GroundWater.Dispose();
		GroundWaterTemperature.Dispose();
		IceTemperature.Dispose();
		IceMass.Dispose();
		CloudMass.Dispose();
		CloudTemperature.Dispose();
		CloudDropletMass.Dispose();
		LavaMass.Dispose();
		LavaTemperature.Dispose();
		MagmaMass.Dispose();
		CrustDepth.Dispose();
		FlowWater.Dispose();
		FlowLava.Dispose();
		for (int i = 0; i < AirTemperaturePotential.Length; i++)
		{
			AirTemperaturePotential[i].Dispose();
			AirVapor[i].Dispose();
			AirCarbon[i].Dispose();
			AirVelocity[i].Dispose();
			Dust[i].Dispose();
		}

		for (int i = 0; i < WaterTemperature.Length; i++)
		{
			WaterTemperature[i].Dispose();
			WaterMass[i].Dispose();
			PlanktonMass[i].Dispose();
			PlanktonGlucose[i].Dispose();
			WaterCarbon[i].Dispose();
			SaltMass[i].Dispose();
			WaterVelocity[i].Dispose();
		}
	}
}

public struct PlanetState {
	public int Ticks;
	public float Gravity;
	public float SpinSpeed;
	public float OrbitSpeed;
	public float GeothermalHeat;
	public float SolarRadiation;
	public float DistanceToSun;
	public float3 Rotation;
	public float3 Position;
	public float AngularSpeed;
	public float Oxygen;
}

