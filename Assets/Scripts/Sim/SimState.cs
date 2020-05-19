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
	public NativeArray<short> Plate;
	public NativeArray<float> Elevation;
	public NativeArray<float> Roughness;
	public NativeArray<float> GroundCarbon;
	public NativeArray<float> GroundTemperature;
	public NativeArray<float> GroundWater;
	public NativeArray<float> GroundWaterTemperature;
	public NativeArray<float> FloraMass;
	public NativeArray<float> FloraWater;
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
	public NativeArray<float> AirTemperaturePotential;
	public NativeArray<float> AirVapor;
	public NativeArray<float> AirCarbon;
	public NativeArray<float3> AirVelocity;
	public NativeArray<float> Dust;
	public NativeArray<float> WaterTemperature;
	public NativeArray<float> WaterMass;
	public NativeArray<float> WaterCarbon;
	public NativeArray<float3> WaterVelocity;
	public NativeArray<float> SaltMass;
	public NativeArray<float> PlanktonMass;
	public NativeArray<float> PlanktonGlucose;

	public void Init(int count, ref WorldData worldData)
	{
		Plate = new NativeArray<short>(count, Allocator.Persistent);
		Roughness = new NativeArray<float>(count, Allocator.Persistent);
		GroundCarbon = new NativeArray<float>(count, Allocator.Persistent);
		Elevation = new NativeArray<float>(count, Allocator.Persistent);
		GroundTemperature = new NativeArray<float>(count, Allocator.Persistent);
		GroundWater = new NativeArray<float>(count, Allocator.Persistent);
		GroundWaterTemperature = new NativeArray<float>(count, Allocator.Persistent);
		FloraMass = new NativeArray<float>(count, Allocator.Persistent);
		FloraWater = new NativeArray<float>(count, Allocator.Persistent);
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

		AirTemperaturePotential = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirVapor = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirCarbon = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirVelocity = new NativeArray<float3>(count * worldData.AirLayers, Allocator.Persistent);
		Dust = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);

		WaterTemperature = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterMass = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterCarbon = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		PlanktonMass = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		PlanktonGlucose = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		SaltMass = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterVelocity = new NativeArray<float3>(count * worldData.WaterLayers, Allocator.Persistent);
	}

	public void CopyFrom(ref SimState from)
	{
		PlanetState = from.PlanetState;
		Plate.CopyFrom(from.Plate);
		Roughness.CopyFrom(from.Roughness);
		GroundCarbon.CopyFrom(from.GroundCarbon);
		Elevation.CopyFrom(from.Elevation);
		GroundTemperature.CopyFrom(from.GroundTemperature);
		GroundWater.CopyFrom(from.GroundWater);
		GroundWaterTemperature.CopyFrom(from.GroundWaterTemperature);
		FloraMass.CopyFrom(from.FloraMass);
		FloraWater.CopyFrom(from.FloraWater);
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

		AirTemperaturePotential.CopyFrom(from.AirTemperaturePotential);
		AirVapor.CopyFrom(from.AirVapor);
		AirCarbon.CopyFrom(from.AirCarbon);
		AirVelocity.CopyFrom(from.AirVelocity);
		Dust.CopyFrom(from.Dust);
		
		WaterTemperature.CopyFrom(from.WaterTemperature);
		WaterMass.CopyFrom(from.WaterMass);
		WaterCarbon.CopyFrom(from.WaterCarbon);
		PlanktonMass.CopyFrom(from.PlanktonMass);
		PlanktonGlucose.CopyFrom(from.PlanktonGlucose);
		SaltMass.CopyFrom(from.SaltMass);
		WaterVelocity.CopyFrom(from.WaterVelocity);
	}

	public void Dispose()
	{
		Plate.Dispose();
		Roughness.Dispose();
		GroundCarbon.Dispose();
		Elevation.Dispose();
		GroundTemperature.Dispose();
		FloraMass.Dispose();
		FloraWater.Dispose();
		FloraGlucose.Dispose();
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

		AirTemperaturePotential.Dispose();
		AirVapor.Dispose();
		AirCarbon.Dispose();
		AirVelocity.Dispose();
		Dust.Dispose();

		WaterTemperature.Dispose();
		WaterMass.Dispose();
		PlanktonMass.Dispose();
		PlanktonGlucose.Dispose();
		WaterCarbon.Dispose();
		SaltMass.Dispose();
		WaterVelocity.Dispose();
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

