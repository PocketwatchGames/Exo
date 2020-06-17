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
	public NativeArray<float> GroundNitrogen;
	public NativeArray<float> GroundGlucose;
	public NativeArray<float> GroundMinerals;
	public NativeArray<float> GroundSalt;
	public NativeArray<float> GroundTemperature;
	public NativeArray<float> GroundWater;
	public NativeArray<float> GroundWaterTemperature;
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
	public NativeArray<float> AirCarbonDioxide;
	public NativeArray<float> AirNitrogen;
	public NativeArray<float> AirMethane;
	public NativeArray<float> AirOxygen;
	public NativeArray<float3> AirVelocity;
	public NativeArray<float> AirDust;
	public NativeArray<float> AirMinerals;
	public NativeArray<float> WaterTemperature;
	public NativeArray<float> WaterMass;
	public NativeArray<float> WaterCarbonDioxide;
	public NativeArray<float> WaterGlucose;
	public NativeArray<float> WaterMinerals;
	public NativeArray<float> WaterOxygen;
	public NativeArray<float> WaterNitrogen;
	public NativeArray<float3> WaterVelocity;
	public NativeArray<float> WaterSaltMass;

	public void Init(int count, ref WorldData worldData)
	{
		Plate = new NativeArray<short>(count, Allocator.Persistent);
		Roughness = new NativeArray<float>(count, Allocator.Persistent);
		Elevation = new NativeArray<float>(count, Allocator.Persistent);
		GroundTemperature = new NativeArray<float>(count, Allocator.Persistent);
		GroundWater = new NativeArray<float>(count, Allocator.Persistent);
		GroundWaterTemperature = new NativeArray<float>(count, Allocator.Persistent);
		GroundNitrogen = new NativeArray<float>(count, Allocator.Persistent);
		GroundGlucose = new NativeArray<float>(count, Allocator.Persistent);
		GroundSalt = new NativeArray<float>(count, Allocator.Persistent);
		GroundMinerals = new NativeArray<float>(count, Allocator.Persistent);
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
		AirCarbonDioxide = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirOxygen = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirNitrogen = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirMethane = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirMinerals = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirDust = new NativeArray<float>(count * worldData.AirLayers, Allocator.Persistent);
		AirVelocity = new NativeArray<float3>(count * worldData.AirLayers, Allocator.Persistent);

		WaterTemperature = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterMass = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterCarbonDioxide = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterSaltMass = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterGlucose = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterMinerals = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterOxygen = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterNitrogen = new NativeArray<float>(count * worldData.WaterLayers, Allocator.Persistent);
		WaterVelocity = new NativeArray<float3>(count * worldData.WaterLayers, Allocator.Persistent);
	}

	public void CopyFrom(ref SimState from)
	{
		PlanetState = from.PlanetState;
		Plate.CopyFrom(from.Plate);
		Roughness.CopyFrom(from.Roughness);
		Elevation.CopyFrom(from.Elevation);
		GroundTemperature.CopyFrom(from.GroundTemperature);
		GroundWater.CopyFrom(from.GroundWater);
		GroundWaterTemperature.CopyFrom(from.GroundWaterTemperature);
		GroundMinerals.CopyFrom(from.GroundMinerals);
		GroundNitrogen.CopyFrom(from.GroundNitrogen);
		GroundGlucose.CopyFrom(from.GroundGlucose);
		GroundSalt.CopyFrom(from.GroundSalt);
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
		AirCarbonDioxide.CopyFrom(from.AirCarbonDioxide);
		AirOxygen.CopyFrom(from.AirOxygen);
		AirNitrogen.CopyFrom(from.AirNitrogen);
		AirMethane.CopyFrom(from.AirMethane);
		AirMinerals.CopyFrom(from.AirMinerals);
		AirDust.CopyFrom(from.AirDust);
		AirVelocity.CopyFrom(from.AirVelocity);

		WaterTemperature.CopyFrom(from.WaterTemperature);
		WaterMass.CopyFrom(from.WaterMass);
		WaterCarbonDioxide.CopyFrom(from.WaterCarbonDioxide);
		WaterSaltMass.CopyFrom(from.WaterSaltMass);
		WaterMinerals.CopyFrom(from.WaterMinerals);
		WaterGlucose.CopyFrom(from.WaterGlucose);
		WaterOxygen.CopyFrom(from.WaterOxygen);
		WaterVelocity.CopyFrom(from.WaterVelocity);
		WaterNitrogen.CopyFrom(from.WaterNitrogen);
	}

	public void Dispose()
	{
		Plate.Dispose();
		Roughness.Dispose();
		Elevation.Dispose();
		GroundTemperature.Dispose();
		GroundWater.Dispose();
		GroundWaterTemperature.Dispose();
		GroundGlucose.Dispose();
		GroundNitrogen.Dispose();
		GroundMinerals.Dispose();
		GroundSalt.Dispose();
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
		AirCarbonDioxide.Dispose();
		AirOxygen.Dispose();
		AirVelocity.Dispose();
		AirDust.Dispose();
		AirMinerals.Dispose();
		AirNitrogen.Dispose();
		AirMethane.Dispose();

		WaterTemperature.Dispose();
		WaterMass.Dispose();
		WaterCarbonDioxide.Dispose();
		WaterSaltMass.Dispose();
		WaterGlucose.Dispose();
		WaterNitrogen.Dispose();
		WaterMinerals.Dispose();
		WaterOxygen.Dispose();
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
}

