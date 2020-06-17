using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

public struct DiffusionAir {
	public float Temperature;
	public float WaterVapor;
	public float Dust;
	public float Minerals;
	public float CarbonDioxide;
	public float Nitrogen;
	public float Oxygen;
	public float Methane;
	public float3 Velocity;
}
public struct DiffusionCloud {
	public float Mass;
	public float Temperature;
	public float DropletMass;
}
public struct DiffusionWater {
	public float WaterMass;
	public float SaltMass;
	public float CarbonMass;
	public float NitrogenMass;
	public float GlucoseMass;
	public float MineralMass;
	public float OxygenMass;
	public float Temperature;
	public float3 Velocity;
}
public struct DiffusionLava {
	public float Mass;
	public float Temperature;
}

[BurstCompile]
public struct UpdateAirVelocityJob : IJobParallelFor {
	public NativeArray<float3> AirVelocity;
	public NativeArray<float> AirMovementVertical;
	[ReadOnly] public NativeArray<float> WindFriction;
	[ReadOnly] public NativeArray<float3> PressureGradientForce;
	[ReadOnly] public NativeArray<float3> Position;
	[ReadOnly] public float WindFrictionMultiplier;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float3 wind = AirVelocity[i];
		wind *= (1.0f - WindFriction[i] * WindFrictionMultiplier);
		wind += PressureGradientForce[i] * SecondsPerTick;
		AirVelocity[i] = wind;

		AirMovementVertical[i] = math.dot(Position[i], AirVelocity[i]) * SecondsPerTick;
	}
}


[BurstCompile]
public struct AdvectionAirJob : IJobParallelFor {
	public NativeSlice<DiffusionAir> Delta;
	[ReadOnly] public NativeSlice<float> Temperature;
	[ReadOnly] public NativeSlice<float> AirMass;
	[ReadOnly] public NativeSlice<float> Vapor;
	[ReadOnly] public NativeSlice<float> CarbonDioxide;
	[ReadOnly] public NativeSlice<float> Dust;
	[ReadOnly] public NativeSlice<float3> Velocity;
	[ReadOnly] public NativeSlice<float3> Positions;
	[ReadOnly] public NativeSlice<float> Destination;
	[ReadOnly] public NativeSlice<int> NeighborsVert;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public int CellsPerLayer;
	public void Execute(int i)
	{
		Debug.Assert(CellsPerLayer > 0);
		int columnIndex = i % CellsPerLayer;
		int fullRangeIndex = i + CellsPerLayer;

		float newTemperature = 0;
		float newWaterVapor = 0;
		float newCO2 = 0;
		float newDust = 0;
		float3 newVelocity = 0;

		// TODO: remove this when we have incompressibility
		float totalMass = 0;

		float massRemaining = AirMass[fullRangeIndex];
		for (int j=0;j<StaticState.MaxNeighborsVert; j++)
		{
			massRemaining -= math.max(0, Destination[fullRangeIndex * StaticState.MaxNeighborsVert + j]);
		}

		float remaining = massRemaining / AirMass[fullRangeIndex]; 
		totalMass += massRemaining;
		newTemperature += Temperature[fullRangeIndex] * massRemaining;
		newVelocity += Velocity[fullRangeIndex] * massRemaining;
		newWaterVapor += Vapor[fullRangeIndex] * remaining;
		newCO2 += CarbonDioxide[fullRangeIndex] * remaining;
		newDust += Dust[fullRangeIndex] * remaining;

		// TODO: subtract vertical motion first before applying horizontal motion (right now it adds up to more than 1)
		// Wait, does it? is this comment outdated?
		for (int j = 0; j < StaticState.MaxNeighborsVert; j++)
		{
			// TODO: this removes all vert motion, fix vert motion and remove!
			if (j == StaticState.NeighborDown || j == StaticState.NeighborUp)
			{
				continue;
			}
			int n = NeighborsVert[fullRangeIndex * StaticState.MaxNeighborsVert + j];
			if (n >= 0)
			{
				float incomingMass = math.max(0, -Destination[fullRangeIndex * StaticState.MaxNeighborsVert + j]);
				float incoming = incomingMass / AirMass[n];
				totalMass += incomingMass;
				newTemperature += Temperature[n] * incomingMass;
				newWaterVapor += Vapor[n] * incoming;
				newCO2 += CarbonDioxide[n] * incoming;
				newDust += Dust[n] * incoming;

				// TODO: this is increasing speed, is that right???  Shouldnt it only rotate?
				var deflectedVelocity = Velocity[n];

				int nColumnIndex = n % CellsPerLayer;
				if (j != StaticState.NeighborDown && j != StaticState.NeighborUp) {
					deflectedVelocity += math.cross(Positions[nColumnIndex], Velocity[n]) * CoriolisMultiplier[nColumnIndex] * CoriolisTerm * SecondsPerTick;

				}
				// TODO: turn velocity along great circle, instead of just erasing the vertical component as we are doing here
				//var deflectedVertical = math.dot(deflectedVelocity, Positions[nColumnIndex]);
				//deflectedVelocity -= Positions[nColumnIndex] * deflectedVertical;
				deflectedVelocity -= Positions[columnIndex] * math.dot(Positions[columnIndex], deflectedVelocity);
				//deflectedVelocity += deflectedVertical * Positions[columnIndex];


				// adjust vertical wind velocity when we hit a mountain or go down a valley
				// TODO: this should be based on the difference of the current velocity's slope and the terrain slope
				// probably apply this when we generate wind forces
				//deflectedVelocity += Positions[i] * (layerMiddle - LayerMiddle[n]) * TicksPerSecond;

				// TODO: this is temp
				// need to deal with centrifugal force/gravity so that as air moves horizontally, it can fly into the air or get pulled to earth
				//newVelocity += incomingMass * (deflectedVelocity - Positions[i] * math.dot(Positions[i], deflectedVelocity));

				newVelocity += incomingMass * deflectedVelocity;
			}
		}


		if (totalMass > 0)
		{
			newTemperature /= totalMass;
			newVelocity /= totalMass;
		}
		else
		{
			// TODO: remove once we have incompressibility
			newTemperature = Temperature[fullRangeIndex];
			newVelocity = Velocity[fullRangeIndex];
		}

		Delta[i] = new DiffusionAir()
		{
			Temperature = newTemperature,
			WaterVapor = newWaterVapor,
			CarbonDioxide = newCO2,
			Dust = newDust,
			Velocity = newVelocity,
		};



	}
}


[BurstCompile]
public struct AdvectionWaterJob : IJobParallelFor {
	public NativeSlice<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<float> Destination;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float> Carbon;
	[ReadOnly] public NativeArray<float> Nitrogen;
	[ReadOnly] public NativeArray<float> Mineral;
	[ReadOnly] public NativeArray<float> Oxygen;
	[ReadOnly] public NativeArray<float> Glucose;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<int> NeighborsVert;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public int CellsPerLayer;
	public void Execute(int i)
	{
		Debug.Assert(CellsPerLayer > 0);
		int columnIndex = i % CellsPerLayer;
		int fullRangeIndex = i + CellsPerLayer;

		float waterMass = Mass[fullRangeIndex];

		float newMass = waterMass;
		float newSaltMass = 0;
		float newCarbon = 0;
		float newGlucose = 0;
		float newNitrogen = 0;
		float newMineral = 0;
		float newOxygen = 0;
		float newTemperature = 0;
		float3 newVelocity = 0;

		if (waterMass > 0)
		{
			for (int j = 0; j < StaticState.MaxNeighborsVert; j++)
			{
				newMass -= math.max(0, Destination[fullRangeIndex * StaticState.MaxNeighborsVert + j]);
			}
			newMass = math.max(0, newMass);

			float percentRemaining = newMass / waterMass;

			newSaltMass = Salt[fullRangeIndex] * percentRemaining;
			newCarbon = Carbon[fullRangeIndex] * percentRemaining;
			newNitrogen = Nitrogen[fullRangeIndex] * percentRemaining;
			newGlucose = Glucose[fullRangeIndex] * percentRemaining;
			newMineral = Mineral[fullRangeIndex] * percentRemaining;
			newOxygen = Oxygen[fullRangeIndex] * percentRemaining;
			newTemperature = Temperature[fullRangeIndex] * newMass;
			newVelocity = Velocity[fullRangeIndex] * newMass;

			for (int j = 0; j < StaticState.MaxNeighborsVert; j++)
			{
				// TODO: this removes all vert motion, fix vert motion and remove!
				if (j == StaticState.NeighborDown || j == StaticState.NeighborUp)
				{
					continue;
				}
				int edgeIndex = fullRangeIndex * StaticState.MaxNeighborsVert + j;
				int n = NeighborsVert[edgeIndex];
				if (n >= 0)
				{
					int nColumnIndex = n % CellsPerLayer;
					float incomingMass = math.max(0, -Destination[edgeIndex]);
					float incoming;
					if (Mass[n] > 0)
					{
						incoming = incomingMass / Mass[n];
					} else
					{
						incoming = 0;
					}

					newMass += incomingMass;
					newSaltMass += Salt[n] * incoming;
					newCarbon += Carbon[n] * incoming;
					newNitrogen += Nitrogen[n] * incoming;
					newMineral += Mineral[n] * incoming;
					newOxygen += Oxygen[n] * incoming;
					newGlucose += Glucose[n] * incoming;
					newTemperature += Temperature[n] * incomingMass;

					// TODO: this is increasing speed, is that right???  Shouldnt it only rotate?
					var deflectedVelocity = Velocity[n];

					if (j != StaticState.NeighborDown && j != StaticState.NeighborUp)
					{
						deflectedVelocity += math.cross(Positions[nColumnIndex], Velocity[n]) * CoriolisMultiplier[nColumnIndex] * CoriolisTerm * SecondsPerTick;
					}

					// TODO: turn velocity along great circle, instead of just erasing the vertical component as we are doing here
					//var deflectedVertical = math.dot(deflectedVelocity, Positions[nColumnIndex]);
					//deflectedVelocity -= Positions[nColumnIndex] * deflectedVertical;
					deflectedVelocity -= Positions[columnIndex] * math.dot(Positions[columnIndex], deflectedVelocity);
					//deflectedVelocity += deflectedVertical * Positions[columnIndex];

					// TODO: adjust vertical wind velocity when we hit a mountain or go down a valley
					// deflectedVelocity += (LayerMiddle[n] - layerMiddle) * NeighborDistInverse[i * 6 + j] * TicksPerSecond;

					// TODO: this is temp
					// need to deal with centrifugal force/gravity so that as air moves horizontally, it can fly into the air or get pulled to earth
					//newVelocity += incomingMass * (deflectedVelocity - Positions[i] * math.dot(Positions[i], deflectedVelocity));

					newVelocity += incomingMass * deflectedVelocity;

				}
			}

			// TODO: remove this when we have incompressibility??
			// OR NOT? water is compressible!
			if (newMass > 0)
			{
				float inverseNewMass = 1.0f / newMass;
				newTemperature *= inverseNewMass;
				newVelocity *= inverseNewMass;
				newSaltMass *= inverseNewMass * waterMass;
				newCarbon *= inverseNewMass * waterMass;
				newNitrogen *= inverseNewMass * waterMass;
				newGlucose *= inverseNewMass * waterMass;
				newMineral *= inverseNewMass * waterMass;
				newOxygen *= inverseNewMass * waterMass;
				newMass = waterMass;
			}
			else
			{
				// TODO: remove once we have incompressibility
				newTemperature = Temperature[fullRangeIndex];
				newVelocity = Velocity[fullRangeIndex];
			}
		}

		Delta[i] = new DiffusionWater()
		{
			WaterMass = newMass,
			SaltMass = newSaltMass,
			CarbonMass = newCarbon,
			NitrogenMass = newNitrogen,
			GlucoseMass = newGlucose,
			MineralMass = newMineral,
			OxygenMass = newOxygen,
			Temperature = newTemperature,
			Velocity = newVelocity
		};

	}
}

[BurstCompile]
public struct AdvectionCloudJob : IJobParallelFor {
	public NativeArray<DiffusionCloud> Delta;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> DropletMass;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<BarycentricValue> Destination;
	public void Execute(int i)
	{
		float newMass = 0;
		float newTemperature = 0;
		float newDropletMass = 0;

		float massPercentRemaining = 0;
		if (Destination[i].indexA == i || Destination[i].indexA < 0)
		{
			massPercentRemaining += Destination[i].valueA;
		}
		if (Destination[i].indexB == i || Destination[i].indexB < 0)
		{
			massPercentRemaining += Destination[i].valueB;
		}
		if (Destination[i].indexC == i || Destination[i].indexC < 0)
		{
			massPercentRemaining += Destination[i].valueC;
		}
		float m = Mass[i] * massPercentRemaining;
		newMass += m;
		newTemperature += Temperature[i] * m;
		newDropletMass += DropletMass[i] * m;

		for (int j = 0; j < StaticState.MaxNeighbors; j++)
		{
			int n = Neighbors[i * StaticState.MaxNeighbors + j];
			if (n >= 0)
			{
				float incoming = 0;
				if (Destination[n].indexA == i)
				{
					incoming += Destination[n].valueA;
				}
				else if (Destination[n].indexB == i)
				{
					incoming += Destination[n].valueB;
				}
				else if (Destination[n].indexC == i)
				{
					incoming += Destination[n].valueC;
				}
				float massIncoming = Mass[n] * incoming;
				newMass += massIncoming;
				newTemperature += Temperature[n] * massIncoming;
				newDropletMass += DropletMass[n] * massIncoming;
			}
		}

		if (newMass > 0)
		{
			float inverseMass = 1.0f / newMass;
			newTemperature *= inverseMass;
			newDropletMass *= inverseMass;
		}

		Delta[i] = new DiffusionCloud()
		{
			Mass = newMass,
			Temperature = newTemperature,
			DropletMass = newDropletMass,
		};
	}
}


[BurstCompile]
public struct ApplyAdvectionAirJob : IJobParallelFor {
	public NativeSlice<float> Temperature;
	public NativeSlice<float> Vapor;
	public NativeSlice<float> Dust;
	public NativeSlice<float> CarbonDioxide;
	public NativeSlice<float> Nitrogen;
	public NativeSlice<float> Oxygen;
	public NativeSlice<float> Minerals;
	public NativeSlice<float> Methane;
	public NativeSlice<float3> AirVelocity;
	[ReadOnly] public NativeSlice<DiffusionAir> Advection;
	public void Execute(int i)
	{
		Temperature[i] = Advection[i].Temperature;
		Vapor[i] = Advection[i].WaterVapor;
		Dust[i] = Advection[i].Dust;
		CarbonDioxide[i] = Advection[i].CarbonDioxide;
		Nitrogen[i] = Advection[i].Nitrogen;
		Oxygen[i] = Advection[i].Oxygen;
		Minerals[i] = Advection[i].Minerals;
		Methane[i] = Advection[i].Methane;
		AirVelocity[i] = Advection[i].Velocity;
	}
}


[BurstCompile]
public struct ApplyAdvectionWaterJob : IJobParallelFor {
	public NativeSlice<float> WaterMass;
	public NativeSlice<float> SaltMass;
	public NativeSlice<float> CarbonMass;
	public NativeSlice<float> NitrogenMass;
	public NativeSlice<float> GlucoseMass;
	public NativeSlice<float> MineralMass;
	public NativeSlice<float> OxygenMass;
	public NativeSlice<float> Temperature;
	public NativeSlice<float3> Velocity;
	[ReadOnly] public NativeSlice<DiffusionWater> Advection;
	public void Execute(int i)
	{
		WaterMass[i] = Advection[i].WaterMass;
		SaltMass[i] = Advection[i].SaltMass;
		CarbonMass[i] = Advection[i].CarbonMass;
		NitrogenMass[i] = Advection[i].NitrogenMass;
		GlucoseMass[i] = Advection[i].GlucoseMass;
		MineralMass[i] = Advection[i].MineralMass;
		OxygenMass[i] = Advection[i].OxygenMass;
		Temperature[i] = Advection[i].Temperature;
		Velocity[i] = Advection[i].Velocity;
	}
}


[BurstCompile]
public struct ApplyAdvectionCloudJob : IJobParallelFor {
	public NativeArray<float> CloudMass;
	public NativeArray<float> DropletMass;
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<DiffusionCloud> Advection;
	public void Execute(int i)
	{
		CloudMass[i] = Advection[i].Mass;
		Temperature[i] = Advection[i].Temperature;
		DropletMass[i] = Advection[i].DropletMass;
	}
}


[BurstCompile]
public struct ApplyAdvectionLavaJob : IJobParallelFor {
	public NativeArray<float> Mass;
	public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<DiffusionLava> Advection;
	public void Execute(int i)
	{
		Mass[i] = Advection[i].Mass;
		Temperature[i] = Advection[i].Temperature;
	}
}


