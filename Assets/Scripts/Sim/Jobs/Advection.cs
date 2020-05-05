using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

public struct DiffusionAir {
	public float Temperature;
	public float WaterVapor;
	public float Dust;
	public float CarbonDioxide;
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
	public float Plankton;
	public float PlanktonGlucose;
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
		for (int j=0;j<StaticState.MaxNeighbors;j++)
		{
			massRemaining -= math.max(0, Destination[fullRangeIndex * StaticState.MaxNeighbors + j]);
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

				if (j != StaticState.NeighborDown && j != StaticState.NeighborUp) {
					int nColumnIndex = n % CellsPerLayer;
					deflectedVelocity += math.cross(Positions[nColumnIndex], Velocity[n]) * CoriolisMultiplier[nColumnIndex] * CoriolisTerm * SecondsPerTick;

					// TODO: turn velocity along great circle, instead of just erasing the vertical component as we are doing here
					var deflectedVertical = math.dot(deflectedVelocity, Positions[nColumnIndex]);
					deflectedVelocity -= Positions[nColumnIndex] * deflectedVertical;
					deflectedVelocity -= Positions[columnIndex] * math.dot(Positions[columnIndex], deflectedVelocity);
					deflectedVelocity += deflectedVertical * Positions[columnIndex];
				}


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
	[ReadOnly] public NativeArray<float> PlanktonMass;
	[ReadOnly] public NativeArray<float> PlanktonGlucose;
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

		float waterMass = Mass[i];

		float newMass = waterMass;
		float newSaltMass = 0;
		float newPlankton = 0;
		float newGlucose = 0;
		float newCarbon = 0;
		float newTemperature = 0;
		float3 newVelocity = 0;

		if (waterMass > 0)
		{
			for (int j = 0; j < StaticState.MaxNeighborsVert; j++)
			{
				newMass -= math.max(0, Destination[fullRangeIndex * StaticState.MaxNeighborsVert + j]);
			}
			float percentRemaining = newMass / waterMass;

			newPlankton = PlanktonMass[fullRangeIndex] * percentRemaining;
			newGlucose = PlanktonGlucose[fullRangeIndex] * percentRemaining;
			newSaltMass += Salt[fullRangeIndex] * percentRemaining;
			newCarbon += Carbon[fullRangeIndex] * percentRemaining;
			newTemperature += Temperature[fullRangeIndex] * newMass;
			newVelocity += Velocity[fullRangeIndex] * newMass;

			for (int j = 0; j < StaticState.MaxNeighborsVert; j++)
			{
				int n = NeighborsVert[fullRangeIndex * StaticState.MaxNeighborsVert + j];
				if (n >= 0)
				{
					int nColumnIndex = n % CellsPerLayer;
					float incomingMass = math.max(0, -Destination[fullRangeIndex * StaticState.MaxNeighborsVert + j]);
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
					newPlankton += PlanktonMass[n] * incoming;
					newGlucose += PlanktonGlucose[n] * incoming;
					newTemperature += Temperature[n] * incomingMass;

					// TODO: this is increasing speed, is that right???  Shouldnt it only rotate?
					var deflectedVelocity = Velocity[n] + math.cross(Positions[nColumnIndex], Velocity[n]) * CoriolisMultiplier[nColumnIndex] * CoriolisTerm * SecondsPerTick;

					// TODO: turn velocity along great circle, instead of just erasing the vertical component as we are doing here
					var deflectedVertical = math.dot(deflectedVelocity, Positions[nColumnIndex]);
					deflectedVelocity -= Positions[nColumnIndex] * deflectedVertical;
					deflectedVelocity -= Positions[columnIndex] * math.dot(Positions[columnIndex], deflectedVelocity);
					deflectedVelocity += deflectedVertical * Positions[columnIndex];



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
				newTemperature /= newMass;
				newVelocity /= newMass;
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
			Plankton = newPlankton,
			PlanktonGlucose = newGlucose,
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
		float newMass;
		float newTemperature;
		float newDropletMass;

#if DISABLE_CLOUD_ADVECTION
		newMass = Mass[i];
		newTemperature = Temperature[i];
		newDropletMass = DropletMass[i];

#else
		newMass = 0;
		newTemperature = 0;
		newDropletMass = 0;
		float totalMass = 0;

		if (Destination[i].indexA == i)
		{
			float v = Destination[i].valueA;
			float m = Mass[i] * v;
			totalMass += m;
			newMass += m;
			newTemperature += Temperature[i] * v * m;
			newDropletMass += DropletMass[i] * v;
		}

		for (int j = 0; j < 6; j++)
		{
			int n = Neighbors[i * 6 + j];
			if (n >= 0)
			{
				float incoming = 0;
				if (Destination[n].indexA == i)
				{
					incoming = Destination[n].valueA;
				}
				else if (Destination[n].indexB == i)
				{
					incoming = Destination[n].valueB;
				}
				else if (Destination[n].indexC == i)
				{
					incoming = Destination[n].valueC;
				}
				float m = Mass[n] * incoming;
				totalMass += m;
				newMass += m;
				newTemperature += Temperature[n] * incoming * m;
				newDropletMass += DropletMass[n] * incoming;
			}
		}

		if (totalMass > 0)
		{
			newTemperature /= totalMass;
		}

#endif
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
	public NativeSlice<float3> AirVelocity;
	[ReadOnly] public NativeSlice<DiffusionAir> Advection;
	public void Execute(int i)
	{
		Temperature[i] = Advection[i].Temperature;
		Vapor[i] = Advection[i].WaterVapor;
		Dust[i] = Advection[i].Dust;
		CarbonDioxide[i] = Advection[i].CarbonDioxide;
		AirVelocity[i] = Advection[i].Velocity;
	}
}


[BurstCompile]
public struct ApplyAdvectionWaterJob : IJobParallelFor {
	public NativeSlice<float> WaterMass;
	public NativeSlice<float> SaltMass;
	public NativeSlice<float> CarbonMass;
	public NativeSlice<float> PlanktonMass;
	public NativeSlice<float> PlanktonGlucose;
	public NativeSlice<float> Temperature;
	public NativeSlice<float3> Velocity;
	[ReadOnly] public NativeSlice<DiffusionWater> Advection;
	public void Execute(int i)
	{
		WaterMass[i] = Advection[i].WaterMass;
		SaltMass[i] = Advection[i].SaltMass;
		CarbonMass[i] = Advection[i].CarbonMass;
		PlanktonMass[i] = Advection[i].Plankton;
		PlanktonGlucose[i] = Advection[i].PlanktonGlucose;
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


