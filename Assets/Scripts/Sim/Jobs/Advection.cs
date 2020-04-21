//#define DISABLE_VERTICAL_AIR_MOVEMENT
//#define DISABLE_AIR_ADVECTION
//#define DISABLE_WATER_ADVECTION
//#define DISABLE_CLOUD_ADVECTION

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
	public NativeArray<DiffusionAir> Delta;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> TemperatureAbove;
	[ReadOnly] public NativeArray<float> TemperatureBelow;
	[ReadOnly] public NativeArray<float> AirMass;
	[ReadOnly] public NativeArray<float> AirMassAbove;
	[ReadOnly] public NativeArray<float> AirMassBelow;
	[ReadOnly] public NativeArray<float> Vapor;
	[ReadOnly] public NativeArray<float> VaporAbove;
	[ReadOnly] public NativeArray<float> VaporBelow;
	[ReadOnly] public NativeArray<float> CarbonDioxide;
	[ReadOnly] public NativeArray<float> CarbonDioxideAbove;
	[ReadOnly] public NativeArray<float> CarbonDioxideBelow;
	[ReadOnly] public NativeArray<float> Dust;
	[ReadOnly] public NativeArray<float> DustAbove;
	[ReadOnly] public NativeArray<float> DustBelow;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> VelocityAbove;
	[ReadOnly] public NativeArray<float3> VelocityBelow;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<BarycentricValueVertical> Destination;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationAbove;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationBelow;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public NativeArray<float> NeighborDistInverse;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<float> LayerMiddle;
	[ReadOnly] public NativeArray<float> LayerMiddleAbove;
	[ReadOnly] public NativeArray<float> LayerMiddleBelow;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float SecondsPerTick;
	[ReadOnly] public float TicksPerSecond;
	public void Execute(int i)
	{
		float newTemperature;
		float newWaterVapor;
		float newCO2;
		float newDust;
		float3 newVelocity;


#if DISABLE_AIR_ADVECTION
		newTemperature = Temperature[i];
		newWaterVapor = Vapor[i];
		newVelocity = Velocity[i];
		newDust = Dust[i];
		newCO2 = CarbonDioxide[i];
#else

		newTemperature = 0;
		newWaterVapor = 0;
		newCO2 = 0;
		newVelocity = 0;
		newDust = 0;

		float layerMiddle = LayerMiddle[i];

		// TODO: remove this when we have incompressibility
		float totalMass = 0;

		if (Destination[i].indexA == i)
		{
			float v = Destination[i].valueA;
			float massRemaining = v * AirMass[i]; 
			totalMass += massRemaining;
			newTemperature += Temperature[i] * massRemaining;
			newVelocity += Velocity[i] * massRemaining;
			newWaterVapor += Vapor[i] * v;
			newCO2 += CarbonDioxide[i] * v;
			newDust += Dust[i] * v;
		}

		// TODO: subtract vertical motion first before applying horizontal motion (right now it adds up to more than 1)
		// Wait, does it? is this comment outdated?
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
				float incomingMass = incoming * AirMass[n];
				totalMass += incomingMass;
				newTemperature += Temperature[n] * incomingMass;
				newWaterVapor += Vapor[n] * incoming;
				newCO2 += CarbonDioxide[n] * incoming;
				newDust += Dust[n] * incoming;

				// TODO: this is increasing speed, is that right???  Shouldnt it only rotate?
				var deflectedVelocity = Velocity[n] + math.cross(Positions[n], Velocity[n]) * CoriolisMultiplier[n] * CoriolisTerm * SecondsPerTick;

				// TODO: turn velocity along great circle, instead of just erasing the vertical component as we are doing here
				var deflectedVertical = Utils.GetVerticalComponent(deflectedVelocity, Positions[n]);
				deflectedVelocity += deflectedVertical * (Positions[i] - Positions[n]);



				// TODO: is this really appropriate? do we have vertical motion due to the atmospheric layer changing height, or is this accounted for in the pressure gradient?
				//deflectedVelocity += (LayerMiddle[n] - layerMiddle) * NeighborDistInverse[i * 6 + j] * TicksPerSecond;

				// TODO: this is temp
				// need to deal with centrifugal force/gravity so that as air moves horizontally, it can fly into the air or get pulled to earth
				//newVelocity += incomingMass * (deflectedVelocity - Positions[i] * math.dot(Positions[i], deflectedVelocity));

				newVelocity += incomingMass * deflectedVelocity;
			}
		}



#if !DISABLE_VERTICAL_AIR_MOVEMENT

		// from top coming down
		{
			var vertMove = DestinationAbove[i];
			if (vertMove.moveVertical < 0)
			{
				float incomingMass = math.max(0, -AirMassAbove[i] * vertMove.moveVertical);
				totalMass += incomingMass;
				newWaterVapor += VaporAbove[i] * -vertMove.moveVertical;
				newCO2 += CarbonDioxideAbove[i] * -vertMove.moveVertical;
				newDust += DustAbove[i] * -vertMove.moveVertical;
				newTemperature += TemperatureAbove[i] * incomingMass;
				newVelocity += VelocityAbove[i] * incomingMass;
			}
		}

		// from bottom going up
		{
			var vertMove = DestinationBelow[i];
			if (vertMove.moveVertical > 0)
			{
				float incomingMass = math.max(0, AirMassBelow[i] * vertMove.moveVertical);
				totalMass += incomingMass;
				newWaterVapor += VaporBelow[i] * vertMove.moveVertical;
				newCO2 += CarbonDioxideBelow[i] * vertMove.moveVertical;
				newDust += DustBelow[i] * vertMove.moveVertical;
				newTemperature += TemperatureBelow[i] * incomingMass;
				newVelocity += VelocityBelow[i] * incomingMass;
			}
		}

#endif
		if (totalMass > 0)
		{
			newTemperature /= totalMass;
			newVelocity /= totalMass;
		}
		else
		{
			// TODO: remove once we have incompressibility
			newTemperature = Temperature[i];
			newVelocity = Velocity[i];
		}

#endif

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
public struct AdvectionWaterJob : IJobParallelFor {
	public NativeArray<DiffusionWater> Delta;
	[ReadOnly] public NativeArray<BarycentricValueVertical> Destination;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationAbove;
	[ReadOnly] public NativeArray<BarycentricValueVertical> DestinationBelow;
	[ReadOnly] public NativeArray<float3> Positions;
	[ReadOnly] public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<float3> VelocityAbove;
	[ReadOnly] public NativeArray<float3> VelocityBelow;
	[ReadOnly] public NativeArray<float> Mass;
	[ReadOnly] public NativeArray<float> MassAbove;
	[ReadOnly] public NativeArray<float> MassBelow;
	[ReadOnly] public NativeArray<float> Temperature;
	[ReadOnly] public NativeArray<float> TemperatureAbove;
	[ReadOnly] public NativeArray<float> TemperatureBelow;
	[ReadOnly] public NativeArray<float> Salt;
	[ReadOnly] public NativeArray<float> SaltAbove;
	[ReadOnly] public NativeArray<float> SaltBelow;
	[ReadOnly] public NativeArray<float> Carbon;
	[ReadOnly] public NativeArray<float> CarbonAbove;
	[ReadOnly] public NativeArray<float> CarbonBelow;
	[ReadOnly] public NativeArray<float> PlanktonMass;
	[ReadOnly] public NativeArray<float> PlanktonGlucose;
	[ReadOnly] public NativeArray<float> CoriolisMultiplier;
	[ReadOnly] public NativeArray<int> Neighbors;
	[ReadOnly] public float CoriolisTerm;
	[ReadOnly] public float SecondsPerTick;
	public void Execute(int i)
	{
		float waterMass = Mass[i];

		float newMass;
		float newSaltMass;
		float newPlankton;
		float newGlucose;
		float newCarbon;
		float newTemperature;
		float3 newVelocity;

#if DISABLE_WATER_ADVECTION
		newMass = waterMass;
		newSaltMass = Salt[i];
		newPlankton = PlanktonMass[i];
		newGlucose = PlanktonGlucose[i];
		newCarbon = Carbon[i];
		newTemperature = Temperature[i];
		newVelocity = Velocity[i];
#else
		newMass = 0;
		newSaltMass = 0;
		newCarbon = 0;
		newPlankton = 0;
		newGlucose = 0;
		newTemperature = 0;
		newVelocity = 0;

		if (waterMass > 0) {

			float valueRemaining = 0;
			float valueRemainingHorizontal = 0;
			int destIndexA = Destination[i].indexA;
			if (destIndexA == i || destIndexA < 0 || Mass[destIndexA] == 0)
			{
				valueRemaining += Destination[i].valueA;
			}
			int destIndexB = Destination[i].indexB;
			if (destIndexB == i || destIndexB < 0 || Mass[destIndexB] == 0)
			{
				valueRemaining += Destination[i].valueB;
			}
			int destIndexC = Destination[i].indexC;
			if (destIndexC == i || destIndexC < 0 || Mass[destIndexC] == 0)
			{
				valueRemaining += Destination[i].valueC;
			}
			valueRemainingHorizontal = valueRemaining + math.abs(Destination[i].moveVertical);
			if (Destination[i].moveVertical > 0)
			{
				if (MassAbove[i] == 0)
				{
					valueRemaining += Destination[i].moveVertical;
				}
			} else
			{
				if (MassBelow[i] == 0)
				{
					valueRemaining -= Destination[i].moveVertical;
				}
			}
			newPlankton = PlanktonMass[i] * valueRemainingHorizontal;
			newGlucose = PlanktonGlucose[i] * valueRemainingHorizontal;
			newMass = Mass[i] * valueRemaining;
			newSaltMass += Salt[i] * valueRemaining;
			newCarbon += Carbon[i] * valueRemaining;
			newTemperature += Temperature[i] * newMass;
			newVelocity += Velocity[i] * newMass;

			for (int j = 0; j < 6; j++)
			{
				int n = Neighbors[i * 6 + j];
				if (n >= 0)
				{
					float nMass = Mass[n];
					if (nMass > 0)
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

						float massIncoming = nMass * incoming;
						newMass += massIncoming;
						newSaltMass += Salt[n] * incoming;
						newCarbon += Carbon[n] * incoming;
						newPlankton += PlanktonMass[n] * incoming;
						newGlucose += PlanktonGlucose[n] * incoming;
						newTemperature += Temperature[n] * massIncoming;

						// TODO: this is increasing speed, is that right???  Shouldnt it only rotate?
						var deflectedVelocity = Velocity[n] + math.cross(Positions[n], Velocity[n]) * CoriolisMultiplier[n] * CoriolisTerm * SecondsPerTick;

						// TODO: turn velocity along great circle, instead of just erasing the vertical component as we are doing here
						var deflectedVertical = Utils.GetVerticalComponent(deflectedVelocity, Positions[n]);
						deflectedVelocity += deflectedVertical * (Positions[i] - Positions[n]);



						// TODO: is this really appropriate? do we have vertical motion due to the atmospheric layer changing height, or is this accounted for in the pressure gradient?
						//deflectedVelocity += (LayerMiddle[n] - layerMiddle) * NeighborDistInverse[i * 6 + j] * TicksPerSecond;

						// TODO: this is temp
						// need to deal with centrifugal force/gravity so that as air moves horizontally, it can fly into the air or get pulled to earth
						//newVelocity += incomingMass * (deflectedVelocity - Positions[i] * math.dot(Positions[i], deflectedVelocity));

						newVelocity += massIncoming * deflectedVelocity;



						//					newVelocity += VelocityDeflected[n] * massIncoming;
					}
				}
			}


			// from top coming down
			{
				var vertMove = DestinationAbove[i].moveVertical;
				float massIncoming = -(MassAbove[i] * vertMove);
				if (massIncoming < 0)
				{
					newMass += massIncoming;
					newTemperature += TemperatureAbove[i] * massIncoming;
					newVelocity += VelocityAbove[i] * massIncoming;
					newSaltMass += SaltAbove[i] * vertMove;
					newCarbon += CarbonAbove[i] * vertMove;
				}
			}

			// from bottom going up
			{
				var vertMove = DestinationBelow[i].moveVertical;
				float massIncoming = MassBelow[i] * vertMove;
				if (massIncoming > 0)
				{
					newMass += massIncoming;
					newTemperature += TemperatureBelow[i] * massIncoming;
					newVelocity += VelocityBelow[i] * massIncoming;
					newSaltMass += SaltBelow[i] * vertMove;
					newCarbon += CarbonBelow[i] * vertMove;
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
				newTemperature = Temperature[i];
				newVelocity = Velocity[i];
			}
		}
#endif

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
public struct ApplyAdvectionAirJob : IJobParallelFor {
	public NativeArray<float> Temperature;
	public NativeArray<float> Vapor;
	public NativeArray<float> Dust;
	public NativeArray<float> CarbonDioxide;
	public NativeArray<float3> AirVelocity;
	[ReadOnly] public NativeArray<DiffusionAir> Advection;
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
	public NativeArray<float> WaterMass;
	public NativeArray<float> SaltMass;
	public NativeArray<float> CarbonMass;
	public NativeArray<float> PlanktonMass;
	public NativeArray<float> PlanktonGlucose;
	public NativeArray<float> Temperature;
	public NativeArray<float3> Velocity;
	[ReadOnly] public NativeArray<DiffusionWater> Advection;
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


