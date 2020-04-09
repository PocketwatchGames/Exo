using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;

public class ToolAddWater : GameTool {

	public float MetersPerSecond = 100;

	override public void OnDragMove(Vector3 worldPos, int cellIndex, Vector2 direction) {
		Gameplay.SetActiveCell(cellIndex, false);
		Gameplay.Sim.Edit((ref SimState last, ref SimState next) => { Activate(ref last, ref next, cellIndex, ref Gameplay.Sim.WorldData, ref Gameplay.Sim.DependentState); });
	}
	public override void OnUpdate(Vector3 worldPos, int cellIndex)
	{
		base.OnUpdate(worldPos, cellIndex);
		Gameplay.SetActiveCell(cellIndex, false);
	}

	private void Activate(ref SimState lastState, ref SimState nextState, int cell, ref WorldData worldData, ref DependentState dependent)
	{
		nextState.CopyFrom(ref lastState);

		if (cell >= 0)
		{
			float waterAdded = MetersPerSecond * WorldData.MassWater * Time.deltaTime;
			float waterMass = lastState.WaterMass[worldData.SurfaceWaterLayer][cell];
			float saltMass = lastState.SaltMass[worldData.SurfaceWaterLayer][cell];
			nextState.WaterMass[worldData.SurfaceWaterLayer][cell] = waterMass + waterAdded;
			nextState.WaterTemperature[worldData.SurfaceWaterLayer][cell] =
				(lastState.WaterTemperature[worldData.SurfaceWaterLayer][cell] * (waterMass + saltMass)
				+ dependent.SurfaceAirTemperatureAbsolute[cell] * waterAdded)
				/ (waterMass + saltMass + waterAdded);
		}
	}
}
