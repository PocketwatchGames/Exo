﻿using System;
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
		Gameplay.Sim.Edit((ref SimState last, ref SimState next) => { Activate(ref last, ref next, cellIndex, ref Gameplay.Sim.WorldData, ref Gameplay.Sim.LastTempState, ref Gameplay.Sim.StaticState); });
	}
	public override void OnUpdate(Vector3 worldPos, int cellIndex)
	{
		base.OnUpdate(worldPos, cellIndex);
		Gameplay.SetActiveCell(cellIndex, false);
	}

	private void Activate(ref SimState lastState, ref SimState nextState, int cell, ref WorldData worldData, ref TempState tempState, ref StaticState staticState)
	{
		nextState.CopyFrom(ref lastState);

		if (cell >= 0)
		{
			int index = staticState.GetWaterIndex(worldData.SurfaceAirLayer, cell);
			float waterAdded = MetersPerSecond * WorldData.MassWater * Time.deltaTime;
			float waterMass = lastState.WaterMass[index];
			float saltMass = lastState.WaterSaltMass[index];
			nextState.WaterMass[index] = waterMass + waterAdded;
			nextState.WaterTemperature[index] =
				(lastState.WaterTemperature[index] * (waterMass + saltMass)
				+ tempState.SurfaceAirTemperatureAbsolute[cell] * waterAdded)
				/ (waterMass + saltMass + waterAdded);
		}
	}
}
