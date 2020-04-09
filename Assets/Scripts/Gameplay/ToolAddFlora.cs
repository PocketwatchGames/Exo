using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;

public class ToolAddFlora : GameTool {

	public float Speed = 1000;

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
			float floraAdded = Speed * Time.deltaTime;
			nextState.FloraMass[cell] += floraAdded;
			nextState.FloraGlucose[cell] += floraAdded;
			nextState.FloraWater[cell] += floraAdded;
		}
	}
}
