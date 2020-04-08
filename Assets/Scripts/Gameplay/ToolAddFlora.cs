using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;

public class ToolAddFlora : GameTool {

	public float Speed = 1000;

	public override void OnSelected()
	{
		base.OnSelected();
		Gameplay.SetActiveCell(-1, false);
	}
	public override void OnDeselected()
	{
		base.OnDeselected();
		Gameplay.SetActiveCell(-1, false);
	}
	override public void OnPointerDown(Vector3 worldPos, int cellIndex) {
		Gameplay.SetActiveCell(cellIndex, false);
	}
	override public void OnPointerUp(Vector3 worldPos, int cellIndex) {
		Gameplay.SetActiveCell(-1, false);
	}
	override public void OnDragMove(Vector3 worldPos, int cellIndex, Vector2 direction) {
		Gameplay.SetActiveCell(cellIndex, false);
		Gameplay.Sim.Edit((ref SimState last, ref SimState next) => { Activate(ref last, ref next, cellIndex, ref Gameplay.Sim.WorldData, ref Gameplay.Sim.DependentState); });
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
