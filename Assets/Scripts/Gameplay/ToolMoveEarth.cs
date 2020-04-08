using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

public class ToolMoveEarth : GameTool {

	public float DigSpeed = 100;

	public override void OnSelected()
	{
		base.OnSelected();
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
		Gameplay.Sim.Edit((ref SimState last, ref SimState next) => { Activate(ref last, ref next, cellIndex, cellIndex); });
	}

	private void Activate(ref SimState lastState, ref SimState nextState, int from, int to)
	{
		nextState.CopyFrom(ref lastState);

		if (to >= 0 && from >= 0)
		{
			float move = DigSpeed * Time.deltaTime;
			nextState.Elevation[from] -= move;
			nextState.Elevation[to] += move;
		}
	}
}
