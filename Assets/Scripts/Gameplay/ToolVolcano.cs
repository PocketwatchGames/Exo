using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;

public class ToolVolcano : GameTool {

	public float AddMagmaDepthSpeed = 100;

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
	}
	override public void OnPointerUp(Vector3 worldPos, int cellIndex) {
		Gameplay.SetActiveCell(-1, false);
	}
	override public void OnDragMove(Vector3 worldPos, int cellIndex, Vector2 direction) {
		Gameplay.SetActiveCell(cellIndex, false);
		Gameplay.Sim.Edit((ref SimState last, ref SimState next) => { Activate(ref last, ref next, cellIndex); });
	}

	private void Activate(ref SimState lastState, ref SimState nextState, int cell)
	{
		nextState.CopyFrom(ref lastState);

		if (cell >= 0)
		{
			float magmaVolume = AddMagmaDepthSpeed * Time.deltaTime;
			nextState.MagmaMass[cell] = lastState.MagmaMass[cell] + magmaVolume * WorldData.MassLava;

			float crustDelta = math.min(lastState.CrustDepth[cell], magmaVolume);
			nextState.CrustDepth[cell] = lastState.CrustDepth[cell] - crustDelta;
			nextState.Elevation[cell] = lastState.Elevation[cell] - crustDelta;
		}
	}
}
