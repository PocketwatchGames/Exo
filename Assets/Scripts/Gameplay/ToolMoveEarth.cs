using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

public class ToolMoveEarth : GameTool {

	public int ActiveCell;
	public bool IsActive { get { return ActiveCell >= 0; } }

	public override void OnSelected()
	{
		base.OnSelected();
		ActiveCell = -1;
		Gameplay.SetActiveCell(-1, false);
	}
	public override void OnDeselected()
	{
		base.OnDeselected();
		ActiveCell = -1;
		Gameplay.SetActiveCell(-1, false);
	}
	override public void OnPointerDown(Vector3 worldPos, int cellIndex) {
		ActiveCell = cellIndex;
		Gameplay.SetActiveCell(cellIndex, false);
	}
	override public void OnPointerUp(Vector3 worldPos, int cellIndex) {
		ActiveCell = -1;
		Gameplay.SetActiveCell(-1, false);
	}
	override public void OnDragMove(Vector3 worldPos, int cellIndex, Vector2 direction) {
	}
}
