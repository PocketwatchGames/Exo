using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

public class ToolSelect : GameTool {

	public int ActiveCell;
	public bool ActiveCellLocked { get; private set; }
	public int ActiveCellIndex { get; private set; }

	public override void OnDeselected()
	{
		base.OnDeselected();
		Gameplay.Sim.SimSettings.CollectGlobalsDebug = false;
	}
	public override void OnSelected()
	{
		base.OnSelected();
		Gameplay.Sim.SimSettings.CollectGlobalsDebug = true;
	}
	override public void OnClick(Vector3 worldPos, int cellIndex)
	{
		Gameplay.SetActiveCell(cellIndex, cellIndex >= 0);
	}
	override public void OnUpdate(Vector3 worldPos, int cellIndex) {
		if (ActiveCell != -1)
		{
			Gameplay.SetActiveCell(ActiveCell, true);
			ActiveCell = -1;
		}
	}
}
