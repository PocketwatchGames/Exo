using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

public class ToolMoveEarth : GameTool {

	public int ActiveCell;
	public bool ActiveCellLocked { get; private set; }
	public int ActiveCellIndex { get; private set; }

	override public void OnClick(Vector3 worldPos, int cellIndex)
	{
	}
	override public void OnUpdate(Vector3 worldPos, int cellIndex) {
	}
}
