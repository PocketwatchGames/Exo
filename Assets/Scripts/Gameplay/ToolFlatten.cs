using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;

public class ToolFlatten : GameTool {

	public float MoveSpeed = 100;

	override public void OnDragMove(Vector3 worldPos, int cellIndex, Vector2 direction) {
		Gameplay.SetActiveCell(cellIndex, false);
		Gameplay.Sim.Edit((ref SimState last, ref SimState next) => { Activate(ref last, ref next, ref Gameplay.Sim.StaticState, cellIndex); });
	}
	public override void OnUpdate(Vector3 worldPos, int cellIndex)
	{
		base.OnUpdate(worldPos, cellIndex);
		Gameplay.SetActiveCell(cellIndex, false);
	}

	private void Activate(ref SimState lastState, ref SimState nextState, ref StaticState staticState, int cell)
	{
		nextState.CopyFrom(ref lastState);

		if (cell >= 0)
		{
			float elevation = lastState.Elevation[cell];
			float avgElevation = lastState.Elevation[cell];
			float neighborElevation = 0;
			int neighborCount = StaticState.GetMaxNeighbors(cell, staticState.Neighbors);
			for (int i=0;i< neighborCount; i++)
			{
				neighborElevation += lastState.Elevation[staticState.Neighbors[cell * StaticState.MaxNeighbors + i]];
			}
			avgElevation = (neighborElevation + elevation) / (neighborCount + 1);

			float moveDir = math.sign(avgElevation - elevation);
			float move = moveDir * math.min(MoveSpeed * Time.deltaTime, math.abs(avgElevation - elevation));

			nextState.Elevation[cell] += move;

			for (int i = 0; i < neighborCount; i++)
			{
				nextState.Elevation[staticState.Neighbors[cell * StaticState.MaxNeighbors + i]] -= move / neighborCount;
			}
		}
	}
}
