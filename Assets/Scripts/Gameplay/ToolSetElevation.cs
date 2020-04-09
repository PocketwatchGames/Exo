using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;

public class ToolSetElevation : GameTool {

	public float MoveSpeed = 100;
	public float Direction = 1;

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
			float move = Direction * MoveSpeed * Time.deltaTime;

			nextState.Elevation[cell] += move;

			int neighborCount = staticState.GetMaxNeighbors(cell);
			for (int i = 0; i < neighborCount; i++)
			{
				nextState.Elevation[staticState.Neighbors[cell * StaticState.MaxNeighbors + i]] -= move / neighborCount;
			}
		}
	}
}
