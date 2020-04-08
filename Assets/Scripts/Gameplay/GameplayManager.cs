using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;
using System.Globalization;
using UnityEngine.EventSystems;

public class GameplayManager : MonoBehaviour {

	public List<GameTool> Tools;

	[Header("References")]
	public WorldSimComponent Sim;

	[HideInInspector]
	public bool ActiveCellLocked { get; private set; }
	[HideInInspector]
	public int ActiveCellIndex { get; private set; }
	[HideInInspector]
	public GameTool ActiveTool;

	public event Action<int> OnSetActiveCell;


	public void SetActiveTool(GameTool tool)
	{
		if (ActiveTool != tool)
		{
			ActiveTool?.OnDeselected();
			ActiveTool = tool;
			ActiveTool.OnSelected();
		}
	}

	public void Update()
	{
	}

	public void OnWorldClicked(Vector3 worldPos, int cell)
	{
		ActiveTool?.OnClick(worldPos, cell);
	}

	public void OnWorldMouseDown(Vector3 worldPos, int cell)
	{
		ActiveTool?.OnPointerDown(worldPos, cell);
	}

	public void OnWorldMouseUp(Vector3 worldPos, int cell)
	{
		ActiveTool?.OnPointerUp(worldPos, cell);
	}

	public void OnWorldHovered(Vector3 worldPos, int cell)
	{
		ActiveTool?.OnUpdate(worldPos, cell);
	}

	public void OnWorldDrag(Vector3 worldPos, int cell, Vector3 direction)
	{
		ActiveTool?.OnDragMove(worldPos, cell, direction);

	}

	public void SetActiveCell(int index, bool locked)
	{
		ActiveCellIndex = index;
		ActiveCellLocked = locked;
		OnSetActiveCell?.Invoke(index);
	}



}
