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

public class GameplayManager : MonoBehaviour {

	public enum GameTool {
		None,
		MoveEarth,
		MovePlate
	}

	public CellInfo.TemperatureUnits ActiveTemperatureUnits = CellInfo.TemperatureUnits.Celsius;
	public int ActiveCell;
	public bool ActiveCellLocked { get; private set; }
	public int ActiveCellIndex { get; private set; }
	public GameTool ActiveTool;
	private float3 DragStart;

	[Header("References")]
	public WorldSimComponent Sim;

	public event Action<int> OnSetActiveCell;

	public void OnGameToolSelected(GameObject button)
	{
		var tool = button.GetComponent<GameToolButton>()?.Tool;
		Debug.Assert(tool.HasValue);
		SetActiveTool(tool.GetValueOrDefault());
	}

	public void SetActiveTool(GameTool tool)
	{
		ActiveTool = tool;
	}

	public void Update()
	{
		if (ActiveCell != -1)
		{
			SetActiveCell(ActiveCell, true);
			ActiveCell = -1;
		}
	}

	public void OnHUDTemperatureUnitsChanged(UnityEngine.UI.Dropdown dropdown)
	{
		ActiveTemperatureUnits = (CellInfo.TemperatureUnits)dropdown.value;
	}

	public string GetCellInfo(CellInfo.CellInfoType cellInfoType)
	{
		switch (cellInfoType)
		{
			case CellInfo.CellInfoType.Global:
				return CellInfo.GetCellInfoGlobal(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Enthalpy:
				return CellInfo.GetCellInfoEnthalpy(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Energy:
				return CellInfo.GetCellInfoEnergy(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Cell:
				return CellInfo.GetCellInfoCell(ActiveTemperatureUnits, ActiveCellIndex, ref Sim.StaticState, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Atmosphere:
				return CellInfo.GetCellInfoAtmosphere(ActiveTemperatureUnits, ActiveCellIndex, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.StaticState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Ground:
				return CellInfo.GetCellInfoGround(ActiveTemperatureUnits, ActiveCellIndex, ref Sim.ActiveSimState, ref Sim.DependentState);
			case CellInfo.CellInfoType.Water:
				return CellInfo.GetCellInfoWater(ActiveTemperatureUnits, ActiveCellIndex, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.StaticState, ref Sim.DisplayState);
		}
		return "";
	}


	public void SetActiveCell(int index, bool locked)
	{
		ActiveCellIndex = index;
		ActiveCellLocked = locked;
		OnSetActiveCell?.Invoke(index);

	}



}
