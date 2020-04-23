using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using UnityEngine.UI;
using UnityEngine.EventSystems;
using System;

public class HUD : MonoBehaviour
{
	public WorldView View;
	public GameplayManager Gameplay;
	public WorldSimComponent Sim;
	public GameToolButton GameToolButtonPrefab;
	public GameObject Toolbar;
	public GameObject OptionsPanel;

	private bool _isDown;
	private Vector3 _dragStart;

	public CellInfo.TemperatureUnits ActiveTemperatureUnits = CellInfo.TemperatureUnits.Celsius;

	private List<List<GameObject>> ModeComponents;

	// Start is called before the first frame update
	void Start()
    {
	}

	public void AddToolbarButton(GameTool tool)
	{
		var b = GameObject.Instantiate(GameToolButtonPrefab, Toolbar.transform);
		b.Tool = tool;
	}

	void Update()
    {
		PointerEventData ped = new PointerEventData(EventSystem.current);
		ped.position = Input.mousePosition;
		List<RaycastResult> results = new List<RaycastResult>();
		GetComponent<GraphicRaycaster>().Raycast(ped, results);
		if (results.Count == 0)
		{
		}
		var p = GetMouseCellIndex();
		if (_isDown)
		{
			Gameplay.OnWorldDrag(p.Item1, p.Item2, _dragStart - Input.mousePosition);
		} else
		{
			Gameplay.OnWorldHovered(p.Item1, p.Item2);
		}
	}

	public void OnOptionsClicked()
	{
		OptionsPanel.SetActive(true);
	}

	public void OnAirLayerChanged(Slider slider)
	{
		View.ActiveMeshLayerAir = (int)slider.value;
	}

	public void OnBackgroundClicked()
	{
		var p = GetMouseCellIndex();
		Gameplay.OnWorldClicked(p.Item1, p.Item2);
	}
	public void OnBackgroundMouseDown()
	{
		_isDown = true;
		_dragStart = Input.mousePosition;
		var p = GetMouseCellIndex();
		Gameplay.OnWorldMouseDown(p.Item1, p.Item2);
	}
	public void OnBackgroundMouseUp()
	{
		_isDown = false;
		var p = GetMouseCellIndex();
		Gameplay.OnWorldMouseUp(p.Item1, p.Item2);
	}

	public string GetCellInfo(CellInfo.CellInfoType cellInfoType)
	{
		switch (cellInfoType)
		{
			case CellInfo.CellInfoType.Global:
				return CellInfo.GetCellInfoGlobal(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.TempState, ref View.DisplayState);
			case CellInfo.CellInfoType.Enthalpy:
				return CellInfo.GetCellInfoEnthalpy(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.TempState, ref View.DisplayState);
			case CellInfo.CellInfoType.Energy:
				return CellInfo.GetCellInfoEnergy(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.TempState, ref View.DisplayState);
			case CellInfo.CellInfoType.Cell:
				return CellInfo.GetCellInfoCell(ActiveTemperatureUnits, Gameplay.ActiveCellIndex, ref Sim.StaticState, ref Sim.LastSimState, ref Sim.TempState, ref View.DisplayState);
			case CellInfo.CellInfoType.Atmosphere:
				return CellInfo.GetCellInfoAtmosphere(ActiveTemperatureUnits, Gameplay.ActiveCellIndex, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.TempState, ref Sim.StaticState, ref View.DisplayState);
			case CellInfo.CellInfoType.Ground:
				return CellInfo.GetCellInfoGround(ActiveTemperatureUnits, Gameplay.ActiveCellIndex, ref Sim.LastSimState, ref Sim.TempState);
			case CellInfo.CellInfoType.Water:
				return CellInfo.GetCellInfoWater(ActiveTemperatureUnits, Gameplay.ActiveCellIndex, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.TempState, ref Sim.StaticState, ref View.DisplayState);
		}
		return "";
	}


	private Tuple<Vector3, int> GetMouseCellIndex()
	{
		var ray = Camera.main.ScreenPointToRay(Input.mousePosition);
		RaycastHit hit;
		int cellIndex = -1;
		Vector3 worldPos = Vector3.zero;
		if (Physics.Raycast(ray, out hit))
		{
			if (hit.collider.gameObject.transform.parent == View.Planet.transform)
			{
				int tIndex;
				if (hit.barycentricCoordinate.x > hit.barycentricCoordinate.y && hit.barycentricCoordinate.x > hit.barycentricCoordinate.z)
				{
					tIndex = 0;
				}
				else if (hit.barycentricCoordinate.y > hit.barycentricCoordinate.z)
				{
					tIndex = 1;
				}
				else
				{
					tIndex = 2;
				}
				cellIndex = View.GetClosestVert(hit.triangleIndex, tIndex);
			}
		}
		return new Tuple<Vector3, int>(worldPos, cellIndex);
	}


}
