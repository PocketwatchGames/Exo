using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using UnityEngine.UI;
using UnityEngine.EventSystems;
using System;
using System.Linq;

public class HUD : MonoBehaviour
{
	public enum FilterType {
		None,
		CarbonDioxide,
		Nitrogen,		
	}

	[Serializable]
	public struct SimSpeed {
		public float Speed;
		public Image Image;
	}

	public WorldView View;
	public GameplayManager Gameplay;
	public WorldSimComponent Sim;
	public GameToolButton GameToolButtonPrefab;
	public GameObject Toolbar;
	public GameObject OptionsPanel;
	public GameObject FiltersPanel;
	public GameObject FilterDetailsPanels;
	public Image PauseButtonPausedImage;
	public Image PauseButtonPlayImage;
	public List<SimSpeed> SimSpeeds;
	public ToggleGroup ModeToggleGroup;
	public List<GameObject> ModePanels;
	public Text EnergyText;
	public Text EnergyDeltaText;
	public Text BiomassText;
	public Text BiomassDeltaText;
	public Text TemperatureText;
	public Text TemperatureDeltaText;
	public Text TimeDateText;
	public ToggleGroup FilterToggleGroup;

	private bool _isDown;
	private Vector3 _dragStart;
	private int _simSpeedIndex;
	private bool _paused = true;
	private int _mode;
	private float _averageTemperature;

	public CellInfo.TemperatureUnits ActiveTemperatureUnits = CellInfo.TemperatureUnits.Celsius;

	private List<List<GameObject>> ModeComponents;

	// Start is called before the first frame update
	void Awake()
    {
		View.OnMeshOverlayChanged += OnFilterChanged;
		Gameplay.Sim.OnTimeScaleChanged += OnTimeScaleChanged;
		Sim.OnTick += OnTick;
		SetSimSpeedIndex(0);
		SetMode(0);
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

		if (Input.GetKeyUp(KeyCode.Period))
		{
			Sim.StepTime();
		}
	}

	public void OnOptionsClicked()
	{
		OptionsPanel.SetActive(!OptionsPanel.activeSelf);
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
				return CellInfo.GetCellInfoGlobal(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.LastTempState, ref View.DisplayState);
			case CellInfo.CellInfoType.Enthalpy:
				return CellInfo.GetCellInfoEnthalpy(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.LastTempState, ref View.DisplayState);
			case CellInfo.CellInfoType.Energy:
				return CellInfo.GetCellInfoEnergy(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.LastTempState, ref View.DisplayState);
			case CellInfo.CellInfoType.Cell:
				return CellInfo.GetCellInfoCell(ActiveTemperatureUnits, Gameplay.ActiveCellIndex, ref Sim.StaticState, ref Sim.LastSimState, ref Sim.LastTempState, ref View.DisplayState);
			case CellInfo.CellInfoType.Atmosphere:
				return CellInfo.GetCellInfoAtmosphere(ActiveTemperatureUnits, Gameplay.ActiveCellIndex, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.LastTempState, ref Sim.StaticState, ref View.DisplayState);
			case CellInfo.CellInfoType.Ground:
				return CellInfo.GetCellInfoGround(ActiveTemperatureUnits, Gameplay.ActiveCellIndex, ref Sim.LastSimState, ref Sim.LastTempState);
			case CellInfo.CellInfoType.Water:
				return CellInfo.GetCellInfoWater(ActiveTemperatureUnits, Gameplay.ActiveCellIndex, ref Sim.WorldData, ref Sim.LastSimState, ref Sim.LastTempState, ref Sim.StaticState, ref View.DisplayState);
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

	public void OnSimSpeedClicked()
	{
		SetSimSpeedIndex((_simSpeedIndex + 1) % SimSpeeds.Count);
		Gameplay.Sim.TimeScale = _paused ? 0 : SimSpeeds[_simSpeedIndex].Speed;
	}

	public void OnPauseClicked()
	{
		_paused = !_paused;
		Gameplay.Sim.TimeScale = _paused ? 0 : SimSpeeds[_simSpeedIndex].Speed;
	}

	private void SetSimSpeedIndex(int index)
	{
		_simSpeedIndex = index;
		for (int i=0;i<SimSpeeds.Count;i++)
		{
			SimSpeeds[i].Image.gameObject.SetActive(_simSpeedIndex == i);
		}
	}

	private void OnTimeScaleChanged(float timeScale)
	{
		PauseButtonPausedImage.gameObject.SetActive(timeScale==0);
		PauseButtonPlayImage.gameObject.SetActive(timeScale > 0);
		if (timeScale > 0)
		{
			for (int i = SimSpeeds.Count - 1; i >= 0; i--)
			{
				if (timeScale >= SimSpeeds[i].Speed)
				{
					SetSimSpeedIndex(i);
					break;
				}
				else
				{
					SimSpeeds[i].Image.gameObject.SetActive(false);
				}
			}
		}
	}

	public void OnFiltersClicked()
	{
		if (FiltersPanel.activeSelf)
		{
			FilterToggleGroup.SetAllTogglesOff();
		}
		FiltersPanel.SetActive(!FiltersPanel.activeSelf);
		if (!FiltersPanel.activeSelf)
		{
			View.SetActiveMeshOverlay(WorldView.MeshOverlay.None);
		}
	}

	public void OnFilterChanged(WorldView.MeshOverlay overlay)
	{
	}

	public void OnFilterClicked(Toggle t)
	{
		for (int i = 0; i < FilterDetailsPanels.transform.childCount; i++)
		{
			FilterDetailsPanels.transform.GetChild(i).gameObject.SetActive(false);
		}
		if (t.isOn)
		{
			View.SetActiveMeshOverlay(t.GetComponent<ToggleOverlayMesh>().Overlay);
			t.GetComponent<ToggleOverlayMesh>().DetailsPanel.SetActive(true);
			t.GetComponent<ToggleOverlayMesh>().DetailsPanel.GetComponent<FilterDetailsPanel>().Show(this);
		} else
		{
			View.SetActiveMeshOverlay(WorldView.MeshOverlay.None);
		}
	}


	public void OnModeClicked(int mode)
	{
		SetMode(mode);
	}

	private void SetMode(int mode)
	{
		_mode = mode;
		for (int i=0;i<ModePanels.Count;i++)
		{
			ModePanels[i].SetActive(mode == i);
		}

	}

	private void OnTick()
	{
		var time = WorldTime.GetTime(Sim.LastSimState.PlanetState.Ticks, Sim.LastSimState.PlanetState.SpinSpeed) / Sim.LastSimState.PlanetState.SpinSpeed;
		float days = WorldTime.GetDays(Sim.LastSimState.PlanetState.Ticks, Sim.LastSimState.PlanetState.SpinSpeed * Sim.WorldData.SecondsPerTick);
		TimeDateText.text = ((int)days).ToString() + "\\" + ((int)(Sim.LastSimState.PlanetState.Ticks * Sim.LastSimState.PlanetState.OrbitSpeed)).ToString("X4");
		//		TimeDateText.text = ;

		int energyDelta = 0;
		EnergyText.text = "0";
		Color deltaColor = Color.gray;
		string deltaString = "+" + energyDelta;
		EnergyDeltaText.text = energyDelta.ToString("0");
		EnergyDeltaText.color = deltaColor;

		float temperatureDelta = ((float)View.DisplayState.GlobalAirTemperaturePotential - _averageTemperature);
		TemperatureText.text = CellInfo.GetTemperatureString((float)View.DisplayState.GlobalSurfaceTemperature, ActiveTemperatureUnits, 1);
		float deltaSign = Math.Sign(temperatureDelta);
		if (deltaSign > 0)
		{
			deltaString = "+" + temperatureDelta;
			deltaColor = new Color(255, 100, 100, 255);
		}
		else if (deltaSign < 0)
		{
			deltaString = "-" + temperatureDelta;
			deltaColor = new Color(100, 100, 255, 255);
		}
		else
		{
			deltaString = "+" + temperatureDelta;
			deltaColor = Color.gray;
		}
		TemperatureDeltaText.text = temperatureDelta.ToString("0.0");
		TemperatureDeltaText.color = deltaColor;

		float biomassDelta = 0;
		deltaSign = Math.Sign(biomassDelta);
		if (deltaSign > 0)
		{
			deltaString = "+" + biomassDelta;
			deltaColor = new Color(0, 255, 0, 255);
		}
		else if (deltaSign < 0)
		{
			deltaString = "-" + biomassDelta;
			deltaColor = new Color(255, 0, 0, 255);
		}
		else
		{
			deltaString = "+" + biomassDelta;
			deltaColor = Color.gray;
		}
		BiomassText.text = "0";
		BiomassDeltaText.text = deltaString;
		BiomassDeltaText.color = deltaColor;
	}


}
