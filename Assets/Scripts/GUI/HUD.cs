using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using UnityEngine.UI;
using UnityEngine.EventSystems;

public class HUD : MonoBehaviour
{
	public WorldView View;
	public GameplayManager Gameplay;

	public enum HUDMode {
		Info,
		Geology,
		Plants,
		Animals,
		Count
	}

	public HUDMode ActiveMode;

	public List<GameObject> ModeComponentsInfo;
	public List<GameObject> ModeComponentsGeo;
	public List<GameObject> ModeComponentsPlants;
	public List<GameObject> ModeComponentsAnimals;

	private List<List<GameObject>> ModeComponents;

	// Start is called before the first frame update
	void Start()
    {
		ModeComponents = new List<List<GameObject>> { ModeComponentsInfo, ModeComponentsGeo, ModeComponentsPlants, ModeComponentsAnimals };
		SetMode(HUDMode.Info);
	}

// Update is called once per frame
	void Update()
    {
		PointerEventData ped = new PointerEventData(EventSystem.current);
		ped.position = Input.mousePosition;
		List<RaycastResult> results = new List<RaycastResult>();
		GetComponent<GraphicRaycaster>().Raycast(ped, results);
		if (results.Count == 0)
		{
			if (!Gameplay.ActiveCellLocked)
			{
				Gameplay.SetActiveCell(GetMouseCellIndex(), false);
			}
		}
	}

	public void OnBackgroundClicked()
	{
		int cellHovered = GetMouseCellIndex();
		Gameplay.SetActiveCell(cellHovered, cellHovered >= 0);
	}

	private int GetMouseCellIndex()
	{
		var ray = Camera.main.ScreenPointToRay(Input.mousePosition);
		RaycastHit hit;
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
				return View.GetClosestVert(hit.triangleIndex, tIndex);
			}
		}
		return -1;
	}

	public void SetMode(HUDMode mode)
	{
		for (int i = 0; i < (int)HUDMode.Count; i++)
		{
			if (i != (int)mode)
			{
				foreach (var c in ModeComponents[i])
				{
					c.SetActive(false);
				}
			}
		}
		foreach (var c in ModeComponents[(int)mode])
		{
			c.SetActive(true);
		}
		ActiveMode = mode;
	}

	public void OnHUDModeClicked(UnityEngine.UI.Button button)
	{
		SetMode(button.GetComponent<HUDModeButton>().Mode);
	}
}
