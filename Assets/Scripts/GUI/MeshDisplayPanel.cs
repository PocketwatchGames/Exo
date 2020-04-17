using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class MeshDisplayPanel : MonoBehaviour {
	public WorldView View;
	public Dropdown OverlayDropdown;
	public Dropdown WindDropdown;
	public GameObject LayerSliderAir;
	public GameObject LayerSliderWater;

	// Start is called before the first frame update
	void Start()
	{
		View.AirLayerChangedEvent += OnAirLayerChanged;
		View.WaterLayerChangedEvent += OnWaterLayerChanged;
		View.MeshOverlayChangedEvent+= OnMeshOverlayChanged;
		View.WindOverlayChangedEvent += OnWindOverlayChanged;

		LayerSliderAir.GetComponentInChildren<Slider>().minValue = 1;
		LayerSliderAir.GetComponentInChildren<Slider>().maxValue = View.Sim.WorldData.AirLayers - 2;
		LayerSliderWater.GetComponentInChildren<Slider>().minValue = 1;
		LayerSliderWater.GetComponentInChildren<Slider>().maxValue = View.Sim.WorldData.WaterLayers - 2;
	}

	// Update is called once per frame
	void Update()
	{

	}

	private enum OverlayType {
		None,
		Air,
		Water
	}

	private Dictionary<WorldView.MeshOverlay, OverlayType> OverlayTypes = new Dictionary<WorldView.MeshOverlay, OverlayType>{
		{ WorldView.MeshOverlay.AbsoluteHumidity, OverlayType.Air },
		{ WorldView.MeshOverlay.CarbonDioxide, OverlayType.Air },
		{ WorldView.MeshOverlay.PotentialTemperature, OverlayType.Air },
		{ WorldView.MeshOverlay.Pressure, OverlayType.Air },
		{ WorldView.MeshOverlay.RelativeHumidity, OverlayType.Air },
		{ WorldView.MeshOverlay.VerticalWind, OverlayType.Air },
		{ WorldView.MeshOverlay.Salinity, OverlayType.Water },
		{ WorldView.MeshOverlay.WaterCarbonDioxide, OverlayType.Water },
		{ WorldView.MeshOverlay.WaterTemperature, OverlayType.Water },
	};
	public void OnOverlayChanged()
	{
		WorldView.MeshOverlay overlay = (WorldView.MeshOverlay)OverlayDropdown.value;
		View.SetActiveMeshOverlay(overlay);
	}
	public void OnWindChanged()
	{
		WorldView.WindOverlay overlay = (WorldView.WindOverlay)OverlayDropdown.value;
		View.SetActiveWindOverlay(overlay);
	}
	public void OnAirSliderChanged()
	{
		View.SetActiveLayerAir((int)LayerSliderAir.GetComponentInChildren<Slider>().value);
	}
	public void OnWaterSliderChanged()
	{
		View.SetActiveLayerWater((int)LayerSliderWater.GetComponentInChildren<Slider>().value);
	}


	private void OnAirLayerChanged()
	{
		LayerSliderAir.GetComponentInChildren<Slider>().SetValueWithoutNotify(View.ActiveMeshLayerAir);
	}
	private void OnWaterLayerChanged()
	{
		LayerSliderWater.GetComponentInChildren<Slider>().SetValueWithoutNotify(View.ActiveMeshLayerWater);
	}
	private void OnMeshOverlayChanged()
	{
		UpdateOverlayType(GetOverlayType(View.ActiveMeshOverlay));
	}
	private void OnWindOverlayChanged()
	{

	}

	private OverlayType GetOverlayType(WorldView.MeshOverlay overlay)
	{
		OverlayType overlayType = OverlayType.None;
		OverlayTypes.TryGetValue(overlay, out overlayType);
		return overlayType;
	}

	private void UpdateOverlayType(OverlayType overlayType)
	{
		LayerSliderAir.SetActive(overlayType == OverlayType.Air);
		LayerSliderWater.SetActive(overlayType == OverlayType.Water);

		switch (overlayType)
		{
			case OverlayType.Air:
				LayerSliderAir.GetComponentInChildren<Slider>().SetValueWithoutNotify(View.ActiveMeshLayerAir);
				break;
			case OverlayType.Water:
				LayerSliderWater.GetComponentInChildren<Slider>().SetValueWithoutNotify(View.ActiveMeshLayerWater);
				break;
			default:
				break;
		}
	}

}
