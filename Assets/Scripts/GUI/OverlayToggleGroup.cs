using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Linq;

public class OverlayToggleGroup : MonoBehaviour
{
	public WorldView View;
	public Slider LayerSliderAir;
	public Slider LayerSliderWater;

	// Start is called before the first frame update
	void Start()
    {
		LayerSliderAir.GetComponentInChildren<Slider>().minValue = 1;
		LayerSliderAir.GetComponentInChildren<Slider>().maxValue = View.Sim.WorldData.AirLayers - 2;
		//LayerSliderWater.GetComponentInChildren<Slider>().minValue = 1;
		//LayerSliderWater.GetComponentInChildren<Slider>().maxValue = View.Sim.WorldData.WaterLayers - 2;

	}

	// Update is called once per frame
	void Update()
    {
        
    }

	public void OnMeshOverlayChanged()
	{
		var group = GetComponent<ToggleGroup>();
		if (!group.AnyTogglesOn())
		{
			View.SetActiveMeshOverlay(WorldView.MeshOverlay.None);
		}
		else
		{
			View.SetActiveMeshOverlay(group.ActiveToggles().FirstOrDefault().GetComponent<ToggleOverlayMesh>().Overlay);
		}
	}

	public void OnAirSliderChanged()
	{
		View.SetActiveLayerAir((int)LayerSliderAir.GetComponentInChildren<Slider>().value);
	}
	public void OnWaterSliderChanged()
	{
		View.SetActiveLayerWater((int)LayerSliderWater.GetComponentInChildren<Slider>().value);
	}

}
