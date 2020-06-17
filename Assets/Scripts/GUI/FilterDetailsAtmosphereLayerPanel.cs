using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class FilterDetailsAtmosphereLayerPanel : FilterDetailsPanel {
	public LegendPanel Legend;

	public override void Show(HUD hud)
	{
		base.Show(hud);
		var o = hud.View.OverlayColors[hud.View.ActiveMeshOverlay];
		Legend.SetValues(o, hud.ActiveTemperatureUnits);
	}
}
