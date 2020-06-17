using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class FilterDetailsPanel : MonoBehaviour {
	public Text Title;

	virtual public void Show(HUD hud)
	{
		var o = hud.View.OverlayColors[hud.View.ActiveMeshOverlay];
		Title.text = o.Title;
	}

	virtual public void OnUpdate()
	{
	}
}
