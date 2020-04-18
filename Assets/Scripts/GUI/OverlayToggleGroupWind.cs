using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Linq;

public class OverlayToggleGroupWind : MonoBehaviour
{
	public WorldView View;

	// Start is called before the first frame update
	void Start()
    {
	}

	// Update is called once per frame
	void Update()
    {
        
    }

	public void OnWindOverlayChanged()
	{
		var group = GetComponent<ToggleGroup>();
		if (!group.AnyTogglesOn())
		{
			View.SetActiveWindOverlay(WorldView.WindOverlay.None);
		}
		else
		{
			View.SetActiveWindOverlay(group.ActiveToggles().FirstOrDefault().GetComponent<ToggleOverlayWind>().Overlay);
		}
	}

}
