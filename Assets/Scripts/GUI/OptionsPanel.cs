using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class OptionsPanel : MonoBehaviour
{
	public HUD HUD;
	public Dropdown TemperatureDropdown;

    // Start is called before the first frame update
    void OnEnable()
    {
		TemperatureDropdown.GetComponent<DropdownOverlay<CellInfo.TemperatureUnits>>().SetValueWithoutNotify(HUD.ActiveTemperatureUnits);
	}

    // Update is called once per frame
    void Update()
    {
        
    }

	public void OnOK()
	{
		gameObject.SetActive(false);
	}

	public void OnHUDTemperatureUnitsChanged()
	{
		HUD.ActiveTemperatureUnits = (CellInfo.TemperatureUnits)TemperatureDropdown.value;
	}

}
