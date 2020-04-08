using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CellInfoButton : MonoBehaviour
{
	public GameObject Panel;
	public void Start()
	{
		Panel.SetActive(GetComponent<ToggleButton>().isOn);
	}
	public void OnCellInfoToggleChanged()
	{
		Panel.SetActive(GetComponent<ToggleButton>().isOn);
	}
}
