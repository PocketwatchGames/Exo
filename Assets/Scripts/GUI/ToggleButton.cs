using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.UI;

public class ToggleButton : Toggle {
	private Text _buttonText;
	private Image _background;

	// Start is called before the first frame update
	public void Awake()
	{
		base.Start();
		_background = GetComponentInChildren<Image>();
		_buttonText = GetComponentInChildren<Text>();
		_buttonText.fontStyle = isOn ? FontStyle.Bold : FontStyle.Normal;
		_background.color = isOn ? colors.selectedColor : colors.normalColor;
	}

	public void OnValueChanged()
	{
		_buttonText.fontStyle = isOn ? FontStyle.Bold : FontStyle.Normal;
		_background.color = isOn ? colors.selectedColor : colors.normalColor;
	}
}
