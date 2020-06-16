using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using Unity.Collections;
using System;

public class LegendPanel : MonoBehaviour
{
	public Text MinValue;
	public Text MaxValue;
	public Text MidValue;
	public RawImage Gradient;

	private static Dictionary<WorldView.LegendType, string> LegendSuffixes = new Dictionary<WorldView.LegendType, string>
	{
		{ WorldView.LegendType.Mass, " kg" },
		{ WorldView.LegendType.Percent, "%" },
		{ WorldView.LegendType.PPM, " ppm" },
		{ WorldView.LegendType.Pressure, " Pa" },
		{ WorldView.LegendType.Volume, " m3" },
		{ WorldView.LegendType.Watts, " W" },
	};
	public void SetValues(WorldView.MeshOverlayColors c, CellInfo.TemperatureUnits temperatureUnits)
	{
		UpdateGradient(c.ColorValuePairs);
		float min = c.Min;
		float max = c.Max;
		string suffix = "";
		System.Globalization.NumberFormatInfo format = new System.Globalization.NumberFormatInfo() { NumberDecimalDigits = c.DecimalPlaces };
		if (c.LegendType == WorldView.LegendType.Temperature)
		{
			min = CellInfo.ConvertTemperature(min, temperatureUnits);
			max = CellInfo.ConvertTemperature(max, temperatureUnits);
			float freezing = CellInfo.ConvertTemperature(WorldData.FreezingTemperature, temperatureUnits);
			if (min < freezing && max > freezing)
			{
				MidValue.gameObject.SetActive(true);
				MidValue.rectTransform.anchoredPosition = new Vector2((freezing - min) / (max - min) * Gradient.rectTransform.rect.width, MidValue.rectTransform.anchoredPosition.y);
				MidValue.text = freezing.ToString("N", format) + suffix;
			}
			else
			{
				MidValue.gameObject.SetActive(false);
			}
			switch (temperatureUnits)
			{
				case CellInfo.TemperatureUnits.Celsius:
					suffix = "°C";
					break;
				case CellInfo.TemperatureUnits.Farenheit:
					suffix = "°F";
					break;
				case CellInfo.TemperatureUnits.Kelvin:
					suffix = "°K";
					break;
			}
		} else
		{
			if (c.LegendType == WorldView.LegendType.Percent)
			{
				min *= 100;
				max *= 100;
			}
			MidValue.gameObject.SetActive(false);
			if (!LegendSuffixes.TryGetValue(c.LegendType, out suffix))
			{
				suffix = "";
			}
		}
		MinValue.text = min.ToString("N", format) + suffix;
		MaxValue.text = max.ToString("N", format) + suffix;
	}
	static string[] MaterialValueNames = new string[]
	{
			"Vector1_765A4984",
			"Vector1_8A3EB1E6",
			"Vector1_40F3F1D3",
			"Vector1_7474FEB7",
			"Vector1_AD62CE61",
			"Vector1_E2411E1",
			"Vector1_B3E621B0",
	};
	static string[] MaterialColorNames = new string[]
	{
			"Color_66F70B0C",
			"Color_257472EE",
			"Color_B9142184",
			"Color_1533FF50",
			"Color_3E1D91D8",
			"Color_5CCB69B5",
			"Color_56BAA9D2",
	};
	private void UpdateGradient(NativeArray<CVP> colors)
	{
		int colorCount = Math.Min(7, colors.Length);
		Gradient.material.SetInt("Vector1_D53C5A56", colorCount);
		float min = colors[0].Value;
		float range = colors[colors.Length - 1].Value - min;
		for (int i = 0; i < colorCount; i++)
		{
			Gradient.material.SetFloat(MaterialValueNames[i], (colors[i].Value - min) / range);
			Gradient.material.SetColor(MaterialColorNames[i], colors[i].Color);
		}
	}
}
