using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class GameToolButton : MonoBehaviour
{
	public GameplayManager Gameplay;
	public GameTool Tool;

	public void Awake()
	{
		var textComponent = GetComponentInChildren<Text>();
		textComponent.text = Tool.Name;
	}


	public void Update()
	{
		
	}

	public void OnGameToolChanged()
	{
		Gameplay.SetActiveTool(Tool);
	}


}
