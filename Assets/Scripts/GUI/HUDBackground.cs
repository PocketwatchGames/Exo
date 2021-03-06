﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.EventSystems;

public class HUDBackground : MonoBehaviour, IPointerClickHandler, IPointerDownHandler, IPointerUpHandler {

	public HUD hud;

    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }

	public void OnPointerDown(PointerEventData ped)
	{
		if (ped.button == PointerEventData.InputButton.Left)
		{
			hud.OnBackgroundMouseDown();
		}
	}
	public void OnPointerUp(PointerEventData ped)
	{
		if (ped.button == PointerEventData.InputButton.Left)
		{
			hud.OnBackgroundMouseUp();
		}
	}
	public void OnPointerClick(PointerEventData ped)
	{
		if (ped.button == PointerEventData.InputButton.Left)
		{
			hud.OnBackgroundClicked();
		}
	}

}
