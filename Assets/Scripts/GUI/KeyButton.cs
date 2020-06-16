using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.EventSystems;

public class KeyButton : MonoBehaviour
{
	public KeyCode Key;

    // Start is called before the first frame update
    void Start()
    {
        
    }

	// Update is called once per frame
	void Update()
    {
        if (Input.GetKeyDown(Key))
		{
			GetComponent<Selectable>().OnPointerDown(new PointerEventData(EventSystem.current));
		}
		if (Input.GetKeyUp(Key))
		{
			GetComponent<Selectable>()?.OnPointerUp(new PointerEventData(EventSystem.current));
			GetComponent<Button>()?.OnPointerClick(new PointerEventData(EventSystem.current));
		}
	}
}
