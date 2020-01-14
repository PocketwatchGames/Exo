using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DropdownOverlay<T> : MonoBehaviour where T : struct, System.IConvertible
{

    // Start is called before the first frame update
    void Start()
    {
		var dd = GetComponent<UnityEngine.UI.Dropdown>();
		dd.options.Clear();
		foreach (var e in System.Enum.GetValues(typeof(T)))
		{
			dd.options.Add(new UnityEngine.UI.Dropdown.OptionData(e.ToString()));
		}
	}

    // Update is called once per frame
    void Update()
    {
        
    }
}

