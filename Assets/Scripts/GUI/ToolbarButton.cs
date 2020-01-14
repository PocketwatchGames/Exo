using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.UI;

public class ToolbarButton : MonoBehaviour {
	public bool Active;
	public Text ButtonText;

	// Start is called before the first frame update
	virtual public void Start()
	{
		ButtonText.fontStyle = Active ? FontStyle.Bold : FontStyle.Normal;
		ButtonText.color = Active ? Color.white : Color.black;
	}

	virtual public void OnClick()
	{
		Active = !Active;
		ButtonText.fontStyle = Active ? FontStyle.Bold : FontStyle.Normal;
		ButtonText.color = Active ? Color.white : Color.black;
	}
}
