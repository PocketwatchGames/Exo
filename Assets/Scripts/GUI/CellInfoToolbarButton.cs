using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class CellInfoToolbarButton : MonoBehaviour
{
	public bool Active;
	public CellInfoPanel Panel;
	public Text ButtonText;

    // Start is called before the first frame update
    void Start()
    {
		Panel.gameObject.SetActive(Active);
		ButtonText.fontStyle = Active ? FontStyle.Bold : FontStyle.Normal;
		ButtonText.color = Active ? Color.white : Color.black;
	}

	// Update is called once per frame
	void Update()
    {
    }

	public void OnClick()
	{
		Active = !Active;
		ButtonText.fontStyle = Active ? FontStyle.Bold : FontStyle.Normal;
		ButtonText.color = Active ? Color.white : Color.black;
		Panel.gameObject.SetActive(Active);
	}
}
