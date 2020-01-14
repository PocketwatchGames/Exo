using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class CellInfoToolbarButton : ToolbarButton
{
	public CellInfoPanel Panel;

    // Start is called before the first frame update
    public override void Start()
    {
		base.Start();
		Panel.gameObject.SetActive(Active);
	}

	// Update is called once per frame
	void Update()
    {
    }

	public override void OnClick()
	{
		base.OnClick();
		Panel.gameObject.SetActive(Active);
	}
}
