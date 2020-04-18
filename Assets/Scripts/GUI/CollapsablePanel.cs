using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
public class CollapsablePanel : MonoBehaviour
{
	public Text Title;
	public Button Button;
	public GameObject Contents;
	public GameObject TitleBar;
	public string OnText;
	public string OffText;
	private bool Active;
	private float _contentSize;


    // Start is called before the first frame update
    void Start()
    {
		_contentSize = GetComponent<RectTransform>().rect.height;
		UpdateSize();
	}

	// Update is called once per frame
	void Update()
    {
	}

	private void UpdateButtonText()
	{
		Button.GetComponentInChildren<Text>().text = Active ? OnText : OffText;
	}

	public void OnCollapseClicked()
	{
		Active = !Active;
		UpdateSize();
	}

	private void UpdateSize()
	{
		Contents.SetActive(Active);
		UpdateButtonText();
		float height;
		if (Active)
		{
			height = _contentSize;
		} else
		{
			height = TitleBar.GetComponent<RectTransform>().rect.height;
		}
		GetComponent<RectTransform>().SetSizeWithCurrentAnchors(RectTransform.Axis.Vertical, height);
	}
}
