using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CellInfoPanel : MonoBehaviour
{
	public UnityEngine.UI.Text TextObject;
	public WorldMesh WorldView;
	public WorldMesh.CellInfoType CellInfo;

	// Start is called before the first frame update
	void Start()
    {
    }

    // Update is called once per frame
    void Update()
    {
		TextObject.text = WorldView.GetCellInfo(CellInfo);
    }



}
