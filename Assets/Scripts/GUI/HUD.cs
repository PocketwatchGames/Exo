using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class HUD : MonoBehaviour
{
	public WorldView View;

    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
		
		if (Input.GetMouseButtonUp(0))
		{
			int cellHovered = GetMouseCellIndex();
			View.SetActiveCell(cellHovered, cellHovered >= 0);
		} else if (!View.ActiveCellLocked)
		{
			View.SetActiveCell(GetMouseCellIndex(), false);
		}
	}

	private int GetMouseCellIndex()
	{
		var ray = Camera.main.ScreenPointToRay(Input.mousePosition);
		RaycastHit hit;
		if (Physics.Raycast(ray, out hit))
		{
			if (hit.collider.gameObject.transform.parent == View.Planet.transform)
			{
				int tIndex;
				if (hit.barycentricCoordinate.x > hit.barycentricCoordinate.y && hit.barycentricCoordinate.x > hit.barycentricCoordinate.z)
				{
					tIndex = 0;
				}
				else if (hit.barycentricCoordinate.y > hit.barycentricCoordinate.z)
				{
					tIndex = 1;
				}
				else
				{
					tIndex = 2;
				}
				return View.GetClosestVert(hit.triangleIndex, tIndex);
			}
		}
		return -1;
	}
}
