using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity;
using UnityEngine;

public class GameTool : MonoBehaviour {

	public string Name;
	public List<GameObject> HUDComponents;
	public GameplayManager Gameplay;
	virtual public void OnSelected()
	{
		foreach (var c in HUDComponents)
		{
			c?.SetActive(true);
		}
	}
	virtual public void OnDeselected()
	{
		foreach (var c in HUDComponents)
		{
			c?.SetActive(false);
		}
	}
	virtual public void OnClick(Vector3 worldPos, int cellIndex) { }
	virtual public void OnPointerDown(Vector3 worldPos, int cellIndex) { }
	virtual public void OnPointerUp(Vector3 worldPos, int cellIndex) { }
	virtual public void OnDragMove(Vector3 worldPos, int cellIndex, Vector2 direction) { }
	virtual public void OnUpdate(Vector3 worldPos, int cellIndex) { }

}
