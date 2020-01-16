using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;
using Unity.Burst;

public class WorldMesh : MonoBehaviour {

	public enum MeshOverlay {
		None,
		TerrainTemperature,
		GroundWater,
		DeepWaterTemperature,
		DeepWaterSalinity,
		ShallowWaterTemperature,
		ShallowWaterSalinity,
		LowerAirTemperature,
		LowerAirPressure,
		UpperAirTemperature,
		UpperAirPressure,
		VerticalWind,
		AbsoluteHumidity,
		RelativeHumidity,
		CloudMass,
		CloudCoalescence,
		Evaporation,
		Rainfall,
		Condensation,
		HeatAbsorbed,
	}

	public enum WindOverlay {
		None,
		DeepWaterCurrent,
		ShallowWaterCurrent,
		LowerAirWind,
		UpperAirWind
	}

	public float DistanceToSun = 100;
	public float TerrainScale = 0.00005f;
	public MeshOverlay ActiveMeshOverlay;
	public WindOverlay ActiveWindOverlay;

	public Material TerrainMaterial;
	public Material WaterMaterial;
	public Material CloudMaterial;

	static Color32 green = new Color32(20, 255, 30, 255);
	static Color32 brown = new Color32(220, 150, 70, 255);
	static Color32 blue = new Color32(0, 0, 255, 255);
	static Color32 white = new Color32(255, 255, 255, 255);
	static Color32 black = new Color32(0, 0, 0, 255);

	public Icosphere Icosphere;

	private Mesh _terrainMesh;
	private Mesh _waterMesh;
	private Mesh _cloudMesh;

	private GameObject _terrainObject;
	private GameObject _waterObject;
	private GameObject _cloudObject;

	bool _indicesInitialized;
	int[] indices;

	public void Init(int subdivisions)
	{
		Icosphere = new Icosphere(subdivisions);

		if (_terrainMesh)
		{
			GameObject.Destroy(_terrainMesh);
		}
		_terrainMesh = new Mesh();
		_cloudMesh = new Mesh();
		_waterMesh = new Mesh();

		_terrainObject = new GameObject("Terrain Mesh");
		_terrainObject.transform.SetParent(gameObject.transform, false);
		var terrainFilter = _terrainObject.AddComponent<MeshFilter>();
		var terrainSurfaceRenderer = _terrainObject.AddComponent<MeshRenderer>();
		terrainSurfaceRenderer.material = TerrainMaterial;
		terrainFilter.mesh = _terrainMesh;

		_waterObject = new GameObject("Water Mesh");
		_waterObject.transform.SetParent(gameObject.transform, false);
		var waterFilter = _waterObject.AddComponent<MeshFilter>();
		var waterSurfaceRenderer = _waterObject.AddComponent<MeshRenderer>();
		waterSurfaceRenderer.material = WaterMaterial;
		waterFilter.mesh = _waterMesh;

		_cloudObject = new GameObject("Cloud Mesh");
		_cloudObject.transform.SetParent(gameObject.transform, false);
		var cloudFilter = _cloudObject.AddComponent<MeshFilter>();
		var cloudSurfaceRenderer = _cloudObject.AddComponent<MeshRenderer>();
		cloudSurfaceRenderer.material = CloudMaterial;
		cloudFilter.mesh = _cloudMesh;


		int indexCount = Icosphere.Polygons.Count * 3;
		indices = new int[indexCount];
		for (int i = 0; i < Icosphere.Polygons.Count; i++)
		{
			var poly = Icosphere.Polygons[i];
			indices[i * 3 + 0] = poly.m_Vertices[0];
			indices[i * 3 + 1] = poly.m_Vertices[1];
			indices[i * 3 + 2] = poly.m_Vertices[2];
		}

	}


	public void LerpRenderState(ref RenderState lastState, ref RenderState nextState, float t, ref RenderState state)
	{
		state.Ticks = (nextState.Ticks - lastState.Ticks) * t + lastState.Ticks;
		state.TiltAngle = Mathf.LerpAngle(lastState.TiltAngle, nextState.TiltAngle, t);
		for (int i = 0; i < state.CloudColor.Length; i++)
		{
			state.TerrainColor[i] = Color32.Lerp(lastState.TerrainColor[i], nextState.TerrainColor[i], t);
			state.TerrainPosition[i] = Vector3.Lerp(lastState.TerrainPosition[i], nextState.TerrainPosition[i], t);
			state.TerrainNormal[i] = Vector3.Lerp(lastState.TerrainNormal[i], nextState.TerrainNormal[i], t);
			state.WaterColor[i] = Color32.Lerp(lastState.WaterColor[i], nextState.WaterColor[i], t);
			state.WaterPosition[i] = Vector3.Lerp(lastState.WaterPosition[i], nextState.WaterPosition[i], t);
			state.WaterNormal[i] = Vector3.Lerp(lastState.WaterNormal[i], nextState.WaterNormal[i], t);
			state.CloudColor[i] = Color32.Lerp(lastState.CloudColor[i], nextState.CloudColor[i], t);
			state.CloudPosition[i] = Vector3.Lerp(lastState.CloudPosition[i], nextState.CloudPosition[i], t);
			state.CloudNormal[i] = Vector3.Lerp(lastState.CloudNormal[i], nextState.CloudNormal[i], t);
		}
	}

	public void BuildRenderState(ref SimState from, ref RenderState to, StaticState staticState)
	{
		to.Ticks = from.Ticks;
		to.OrbitSpeed = from.OrbitSpeed;
		to.SpinAngle = from.SpinAngle;
		to.SpinSpeed = from.SpinSpeed;
		to.TiltAngle = from.TiltAngle;
		for (int i = 0; i < from.Count; i++)
		{
			ref var fromCell = ref from.Cells[i];
			to.TerrainColor[i] = GetTerrainColor(staticState, fromCell);
			to.WaterColor[i] = GetWaterColor(staticState, fromCell);
			to.CloudColor[i] = GetCloudColor(staticState, fromCell);
			to.TerrainNormal[i] = Icosphere.Vertices[i];
			to.WaterNormal[i] = Icosphere.Vertices[i];
			to.CloudNormal[i] = Icosphere.Vertices[i];
			to.TerrainPosition[i] = Icosphere.Vertices[i] * (fromCell.Elevation * TerrainScale + staticState.Radius) / staticState.Radius;
			to.WaterPosition[i] = Icosphere.Vertices[i] * (fromCell.WaterElevation * TerrainScale + staticState.Radius) / staticState.Radius;
			to.CloudPosition[i] = Icosphere.Vertices[i] * (fromCell.CloudElevation * TerrainScale + staticState.Radius) / staticState.Radius;
		}

	}

	public void UpdateMesh(ref RenderState state)
	{
		_terrainMesh.vertices = state.TerrainPosition;
		_terrainMesh.normals = state.TerrainNormal;
		_terrainMesh.colors32 = state.TerrainColor;

		_waterMesh.vertices = state.WaterPosition;
		_waterMesh.normals = state.WaterNormal;
		_waterMesh.colors32 = state.WaterColor;

		_cloudMesh.vertices = state.CloudPosition;
		_cloudMesh.normals = state.CloudNormal;
		_cloudMesh.colors32 = state.CloudColor;

		if (!_indicesInitialized)
		{
			_terrainMesh.SetTriangles(indices, 0);
			_waterMesh.SetTriangles(indices, 0);
			_cloudMesh.SetTriangles(indices, 0);
			_indicesInitialized = true;
		}

		_terrainMesh.RecalculateBounds();
		_waterMesh.RecalculateBounds();
		_cloudMesh.RecalculateBounds();

		_terrainMesh.RecalculateNormals();
		_waterMesh.RecalculateNormals();
		_cloudMesh.RecalculateNormals();

		var time = WorldTime.GetTime(state.Ticks, state.SpinSpeed);
		var spin = Quaternion.Euler(state.TiltAngle, time * 360, 0);
		float orbitAngle = state.Ticks * state.OrbitSpeed * 360;
		transform.SetPositionAndRotation(new Vector3(math.cos(orbitAngle), 0, math.sin(orbitAngle)) * DistanceToSun, spin);
	}

	public void OnWaterDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		_waterObject.SetActive(toggle.isOn);
	}
	public void OnCloudDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		_cloudObject.SetActive(toggle.isOn);
	}
	public void OnHUDOverlayChanged(UnityEngine.UI.Dropdown dropdown)
	{
		ActiveMeshOverlay = (WorldMesh.MeshOverlay)dropdown.value;
	}
	public void OnHUDWindChanged(UnityEngine.UI.Dropdown dropdown)
	{
		ActiveWindOverlay = (WorldMesh.WindOverlay)dropdown.value;
	}


	#region private functions

	private Color32 GetTerrainColor(StaticState staticState, SimStateCell cell)
	{
		return Color32.Lerp(brown, green, cell.Vegetation);
	}

	private Color32 GetWaterColor(StaticState staticState, SimStateCell cell)
	{
		return Color32.Lerp(blue, white, cell.Ice);
	}

	private Color32 GetCloudColor(StaticState staticState, SimStateCell cell)
	{
		var c = Color32.Lerp(white, black, cell.RelativeHumidity);
		float opacity = -math.cos(cell.CloudCoverage * math.PI)/ 2 + 0.5f;
		c.a = (byte)(255 * opacity);
		return c;
	}

	#endregion
}
