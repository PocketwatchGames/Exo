using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class WorldComponent : MonoBehaviour
{
	public enum CellInfoType {
		Global,
		Cell,
		Energy,
		Atmosphere,
		Water,
		Terrain
	}


	public int Seed;
	public float Scale = 0.00005f;
	public int Subdivisions = 5;
	public TextAsset WorldGenAsset;
	public TextAsset WorldDataAsset;
	public WorldData WorldData;

	public Material m_TerrainMaterial;
	public Material m_WaterMaterial;
	public Material m_CloudMaterial;

	private MeshBuilder meshBuilder;
	GameObject m_TerrainMesh;
	GameObject m_WaterMesh;
	GameObject m_CloudMesh;


	private WorldGenData _worldGenData = new WorldGenData();
	

	private int _simVertCount;
	public int SimVertCount {  get { return _simVertCount; } }

	private StaticState _staticState;
	private SimState _state;
	private SimState _lastRenderState;
	private SimState _nextRenderState;
	private float _renderStateLerp;
	private SimState _activeSimState;


	public void Start()
    {
		meshBuilder.Init(Subdivisions);

		_simVertCount = meshBuilder.GetVertices().Count;

		_staticState = new StaticState();
		_staticState.Init(_simVertCount);

		_state = new SimState();
		_state.Init(_simVertCount);

		if (m_TerrainMesh)
		{
			GameObject.Destroy(m_TerrainMesh);
		}
		m_TerrainMesh = new GameObject("Terrain Mesh");
		m_TerrainMesh.transform.parent = gameObject.transform;
		var terrainFilter = m_TerrainMesh.AddComponent<MeshFilter>();
		var terrainSurfaceRenderer = m_TerrainMesh.AddComponent<MeshRenderer>();
		terrainSurfaceRenderer.material = m_TerrainMaterial;

		m_WaterMesh = new GameObject("Water Mesh");
		m_WaterMesh.transform.parent = gameObject.transform;
		var waterFilter = m_WaterMesh.AddComponent<MeshFilter>();
		var waterSurfaceRenderer = m_WaterMesh.AddComponent<MeshRenderer>();
		waterSurfaceRenderer.material = m_WaterMaterial;

		m_CloudMesh = new GameObject("Cloud Mesh");
		m_CloudMesh.transform.parent = gameObject.transform;
		var cloudFilter = m_CloudMesh.AddComponent<MeshFilter>();
		var cloudSurfaceRenderer = m_CloudMesh.AddComponent<MeshRenderer>();
		cloudSurfaceRenderer.material = m_CloudMaterial;


		var terrainMesh = new Mesh();
		var waterMesh = new Mesh();
		var cloudMesh = new Mesh();
		meshBuilder.InitMesh(terrainMesh, waterMesh, cloudMesh);
		terrainFilter.mesh = terrainMesh;
		waterFilter.mesh = waterMesh;
		cloudFilter.mesh = cloudMesh;

		_staticState.ExtractCoordinates(meshBuilder.GetVertices());
		_worldGenData = JsonUtility.FromJson<WorldGenData>(WorldGenAsset.text);
		WorldData = JsonUtility.FromJson<WorldData>(WorldDataAsset.text);
		WorldGen.Generate(_state, _staticState, Seed, _worldGenData, WorldData);

		_activeSimState = _state;
		_lastRenderState = _state;
		_nextRenderState = _state;
		meshBuilder.UpdateMesh(terrainMesh, waterMesh, cloudMesh, _staticState, _state, _state, 0, Scale);
		

	}


	public void Update()
	{
		var terrainMeshFilter = m_TerrainMesh.GetComponent<MeshFilter>();
		var waterMeshFilter = m_WaterMesh.GetComponent<MeshFilter>();
		var cloudMeshFilter = m_CloudMesh.GetComponent<MeshFilter>();
		if (terrainMeshFilter != null && waterMeshFilter != null && cloudMeshFilter != null)
		{
			meshBuilder.UpdateMesh(terrainMeshFilter.mesh, waterMeshFilter.mesh, cloudMeshFilter.mesh, _staticState, _lastRenderState, _nextRenderState, _renderStateLerp, Scale);
		}

		var quat = Quaternion.Euler(Mathf.LerpAngle(_lastRenderState.TiltAngle, _nextRenderState.TiltAngle, _renderStateLerp), Mathf.LerpAngle(_lastRenderState.SpinAngle, _nextRenderState.SpinAngle, _renderStateLerp), 0);
		transform.SetPositionAndRotation(Vector3.zero, quat);

		
	}

	public string GetCellInfo(CellInfoType cellInfoType)
	{
		return "AFA";
	}

	public void OnWaterDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		m_WaterMesh.SetActive(toggle.isOn);
	}
	public void OnCloudDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		m_CloudMesh.SetActive(toggle.isOn);
	}
	public void OnHUDOverlayChanged(UnityEngine.UI.Dropdown dropdown)
	{
		meshBuilder.ActiveMeshOverlay = (MeshBuilder.MeshOverlay)dropdown.value;
	}
	public void OnHUDWindChanged(UnityEngine.UI.Dropdown dropdown)
	{
		meshBuilder.ActiveWindOverlay = (MeshBuilder.WindOverlay)dropdown.value;

	}
}
