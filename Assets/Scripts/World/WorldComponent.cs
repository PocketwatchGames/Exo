using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class WorldComponent : MonoBehaviour
{
	public WorldMesh mesh;
	public Material m_TerrainMaterial;
	public Material m_WaterMaterial;
	public Material m_CloudMaterial;
	GameObject m_TerrainMesh;
	GameObject m_WaterMesh;
	GameObject m_CloudMesh;


	public WorldGenData WorldGenData = new WorldGenData();
	public int Seed;

	private float _scale = 0.00005f;

	public int Subdivisions = 5;
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
		mesh.Init(Subdivisions);

		_simVertCount = mesh.GetVertices().Count;

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
		mesh.InitMesh(terrainMesh, waterMesh, cloudMesh);
		terrainFilter.mesh = terrainMesh;
		waterFilter.mesh = waterMesh;
		cloudFilter.mesh = cloudMesh;

		_staticState.ExtractCoordinates(mesh.GetVertices());
		WorldGen.Generate(_state, _staticState, Seed, WorldGenData);

		_activeSimState = _state;
		_lastRenderState = _state;
		_nextRenderState = _state;
		mesh.UpdateMesh(terrainMesh, waterMesh, cloudMesh, _staticState, _state, _state, 0, _scale);
		

	}


	public void Update()
	{
		var terrainMeshFilter = m_TerrainMesh.GetComponent<MeshFilter>();
		var waterMeshFilter = m_WaterMesh.GetComponent<MeshFilter>();
		var cloudMeshFilter = m_CloudMesh.GetComponent<MeshFilter>();
		if (terrainMeshFilter != null && waterMeshFilter != null && cloudMeshFilter != null)
		{
			mesh.UpdateMesh(terrainMeshFilter.mesh, waterMeshFilter.mesh, cloudMeshFilter.mesh, _staticState, _lastRenderState, _nextRenderState, _renderStateLerp, _scale);
		}
	}

}
