using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class WorldComponent : MonoBehaviour
{
	public WorldMesh mesh;
	public Material m_Material;
	Mesh terrainMesh;
	GameObject m_TerrainMesh;


	public WorldGenData WorldGenData = new WorldGenData();
	public int Seed;

	private float _scale = 0.00001f;

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

		_simVertCount = Utils.Sqr(Utils.Sqr(Subdivisions)) * 10 + 2;

		_staticState = new StaticState();
		_staticState.Init(_simVertCount);

		_state = new SimState();
		_state.Init(_simVertCount);

		WorldGen.Generate(_state, _staticState, Seed, WorldGenData);

		_activeSimState = _state;
		_lastRenderState = _state;
		_nextRenderState = _state;

		if (m_TerrainMesh)
		{
			GameObject.Destroy(m_TerrainMesh);
		}
		m_TerrainMesh = new GameObject("Terrain Mesh");
		m_TerrainMesh.transform.parent = gameObject.transform;
		MeshFilter terrainFilter = m_TerrainMesh.AddComponent<MeshFilter>();
		MeshRenderer surfaceRenderer = m_TerrainMesh.AddComponent<MeshRenderer>();
		surfaceRenderer.material = m_Material;

		Mesh terrainMesh = new Mesh();
		mesh.InitMesh(terrainMesh);
		terrainFilter.mesh = terrainMesh;

		mesh.UpdateMesh(terrainFilter, _staticState, _state, _state, 0, _scale);

		_staticState.ExtractCoordinates(mesh.GetVertices());

	}


	public void Update()
	{
		var terrainMeshFilter = m_TerrainMesh.GetComponent<MeshFilter>();
		if (terrainMeshFilter != null)
		{
			mesh.UpdateMesh(terrainMeshFilter, _staticState, _lastRenderState, _nextRenderState, _renderStateLerp, _scale);
		}
	}

}
