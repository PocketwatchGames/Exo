using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Entities;

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

	public Material TerrainMaterial;
	public Material WaterMaterial;
	public Material CloudMaterial;

	public float TimeScale;



	public ref SimState LastSimState {  get { return ref _simStates[_activeSimState]; } }

	GameObject _terrainMesh;
	GameObject _waterMesh;
	GameObject _cloudMesh;
	private WorldRenderer _meshBuilder;
	private WorldGenData _worldGenData = new WorldGenData();
	private const int _simStateCount = 2;
	private const int _renderStateCount = 3;
	private StaticState _staticState;
	private SimState[] _simStates;
	private RenderState[] _renderStates;
	private int _curRenderState;
	private int _lastRenderState;
	private int _nextRenderState;
	private float _renderStateLerp;
	private int _activeSimState;
	private float _timeTillTick = 0.00001f;
	private float _tickLerpTime;
	private float _ticksPerSecond = 1;


	public void Start()
    {
		_meshBuilder = new WorldRenderer();
		_meshBuilder.Init(Subdivisions);

		int cellCount = _meshBuilder.Icosphere.Vertices.Count;

		_staticState = new StaticState();
		_staticState.Init(cellCount);

		_activeSimState = 0;
		_curRenderState = 1;
		_nextRenderState = 0;
		_lastRenderState = 0;
		_simStates = new SimState[_simStateCount];
		for (int i = 0; i < _simStateCount; i++)
		{
			_simStates[i] = new SimState(cellCount, "Sim " + i);
			_simStates[i].World.AddSystem(_simStates[i].World.CreateSystem(typeof(WorldSimSystem)));
			_simStates[i].World.AddSystem(_simStates[i].World.CreateSystem(typeof(SimToRenderCellSystem)));
		}
		_renderStates = new RenderState[_renderStateCount];
		for (int i=0;i<_renderStateCount;i++)
		{
			_renderStates[i] = new RenderState(cellCount);
		}

		if (_terrainMesh)
		{
			GameObject.Destroy(_terrainMesh);
		}
		_terrainMesh = new GameObject("Terrain Mesh");
		_terrainMesh.transform.parent = gameObject.transform;
		var terrainFilter = _terrainMesh.AddComponent<MeshFilter>();
		var terrainSurfaceRenderer = _terrainMesh.AddComponent<MeshRenderer>();
		terrainSurfaceRenderer.material = TerrainMaterial;

		_waterMesh = new GameObject("Water Mesh");
		_waterMesh.transform.parent = gameObject.transform;
		var waterFilter = _waterMesh.AddComponent<MeshFilter>();
		var waterSurfaceRenderer = _waterMesh.AddComponent<MeshRenderer>();
		waterSurfaceRenderer.material = WaterMaterial;

		_cloudMesh = new GameObject("Cloud Mesh");
		_cloudMesh.transform.parent = gameObject.transform;
		var cloudFilter = _cloudMesh.AddComponent<MeshFilter>();
		var cloudSurfaceRenderer = _cloudMesh.AddComponent<MeshRenderer>();
		cloudSurfaceRenderer.material = CloudMaterial;


		_staticState.ExtractCoordinates(_meshBuilder.Icosphere.Vertices);
		_worldGenData = JsonUtility.FromJson<WorldGenData>(WorldGenAsset.text);
		WorldData = JsonUtility.FromJson<WorldData>(WorldDataAsset.text);
		WorldGen.Generate(_staticState, Seed, _worldGenData, WorldData, ref _simStates[0]);

		_meshBuilder.BuildRenderState(_staticState, _simStates[_activeSimState], _renderStates[_nextRenderState], Scale);


		var terrainMesh = new Mesh();
		var waterMesh = new Mesh();
		var cloudMesh = new Mesh();
		_meshBuilder.UpdateMesh(terrainMesh, waterMesh, cloudMesh, _staticState, _renderStates[_lastRenderState], Scale);
		_meshBuilder.InitMesh(terrainMesh, waterMesh, cloudMesh);
		terrainFilter.mesh = terrainMesh;
		waterFilter.mesh = waterMesh;
		cloudFilter.mesh = cloudMesh;


	}



	public void Update()
	{
		if (_timeTillTick > -1)
		{
			_timeTillTick -= Time.deltaTime * TimeScale;
		}
		if (_timeTillTick <= 0)
		{
			while (_timeTillTick <= 0)
			{
				DoSimTick();
			}
			_lastRenderState = _curRenderState;
			_nextRenderState = (_curRenderState + 1) % _renderStateCount;
			_curRenderState = (_nextRenderState + 1) % _renderStateCount;
			_meshBuilder.BuildRenderState(_staticState, _simStates[_activeSimState], _renderStates[_nextRenderState], Scale);
			_renderStateLerp = 0;
			_tickLerpTime = _timeTillTick;
		}
		ref var renderState = ref _renderStates[_curRenderState];
		if (_tickLerpTime > 0)
		{
			_renderStateLerp = Mathf.Clamp01(1.0f - _timeTillTick / _tickLerpTime);
		}
		_meshBuilder.LerpRenderState(_renderStates[_lastRenderState], _renderStates[_nextRenderState], _renderStateLerp, renderState);


		var terrainMeshFilter = _terrainMesh.GetComponent<MeshFilter>();
		var waterMeshFilter = _waterMesh.GetComponent<MeshFilter>();
		var cloudMeshFilter = _cloudMesh.GetComponent<MeshFilter>();
		if (terrainMeshFilter != null && waterMeshFilter != null && cloudMeshFilter != null)
		{
			_meshBuilder.UpdateMesh(terrainMeshFilter.mesh, waterMeshFilter.mesh, cloudMeshFilter.mesh, _staticState, _renderStates[_curRenderState], Scale);
		}

		var quat = Quaternion.Euler(renderState.TiltAngle, renderState.SpinAngle, 0);
		transform.SetPositionAndRotation(Vector3.zero, quat);

	}
	private void DoSimTick()
	{
		if (_timeTillTick <= 0)
		{
			_timeTillTick += _ticksPerSecond;

			int nextStateIndex = (_activeSimState + 1) % _simStateCount;
			//WorldSim.Tick(ref _simStates[_activeSimState], ref _simStates[nextStateIndex]);
			_simStates[nextStateIndex].Update(_simStates[_activeSimState]);
			_activeSimState = nextStateIndex;
		}

	}

	public string GetCellInfo(CellInfoType cellInfoType)
	{
		return "AAAAA";
	}

	public void OnWaterDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		_waterMesh.SetActive(toggle.isOn);
	}
	public void OnCloudDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		_cloudMesh.SetActive(toggle.isOn);
	}
	public void OnHUDOverlayChanged(UnityEngine.UI.Dropdown dropdown)
	{
		_meshBuilder.ActiveMeshOverlay = (WorldRenderer.MeshOverlay)dropdown.value;
	}
	public void OnHUDWindChanged(UnityEngine.UI.Dropdown dropdown)
	{
		_meshBuilder.ActiveWindOverlay = (WorldRenderer.WindOverlay)dropdown.value;
	}


	public void StepTime()
	{
		TimeScale = 0;
		_timeTillTick = 0;
	}

}
