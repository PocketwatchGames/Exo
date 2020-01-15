using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Jobs;
using Unity.Collections;

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



	private int _simVertCount;
	public int SimVertCount {  get { return _simVertCount; } }
	public ref SimState ActiveSimState {  get { return ref _simStates[_activeSimState]; } }

	GameObject _terrainMesh;
	GameObject _waterMesh;
	GameObject _cloudMesh;
	private MeshBuilder _meshBuilder;
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
		if (_terrainMesh)
		{
			GameObject.Destroy(_terrainMesh);
		}
		_terrainMesh = new GameObject("Terrain Mesh");
		_terrainMesh.transform.parent = gameObject.transform;
		var terrainFilter = _terrainMesh.AddComponent<MeshFilter>();
		var terrainSurfaceRenderer = _terrainMesh.AddComponent<MeshRenderer>();
		terrainSurfaceRenderer.material = TerrainMaterial;
		terrainFilter.mesh = new Mesh();

		_waterMesh = new GameObject("Water Mesh");
		_waterMesh.transform.parent = gameObject.transform;
		var waterFilter = _waterMesh.AddComponent<MeshFilter>();
		var waterSurfaceRenderer = _waterMesh.AddComponent<MeshRenderer>();
		waterSurfaceRenderer.material = WaterMaterial;
		waterFilter.mesh = new Mesh();

		_cloudMesh = new GameObject("Cloud Mesh");
		_cloudMesh.transform.parent = gameObject.transform;
		var cloudFilter = _cloudMesh.AddComponent<MeshFilter>();
		var cloudSurfaceRenderer = _cloudMesh.AddComponent<MeshRenderer>();
		cloudSurfaceRenderer.material = CloudMaterial;
		cloudFilter.mesh = new Mesh();

		_meshBuilder = new MeshBuilder(Subdivisions, terrainFilter.mesh, waterFilter.mesh, cloudFilter.mesh);

		_simVertCount = _meshBuilder.GetVertices().Count;

		_staticState = new StaticState();
		_staticState.Init(_simVertCount);

		_activeSimState = 0;
		_simStates = new SimState[_simStateCount];
		for (int i = 0; i < _simStateCount; i++)
		{
			_simStates[i] = new SimState();
			_simStates[i].Init(_simVertCount);
		}

		_nextRenderState = 0;
		_lastRenderState = 0;
		_curRenderState = 1;
		_renderStates = new RenderState[_renderStateCount];
		for (int i = 0; i < _renderStateCount; i++)
		{
			_renderStates[i] = new RenderState();
			_renderStates[i].Init(_simVertCount);
		}

		_staticState.ExtractCoordinates(_meshBuilder.GetVertices());
		_worldGenData = JsonUtility.FromJson<WorldGenData>(WorldGenAsset.text);
		WorldData = JsonUtility.FromJson<WorldData>(WorldDataAsset.text);
		WorldGen.Generate(_staticState, Seed, _worldGenData, WorldData, ref _simStates[0]);

		_meshBuilder.BuildRenderState(ref _simStates[0], ref _renderStates[0], _staticState, Scale);
		_meshBuilder.UpdateMesh(ref _renderStates[0]);

	}

	bool _simulating;
	public void Update()
	{
		var terrainMeshFilter = _terrainMesh.GetComponent<MeshFilter>()?.mesh;
		var waterMeshFilter = _waterMesh.GetComponent<MeshFilter>()?.mesh;
		var cloudMeshFilter = _cloudMesh.GetComponent<MeshFilter>()?.mesh;
		if (_timeTillTick > -1)
		{
			_timeTillTick -= Time.deltaTime * TimeScale;
		}
		if (_timeTillTick <= 0)
		{
			int iterations = 0;
			while (_timeTillTick <= 0)
			{
				_timeTillTick += _ticksPerSecond;
				iterations++;
			}
			Tick(ref _simStates[_activeSimState], iterations);
			_meshBuilder.UpdateMesh(ref _renderStates[_lastRenderState]);
		}
		if (_tickLerpTime > 0)
		{
			_renderStateLerp = Mathf.Clamp01(1.0f - _timeTillTick / _tickLerpTime);
			_meshBuilder.LerpRenderState(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], _renderStateLerp, ref _renderStates[_curRenderState]);
			_meshBuilder.UpdateMesh(ref _renderStates[_curRenderState]);
		}
		var quat = Quaternion.Euler(_renderStates[_curRenderState].TiltAngle, _renderStates[_curRenderState].SpinAngle, 0);
		transform.SetPositionAndRotation(Vector3.zero, quat);

	}

	private void Tick(ref SimState state, int ticksToAdvance)
	{
		JobHandle lastJobHandle = default(JobHandle);
		var cells = new NativeArray<SimStateCell>(state.Cells, Allocator.TempJob);
		int ticks = state.Ticks;		
		for (int i=0;i<ticksToAdvance;i++)
		{
			ticks++;
			var tickJob = new WorldSim.TickCellJob();
			tickJob.Cells = cells;
			tickJob.Ticks = ticks;
			lastJobHandle = tickJob.Schedule(state.Count, 100, lastJobHandle);
		}
		lastJobHandle.Complete();
		
		int nextStateIndex = (_activeSimState + 1) % _simStateCount;
		_activeSimState = nextStateIndex;

		ref var nextState = ref _simStates[_activeSimState];
		nextState.Ticks = ticks;
		nextState.Gravity = state.Gravity;
		nextState.OrbitSpeed = state.OrbitSpeed;
		nextState.SpinAngle = state.SpinAngle;
		nextState.SpinSpeed = state.SpinSpeed;
		nextState.TiltAngle = state.TiltAngle;
		for (int i = 0; i < nextState.Count; i++)
		{
			nextState.Cells[i] = cells[i];
		}
		cells.Dispose();


		_lastRenderState = _curRenderState;
		_nextRenderState = (_curRenderState + 1) % _renderStateCount;
		_curRenderState = (_nextRenderState + 1) % _renderStateCount;
		_meshBuilder.BuildRenderState(ref _simStates[_activeSimState], ref _renderStates[_nextRenderState], _staticState, Scale);
		_renderStateLerp = 0;
		_tickLerpTime = _timeTillTick;


	}

	public string GetCellInfo(CellInfoType cellInfoType)
	{
		return "AFA";
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
		_meshBuilder.ActiveMeshOverlay = (MeshBuilder.MeshOverlay)dropdown.value;
	}
	public void OnHUDWindChanged(UnityEngine.UI.Dropdown dropdown)
	{
		_meshBuilder.ActiveWindOverlay = (MeshBuilder.WindOverlay)dropdown.value;
	}


	public void StepTime()
	{
		TimeScale = 0;
		_timeTillTick = 0;
	}

}
