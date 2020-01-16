using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Jobs;
using Unity.Collections;

public class WorldSim : MonoBehaviour
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
	public int Subdivisions = 5;
	public TextAsset WorldGenAsset;
	public TextAsset WorldDataAsset;
	public WorldData WorldData;
	public WorldMesh Mesh;


	public float TimeScale;



	public int SimVertCount { get; private set; }
	public ref SimState ActiveSimState {  get { return ref _simStates[_activeSimState]; } }

	private WorldGenData _worldGenData = new WorldGenData();
	private const int _simStateCount = 2;
	private const int _renderStateCount = 3;
	private StaticState _staticState;
	private SimState[] _simStates;
	private RenderState[] _renderStates;
	private int _curRenderState;
	private int _lastRenderState;
	private int _nextRenderState;
	private int _activeSimState;
	private float _timeTillTick = 0.00001f;
	private float _tickLerpTime;
	private float _ticksPerSecond = 1;


	public void Start()
    {
		Mesh.Init(Subdivisions);

		SimVertCount = Mesh.Icosphere.Vertices.Count;

		_staticState = new StaticState();
		_staticState.Init(SimVertCount, Mesh.Icosphere);

		_activeSimState = 0;
		_simStates = new SimState[_simStateCount];
		for (int i = 0; i < _simStateCount; i++)
		{
			_simStates[i] = new SimState();
			_simStates[i].Init(SimVertCount);
		}

		_nextRenderState = 0;
		_lastRenderState = 0;
		_curRenderState = 0;
		_renderStates = new RenderState[_renderStateCount];
		for (int i = 0; i < _renderStateCount; i++)
		{
			_renderStates[i] = new RenderState();
			_renderStates[i].Init(SimVertCount);
		}

		_worldGenData = JsonUtility.FromJson<WorldGenData>(WorldGenAsset.text);
		WorldData = JsonUtility.FromJson<WorldData>(WorldDataAsset.text);
		WorldGen.Generate(_staticState, Seed, _worldGenData, WorldData, ref _simStates[0]);

		Mesh.BuildRenderState(ref _simStates[0], ref _renderStates[0], _staticState);
		Mesh.UpdateMesh(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], 0, ref _renderStates[_curRenderState]);

	}

	public void OnDestroy()
	{
		_staticState.Dispose();
	}

	bool _simulating;
	public void Update()
	{
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
		}
		if (_tickLerpTime > 0)
		{

			float renderStateLerp = Mathf.Clamp01(1.0f - _timeTillTick / _tickLerpTime);
			Mesh.UpdateMesh(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], renderStateLerp, ref _renderStates[_curRenderState]);
		}

	}

	private void Tick(ref SimState state, int ticksToAdvance)
	{
		int nextStateIndex = (_activeSimState + 1) % _simStateCount;
		_activeSimState = nextStateIndex;
		ref var nextState = ref _simStates[_activeSimState];

		WorldTick.Tick(ref state, ref nextState, ticksToAdvance, _staticState, WorldData);

		_lastRenderState = _curRenderState;
		_nextRenderState = (_curRenderState + 1) % _renderStateCount;
		_curRenderState = (_nextRenderState + 1) % _renderStateCount;
		Mesh.BuildRenderState(ref _simStates[_activeSimState], ref _renderStates[_nextRenderState], _staticState);
		_tickLerpTime = _timeTillTick;

	}

	public string GetCellInfo(CellInfoType cellInfoType)
	{
		return "AFA";
	}


	public void StepTime()
	{
		TimeScale = 0;
		_timeTillTick = 0;
	}

}
