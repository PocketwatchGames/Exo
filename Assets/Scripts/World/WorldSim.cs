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
	public StaticState StaticState;
	public Icosphere Icosphere;
	public float TimeScale;

	public int CellCount { get; private set; }
	public ref SimState ActiveSimState {  get { return ref _simStates[_activeSimState]; } }
	public float TimeTillTick { get { return _timeTillTick; } private set { _timeTillTick = value; } }

	private WorldGenData _worldGenData = new WorldGenData();
	private const int _simStateCount = 2;
	private SimState[] _simStates;
	private int _activeSimState;
	public float _timeTillTick = 0.00001f;
	private float _ticksPerSecond = 1;

	public delegate void TickEventHandler();
	public event TickEventHandler OnTick;

	public void Start()
    {
		Icosphere = new Icosphere(Subdivisions);

		CellCount = Icosphere.Vertices.Count;

		StaticState = new StaticState();
		StaticState.Init(CellCount, Icosphere);

		_activeSimState = 0;
		_simStates = new SimState[_simStateCount];
		for (int i = 0; i < _simStateCount; i++)
		{
			_simStates[i] = new SimState();
			_simStates[i].Init(CellCount);
		}

		_worldGenData = JsonUtility.FromJson<WorldGenData>(WorldGenAsset.text);
		WorldData = JsonUtility.FromJson<WorldData>(WorldDataAsset.text);
		WorldGen.Generate(StaticState, Seed, _worldGenData, WorldData, ref _simStates[0]);

	}

	public void OnDestroy()
	{
		StaticState.Dispose();
	}

	bool _simulating;
	public void Update()
	{
		if (TimeTillTick > -1)
		{
			TimeTillTick -= Time.deltaTime * TimeScale;
		}
		if (TimeTillTick <= 0)
		{
			int iterations = 0;
			while (TimeTillTick <= 0)
			{
				TimeTillTick += _ticksPerSecond;
				iterations++;
			}
			Tick(ref _simStates[_activeSimState], iterations);

		}
	}

	private void Tick(ref SimState state, int ticksToAdvance)
	{
		int nextStateIndex = (_activeSimState + 1) % _simStateCount;
		_activeSimState = nextStateIndex;
		ref var nextState = ref _simStates[_activeSimState];

		WorldTick.Tick(ref state, ref nextState, ticksToAdvance, StaticState, WorldData);

		OnTick?.Invoke();
	}

	public string GetCellInfo(CellInfoType cellInfoType)
	{
		return "AFA";
	}


	public void StepTime()
	{
		TimeScale = 0;
		TimeTillTick = 0;
	}


}
