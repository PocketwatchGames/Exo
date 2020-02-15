using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;

public class WorldView : MonoBehaviour {

	public enum MeshOverlay {
		None,
		GroundTemperature,
		DeepWaterTemperature,
		DeepWaterSalinity,
		ShallowWaterTemperature,
		ShallowWaterSalinity,
		LowerAirTemperature,
		UpperAirTemperature,
		Pressure,
		VerticalWind,
		AbsoluteHumidityLower,
		RelativeHumidityLower,
		AbsoluteHumidityUpper,
		CloudMass,
		CloudTemperature,
		CloudDropletMass,
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

	public enum TemperatureUnits {
		Celsius,
		Farenheit,
		Kelvin,
	}

	public enum CellInfoType {
		Global,
		Cell,
		Energy,
		Atmosphere,
		Water,
		Ground
	}

	public bool ActiveCellLocked { get; private set; }
	public int ActiveCellIndex { get; private set; }
	public TemperatureUnits ActiveTemperatureUnits = TemperatureUnits.Celsius;
	public bool LerpStates = true;


	[Header("Display")]
	public float TerrainScale = 100f;
//	public TemperatureDisplayType TemperatureDisplay;
	public float MinElevation = -11000;
	public float MaxElevation = 10000;
	public float MaxDepth = 11000;
	public float maxCloudColor = 300.0f;
	public float WaterDepthThreshold = 10;

	public float DisplayWindMax = 100;
	public float DisplayCurrentMax = 10;
	public float DisplayCurrentMinDepth = 50;
	public float DisplayEnergyAborsobedMax = 300;
	public float DisplayRainfallMax = 5.0f;
	public float DisplaySalinityMin = 0;
	public float DisplaySalinityMax = 50;
	public float DisplayWindSpeedLowerAirMax = 50;
	public float DisplayWindSpeedUpperAirMax = 250;
	public float DisplayWindSpeedSurfaceWaterMax = 5;
	public float DisplayWindSpeedDeepWaterMax = 0.5f;
	public float DisplayVerticalWindSpeedMax = 1.0f;
	public float DisplayEvaporationMax = 5.0f;
	public float DisplayVegetationMax = 1000;
	public float DisplayTemperatureMin = 223;
	public float DisplayTemperatureMax = 323;
	public float DisplayAbsoluteHumidityMax = 0.05f;
	public float DisplayAirPressureMin = 97000;
	public float DisplayAirPressureMax = 110000;
	public float DisplayHeatAbsorbedMax = 1000;

	[Header("References")]
	public WorldSimComponent Sim;
	public GameObject Planet;
	public GameObject Sun;
	public GameObject Moon;
	public MeshOverlay ActiveMeshOverlay;
	public WindOverlay ActiveWindOverlay;
	public Material TerrainMaterial;
	public Material WaterMaterial;
	public Material CloudMaterial;
	public GameObject SelectionCirclePrefab;
	public GameObject WindArrowPrefab;


	private RenderState[] _renderStates;
	private int _curRenderState;
	private int _lastRenderState;
	private int _nextRenderState;
	private const int _renderStateCount = 3;

	private Vector3[] _terrainVertices;
	private Vector3[] _terrainNormals;
	private Color32[] _terrainColors;
	private Vector3[] _waterVertices;
	private Vector3[] _waterNormals;
	private Color32[] _waterColors;
	private Vector3[] _cloudVertices;
	private Vector3[] _cloudNormals;
	private Color32[] _cloudColors;

	private Mesh _terrainMesh;
	private Mesh _waterMesh;
	private Mesh _cloudMesh;

	private GameObject _terrainObject;
	private GameObject _waterObject;
	private GameObject _cloudObject;
	private GameObject _selectionCircle;
	private List<GameObject> _selectionCircleNeighbors;

	private GameObject[] _windArrows;

	private float _tickLerpTime = 0;
	private float _tickLerpTimeTotal = 1;

	private bool _indicesInitialized;
	private int[] indices;

	private NativeArray<CVP> _normalizedRainbow;
	private NativeArray<CVP> _normalizedBlueBlackRed;

	private int _batchCount = 128;

	public void Start()
	{
		Sim.OnTick += OnSimTick;

		_normalizedRainbow = new NativeArray<CVP>(new CVP[] {
											new CVP(Color.black, 0),
											new CVP(Color.white, 0.1667f),
											new CVP(Color.blue, 0.3333f),
											new CVP(Color.green, 0.5f),
											new CVP(Color.yellow, 0.6667f),
											new CVP(Color.red, 0.8333f),
											new CVP(Color.magenta, 1) },
											Allocator.Persistent);
		_normalizedBlueBlackRed = new NativeArray<CVP>( new CVP[] {
											new CVP(Color.blue, 0),
											new CVP(Color.black, 0.5f),
											new CVP(Color.red, 1) },
											Allocator.Persistent);


		if (_terrainMesh)
		{
			GameObject.Destroy(_terrainMesh);
		}
		_terrainMesh = new Mesh();
		_cloudMesh = new Mesh();
		_waterMesh = new Mesh();

		_terrainObject = new GameObject("Terrain Mesh");
		_terrainObject.transform.SetParent(Planet.transform, false);
		var terrainFilter = _terrainObject.AddComponent<MeshFilter>();
		var terrainSurfaceRenderer = _terrainObject.AddComponent<MeshRenderer>();
		var terrainCollider = _terrainObject.AddComponent<MeshCollider>();
		terrainSurfaceRenderer.material = TerrainMaterial;
		terrainFilter.mesh = terrainCollider.sharedMesh = _terrainMesh;

		_waterObject = new GameObject("Water Mesh");
		_waterObject.transform.SetParent(Planet.transform, false);
		var waterFilter = _waterObject.AddComponent<MeshFilter>();
		var waterSurfaceRenderer = _waterObject.AddComponent<MeshRenderer>();
		var waterCollider = _waterObject.AddComponent<MeshCollider>();
		waterSurfaceRenderer.material = WaterMaterial;
		waterFilter.mesh = waterCollider.sharedMesh = _waterMesh;

		_cloudObject = new GameObject("Cloud Mesh");
		_cloudObject.transform.SetParent(Planet.transform, false);
		var cloudFilter = _cloudObject.AddComponent<MeshFilter>();
		var cloudSurfaceRenderer = _cloudObject.AddComponent<MeshRenderer>();
		cloudSurfaceRenderer.material = CloudMaterial;
		cloudFilter.mesh = _cloudMesh;

		_selectionCircle = GameObject.Instantiate(SelectionCirclePrefab, Planet.transform);
		_selectionCircle.transform.localScale *= 0.02f;

		_selectionCircleNeighbors = new List<GameObject>();
		for (int i=0;i<6;i++)
		{
			var s = GameObject.Instantiate(SelectionCirclePrefab, Planet.transform);
			s.transform.localScale *= 0.01f;
			_selectionCircleNeighbors.Add(s);
		}

		_windArrows = new GameObject[Sim.CellCount];
		for (int i = 0; i < Sim.CellCount;i++)
		{
			_windArrows[i] = GameObject.Instantiate(WindArrowPrefab, Planet.transform);
			_windArrows[i].SetActive(false);
			_windArrows[i].hideFlags |= HideFlags.HideInHierarchy;
		}

		_nextRenderState = 0;
		_lastRenderState = 0;
		_curRenderState = 1;
		_renderStates = new RenderState[_renderStateCount];
		for (int i = 0; i < _renderStateCount; i++)
		{
			_renderStates[i] = new RenderState();
			_renderStates[i].Init(Sim.CellCount);
		}

		int indexCount = Sim.Icosphere.Polygons.Count * 3;
		indices = new int[indexCount];
		for (int i = 0; i < Sim.Icosphere.Polygons.Count; i++)
		{
			var poly = Sim.Icosphere.Polygons[i];
			indices[i * 3 + 0] = poly.m_Vertices[0];
			indices[i * 3 + 1] = poly.m_Vertices[1];
			indices[i * 3 + 2] = poly.m_Vertices[2];
		}

		_terrainVertices = new Vector3[Sim.CellCount];
		_terrainNormals = new Vector3[Sim.CellCount];
		_terrainColors = new Color32[Sim.CellCount];
		_waterVertices = new Vector3[Sim.CellCount];
		_waterNormals = new Vector3[Sim.CellCount];
		_waterColors = new Color32[Sim.CellCount];
		_cloudVertices = new Vector3[Sim.CellCount];
		_cloudNormals = new Vector3[Sim.CellCount];
		_cloudColors = new Color32[Sim.CellCount];

		BuildRenderState(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState, ref _renderStates[0], ref Sim.WorldData, ref Sim.StaticState);
		UpdateMesh(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], ref _renderStates[_curRenderState]);
	}
	public void OnDestroy()
	{
		Sim.OnTick -= OnSimTick;

		_normalizedRainbow.Dispose();
		_normalizedBlueBlackRed.Dispose();

		for (int i=0;i<_renderStateCount;i++)
		{
			_renderStates[i].Dispose();
		}

	}

	public void Update()
	{
		if (Sim.TimeScale == 0)
		{
			_tickLerpTime -= Time.deltaTime * 20;
		} else
		{
			_tickLerpTime -= Time.deltaTime * Sim.TimeScale;
		}
		UpdateMesh(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], ref _renderStates[_curRenderState]);

	}
	public void StartLerp(float lerpTime)
	{
		_tickLerpTime = lerpTime;
		_tickLerpTimeTotal = lerpTime;
	}

	private void OnSimTick()
	{
		StartLerp(Sim.TimeTillTick);
		_lastRenderState = _curRenderState;
		_nextRenderState = (_curRenderState + 1) % _renderStateCount;
		_curRenderState = (_nextRenderState + 1) % _renderStateCount;
		BuildRenderState(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState);
	}


	public void BuildRenderState(ref SimState from, ref DependentState dependent, ref DisplayState display, ref RenderState to, ref WorldData worldData, ref StaticState staticState)
	{
		to.Ticks = from.PlanetState.Ticks;
		to.Position = from.PlanetState.Position;
		to.Rotation = math.degrees(from.PlanetState.Rotation);

		MeshOverlayData meshOverlay;
		bool useMeshOverlay = GetMeshOverlayData(ActiveMeshOverlay, ref from, ref dependent, ref display, out meshOverlay);

		WindOverlayData windOverlayData;
		bool useWindOverlay = GetWindOverlayData(ActiveWindOverlay, ref from, ref dependent, out windOverlayData);

		var buildRenderStateJob = new BuildRenderStateJob()
		{ 
			TerrainColor = to.TerrainColor,
			TerrainNormal = to.TerrainNormal,
			TerrainPosition = to.TerrainPosition,
			WaterColor = to.WaterColor,
			WaterNormal = to.WaterNormal,
			WaterPosition = to.WaterPosition,
			CloudColor = to.CloudColor,
			CloudNormal = to.CloudNormal,
			CloudPosition = to.CloudPosition,
			VelocityArrow = to.VelocityArrow,
			SurfacePosition = to.SurfacePosition,

			TerrainScale = TerrainScale,
			PlanetRadius = staticState.PlanetRadius,
			CloudDropletSizeMin = worldData.rainDropMinSize,
			InverseCloudDropletSizeRange = 1.0f / (worldData.rainDropMaxSize - worldData.rainDropMinSize),
			MeshOverlayActive = useMeshOverlay,
			WindOverlayActive = useWindOverlay,
			WindVelocityMax = windOverlayData.MaxVelocity,
			WindMaskedByLand = windOverlayData.MaskLand,
			MeshOverlayMin = meshOverlay.Min,
			MeshOverlayInverseRange = meshOverlay.InverseRange,
			CloudElevation = from.CloudElevation,
			Icosphere = Sim.Icosphere.Vertices,
			Terrain = from.Terrain,
			CloudDropletMass = from.CloudDropletMass,
			CloudCoverage = dependent.CloudCoverage,
			IceCoverage = dependent.IceCoverage,
			VegetationCoverage = dependent.VegetationCoverage,
			WaterCoverage = dependent.WaterCoverage[Sim.WorldData.WaterLayers-1],
			WaterDepth = dependent.WaterDepth,
			SurfaceElevation = dependent.SurfaceElevation,
			MeshOverlayData = meshOverlay.Values,
			MeshOverlayColors = meshOverlay.ColorValuePairs,
			WindOverlayData = windOverlayData.Values
		};

		var buildRenderStateJobHandle = buildRenderStateJob.Schedule(Sim.CellCount, 100);
		buildRenderStateJobHandle.Complete();

	}

	public void UpdateMesh(ref RenderState lastState, ref RenderState nextState, ref RenderState state)
	{
		float t = Mathf.Clamp01(1.0f - _tickLerpTime / _tickLerpTimeTotal);
		if (!LerpStates)
		{
			t = 1;
		}
		state.Ticks = (nextState.Ticks - lastState.Ticks) * t + lastState.Ticks;
		state.Rotation = new Vector3(Mathf.LerpAngle(lastState.Rotation.x, nextState.Rotation.x, t), Mathf.LerpAngle(lastState.Rotation.y, nextState.Rotation.y, t), Mathf.LerpAngle(lastState.Rotation.z, nextState.Rotation.z, t));
		state.Position = math.lerp(lastState.Position, nextState.Position, t);


		NativeList<JobHandle> dependencies = new NativeList<JobHandle>(Allocator.Temp);
		dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.TerrainPosition, Start = lastState.TerrainPosition, End = nextState.TerrainPosition }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.TerrainNormal, Start = lastState.TerrainNormal, End = nextState.TerrainNormal }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.TerrainColor, Start = lastState.TerrainColor, End = nextState.TerrainColor }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.WaterPosition, Start = lastState.WaterPosition, End = nextState.WaterPosition }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.WaterNormal, Start = lastState.WaterNormal, End = nextState.WaterNormal }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.WaterColor, Start = lastState.WaterColor, End = nextState.WaterColor }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat3 { Progress = t, Out = state.SurfacePosition, Start = lastState.SurfacePosition, End = nextState.SurfacePosition }).Schedule(Sim.CellCount, _batchCount));

		if (true /* cloudsVisible*/)
		{
			dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.CloudPosition, Start = lastState.CloudPosition, End = nextState.CloudPosition }).Schedule(Sim.CellCount, _batchCount));
			dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.CloudNormal, Start = lastState.CloudNormal, End = nextState.CloudNormal }).Schedule(Sim.CellCount, _batchCount));
			dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.CloudColor, Start = lastState.CloudColor, End = nextState.CloudColor }).Schedule(Sim.CellCount, _batchCount));
		}
		if (ActiveWindOverlay != WindOverlay.None)
		{
			dependencies.Add((new LerpJobfloat2 { Progress = t, Out = state.VelocityArrow, Start = lastState.VelocityArrow, End = nextState.VelocityArrow }).Schedule(Sim.CellCount, _batchCount));
		}

		JobHandle.CompleteAll(dependencies);
		dependencies.Dispose();

		state.TerrainPosition.CopyTo(_terrainVertices);
		state.TerrainNormal.CopyTo(_terrainNormals);
		state.TerrainColor.CopyTo(_terrainColors);
		state.WaterPosition.CopyTo(_waterVertices);
		state.WaterNormal.CopyTo(_waterNormals);
		state.WaterColor.CopyTo(_waterColors);
		state.CloudPosition.CopyTo(_cloudVertices);
		state.CloudNormal.CopyTo(_cloudNormals);
		state.CloudColor.CopyTo(_cloudColors);

		_terrainMesh.vertices = _terrainVertices;
		_terrainMesh.normals = _terrainNormals;
		_terrainMesh.colors32 = _terrainColors;

		_waterMesh.vertices = _waterVertices;
		_waterMesh.normals = _waterNormals;
		_waterMesh.colors32 = _waterColors;

		_cloudMesh.vertices = _cloudVertices;
		_cloudMesh.normals = _cloudNormals;
		_cloudMesh.colors32 = _cloudColors;

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

		Planet.transform.SetPositionAndRotation(state.Position, Quaternion.Euler(state.Rotation));

		if (ActiveWindOverlay != WindOverlay.None)
		{
			for (int i = 0; i < _windArrows.Length; i++)
			{
				var wind = _renderStates[_curRenderState].VelocityArrow[i];
				bool visible = (wind.x == 0 && wind.y  == 0);
				_windArrows[i].SetActive(visible);
				if (visible)
				{
					var pos = _renderStates[_curRenderState].SurfacePosition[i];
					_windArrows[i].transform.localPosition = pos;
					_windArrows[i].transform.localRotation = Quaternion.LookRotation(-pos, Vector3.Cross(new Vector3(wind.x, wind.y, 0), pos));
					_windArrows[i].transform.GetChild(1).localScale = Vector3.one * math.length(wind);
				}
			}
		}

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
		ActiveMeshOverlay = (WorldView.MeshOverlay)dropdown.value;
		_tickLerpTime = 0;
		BuildRenderState(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState);
	}
	public void OnHUDWindChanged(UnityEngine.UI.Dropdown dropdown)
	{
		ActiveWindOverlay = (WorldView.WindOverlay)dropdown.value;
		_tickLerpTime = 0;
		BuildRenderState(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState);
	}

	public void OnHUDTemperatureUnitsChanged(UnityEngine.UI.Dropdown dropdown)
	{
		ActiveTemperatureUnits = (WorldView.TemperatureUnits)dropdown.value;
	}

	public static float ConvertTemperature(float kelvin, TemperatureUnits units)
	{
		switch (units)
		{
			case TemperatureUnits.Celsius:
				return kelvin - WorldData.FreezingTemperature;
			case TemperatureUnits.Farenheit:
				return (kelvin - WorldData.FreezingTemperature) * 9 / 5 + 32;
			case TemperatureUnits.Kelvin:
			default:
				return kelvin;
		}
	}

	public static string GetTemperatureString(float kelvin, TemperatureUnits units, int decimals)
	{
		string tFormat = "0";
		if (decimals > 0)
		{
			tFormat += ".";
		}
		for (int i=0;i<decimals;i++)
		{
			tFormat += "0";
		}
		string t = ConvertTemperature(kelvin, units).ToString(tFormat);
		switch (units)
		{
			case TemperatureUnits.Celsius:
				return t + " C";
			case TemperatureUnits.Farenheit:
				return t + " F";
			case TemperatureUnits.Kelvin:
			default:
				return t + " K";
		}
	}

	public string GetCellInfo(CellInfoType cellInfoType)
	{
		switch (cellInfoType)
		{
			case CellInfoType.Global:
				return GetCellInfoGlobal(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfoType.Energy:
				return GetCellInfoEnergy(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfoType.Cell:
				return GetCellInfoCell(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfoType.Atmosphere:
				return GetCellInfoAtmosphere(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfoType.Ground:
				return GetCellInfoGround(ref Sim.ActiveSimState, ref Sim.DependentState);
			case CellInfoType.Water:
				return GetCellInfoWater(ref Sim.ActiveSimState, ref Sim.DependentState);
		}
		return "";
	}

	public int GetClosestVert(int triangleIndex, int vIndex)
	{
		return indices[triangleIndex * 3 + vIndex];
	}

	public void SetActiveCell(int index, bool locked)
	{
		ActiveCellIndex = index;
		ActiveCellLocked = locked;
		_selectionCircle.SetActive(index >= 0);
		if (index >= 0)
		{
			//			var p = Sim.Icosphere.Vertices[index];
			var pos = _renderStates[_curRenderState].SurfacePosition[index];
			_selectionCircle.transform.localPosition = pos;
			_selectionCircle.transform.localRotation = Quaternion.LookRotation(-pos);
		}

		for (int i=0;i<_selectionCircleNeighbors.Count;i++)
		{
			_selectionCircleNeighbors[i].SetActive(false);
			if (index >= 0)
			{
				int n = Sim.StaticState.Neighbors[index * 6 + i];
				if (n >= 0)
				{
					_selectionCircleNeighbors[i].SetActive(true);
					var pos = _renderStates[_curRenderState].SurfacePosition[n];
					_selectionCircleNeighbors[i].transform.localPosition = pos;
					_selectionCircleNeighbors[i].transform.localRotation = Quaternion.LookRotation(-pos);
				}
			}
		}
	}

	public float ConvertTileEnergyToWatts(float energy)
	{
		return energy * 1000 / Sim.WorldData.SecondsPerTick;
	}


	public struct MeshOverlayData {
		public MeshOverlayData(float min, float max, NativeArray<CVP> colors, NativeArray<float> values)
		{
			Values = values;
			Min = min;
			Max = max;
			ColorValuePairs = colors;
			InverseRange = 1.0f / (Max - Min);
		}
		public float Min { get; private set; }
		public float Max { get; private set; }
		public NativeArray<float> Values { get; private set; }
		public NativeArray<CVP> ColorValuePairs { get; private set; }
		public float InverseRange { get; private set; }
	}

	public struct WindOverlayData {
		public WindOverlayData(float maxVelocity, bool maskLand, NativeArray<float2> values)
		{
			MaxVelocity = maxVelocity;
			MaskLand = maskLand;
			Values = values;
		}
		public float MaxVelocity { get; private set; }
		public bool MaskLand { get; private set; }
		public NativeArray<float2> Values { get; private set; }
	}



	#region private functions

	private bool GetMeshOverlayData(MeshOverlay activeOverlay, ref SimState simState, ref DependentState dependentState, ref DisplayState display, out MeshOverlayData overlay)
	{
		float ticksPerYear = Sim.WorldData.TicksPerSecond * 60 * 60 * 24 * 365;
		switch (activeOverlay)
		{
			case MeshOverlay.AbsoluteHumidityLower:
				overlay = new MeshOverlayData(0, DisplayAbsoluteHumidityMax, _normalizedRainbow, dependentState.AirHumidityAbsolute[0]);
				return true;
			case MeshOverlay.AbsoluteHumidityUpper:
				overlay = new MeshOverlayData(0, 1.0f, _normalizedRainbow, dependentState.AirHumidityAbsolute[1]);
				return true;
			case MeshOverlay.RelativeHumidityLower:
				overlay = new MeshOverlayData(0, 1.0f, _normalizedRainbow, dependentState.AirHumidityRelative[0]);
				return true;
			case MeshOverlay.LowerAirTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.AirTemperature[0]);
				return true;
			case MeshOverlay.UpperAirTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.AirTemperature[1]);
				return true;
			case MeshOverlay.Pressure:
				overlay = new MeshOverlayData(DisplayAirPressureMin, DisplayAirPressureMax, _normalizedRainbow, display.Pressure);
				return true;
			case MeshOverlay.ShallowWaterTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.WaterTemperature[Sim.WorldData.WaterLayers - 1]);
				return true;
			case MeshOverlay.ShallowWaterSalinity:
				overlay = new MeshOverlayData(DisplaySalinityMin, DisplaySalinityMax, _normalizedRainbow, dependentState.WaterSalinity[Sim.WorldData.WaterLayers - 1]);
				return true;
			case MeshOverlay.DeepWaterTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.WaterTemperature[0]);
				return true;
			case MeshOverlay.DeepWaterSalinity:
				overlay = new MeshOverlayData(DisplaySalinityMin, DisplaySalinityMax, _normalizedRainbow, dependentState.WaterSalinity[0]);
				return true;
			case MeshOverlay.GroundTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.TerrainTemperature);
				return true;
			case MeshOverlay.CloudTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.CloudTemperature);
				return true;
			case MeshOverlay.HeatAbsorbed:
				overlay = new MeshOverlayData(0, DisplayHeatAbsorbedMax, _normalizedRainbow, display.SolarRadiationAbsorbedSurface);
				return true;
			case MeshOverlay.Rainfall:
				overlay = new MeshOverlayData(0, DisplayRainfallMax, _normalizedRainbow, display.Rainfall);
				return true;
			case MeshOverlay.Evaporation:
				overlay = new MeshOverlayData(0, DisplayEvaporationMax, _normalizedRainbow, display.Evaporation);
				return true;
		}
		overlay = new MeshOverlayData(0, DisplayEvaporationMax, _normalizedRainbow, display.Evaporation);
		return false;

	}
	private bool GetWindOverlayData(WindOverlay activeOverlay, ref SimState simState, ref DependentState dependentState, out WindOverlayData overlay)
	{
		switch (activeOverlay)
		{
			case WindOverlay.LowerAirWind:
				overlay = new WindOverlayData(DisplayWindSpeedLowerAirMax, false, simState.Wind[0]);
				return true;
			case WindOverlay.UpperAirWind:
				overlay = new WindOverlayData(DisplayWindSpeedUpperAirMax, false, simState.Wind[Sim.WorldData.AirLayers-1]);
				return true;
			case WindOverlay.ShallowWaterCurrent:
				overlay = new WindOverlayData(DisplayWindSpeedSurfaceWaterMax, false, simState.WaterVelocity[Sim.WorldData.WaterLayers-1]);
				return true;
			case WindOverlay.DeepWaterCurrent:
				overlay = new WindOverlayData(DisplayWindSpeedDeepWaterMax, false, simState.WaterVelocity[0]);
				return true;
		}
		overlay = new WindOverlayData(DisplayWindSpeedDeepWaterMax, false, simState.WaterVelocity[0]);
		return false;
	}


	private string GetCellInfoGlobal(ref SimState state, ref DependentState dependent, ref DisplayState display)
	{
		string s = "";
		s += "CO2: " + state.PlanetState.CarbonDioxide;
		s += "\nCloud Coverage: " + (display.GlobalCloudCoverage * 100 * Sim.InverseCellCount).ToString("0.0") + "%";
		s += "\nGlobal Sea Level: " + (display.GlobalSeaLevel * Sim.InverseCellCount).ToString("0.00");
		s += "\nOcean Coverage: " + (display.GlobalOceanCoverage * 100 * Sim.InverseCellCount).ToString("0.0") + "%";
		s += "\nOcean Volume: " + (display.GlobalOceanVolume / 1000000000 * Sim.InverseCellCount).ToString("0.00") + " B";
		s += "\nTemperature: " + GetTemperatureString(display.GlobalTemperature * Sim.InverseCellCount, ActiveTemperatureUnits, 2);
		s += "\nCloud Mass: " + (display.GlobalCloudMass).ToString("0.00");
		s += "\nWater Vapor: " + (display.GlobalWaterVapor).ToString("0.00");
		s += "\nRainfall: " + (display.GlobalRainfall * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater).ToString("0.00");
		s += "\nEvaporation: " + (display.GlobalEvaporation * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater).ToString("0.00");

		return s;
	}
	private string GetCellInfoEnergy(ref SimState state, ref DependentState dependent, ref DisplayState display)
	{
		string s = "";

		var totalReflected = display.EnergySolarReflectedCloud + display.EnergySolarReflectedAtmosphere + display.EnergySolarReflectedSurface;
		var totalOutgoing = display.EnergyThermalSurfaceOutAtmosphericWindow + display.EnergyThermalOutAtmosphere;
		s += "Delta: " + ConvertTileEnergyToWatts((display.SolarRadiation - totalReflected - totalOutgoing) * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Incoming: " + ConvertTileEnergyToWatts(display.SolarRadiation * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Reflected: " + ConvertTileEnergyToWatts((totalReflected) * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Reflected Cloud: " + ConvertTileEnergyToWatts(display.EnergySolarReflectedCloud * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Reflected Atmos: " + ConvertTileEnergyToWatts(display.EnergySolarReflectedAtmosphere * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Reflected Surf: " + ConvertTileEnergyToWatts(display.EnergySolarReflectedSurface * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Abs Atm Total: " + ConvertTileEnergyToWatts(display.EnergySolarAbsorbedAtmosphere * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Abs Clouds: " + ConvertTileEnergyToWatts(display.EnergySolarAbsorbedCloud * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Abs Surface Total: " + ConvertTileEnergyToWatts(display.EnergySolarAbsorbedSurface * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Abs Ocean: " + ConvertTileEnergyToWatts(display.EnergySolarAbsorbedOcean * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Outgoing: " + ConvertTileEnergyToWatts(totalOutgoing * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Out Atm Window: " + ConvertTileEnergyToWatts(display.EnergyThermalSurfaceOutAtmosphericWindow * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Out Radiation: " + ConvertTileEnergyToWatts(display.EnergyThermalOutAtmosphere * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Surface Radiation: " + ConvertTileEnergyToWatts(display.EnergyThermalSurfaceRadiation * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Atm Absorbed: " + ConvertTileEnergyToWatts(display.EnergyThermalAbsorbedAtmosphere * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Back Radiation: " + ConvertTileEnergyToWatts(display.EnergyThermalBackRadiation * Sim.InverseCellCount).ToString("0.0");
		s += "\nEvapotranspiration: " + ConvertTileEnergyToWatts(display.EnergyEvapotranspiration * Sim.InverseCellCount).ToString("0.0");
		s += "\nSurface Conduction: " + ConvertTileEnergyToWatts(display.EnergySurfaceConduction * Sim.InverseCellCount).ToString("0.0");
		s += "\nOcean Radiation: " + ConvertTileEnergyToWatts(display.EnergyThermalOceanRadiation * Sim.InverseCellCount).ToString("0.0");
		s += "\nOcean Conduction: " + ConvertTileEnergyToWatts(display.EnergyOceanConduction * Sim.InverseCellCount).ToString("0.0");



		return s;
	}
	private string GetCellInfoCell(ref SimState state, ref DependentState dependent, ref DisplayState display)
	{
		if (ActiveCellIndex < 0)
			return "";

		var terrain = state.Terrain[ActiveCellIndex];
		var coord = Sim.StaticState.Coordinate[ActiveCellIndex];
		string s = "";
		s += "Index: " + ActiveCellIndex + "\n";
		s += "Coord: (" + math.degrees(coord.x).ToString("0.0") + ", " + math.degrees(coord.y).ToString("0.0") + ")\n";
		s += "Surface: " + (dependent.SurfaceElevation[ActiveCellIndex]).ToString("0.000") + " m\n";
		s += "Elevation: " + terrain.Elevation.ToString("0.000") + " m\n";
		s += "H2O Depth: " + dependent.WaterDepth[ActiveCellIndex].ToString("0.000") + " m\n";
		s += "Ice Depth: " + (state.IceMass[ActiveCellIndex] / WorldData.MassIce).ToString("0.000") + " m\n";
		s += "Ice Temperature: " + GetTemperatureString(state.IceTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1) + " m\n";
		s += "Rainfall: " + (display.Rainfall[ActiveCellIndex] / WorldData.MassWater * 1000).ToString("0.000") + " mm\n";
		s += "Evaporation: " + (display.Evaporation[ActiveCellIndex] / WorldData.MassWater * 1000).ToString("0.000") + " mm\n";
		//		s += "Condensation: " + (display.Condensation / WorldData.MassWater * 1000000).ToString("0.000") + " nm3\n";

		return s;
	}
	private string GetCellInfoAtmosphere(ref SimState state, ref DependentState dependent, ref DisplayState display)
	{
		if (ActiveCellIndex < 0)
			return "";

		float cloudMass = state.CloudMass[ActiveCellIndex];
		int upperAtmosphereLayerIndex = 1;
		string s = "";
		s += "Surface Temp: " + GetTemperatureString(state.AirTemperature[0][ActiveCellIndex], ActiveTemperatureUnits, 1) + "\n";
		s += "Surface Pressure: " + display.Pressure[ActiveCellIndex].ToString("0") + " Pa\n";
		s += "Surface Wind Horz: (" + state.Wind[0][ActiveCellIndex].x.ToString("0.0") + ", " + state.Wind[0][ActiveCellIndex].y.ToString("0.0") + ") m/s\n";
		s += "Surface Humidity: " + (dependent.AirHumidityRelative[0][ActiveCellIndex] * 100).ToString("0.0") + "%\n";
//		s += "Mid Temp: " + GetTemperatureString(state.AirTemperature[1][ActiveCellIndex], ActiveTemperatureUnits, 0) + "\n";
		s += "Upper Temp: " + GetTemperatureString(state.AirTemperature[upperAtmosphereLayerIndex][ActiveCellIndex], ActiveTemperatureUnits, 0) + "\n";
		s += "Upper Humidity: " + (dependent.AirHumidityRelative[upperAtmosphereLayerIndex][ActiveCellIndex] * 100).ToString("0.0") + "%\n";
		s += "Cloud Mass: " + (state.CloudMass[ActiveCellIndex]).ToString("0.000") + " kg\n";
		if (cloudMass > 0) {
			s += "Cloud Temp: " + GetTemperatureString(state.CloudTemperature[ActiveCellIndex], ActiveTemperatureUnits, 0) + "\n";
			s += "Cloud Elevation: " + (state.CloudElevation[ActiveCellIndex]).ToString("0") + " m\n";
			s += "Droplets: " + (state.CloudDropletMass[ActiveCellIndex] * 1000000 / (WorldData.MassWater* cloudMass)).ToString("0.000") + " nm3\n";
		}
		s += "Air Mass 0: " + dependent.AirMass[0][ActiveCellIndex].ToString("0") + " kg\n";
		s += "Air Mass 1: " + dependent.AirMass[1][ActiveCellIndex].ToString("0") + " kg\n";
//		s += "Air Mass 2: " + dependent.AirMass[2][ActiveCellIndex].ToString("0") + " kg\n";
		return s;
	}
	private string GetCellInfoGround(ref SimState state, ref DependentState dependent)
	{
		if (ActiveCellIndex < 0)
			return "";

		var terrain = state.Terrain[ActiveCellIndex];
		string s = "";
		s += "Fertility: " + terrain.SoilFertility + "\n";
		s += "Temp: " + GetTemperatureString(state.TerrainTemperature[ActiveCellIndex], ActiveTemperatureUnits, 0) + "\n";
		s += "Roughness: " + terrain.Roughness.ToString("0") + " m\n";
		return s;
	}
	private string GetCellInfoWater(ref SimState state, ref DependentState dependent)
	{
		if (ActiveCellIndex < 0)
			return "";

		string s = "";
		int waterSurfaceLayer = Sim.WorldData.WaterLayers - 1;
		s += "Surface Temp: " + GetTemperatureString(state.WaterTemperature[Sim.WorldData.WaterLayers - 1][ActiveCellIndex], ActiveTemperatureUnits, 0) + "\n";
		s += "Surface Salinity: " + (1000000 * dependent.WaterSalinity[waterSurfaceLayer][ActiveCellIndex]).ToString("0.0") + " ppm\n";
		return s;
	}

	#endregion
}
