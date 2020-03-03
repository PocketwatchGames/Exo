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
using System.Globalization;

public class WorldView : MonoBehaviour {

	public enum MeshOverlay {
		None,
		GroundTemperature,
		TemperatureSurface,
		PotentialTemperature0,
		PotentialTemperature1,
		PotentialTemperature2,
		Pressure0,
		Pressure1,
		Pressure2,
		WaterTemperature0,
		WaterTemperature1,
		WaterTemperature2,
		Salinity0,
		Salinity1,
		Salinity2,
		VerticalWind,
		AbsoluteHumidity0,
		AbsoluteHumidity1,
		AbsoluteHumidity2,
		RelativeHumidity0,
		Evaporation,
		Rainfall,
		Condensation,
		HeatAbsorbed,
	}

	public enum WindOverlay {
		None,
		Current0,
		Current1,
		Current2,
		Wind0,
		Wind1,
		Wind2,
		PGF0,
		PGF1,
		PGF2,
		WindCloud,

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

	public int ActiveCell;
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
	public float DisplayPressureGradientForceMax = 0.01f;
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


		if (ActiveCell != -1)
		{
			SetActiveCell(ActiveCell, true);
			ActiveCell = -1;
		}
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
		bool useWindOverlay = GetWindOverlayData(ActiveWindOverlay, ref from, ref dependent, ref display, out windOverlayData);

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
			CloudElevation = dependent.CloudElevation,
			Icosphere = Sim.Icosphere.Vertices,
			Terrain = from.Terrain,
			CloudDropletMass = from.CloudDropletMass,
			CloudCoverage = dependent.CloudCoverage,
			IceCoverage = dependent.IceCoverage,
			VegetationCoverage = dependent.VegetationCoverage,
			WaterCoverage = dependent.WaterCoverage[Sim.WorldData.WaterLayers-2],
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
			dependencies.Add((new LerpJobfloat3 { Progress = t, Out = state.VelocityArrow, Start = lastState.VelocityArrow, End = nextState.VelocityArrow }).Schedule(Sim.CellCount, _batchCount));
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
				var windHorizontal = math.cross(math.cross(Sim.StaticState.SphericalPosition[i], wind), Sim.StaticState.SphericalPosition[i]);
				float windSpeed = math.length(windHorizontal);
				var pos = _renderStates[_curRenderState].SurfacePosition[i];
				bool visible = windSpeed > 0;
				_windArrows[i].SetActive(visible);
				if (visible)
				{
					_windArrows[i].transform.localPosition = pos;
					_windArrows[i].transform.localRotation = Quaternion.LookRotation(windHorizontal / windSpeed, pos);
					_windArrows[i].transform.GetChild(1).localScale = Vector3.one * math.min(1, windSpeed);
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

		for (int i = 0; i < _windArrows.Length; i++)
		{
			var wind = _renderStates[_curRenderState].VelocityArrow[i];
			bool visible = ActiveWindOverlay != WindOverlay.None && !(wind.x == 0 && wind.y == 0);
			_windArrows[i].SetActive(visible);

		}
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
				return GetCellInfoAtmosphere(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.StaticState, ref Sim.DisplayState);
			case CellInfoType.Ground:
				return GetCellInfoGround(ref Sim.ActiveSimState, ref Sim.DependentState);
			case CellInfoType.Water:
				return GetCellInfoWater(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.StaticState);
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
		public WindOverlayData(float maxVelocity, bool maskLand, NativeArray<float3> values)
		{
			MaxVelocity = maxVelocity;
			MaskLand = maskLand;
			Values = values;
		}
		public float MaxVelocity { get; private set; }
		public bool MaskLand { get; private set; }
		public NativeArray<float3> Values { get; private set; }
	}



	#region private functions

	private bool GetMeshOverlayData(MeshOverlay activeOverlay, ref SimState simState, ref DependentState dependentState, ref DisplayState display, out MeshOverlayData overlay)
	{
		float ticksPerYear = Sim.WorldData.TicksPerSecond * 60 * 60 * 24 * 365;
		switch (activeOverlay)
		{
			case MeshOverlay.AbsoluteHumidity0:
				overlay = new MeshOverlayData(0, DisplayAbsoluteHumidityMax, _normalizedRainbow, dependentState.AirHumidityAbsolute[1]);
				return true;
			case MeshOverlay.AbsoluteHumidity1:
				overlay = new MeshOverlayData(0, 1.0f, _normalizedRainbow, dependentState.AirHumidityAbsolute[2]);
				return true;
			case MeshOverlay.AbsoluteHumidity2:
				overlay = new MeshOverlayData(0, 1.0f, _normalizedRainbow, dependentState.AirHumidityAbsolute[3]);
				return true;
			case MeshOverlay.RelativeHumidity0:
				overlay = new MeshOverlayData(0, 1.0f, _normalizedRainbow, dependentState.AirHumidityRelative[1]);
				return true;
			case MeshOverlay.TemperatureSurface:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, dependentState.SurfaceAirTemperature);
				return true;
			case MeshOverlay.PotentialTemperature0:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, display.PotentialTemperature[1]);
				return true;
			case MeshOverlay.PotentialTemperature1:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, display.PotentialTemperature[2]);
				return true;
			case MeshOverlay.PotentialTemperature2:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, display.PotentialTemperature[3]);
				return true;
			case MeshOverlay.Pressure0:
				overlay = new MeshOverlayData(DisplayAirPressureMin, DisplayAirPressureMax, _normalizedRainbow, display.Pressure[1]);
				return true;
			case MeshOverlay.Pressure1:
				overlay = new MeshOverlayData(DisplayAirPressureMin, DisplayAirPressureMax, _normalizedRainbow, display.Pressure[2]);
				return true;
			case MeshOverlay.Pressure2:
				overlay = new MeshOverlayData(DisplayAirPressureMin, DisplayAirPressureMax, _normalizedRainbow, display.Pressure[3]);
				return true;
			case MeshOverlay.WaterTemperature0:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.WaterTemperature[Sim.WorldData.WaterLayers - 2]);
				return true;
			case MeshOverlay.WaterTemperature1:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.WaterTemperature[Sim.WorldData.WaterLayers - 3]);
				return true;
			case MeshOverlay.WaterTemperature2:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.WaterTemperature[Sim.WorldData.WaterLayers - 4]);
				return true;
			case MeshOverlay.Salinity0:
				overlay = new MeshOverlayData(DisplaySalinityMin, DisplaySalinityMax, _normalizedRainbow, dependentState.WaterSalinity[Sim.WorldData.WaterLayers - 2]);
				return true;
			case MeshOverlay.Salinity1:
				overlay = new MeshOverlayData(DisplaySalinityMin, DisplaySalinityMax, _normalizedRainbow, dependentState.WaterSalinity[Sim.WorldData.WaterLayers - 3]);
				return true;
			case MeshOverlay.Salinity2:
				overlay = new MeshOverlayData(DisplaySalinityMin, DisplaySalinityMax, _normalizedRainbow, dependentState.WaterSalinity[Sim.WorldData.WaterLayers - 4]);
				return true;
			case MeshOverlay.GroundTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.TerrainTemperature);
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
			//case MeshOverlay.VerticalWind:
			//	overlay = new MeshOverlayData(-DisplayVerticalWindSpeedMax, DisplayVerticalWindSpeedMax, _normalizedBlueBlackRed, simState.WindVertical[1]);
			//	return true;
		}
		overlay = new MeshOverlayData(0, DisplayEvaporationMax, _normalizedRainbow, display.Evaporation);
		return false;

	}
	private bool GetWindOverlayData(WindOverlay activeOverlay, ref SimState simState, ref DependentState dependentState, ref DisplayState displayState, out WindOverlayData overlay)
	{
		switch (activeOverlay)
		{
			case WindOverlay.Wind0:
				overlay = new WindOverlayData(DisplayWindSpeedLowerAirMax, false, simState.Wind[1]);
				return true;
			case WindOverlay.Wind1:
				overlay = new WindOverlayData(DisplayWindSpeedUpperAirMax, false, simState.Wind[2]);
				return true;
			case WindOverlay.Wind2:
				overlay = new WindOverlayData(DisplayWindSpeedUpperAirMax, false, simState.Wind[3]);
				return true;
			case WindOverlay.WindCloud:
				overlay = new WindOverlayData(DisplayWindSpeedUpperAirMax, false, simState.CloudVelocity);
				return true;
			case WindOverlay.PGF0:
				overlay = new WindOverlayData(DisplayPressureGradientForceMax, false, displayState.PressureGradientForce[1]);
				return true;
			case WindOverlay.PGF1:
				overlay = new WindOverlayData(DisplayPressureGradientForceMax, false, displayState.PressureGradientForce[2]);
				return true;
			case WindOverlay.PGF2:
				overlay = new WindOverlayData(DisplayPressureGradientForceMax, false, displayState.PressureGradientForce[3]);
				return true;
			case WindOverlay.Current0:
				overlay = new WindOverlayData(DisplayWindSpeedSurfaceWaterMax, true, simState.WaterVelocity[Sim.WorldData.WaterLayers-2]);
				return true;
			case WindOverlay.Current1:
				overlay = new WindOverlayData(DisplayWindSpeedDeepWaterMax, true, simState.WaterVelocity[Sim.WorldData.WaterLayers - 3]);
				return true;
			case WindOverlay.Current2:
				overlay = new WindOverlayData(DisplayWindSpeedDeepWaterMax, true, simState.WaterVelocity[Sim.WorldData.WaterLayers - 4]);
				return true;
		}
		overlay = new WindOverlayData(DisplayWindSpeedDeepWaterMax, false, simState.WaterVelocity[1]);
		return false;
	}


	private string GetCellInfoGlobal(ref SimState state, ref DependentState dependent, ref DisplayState display)
	{
		StringBuilder s = new StringBuilder();
		NumberFormatInfo nfi1 = new NumberFormatInfo() { NumberDecimalDigits = 1 };
		NumberFormatInfo nfi2 = new NumberFormatInfo() { NumberDecimalDigits = 2 };
		s.AppendFormat("CO2: {0}", state.PlanetState.CarbonDioxide);
		s.AppendFormat("\nCloud Coverage: {0:N1}%", display.GlobalCloudCoverage * 100 * Sim.InverseCellCount);
		s.AppendFormat("\nGlobal Sea Level: {0:N2}", display.GlobalSeaLevel * Sim.InverseCellCount);
		s.AppendFormat("\nOcean Coverage: {0:N1}%", display.GlobalOceanCoverage * 100 * Sim.InverseCellCount);
		s.AppendFormat("\nOcean Volume: {0:N2} B", display.GlobalOceanVolume / 1000000000 * Sim.InverseCellCount);
		s.AppendFormat("\nTemperature: {0}", GetTemperatureString(display.GlobalTemperature * Sim.InverseCellCount, ActiveTemperatureUnits, 2));
		s.AppendFormat("\nCloud Mass: {0:N2}", display.GlobalCloudMass);
		s.AppendFormat("\nWater Vapor: {0:N0}", display.GlobalWaterVapor);
		s.AppendFormat("\nRainfall: {0:N2}", display.GlobalRainfall * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater);
		s.AppendFormat("\nEvaporation: {0:N2}", display.GlobalEvaporation * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater);

		return s.ToString();
	}
	private string GetCellInfoEnergy(ref SimState state, ref DependentState dependent, ref DisplayState display)
	{
		StringBuilder s = new StringBuilder();
		NumberFormatInfo nfi1 = new NumberFormatInfo() { NumberDecimalDigits = 1 };
		NumberFormatInfo nfi2 = new NumberFormatInfo() { NumberDecimalDigits = 2 };

		var totalReflected = display.EnergySolarReflectedAtmosphere + display.EnergySolarReflectedSurface;
		var totalOutgoing = display.EnergyThermalSurfaceOutAtmosphericWindow + display.EnergyThermalOutAtmosphere;
		s.AppendFormat("Delta: {0:N1}", ConvertTileEnergyToWatts((display.SolarRadiation - totalReflected - totalOutgoing) * Sim.InverseCellCount));
		s.AppendFormat("\nS Incoming: {0:N1}", ConvertTileEnergyToWatts(display.SolarRadiation * Sim.InverseCellCount));
		s.AppendFormat("\nS Reflected: {0:N1}", ConvertTileEnergyToWatts((totalReflected) * Sim.InverseCellCount));
		s.AppendFormat("\nS Reflected Atmos: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarReflectedAtmosphere * Sim.InverseCellCount));
		s.AppendFormat("\nS Reflected Surf: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarReflectedSurface * Sim.InverseCellCount));
		s.AppendFormat("\nS Abs Atm: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarAbsorbedAtmosphere * Sim.InverseCellCount));
		s.AppendFormat("\nS Abs Surface Total: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarAbsorbedSurface * Sim.InverseCellCount));
		s.AppendFormat("\nS Abs Ocean: {0:N1}", ConvertTileEnergyToWatts(display.EnergySolarAbsorbedOcean * Sim.InverseCellCount));
		s.AppendFormat("\nT Outgoing: {0:N1}", ConvertTileEnergyToWatts(totalOutgoing * Sim.InverseCellCount));
		s.AppendFormat("\nT Out Atm Window: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalSurfaceOutAtmosphericWindow * Sim.InverseCellCount));
		s.AppendFormat("\nT Out Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalOutAtmosphere * Sim.InverseCellCount));
		s.AppendFormat("\nT Surface Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalSurfaceRadiation * Sim.InverseCellCount));
		s.AppendFormat("\nT Atm Absorbed: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalAbsorbedAtmosphere * Sim.InverseCellCount));
		s.AppendFormat("\nT Back Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalBackRadiation * Sim.InverseCellCount));
		s.AppendFormat("\nEvapotranspiration: {0:N1}", ConvertTileEnergyToWatts(display.EnergyEvapotranspiration * Sim.InverseCellCount));
		s.AppendFormat("\nSurface Conduction: {0:N1}", ConvertTileEnergyToWatts(display.EnergySurfaceConduction * Sim.InverseCellCount));
		s.AppendFormat("\nOcean Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalOceanRadiation * Sim.InverseCellCount));
		s.AppendFormat("\nOcean Conduction: {0:N1}", ConvertTileEnergyToWatts(display.EnergyOceanConduction * Sim.InverseCellCount));

		return s.ToString();
	}
	private string GetCellInfoCell(ref SimState state, ref DependentState dependent, ref DisplayState display)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		var terrain = state.Terrain[ActiveCellIndex];
		var coord = Sim.StaticState.Coordinate[ActiveCellIndex];
		var pos = Sim.StaticState.SphericalPosition[ActiveCellIndex];
		s.AppendFormat("INDEX: {0}\n", ActiveCellIndex);
		s.AppendFormat("COORD: ({0:N1}, {1:N1})\n", math.degrees(coord.x), math.degrees(coord.y));
		s.AppendFormat("POS: ({0:N2}, {1:N2}, {2:N2})\n", pos.x, pos.y, pos.z);

		return s.ToString();
	}
	private string GetCellInfoAtmosphere(ref SimState state, ref DependentState dependent, ref StaticState staticState, ref DisplayState display)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		float cloudMass = state.CloudMass[ActiveCellIndex];
		int upperAtmosphereLayerIndex = Sim.WorldData.AirLayers - 2;
		s.AppendFormat("RAIN: {0:N3} kg\n", display.Rainfall[ActiveCellIndex]);
		s.AppendFormat("EVAP: {0:N3} kg\n", display.Evaporation[ActiveCellIndex]);

		if (cloudMass > 0)
		{
			float dropletSize = 1000 * Atmosphere.GetDropletRadius(state.CloudDropletMass[ActiveCellIndex], Atmosphere.GetWaterDensityAtElevation(dependent.DewPoint[ActiveCellIndex], dependent.CloudElevation[ActiveCellIndex]));
			float3 vel = Utils.GetPolarCoordinates(staticState.SphericalPosition[ActiveCellIndex],state.CloudVelocity[ActiveCellIndex]);
			s.AppendFormat("CLOUD: {0:N3} kg ELE: {1:N0} m R: {2:N3} mm VEL: ({3:N1}, {4:N1}, {5:N1})\n", 
				state.CloudMass[ActiveCellIndex],
				dependent.CloudElevation[ActiveCellIndex],
				state.CloudDropletMass[ActiveCellIndex],
				vel.x, vel.y, vel.z);
		}
		else
		{
			s.AppendFormat("CLOUD: 0 kg\n");
		}
		s.AppendLine();

		for (int i = 1; i < Sim.WorldData.AirLayers - 1; i++)
		{
			var wind = Utils.GetPolarCoordinates(staticState.SphericalPosition[ActiveCellIndex], state.Wind[i][ActiveCellIndex]);
			s.AppendFormat("LAYER {0} | TEMP: {1} RH: {2:P1}\n",
				i, 
				GetTemperatureString(state.AirTemperature[i][ActiveCellIndex], ActiveTemperatureUnits, 1),
				(dependent.AirHumidityRelative[i][ActiveCellIndex]));
			s.AppendFormat("ELE: {0:N0} m P: {1:N0} Pa WIND: ({2:N1}, {3:N1}, {4:N1})\n",
				dependent.LayerElevation[i][ActiveCellIndex],
				display.Pressure[i][ActiveCellIndex],
				wind.x, wind.y, wind.z);
			s.AppendFormat("MASS: {0:N0} kg VAPOR: {1:N0} kg\n", dependent.AirMass[i][ActiveCellIndex], state.AirVapor[i][ActiveCellIndex]);
			s.AppendLine();
		}
		return s.ToString();
	}
	private string GetCellInfoGround(ref SimState state, ref DependentState dependent)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		var terrain = state.Terrain[ActiveCellIndex];
		s.AppendFormat("ELE: {0:N0} m\n", terrain.Elevation);
		s.AppendFormat("ROUGH: {0:N0} m\n", terrain.Roughness);
		s.AppendFormat("TEMP: {0}\n", GetTemperatureString(state.TerrainTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		s.AppendFormat("FERT: {0:N2}\n", terrain.SoilFertility);
		s.AppendFormat("VEG: {0:N2}\n", terrain.Vegetation);
		return s.ToString();
	}
	private string GetCellInfoWater(ref SimState state, ref DependentState dependent, ref StaticState staticState)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		if (state.IceMass[ActiveCellIndex] > 0)
		{
			s.AppendFormat("ICE: {0:N3} m TEMP: {1}\n", 
				(state.IceMass[ActiveCellIndex] / WorldData.MassIce), 
				GetTemperatureString(state.IceTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		} else
		{
			s.AppendFormat("ICE: 0 m\n");
		}
		s.AppendFormat("DEPTH: {0:N3} m\n", dependent.WaterDepth[ActiveCellIndex]);
		s.AppendLine();

		for (int i = Sim.WorldData.WaterLayers - 2; i >= 1; i--) {
			int layerIndex = (Sim.WorldData.WaterLayers - 2 - i);
			if (state.WaterMass[i][ActiveCellIndex] > 0)
			{
				var current = Utils.GetPolarCoordinates(staticState.SphericalPosition[ActiveCellIndex], state.WaterVelocity[i][ActiveCellIndex]);
				s.AppendFormat("LAYER {0} | TEMP: {1} SALT: {2:P3}\n",
					layerIndex,
					GetTemperatureString(state.WaterTemperature[i][ActiveCellIndex], ActiveTemperatureUnits, 0),
					dependent.WaterSalinity[i][ActiveCellIndex]);
				s.AppendFormat("MASS: {0:N1} kg VEL: ({1:N3}, {2:N3}, {3:N3}) P: {4}\n",
					state.SaltMass[i][ActiveCellIndex],
					current.x, current.y, current.z,
					dependent.WaterPressure[i][ActiveCellIndex]);
				s.AppendLine();
			}
		}
		return s.ToString();
	}

	#endregion
}
