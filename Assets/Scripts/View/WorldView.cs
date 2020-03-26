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
		TemperatureSurface,
		Evaporation,
		Rainfall,
		Condensation,
		HeatAbsorbed,
		GroundTemperature,
		GroundWater,
		GroundWaterTemperature,
		FloraTemperature,
		FloraWater,
		CrustDepth,
		MagmaMass,
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
		VerticalWind0,
		VerticalWind1,
		VerticalWind2,
		AbsoluteHumidity0,
		AbsoluteHumidity1,
		AbsoluteHumidity2,
		RelativeHumidity0,
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
	public float AtmosphereScale = 1000f;
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
	public float DisplayFloraMax = 1000;
	public float DisplayTemperatureMin = 223;
	public float DisplayTemperatureMax = 323;
	public float DisplayAbsoluteHumidityMax = 0.05f;
	public float DisplayAirPressureMin = 97000;
	public float DisplayAirPressureMax = 110000;
	public float DisplayHeatAbsorbedMax = 1000;
	public float DisplayCrustDepthMax = 10000;
	public float DisplayMagmaMassMax = 1000000;
	public float DisplayDustHeight = 1000;
	public float DisplayDustMax = 100;
	public float DisplayLavaTemperatureMax = 1200;

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
	public Material DustMaterial;
	public Material LavaMaterial;
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
	private Vector3[] _lavaVertices;
	private Vector3[] _lavaNormals;
	private Color32[] _lavaColors;
	private Vector3[] _dustVertices;
	private Vector3[] _dustNormals;
	private Color32[] _dustColors;

	private Mesh _terrainMesh;
	private Mesh _waterMesh;
	private Mesh _cloudMesh;
	private Mesh _dustMesh;
	private Mesh _lavaMesh;

	private GameObject _terrainObject;
	private GameObject _waterObject;
	private GameObject _cloudObject;
	private GameObject _dustObject;
	private GameObject _lavaObject;
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
		_lavaMesh = new Mesh();
		_dustMesh = new Mesh();

		_terrainObject = new GameObject("Terrain Mesh");
		_terrainObject.transform.SetParent(Planet.transform, false);
		var terrainFilter = _terrainObject.AddComponent<MeshFilter>();
		var terrainSurfaceRenderer = _terrainObject.AddComponent<MeshRenderer>();
		var terrainCollider = _terrainObject.AddComponent<MeshCollider>();
		terrainSurfaceRenderer.material = TerrainMaterial;
		terrainFilter.mesh = terrainCollider.sharedMesh = _terrainMesh;

		_lavaObject = new GameObject("Lava Mesh");
		_lavaObject.transform.SetParent(Planet.transform, false);
		var lavaFilter = _lavaObject.AddComponent<MeshFilter>();
		var lavaSurfaceRenderer = _lavaObject.AddComponent<MeshRenderer>();
		var lavaCollider = _lavaObject.AddComponent<MeshCollider>();
		lavaSurfaceRenderer.material = LavaMaterial;
		lavaFilter.mesh = lavaCollider.sharedMesh = _lavaMesh;


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

		_dustObject = new GameObject("Dust Mesh");
		_dustObject.transform.SetParent(Planet.transform, false);
		var dustFilter = _dustObject.AddComponent<MeshFilter>();
		var dustSurfaceRenderer = _dustObject.AddComponent<MeshRenderer>();
		dustSurfaceRenderer.material = DustMaterial;
		dustFilter.mesh = _dustMesh;

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
		_dustVertices = new Vector3[Sim.CellCount];
		_dustNormals = new Vector3[Sim.CellCount];
		_dustColors = new Color32[Sim.CellCount];
		_lavaVertices = new Vector3[Sim.CellCount];
		_lavaNormals = new Vector3[Sim.CellCount];
		_lavaColors = new Color32[Sim.CellCount];

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
			_tickLerpTime -= Time.deltaTime;
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
			LavaColor = to.LavaColor,
			LavaPosition = to.LavaPosition,
			LavaNormal = to.LavaNormal,
			WaterColor = to.WaterColor,
			WaterNormal = to.WaterNormal,
			WaterPosition = to.WaterPosition,
			CloudColor = to.CloudColor,
			CloudNormal = to.CloudNormal,
			CloudPosition = to.CloudPosition,
			DustColor = to.DustColor,
			DustNormal = to.DustNormal,
			DustPosition = to.DustPosition,
			VelocityArrow = to.VelocityArrow,
			SurfacePosition = to.SurfacePosition,

			TerrainScale = TerrainScale,
			AtmosphereScale = AtmosphereScale,
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
			SoilFertility = from.SoilFertility,
			Roughness = from.Roughness,
			Elevation = from.Elevation,
			CloudDropletMass = from.CloudDropletMass,
			CloudAbsorption = dependent.CloudAbsorptivity,
			IceCoverage = dependent.IceCoverage,
			FloraCoverage = dependent.FloraCoverage,
			WaterCoverage = dependent.WaterCoverage[Sim.WorldData.WaterLayers - 2],
			WaterDepth = dependent.WaterLayerDepth[1],
			WaterTemperature = from.WaterTemperature[Sim.WorldData.WaterLayers - 2],
			LavaMass = from.LavaMass,
			LavaTemperature = from.LavaTemperature,
			SurfaceElevation = dependent.LayerElevation[1],
			GroundWater = from.GroundWater,
			MeshOverlayData = meshOverlay.Values,
			MeshOverlayColors = meshOverlay.ColorValuePairs,
			WindOverlayData = windOverlayData.Values,
			GroundWaterMax = worldData.GroundWaterMax,
			DustHeight = DisplayDustHeight,
			LavaCrystalizationTemperature = worldData.LavaCrystalizationTemperature,
			LavaTemperatureRangeInverse = 1.0f / DisplayLavaTemperatureMax,
			DustCoverage = display.DustMass,
			DustMaxInverse = 1.0f / DisplayDustMax,
			LavaDensityAdjustment = worldData.LavaDensityAdjustment
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
		dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.DustPosition, Start = lastState.DustPosition, End = nextState.DustPosition }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.DustNormal, Start = lastState.DustNormal, End = nextState.DustNormal }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.DustColor, Start = lastState.DustColor, End = nextState.DustColor }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.LavaPosition, Start = lastState.LavaPosition, End = nextState.LavaPosition }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobVector3 { Progress = t, Out = state.LavaNormal, Start = lastState.LavaNormal, End = nextState.LavaNormal }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.LavaColor, Start = lastState.LavaColor, End = nextState.LavaColor }).Schedule(Sim.CellCount, _batchCount));
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
		state.LavaPosition.CopyTo(_lavaVertices);
		state.LavaNormal.CopyTo(_lavaNormals);
		state.LavaColor.CopyTo(_lavaColors);
		state.DustPosition.CopyTo(_dustVertices);
		state.DustNormal.CopyTo(_dustNormals);
		state.DustColor.CopyTo(_dustColors);

		_terrainMesh.vertices = _terrainVertices;
		_terrainMesh.normals = _terrainNormals;
		_terrainMesh.colors32 = _terrainColors;

		_waterMesh.vertices = _waterVertices;
		_waterMesh.normals = _waterNormals;
		_waterMesh.colors32 = _waterColors;

		_cloudMesh.vertices = _cloudVertices;
		_cloudMesh.normals = _cloudNormals;
		_cloudMesh.colors32 = _cloudColors;

		_lavaMesh.vertices = _lavaVertices;
		_lavaMesh.normals = _lavaNormals;
		_lavaMesh.colors32 = _lavaColors;

		_dustMesh.vertices = _dustVertices;
		_dustMesh.normals = _dustNormals;
		_dustMesh.colors32 = _dustColors;

		if (!_indicesInitialized)
		{
			_terrainMesh.SetTriangles(indices, 0);
			_waterMesh.SetTriangles(indices, 0);
			_cloudMesh.SetTriangles(indices, 0);
			_lavaMesh.SetTriangles(indices, 0);
			_dustMesh.SetTriangles(indices, 0);
			_indicesInitialized = true;
		}

		_terrainMesh.RecalculateBounds();
		_waterMesh.RecalculateBounds();
		_cloudMesh.RecalculateBounds();
		_lavaMesh.RecalculateBounds();
		_dustMesh.RecalculateBounds();

		_terrainMesh.RecalculateNormals();
		_waterMesh.RecalculateNormals();
		_cloudMesh.RecalculateNormals();
		_lavaMesh.RecalculateNormals();
		_dustMesh.RecalculateNormals();

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
		_dustObject.SetActive(toggle.isOn);
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
				return GetCellInfoWater(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.StaticState, ref Sim.DisplayState);
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
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, dependentState.SurfaceAirTemperatureAbsolute);
				return true;
			case MeshOverlay.PotentialTemperature0:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.AirTemperaturePotential[1]);
				return true;
			case MeshOverlay.PotentialTemperature1:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.AirTemperaturePotential[2]);
				return true;
			case MeshOverlay.PotentialTemperature2:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.AirTemperaturePotential[3]);
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
				overlay = new MeshOverlayData(DisplaySalinityMin, DisplaySalinityMax, _normalizedRainbow, display.Salinity[Sim.WorldData.WaterLayers - 2]);
				return true;
			case MeshOverlay.Salinity1:
				overlay = new MeshOverlayData(DisplaySalinityMin, DisplaySalinityMax, _normalizedRainbow, display.Salinity[Sim.WorldData.WaterLayers - 3]);
				return true;
			case MeshOverlay.Salinity2:
				overlay = new MeshOverlayData(DisplaySalinityMin, DisplaySalinityMax, _normalizedRainbow, display.Salinity[Sim.WorldData.WaterLayers - 4]);
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
			case MeshOverlay.GroundWater:
				overlay = new MeshOverlayData(0, Sim.WorldData.GroundWaterMax, _normalizedRainbow, simState.GroundWater);
				return true;
			case MeshOverlay.GroundWaterTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.GroundWaterTemperature);
				return true;
			case MeshOverlay.FloraTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.FloraTemperature);
				return true;
			case MeshOverlay.FloraWater:
				overlay = new MeshOverlayData(0, Sim.WorldData.FullCoverageFlora, _normalizedRainbow, simState.FloraWater);
				return true;
			case MeshOverlay.VerticalWind0:
				overlay = new MeshOverlayData(-1, 1, _normalizedBlueBlackRed, display.WindVertical[1]);
				return true;
			case MeshOverlay.VerticalWind1:
				overlay = new MeshOverlayData(-1, 1, _normalizedBlueBlackRed, display.WindVertical[2]);
				return true;
			case MeshOverlay.VerticalWind2:
				overlay = new MeshOverlayData(-1, 1, _normalizedBlueBlackRed, display.WindVertical[3]);
				return true;
			case MeshOverlay.CrustDepth:
				overlay = new MeshOverlayData(DisplayCrustDepthMax, 0, _normalizedRainbow, simState.CrustDepth);
				return true;
			case MeshOverlay.MagmaMass:
				overlay = new MeshOverlayData(0, DisplayMagmaMassMax, _normalizedRainbow, simState.MagmaMass);
				return true;
		}
		overlay = new MeshOverlayData(0, DisplayEvaporationMax, _normalizedRainbow, display.Evaporation);
		return false;

	}
	private bool GetWindOverlayData(WindOverlay activeOverlay, ref SimState simState, ref DependentState dependentState, ref DisplayState displayState, out WindOverlayData overlay)
	{
		switch (activeOverlay)
		{
			case WindOverlay.Wind0:
				overlay = new WindOverlayData(DisplayWindSpeedLowerAirMax, false, simState.AirVelocity[1]);
				return true;
			case WindOverlay.Wind1:
				overlay = new WindOverlayData(DisplayWindSpeedUpperAirMax, false, simState.AirVelocity[2]);
				return true;
			case WindOverlay.Wind2:
				overlay = new WindOverlayData(DisplayWindSpeedUpperAirMax, false, simState.AirVelocity[3]);
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
		s.AppendFormat("\nSurface Temp Air: {0}", GetTemperatureString(display.GlobalSurfaceTemperature * Sim.InverseCellCount, ActiveTemperatureUnits, 2));
		s.AppendFormat("\nSurface Temp Ocean: {0:N0}", GetTemperatureString(display.GlobalOceanSurfaceTemperature, ActiveTemperatureUnits, 2));
		s.AppendFormat("\nTemperature Air: {0}", GetTemperatureString((float)(display.GlobalAirTemperaturePotential), ActiveTemperatureUnits, 2));
		s.AppendFormat("\nTemperature Ocean: {0:N0}", GetTemperatureString(display.GlobalOceanTemperature, ActiveTemperatureUnits, 4));
		s.AppendFormat("\nTemperature Terrain: {0:N0}", GetTemperatureString((float)display.GlobalTerrainTemperature, ActiveTemperatureUnits, 4));
		s.AppendFormat("\nWater Vapor: {0:N0}", display.GlobalWaterVapor);
		s.AppendFormat("\nRainfall: {0:N3}", display.GlobalRainfall * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater);
		s.AppendFormat("\nCondensationCloud: {0:N3}", display.GlobalCondensationCloud * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater);
		s.AppendFormat("\nCondensationGround: {0:N3}", display.GlobalCondensationGround * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater);
		s.AppendFormat("\nEvaporation: {0:N3}", display.GlobalEvaporation * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater);
		s.AppendFormat("\nCloud Mass: {0:N2}", display.GlobalCloudMass);
		s.AppendFormat("\nGlobal Sea Level: {0:N2}", display.GlobalSeaLevel * Sim.InverseCellCount);
		s.AppendFormat("\nOcean Coverage: {0:N1}%", display.GlobalOceanCoverage * 100 * Sim.InverseCellCount);
		s.AppendFormat("\nOcean Mass: {0:N} M", display.GlobalOceanMass / 1000000);

		return s.ToString();
	}
	private string GetCellInfoEnergy(ref SimState state, ref DependentState dependent, ref DisplayState display)
	{
		StringBuilder s = new StringBuilder();
		NumberFormatInfo nfi1 = new NumberFormatInfo() { NumberDecimalDigits = 1 };
		NumberFormatInfo nfi2 = new NumberFormatInfo() { NumberDecimalDigits = 2 };

		var totalReflected = display.EnergySolarReflectedAtmosphere + display.EnergySolarReflectedSurface;
		var totalOutgoing = display.EnergyThermalSurfaceOutAtmosphericWindow + display.EnergyThermalOutAtmosphere;
		var emittedByAtmosphere = display.EnergyThermalOutAtmosphere - display.EnergyThermalSurfaceOutAtmosphericWindow;
		s.AppendFormat("Delta: {0:N1}", ConvertTileEnergyToWatts((display.SolarRadiation + display.GeothermalRadiation - totalReflected - totalOutgoing) * Sim.InverseCellCount));
		s.AppendFormat("\nEnthalpy Delta: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDelta * Sim.InverseCellCount)));
		s.AppendFormat("\nEnthalpy Delta Terrain: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaTerrain * Sim.InverseCellCount)));
		s.AppendFormat("\nEnthalpy Delta Water: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaWater * Sim.InverseCellCount)));
		s.AppendFormat("\nEnthalpy Delta Air: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaAir * Sim.InverseCellCount)));
		s.AppendFormat("\nEnthalpy Delta Ice: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaIce * Sim.InverseCellCount)));
		s.AppendFormat("\nEnthalpy Delta Cloud: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaCloud * Sim.InverseCellCount)));
		s.AppendFormat("\nEnthalpy Delta Terrain: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaTerrain * Sim.InverseCellCount)));
		s.AppendFormat("\nEnthalpy Delta Flora: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaFlora * Sim.InverseCellCount)));
		s.AppendFormat("\nEnthalpy Delta Ground Water: {0:N1}", ConvertTileEnergyToWatts((float)(display.GlobalEnthalpyDeltaGroundWater * Sim.InverseCellCount)));
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
		s.AppendFormat("\nT Atm Emitted: {0:N1}", ConvertTileEnergyToWatts(emittedByAtmosphere * Sim.InverseCellCount));
		s.AppendFormat("\nT Back Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalBackRadiation * Sim.InverseCellCount));
		s.AppendFormat("\nGeothermal Incoming: {0:N1}", ConvertTileEnergyToWatts(display.GeothermalRadiation * Sim.InverseCellCount));
		s.AppendFormat("\nEvapotranspiration: {0:N1}", ConvertTileEnergyToWatts(display.EnergyEvapotranspiration * Sim.InverseCellCount));
		s.AppendFormat("\nSurface Conduction: {0:N1}", ConvertTileEnergyToWatts(display.EnergySurfaceConduction * Sim.InverseCellCount));
		s.AppendFormat("\nOcean Radiation: {0:N1}", ConvertTileEnergyToWatts(display.EnergyThermalOceanRadiation * Sim.InverseCellCount / display.GlobalOceanCoverage));
		s.AppendFormat("\nOcean Conduction: {0:N1}", ConvertTileEnergyToWatts(display.EnergyOceanConduction * Sim.InverseCellCount / display.GlobalOceanCoverage));

		return s.ToString();
	}
	private string GetCellInfoCell(ref SimState state, ref DependentState dependent, ref DisplayState display)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		var coord = Sim.StaticState.Coordinate[ActiveCellIndex];
		var pos = Sim.StaticState.SphericalPosition[ActiveCellIndex];
		s.AppendFormat("INDEX: {0}\n", ActiveCellIndex);
		s.AppendFormat("COORD: ({0:N1}, {1:N1})\n", math.degrees(coord.x), math.degrees(coord.y));
		s.AppendFormat("POS: ({0:N2}, {1:N2}, {2:N2})\n", pos.x, pos.y, pos.z);
		//s.AppendFormat("\nSolar Terrain:   {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[0][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Terrain: {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[0][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Water0:    {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[2][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Water0:  {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[2][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Water1:    {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[3][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Water1:  {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[3][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Water2:    {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[4][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Water2:  {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[4][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Ice:       {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[6][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Ice:     {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[6][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Air0:       {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[8][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Air0:     {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[8][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Air1:       {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[1][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Air1:     {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[1][ActiveCellIndex]));
		//s.AppendFormat("\nSolar Air2:       {0:N1}", ConvertTileEnergyToWatts(display.SolarDelta[2][ActiveCellIndex]));
		//s.AppendFormat("\nThermal Air2:     {0:N1}", ConvertTileEnergyToWatts(display.ThermalDelta[2][ActiveCellIndex]));

		return s.ToString();
	}
	private string GetCellInfoAtmosphere(ref SimState state, ref DependentState dependent, ref StaticState staticState, ref DisplayState display)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		float cloudMass = state.CloudMass[ActiveCellIndex];
		int upperAtmosphereLayerIndex = Sim.WorldData.AirLayers - 2;
		s.AppendFormat("RAIN: {0:N3} kg", display.Rainfall[ActiveCellIndex]);
		s.AppendFormat("\nEVAP: {0:N3} kg", display.Evaporation[ActiveCellIndex]);
		s.AppendFormat("\nSURFACE TEMP: {0:N3}", GetTemperatureString(dependent.SurfaceAirTemperatureAbsolute[ActiveCellIndex], ActiveTemperatureUnits, 1));
		s.AppendFormat("\nTROPOPAUSE ELE: {0:N0}m", dependent.LayerElevation[Sim.WorldData.AirLayers-1][ActiveCellIndex]);
		s.AppendFormat("\nDUST: {0:N3} kg", display.DustMass[ActiveCellIndex]);

		if (cloudMass > 0)
		{
			float dropletSize = 1000 * Atmosphere.GetDropletRadius(state.CloudDropletMass[ActiveCellIndex], Atmosphere.GetWaterDensityAtElevation(dependent.DewPoint[ActiveCellIndex], dependent.CloudElevation[ActiveCellIndex]));
			float3 vel = Utils.GetPolarCoordinates(staticState.SphericalPosition[ActiveCellIndex],state.CloudVelocity[ActiveCellIndex]);
			s.AppendFormat("\nCLOUD: {0:N3} kg ELE: {1:N0} m R: {2:N3} mm VEL: ({3:N1}, {4:N1}, {5:N1})", 
				state.CloudMass[ActiveCellIndex],
				dependent.CloudElevation[ActiveCellIndex],
				state.CloudDropletMass[ActiveCellIndex],
				vel.x, vel.y, vel.z);
		}
		else
		{
			s.AppendFormat("\nCLOUD: 0 kg");
		}
		s.AppendLine();

		for (int i = 1; i < Sim.WorldData.AirLayers - 1; i++)
		{
			var wind = Utils.GetPolarCoordinates(staticState.SphericalPosition[ActiveCellIndex], state.AirVelocity[i][ActiveCellIndex]);
			s.AppendFormat("\nLAYER {0} | TEMP: {1} RH: {2:P1}",
				i, 
				GetTemperatureString(Atmosphere.GetAbsoluteTemperature(state.AirTemperaturePotential[i][ActiveCellIndex], dependent.LayerMiddle[i][ActiveCellIndex]), ActiveTemperatureUnits, 1),
				(dependent.AirHumidityRelative[i][ActiveCellIndex]));
			s.AppendFormat("\nELE: {0:N0} m P: {1:N0} Pa WIND: ({2:N1}, {3:N1}, {4:N1})",
				dependent.LayerElevation[i][ActiveCellIndex],
				display.Pressure[i][ActiveCellIndex],
				wind.x, wind.y, wind.z);
			s.AppendFormat("\nMASS: {0:N0} kg VAPOR: {1:N0} kg", dependent.AirMass[i][ActiveCellIndex], state.AirVapor[i][ActiveCellIndex]);
			s.AppendFormat("\nSOLAR ABSORB: {0:P0} REFLECT: {1:P0}", display.AbsorptionSolar[i][ActiveCellIndex].AbsorptivityAirAbove + (1.0f - display.AbsorptionSolar[i][ActiveCellIndex].AbsorptivityAirAbove) * display.AbsorptionSolar[i][ActiveCellIndex].AbsorptivityAirBelow, display.AbsorptionSolar[i][ActiveCellIndex].ReflectivityAirAbove + (1.0f - display.AbsorptionSolar[i][ActiveCellIndex].ReflectivityAirAbove) * display.AbsorptionSolar[i][ActiveCellIndex].ReflectivityAirBelow);
			s.AppendFormat("\nCLOUD ABSORB: {0:P0} REFLECT: {1:P0}", display.AbsorptionSolar[i][ActiveCellIndex].AbsorptivityCloud, display.AbsorptionSolar[i][ActiveCellIndex].ReflectivityCloud);
			s.AppendFormat("\nTHERMAL ABSORB: {0:P0} CLOUD: {1:P0}", display.AbsorptionThermal[i][ActiveCellIndex].AbsorptivityAirAbove + (1.0f - display.AbsorptionThermal[i][ActiveCellIndex].AbsorptivityAirAbove) * display.AbsorptionThermal[i][ActiveCellIndex].AbsorptivityAirBelow, display.AbsorptionThermal[i][ActiveCellIndex].AbsorptivityCloud);
			s.AppendLine();
		}
		return s.ToString();
	}
	private string GetCellInfoGround(ref SimState state, ref DependentState dependent)
	{
		if (ActiveCellIndex < 0)
			return "";

		StringBuilder s = new StringBuilder();

		s.AppendFormat("ELE: {0:N0} m\n", state.Elevation[ActiveCellIndex]);
		s.AppendFormat("ROUGH: {0:N0} m\n", state.Roughness[ActiveCellIndex]);
		s.AppendFormat("TEMP: {0}\n", GetTemperatureString(state.TerrainTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		s.AppendFormat("FERT: {0:N2}\n", state.SoilFertility[ActiveCellIndex]);
		s.AppendFormat("FLORA: {0:N2} WATER: {1:N2} kg TEMP: {2:N2}\n", 
			state.FloraMass[ActiveCellIndex],
			state.FloraWater[ActiveCellIndex],
			GetTemperatureString(state.FloraTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		s.AppendFormat("GWATER: {0:N2} TEMP: {1:N2}\n",
			state.GroundWater[ActiveCellIndex],
			GetTemperatureString(state.GroundWaterTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		s.AppendFormat("LAVA: {0:N2} TEMP: {1:N0}\n",
			state.LavaMass[ActiveCellIndex],
			GetTemperatureString(state.LavaTemperature[ActiveCellIndex], ActiveTemperatureUnits, 1));
		return s.ToString();
	}
	private string GetCellInfoWater(ref SimState state, ref DependentState dependent, ref StaticState staticState, ref DisplayState display)
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
		float depth = dependent.WaterLayerDepth[1][ActiveCellIndex];
		var nfi = new NumberFormatInfo() { NumberDecimalDigits = (depth >= 1) ? 0 : 3 };
		s.AppendFormat(nfi, "DEPTH: {0:N} m\n", depth);
		s.AppendLine();

		for (int i = Sim.WorldData.WaterLayers - 2; i >= 1; i--) {
			int layerIndex = (Sim.WorldData.WaterLayers - 2 - i);
			if (state.WaterMass[i][ActiveCellIndex] > 0)
			{
				var current = Utils.GetPolarCoordinates(staticState.SphericalPosition[ActiveCellIndex], state.WaterVelocity[i][ActiveCellIndex]);
				s.AppendFormat("LAYER {0} | TEMP: {1} SALT: {2:P4}\n",
					layerIndex,
					GetTemperatureString(state.WaterTemperature[i][ActiveCellIndex], ActiveTemperatureUnits, 2),
					display.Salinity[i][ActiveCellIndex]);
				s.AppendFormat("VEL: ({0:N3}, {1:N3}, {2:N3}) P: {3} D: {4}\n",
					current.x, current.y, current.z,
					dependent.WaterPressure[i][ActiveCellIndex],
					Atmosphere.GetWaterDensity(display.Salinity[i][ActiveCellIndex], state.WaterTemperature[i][ActiveCellIndex]));
				s.AppendFormat(nfi, "DEPTH: {0:N} m HEIGHT: {1:N} m\n",
					dependent.WaterLayerDepth[i][ActiveCellIndex],
					dependent.WaterLayerHeight[i][ActiveCellIndex]
					);
				s.AppendLine();
			}
		}
		return s.ToString();
	}

	#endregion
}
