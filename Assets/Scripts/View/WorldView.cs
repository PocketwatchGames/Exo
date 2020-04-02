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
		WaterCarbonDioxide0,
		WaterCarbonDioxide1,
		WaterCarbonDioxide2,
		Salinity0,
		Salinity1,
		Salinity2,
		VerticalWind0,
		VerticalWind1,
		VerticalWind2,
		CarbonDioxide0,
		CarbonDioxide1,
		CarbonDioxide2,
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
	public const int VertsPerCell = 25;
	public const int MaxNeighbors = 6;

	public int ActiveCell;
	public bool ActiveCellLocked { get; private set; }
	public int ActiveCellIndex { get; private set; }
	public CellInfo.TemperatureUnits ActiveTemperatureUnits = CellInfo.TemperatureUnits.Celsius;
	public bool LerpStates = true;


	[Header("Display")]
	public float SlopeAmount = 1.5f;
	public float TerrainScale = 100f;
	public float AtmosphereScale = 1000f;
//	public TemperatureDisplayType TemperatureDisplay;
	public float MinElevation = -11000;
	public float MaxElevation = 10000;
	public float MaxDepth = 11000;
	public float maxCloudColor = 300.0f;
	public float WaterDepthThreshold = 10;
	public Color32 WallColor = new Color32(50,50,50,255);

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
	public float DisplayFloraMax = 100;
	public float DisplayTemperatureMin = 223;
	public float DisplayTemperatureMax = 323;
	public float DisplayWaterCarbonMax = 0.001f;
	public float DisplayAbsoluteHumidityMax = 0.05f;
	public float DisplayAirPressureMin = 97000;
	public float DisplayAirPressureMax = 110000;
	public float DisplayHeatAbsorbedMax = 1000;
	public float DisplayCrustDepthMax = 10000;
	public float DisplayMagmaMassMax = 1000000;
	public float DisplayDustHeight = 1000;
	public float DisplayDustMax = 100;
	public float DisplayLavaTemperatureMax = 1200;
	public float DisplayCarbonDioxideMax = 0.002f;
	public float DisplayOxygenMax = 0.35f;
	public float DisplaySoilFertilityMax = 5.0f;

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

	private JobHelper _renderJobHelper;
	private JobHelper _renderJobVertsHelper;

	private RenderState[] _renderStates;
	private int _curRenderState;
	private int _lastRenderState;
	private int _nextRenderState;
	private const int _renderStateCount = 3;

	private NativeArray<Vector3> _terrainVerticesArray;
	private NativeArray<Color32> _terrainColorsArray;
	private NativeArray<Vector3> _waterVerticesArray;
	private NativeArray<Vector3> _waterNormalsArray;
	private NativeArray<Color32> _waterColorsArray;
	private NativeArray<Vector3> _cloudVerticesArray;
	private NativeArray<Color32> _cloudColorsArray;
	private NativeArray<Vector3> _lavaVerticesArray;
	private NativeArray<Color32> _lavaColorsArray;
	private NativeArray<Vector3> _dustVerticesArray;
	private NativeArray<Color32> _dustColorsArray;

	private Vector3[] _terrainVertices;
	private Color32[] _terrainColors;
	private Vector3[] _waterVertices;
	private Vector3[] _waterNormals;
	private Color32[] _waterColors;
	private Vector3[] _cloudVertices;
	private Color32[] _cloudColors;
	private Vector3[] _lavaVertices;
	private Color32[] _lavaColors;
	private Vector3[] _dustVertices;
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
	private int[] _indices;
	private int[] _indicesTerrain;

	private NativeArray<CVP> _normalizedRainbow;
	private NativeArray<CVP> _normalizedBlueBlackRed;
	private NativeArray<float3> _hexVerts;

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

		InitVerts(Sim.Icosphere);

		_renderJobHelper = new JobHelper(Sim.CellCount);
		_renderJobVertsHelper = new JobHelper(Sim.CellCount* VertsPerCell);
		_terrainVerticesArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_terrainColorsArray = new NativeArray<Color32>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_waterVerticesArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_waterNormalsArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_waterColorsArray = new NativeArray<Color32>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_cloudVerticesArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_cloudColorsArray = new NativeArray<Color32>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_lavaVerticesArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_lavaColorsArray = new NativeArray<Color32>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_dustVerticesArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_dustColorsArray = new NativeArray<Color32>(Sim.CellCount * VertsPerCell, Allocator.Persistent);

		BuildRenderState(ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState, ref _renderStates[0], ref Sim.WorldData, ref Sim.StaticState);
		UpdateMesh(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], ref _renderStates[_curRenderState]);

		const int boundsSize = 1000;
		_terrainMesh.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));
		_waterMesh.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));
		_cloudMesh.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));
		_dustMesh.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));
		_lavaMesh.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));

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

		_hexVerts.Dispose();
		_terrainVerticesArray.Dispose();
		_terrainColorsArray.Dispose();
		_waterVerticesArray.Dispose();
		_waterNormalsArray.Dispose();
		_waterColorsArray.Dispose();
		_cloudVerticesArray.Dispose();
		_cloudColorsArray.Dispose();
		_lavaVerticesArray.Dispose();
		_lavaColorsArray.Dispose();
		_dustVerticesArray.Dispose();
		_dustColorsArray.Dispose();

	}

	public void Update()
	{
		if (Sim.TimeScale == 0)
		{
			_tickLerpTime -= Time.deltaTime * 0.5f;
		} else
		{
			_tickLerpTime -= Time.deltaTime * 0.5f;
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

		var buildRenderStateJobHandle = _renderJobHelper.Schedule(new BuildRenderStateCellJob()
		{
			TerrainColor = to.TerrainColor,
			TerrainElevation = to.TerrainElevation,
			LavaColor = to.LavaColor,
			LavaElevation = to.LavaElevation,
			WaterColor = to.WaterColor,
			WaterElevation = to.WaterElevation,
			CloudColor = to.CloudColor,
			CloudElevation = to.CloudElevation,
			DustColor = to.DustColor,
			DustElevation = to.DustElevation,
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
			CloudElevationSim = dependent.CloudElevation,
			Icosphere = Sim.Icosphere.Vertices,
			SoilFertility = from.GroundCarbon,
			Roughness = from.Roughness,
			Elevation = from.Elevation,
			CloudDropletMass = from.CloudDropletMass,
			CloudAbsorption = dependent.CloudAbsorptivity,
			IceCoverage = dependent.IceCoverage,
			FloraCoverage = dependent.FloraCoverage,
			WaterCoverage = dependent.WaterCoverage[Sim.WorldData.WaterLayers - 2],
			WaterDepth = dependent.WaterLayerDepth[1],
			WaterTemperature = from.WaterTemperature[Sim.WorldData.WaterLayers - 2],
			PlanktonMass = from.PlanktonMass[Sim.WorldData.WaterLayers - 2],
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
		});

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
		dependencies.Add((new LerpJobfloat { Progress = t, Out = state.TerrainElevation, Start = lastState.TerrainElevation, End = nextState.TerrainElevation }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.TerrainColor, Start = lastState.TerrainColor, End = nextState.TerrainColor }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat { Progress = t, Out = state.WaterElevation, Start = lastState.WaterElevation, End = nextState.WaterElevation }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.WaterColor, Start = lastState.WaterColor, End = nextState.WaterColor }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat { Progress = t, Out = state.DustElevation, Start = lastState.DustElevation, End = nextState.DustElevation }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.DustColor, Start = lastState.DustColor, End = nextState.DustColor }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat { Progress = t, Out = state.LavaElevation, Start = lastState.LavaElevation, End = nextState.LavaElevation }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.LavaColor, Start = lastState.LavaColor, End = nextState.LavaColor }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat3 { Progress = t, Out = state.SurfacePosition, Start = lastState.SurfacePosition, End = nextState.SurfacePosition }).Schedule(Sim.CellCount, _batchCount));

		if (true /* cloudsVisible*/)
		{
			dependencies.Add((new LerpJobfloat { Progress = t, Out = state.CloudElevation, Start = lastState.CloudElevation, End = nextState.CloudElevation }).Schedule(Sim.CellCount, _batchCount));
			dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.CloudColor, Start = lastState.CloudColor, End = nextState.CloudColor }).Schedule(Sim.CellCount, _batchCount));
		}
		if (ActiveWindOverlay != WindOverlay.None)
		{
			dependencies.Add((new LerpJobfloat3 { Progress = t, Out = state.VelocityArrow, Start = lastState.VelocityArrow, End = nextState.VelocityArrow }).Schedule(Sim.CellCount, _batchCount));
		}

		var getVertsHandle = _renderJobVertsHelper.Schedule(new BuildHexVertsJob()
		{
			VTerrainPosition = _terrainVerticesArray,
			VTerrainColor = _terrainColorsArray,
			VWaterPosition = _waterVerticesArray,
			VWaterNormal = _waterNormalsArray,
			VWaterColor = _waterColorsArray,
			VCloudPosition = _cloudVerticesArray,
			VCloudColor = _cloudColorsArray,
			VLavaPosition = _lavaVerticesArray,
			VLavaColor = _lavaColorsArray,
			VDustPosition = _dustVerticesArray,
			VDustColor = _dustColorsArray,

			TerrainElevation = state.TerrainElevation,
			TerrainColor = state.TerrainColor,
			WaterElevation = state.WaterElevation,
			WaterColor = state.WaterColor,
			CloudElevation = state.CloudElevation,
			CloudColor = state.CloudColor,
			LavaElevation = state.LavaElevation,
			LavaColor = state.LavaColor,
			DustElevation = state.DustElevation,
			DustColor = state.DustColor,
			HexVerts = _hexVerts,
			IcosphereVerts = Sim.Icosphere.Vertices,
			WallColor = WallColor
		}, JobHandle.CombineDependencies(dependencies));

		getVertsHandle.Complete();
		dependencies.Dispose();

		_terrainVerticesArray.CopyTo(_terrainVertices);
		_terrainColorsArray.CopyTo(_terrainColors);
		_waterVerticesArray.CopyTo(_waterVertices);
		_waterNormalsArray.CopyTo(_waterNormals);
		_waterColorsArray.CopyTo(_waterColors);
		_cloudVerticesArray.CopyTo(_cloudVertices);
		_cloudColorsArray.CopyTo(_cloudColors);
		_lavaVerticesArray.CopyTo(_lavaVertices);
		_lavaColorsArray.CopyTo(_lavaColors);
		_dustVerticesArray.CopyTo(_dustVertices);
		_dustColorsArray.CopyTo(_dustColors);

		_terrainMesh.vertices = _terrainVertices;
		_terrainMesh.colors32 = _terrainColors;

		_waterMesh.vertices = _waterVertices;
		_waterMesh.normals = _waterNormals;
		_waterMesh.colors32 = _waterColors;

		_cloudMesh.vertices = _cloudVertices;
		_cloudMesh.colors32 = _cloudColors;

		_lavaMesh.vertices = _lavaVertices;
		_lavaMesh.colors32 = _lavaColors;

		_dustMesh.vertices = _dustVertices;
		_dustMesh.colors32 = _dustColors;

		if (!_indicesInitialized)
		{
			_terrainMesh.SetTriangles(_indicesTerrain, 0);
			_waterMesh.SetTriangles(_indicesTerrain, 0);
			_lavaMesh.SetTriangles(_indicesTerrain, 0);
			_cloudMesh.SetTriangles(_indices, 0);
			_dustMesh.SetTriangles(_indices, 0);
			_indicesInitialized = true;
		}

		Planet.transform.SetPositionAndRotation(state.Position, Quaternion.Euler(state.Rotation));

		//_terrainMesh.RecalculateBounds();
		//_waterMesh.RecalculateBounds();
		//_cloudMesh.RecalculateBounds();
		//_lavaMesh.RecalculateBounds();
		//_dustMesh.RecalculateBounds();

		_terrainMesh.RecalculateNormals();
	//	_waterMesh.RecalculateNormals();
		_cloudMesh.RecalculateNormals();
		_lavaMesh.RecalculateNormals();
		_dustMesh.RecalculateNormals();
		_terrainMesh.RecalculateTangents();
		_waterMesh.RecalculateTangents();

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
		ActiveTemperatureUnits = (CellInfo.TemperatureUnits)dropdown.value;
	}

	public string GetCellInfo(CellInfo.CellInfoType cellInfoType)
	{
		switch (cellInfoType)
		{
			case CellInfo.CellInfoType.Global:
				return CellInfo.GetCellInfoGlobal(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Enthalpy:
				return CellInfo.GetCellInfoEnthalpy(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Energy:
				return CellInfo.GetCellInfoEnergy(ActiveTemperatureUnits, Sim.InverseCellCount, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Cell:
				return CellInfo.GetCellInfoCell(ActiveTemperatureUnits, ActiveCellIndex, ref Sim.StaticState, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Atmosphere:
				return CellInfo.GetCellInfoAtmosphere(ActiveTemperatureUnits, ActiveCellIndex, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.StaticState, ref Sim.DisplayState);
			case CellInfo.CellInfoType.Ground:
				return CellInfo.GetCellInfoGround(ActiveTemperatureUnits, ActiveCellIndex, ref Sim.ActiveSimState, ref Sim.DependentState);
			case CellInfo.CellInfoType.Water:
				return CellInfo.GetCellInfoWater(ActiveTemperatureUnits, ActiveCellIndex, ref Sim.WorldData, ref Sim.ActiveSimState, ref Sim.DependentState, ref Sim.StaticState, ref Sim.DisplayState);
		}
		return "";
	}

	public int GetClosestVert(int triangleIndex, int vIndex)
	{
		return _indicesTerrain[triangleIndex * 3 + vIndex] / 8;
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

	private void InitVerts(Icosphere icosphere)
	{
		_terrainVertices = new Vector3[Sim.CellCount * VertsPerCell];
		_terrainColors = new Color32[Sim.CellCount * VertsPerCell];
		_waterVertices = new Vector3[Sim.CellCount * VertsPerCell];
		_waterNormals = new Vector3[Sim.CellCount * VertsPerCell];
		_waterColors = new Color32[Sim.CellCount * VertsPerCell];
		_cloudVertices = new Vector3[Sim.CellCount * VertsPerCell];
		_cloudColors = new Color32[Sim.CellCount * VertsPerCell];
		_dustVertices = new Vector3[Sim.CellCount * VertsPerCell];
		_dustColors = new Color32[Sim.CellCount * VertsPerCell];
		_lavaVertices = new Vector3[Sim.CellCount * VertsPerCell];
		_lavaColors = new Color32[Sim.CellCount * VertsPerCell];

		_hexVerts = new NativeArray<float3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		List<int> indices = new List<int>();
		List<int> indicesTerrain = new List<int>();
		for (int i=0;i< icosphere.Vertices.Length;i++)
		{
			float3 pos = icosphere.Vertices[i];
			_hexVerts[i * VertsPerCell] = pos;
			int neighborCount = (icosphere.Neighbors[(i + 1) * MaxNeighbors - 1] >= 0) ? MaxNeighbors : (MaxNeighbors - 1);
			for (int j=0;j< neighborCount; j++)
			{
				int neighborIndex1 = icosphere.Neighbors[i * MaxNeighbors + j];
				int neighborIndex2 = icosphere.Neighbors[i * MaxNeighbors + (j + 1) % neighborCount];

				float3 midPoint = (icosphere.Vertices[neighborIndex1] + icosphere.Vertices[neighborIndex2] + pos * (1 + SlopeAmount)) / (3 + SlopeAmount);
				float midPointLength = math.length(midPoint);
				float3 extendedMidPoint = midPoint / (midPointLength * midPointLength);
				_hexVerts[i * VertsPerCell + 1 + j] = extendedMidPoint; // surface
				_hexVerts[i * VertsPerCell + 1 + j + MaxNeighbors] = extendedMidPoint; // wall
				_hexVerts[i * VertsPerCell + 1 + j + MaxNeighbors * 2] = extendedMidPoint; // wall
				_hexVerts[i * VertsPerCell + 1 + j + MaxNeighbors * 3] = extendedMidPoint; // corner

				indices.Add(i * VertsPerCell);
				indices.Add(i * VertsPerCell + 1 + ((j + 1) % neighborCount));
				indices.Add(i * VertsPerCell + 1 + j);

				indicesTerrain.Add(i * VertsPerCell);
				indicesTerrain.Add(i * VertsPerCell + 1 + ((j + 1) % neighborCount));
				indicesTerrain.Add(i * VertsPerCell + 1 + j);


				int neighbor1 = -1;
				{
					int neighborNeighborCount = (icosphere.Neighbors[(neighborIndex1 + 1) * MaxNeighbors - 1] >= 0) ? MaxNeighbors : (MaxNeighbors - 1);
					for (int k = 0; k < neighborNeighborCount; k++)
					{
						if (icosphere.Neighbors[neighborIndex1 * MaxNeighbors + k] == i)
						{
							neighbor1 = (k - 1 + neighborNeighborCount) % neighborNeighborCount;
							indicesTerrain.Add(i * VertsPerCell + 1 + 2 * MaxNeighbors + ((j - 1 + neighborCount) % neighborCount));
							indicesTerrain.Add(i * VertsPerCell + 1 + MaxNeighbors + j);
							indicesTerrain.Add(neighborIndex1 * VertsPerCell + 1 + MaxNeighbors + k);

							break;
						}
					}
				}
				if (neighbor1 >= 0 && i < neighborIndex1 && i < neighborIndex2)
				{
					int neighborNeighborCount = (icosphere.Neighbors[(neighborIndex2 + 1) * MaxNeighbors - 1] >= 0) ? MaxNeighbors : (MaxNeighbors - 1);
					for (int k = 0; k < neighborNeighborCount; k++)
					{
						if (icosphere.Neighbors[neighborIndex2 * MaxNeighbors + k] == i)
						{
							indicesTerrain.Add(i * VertsPerCell + 1 + 3 * MaxNeighbors + j);
							indicesTerrain.Add(neighborIndex2 * VertsPerCell + 1 + 3 * MaxNeighbors + k);
							indicesTerrain.Add(neighborIndex1 * VertsPerCell + 1 + 3 * MaxNeighbors + neighbor1);

							break;
						}
					}
				}

			}
		}
		_indices = indices.ToArray();
		_indicesTerrain = indicesTerrain.ToArray();

	}

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
			case MeshOverlay.WaterCarbonDioxide0:
				overlay = new MeshOverlayData(0, DisplayWaterCarbonMax, _normalizedRainbow, display.WaterCarbonDioxidePercent[Sim.WorldData.WaterLayers - 2]);
				return true;
			case MeshOverlay.WaterCarbonDioxide1:
				overlay = new MeshOverlayData(0, DisplayWaterCarbonMax, _normalizedRainbow, display.WaterCarbonDioxidePercent[Sim.WorldData.WaterLayers - 3]);
				return true;
			case MeshOverlay.WaterCarbonDioxide2:
				overlay = new MeshOverlayData(0, DisplayWaterCarbonMax, _normalizedRainbow, display.WaterCarbonDioxidePercent[Sim.WorldData.WaterLayers - 4]);
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
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.GroundTemperature);
				return true;
			case MeshOverlay.CarbonDioxide0:
				overlay = new MeshOverlayData(0, DisplayCarbonDioxideMax, _normalizedRainbow, display.CarbonDioxidePercent[1]);
				return true;
			case MeshOverlay.CarbonDioxide1:
				overlay = new MeshOverlayData(0, DisplayCarbonDioxideMax, _normalizedRainbow, display.CarbonDioxidePercent[2]);
				return true;
			case MeshOverlay.CarbonDioxide2:
				overlay = new MeshOverlayData(0, DisplayCarbonDioxideMax, _normalizedRainbow, display.CarbonDioxidePercent[3]);
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
				overlay = new WindOverlayData(DisplayWindSpeedUpperAirMax, false, dependentState.CloudVelocity);
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


	#endregion
}
