﻿using System;
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
using UnityEngine.Profiling;
using UnityEngine.VFX;
using UnityEngine.Rendering;

public class WorldView : MonoBehaviour {

	
	public enum MeshOverlay {
		None,
		TemperatureSurface,
		PotentialTemperature,
		Pressure,
		AbsoluteHumidity,
		RelativeHumidity,
		CarbonDioxide,
		Evaporation,
		Rainfall,
		Condensation,
		HeatAbsorbed,
		WaterTemperature,
		WaterCarbonDioxide,
		Salinity,
		GroundTemperature,
		GroundWater,
		GroundWaterTemperature,
		CrustDepth,
		MagmaMass,
		DivergenceAir,
		DivergenceWater,
		PlateTectonics
	}

	public enum WindOverlay {
		None,
		Current,
		Wind,
		PGF,
		WindCloud,

	}

	public const int VertsPerCell = 25;
	public const int VertsPerCloud = 25;
	public const int MaxNeighbors = 6;

	public bool LerpStates = true;


	[Header("Clouds")]
	public int CloudLayers = 10;
	public float CloudVerticalScale = 0.1f;
	public float CloudLayerOpacityPower = 0.5f;
	public float CloudHeight = 0.1f;
	public float DisplayCloudMin = 0.1f;
	public float DisplayCloudPower = 0.5f;
	public float AtmosphereScale = 1000f;
	public float SlopeAmountCloud = 10;

	[Header("Terrain")]
	public float SlopeAmountTerrain = 3;
	public float TerrainScale = 100f;
	public List<Color> PlateColors;

	[Header("Display")]
	public float DisplayFloraWeight = 1;
	public float DisplaySandWeight = 1;
	public float DisplaySoilWeight = 1;
	public float DisplayDustMax = 100;
	public float DisplayLavaTemperatureMax = 1200;
	public float DisplaySoilFertilityMax = 5.0f;
	public float DisplayPlanktonMax = 10;
	public float DisplayPlanktonPower = 0.5f;
	public int DisplayPlanktonLevels = 5;
	public int DisplayIceLevels = 5;
	public int DisplayIcePower = 5;
	public float DisplayWaterTemperatureMax = 50;
	public int DisplayWaterTemperatureLevels = 10;

	[Header("Wind Overlay")]
	public float DisplayWindSpeedLowerAirMax = 50;
	public float DisplayWindSpeedUpperAirMax = 250;
	public float DisplayPressureGradientForceMax = 0.01f;
	public float DisplayWindSpeedSurfaceWaterMax = 5;
	public float DisplayWindSpeedDeepWaterMax = 0.5f;
	public float DisplayVerticalWindSpeedMax = 0.5f;
	public float DisplayVerticalWindSpeedMin = 0.05f;

	[Header("Overlays")]
	public float DisplayRainfallMax = 5.0f;
	public float DisplaySalinityMin = 0;
	public float DisplaySalinityMax = 50;
	public float DisplayEvaporationMax = 5.0f;
	public float DisplayTemperatureMin = 223;
	public float DisplayTemperatureMax = 323;
	public float DisplayWaterCarbonMax = 0.001f;
	public float DisplayAbsoluteHumidityMax = 0.05f;
	public float DisplayAirPressureMin = 97000;
	public float DisplayAirPressureMax = 110000;
	public float DisplayHeatAbsorbedMax = 1000;
	public float DisplayCrustDepthMax = 10000;
	public float DisplayMagmaMassMax = 1000000;
	public float DisplayCarbonDioxideMax = 0.002f;
	public float DisplayOxygenMax = 0.35f;
	public float DisplayDivergenceMax = 1000;

	[Header("References")]
	public WorldSimComponent Sim;
	public GameplayManager GameplayManager;
	public FoliageManager Foliage;
	public GameObject Planet;
	public GameObject Sun;
	public GameObject Moon;
	public GameObject SunLight;
	public Material TerrainMaterial;
	public Material WaterMaterialFront, WaterMaterialBack;
	public Material CloudMaterialFront, CloudMaterialBack;
	public Material OverlayMaterial;
	public GameObject SelectionCirclePrefab;
	public GameObject WindArrowPrefab;
	public Volume StandardPostProcessingVolume;
	public Volume OverlayPostProcessingVolume;

	[Header("Debug")]
	public bool RunSynchronously;

	[HideInInspector] public WindOverlay ActiveWindOverlay;
	private MeshOverlay _activeMeshOverlay;
	[HideInInspector] public MeshOverlay ActiveMeshOverlay { get { return _activeMeshOverlay; } set { _activeMeshOverlay = value; OnMeshOverlayChanged?.Invoke(value); } }
	[HideInInspector] public int ActiveMeshLayerWater = 1;
	[HideInInspector] public int ActiveMeshLayerAir = 1;
	[HideInInspector] public DisplayState DisplayState;
	[HideInInspector] public RenderState CurRenderState { get { return _renderStates[_curRenderState]; } }
	public Dictionary<MeshOverlay, MeshOverlayColors> OverlayColors;

	private JobHelper _renderJobHelper;
	private JobHelper _renderTerrainHelper;
	private JobHelper _renderCloudHelper;

	private RenderState[] _renderStates;
	private int _curRenderState;
	private int _lastRenderState;
	private int _nextRenderState;
	private const int _renderStateCount = 3;

	private NativeArray<Vector3> _overlayVerticesArray;
	private NativeArray<Vector4> _overlayUVArray;
	private NativeArray<Vector3> _terrainVerticesArray;
	private NativeArray<Vector4> _terrainColorsArray1;
	private NativeArray<Vector4> _terrainColorsArray2;
	private NativeArray<Vector3> _waterVerticesArray;
	private NativeArray<Vector3> _waterNormalsArray;
	private NativeArray<Vector4> _waterColorsArray;
	private NativeArray<Vector4> _waterCurrentArray;
	private NativeArray<Vector3> _cloudVerticesArray;
	private NativeArray<Color32> _cloudColorsArray;
	private NativeArray<Vector3> _cloudNormalsArray;
	private NativeArray<Color32> _overlayColorsArray;

	private Vector3[] _terrainVertices;
	private Vector3[] _waterVertices;
	private Vector3[] _waterNormals;
	private Color32[] _waterColors;
	private Vector3[] _cloudVertices;
	private Color32[] _cloudColors;
	private Vector3[] _cloudNormals;

	private Mesh _overlayMesh;
	private Mesh _terrainMesh;
	private Mesh _waterMeshFront, _waterMeshBack;
	private Mesh _cloudMeshFront, _cloudMeshBack;

	private GameObject _terrainObject, _overlayObject;
	private GameObject _waterObjectFront, _waterObjectBack;
	private List<GameObject> _cloudObjectFront = new List<GameObject>();
	private List<GameObject> _cloudObjectBack = new List<GameObject>();
	private GameObject _selectionCircle;
	private List<GameObject> _selectionCircleNeighbors;

	private GameObject[] _windArrows;

	private float _tickLerpTime = 0;
	private float _tickLerpTimeTotal = 1;

	private bool _indicesInitialized;
	private int[] _indices;
	private int[] _indicesTerrain;
	private int[] _indicesTerrainBack;
	private int[] _indicesCloudFront;
	private int[] _indicesCloudBack;

	private NativeArray<CVP> _plateColors;
	private NativeArray<CVP> _normalizedRainbow;
	private NativeArray<CVP> _normalizedBlueBlackRed;
	private NativeArray<float3> _terrainVerts;
	private NativeArray<float3> _cloudVerts;

	private const int _batchCount = 128;

	private bool _skyboxActive = false;
	private float _skyboxExposure = 1;
	private float _skyboxExposureDest = 1;
	private float _skyboxExposureStart = 1;
	private float _skyboxExposureTime = 0;

	public event Action<MeshOverlay> OnMeshOverlayChanged;

	public void Awake()
	{
		Sim.OnTick += OnSimTick;
		Sim.OnTimeScaleChanged += OnTimeScaleChanged;
		GameplayManager.OnSetActiveCell += OnSetActiveCell;
		GameplayManager.OnSetActiveTool += OnSetActiveTool;

		_normalizedRainbow = new NativeArray<CVP>(new CVP[] {
											new CVP(Color.black, 0),
											new CVP(new Color(0.25f,0,0.5f,1), 1.0f / 7),
											new CVP(Color.blue, 2.0f / 7),
											new CVP(Color.green, 3.0f / 7),
											new CVP(Color.yellow, 4.0f / 7),
											new CVP(Color.red, 5.0f / 7),
											new CVP(Color.magenta, 6.0f / 7),
											new CVP(Color.white, 1),
											},
											Allocator.Persistent);
		_plateColors = new NativeArray<CVP>(new CVP[] {
											new CVP(Color.black, 0),
											new CVP(new Color(0.25f,0,0.5f,1), 1.0f / 7),
											new CVP(Color.blue, 2.0f / 7),
											new CVP(Color.green, 3.0f / 7),
											new CVP(Color.yellow, 4.0f / 7),
											new CVP(Color.red, 5.0f / 7),
											new CVP(Color.magenta, 6.0f / 7),
											new CVP(Color.white, 1),
											},
											Allocator.Persistent);
		_normalizedBlueBlackRed = new NativeArray<CVP>(new CVP[] {
											new CVP(Color.blue, 0),
											new CVP(Color.black, 0.5f),
											new CVP(Color.red, 1) },
											Allocator.Persistent);

		OverlayColors = new Dictionary<MeshOverlay, MeshOverlayColors>
		{
			{ MeshOverlay.AbsoluteHumidity, new MeshOverlayColors{ Title="Absolute Humidity", Min=0,Max=DisplayAbsoluteHumidityMax,ColorValuePairs=_normalizedRainbow, LegendType=LegendType.Percent, DecimalPlaces=2 } },
			{ MeshOverlay.RelativeHumidity, new MeshOverlayColors{ Title="Relative Humidity", Min=0,Max=1,ColorValuePairs=_normalizedRainbow, LegendType=LegendType.Percent, DecimalPlaces=0 } },
			{ MeshOverlay.TemperatureSurface, new MeshOverlayColors{ Title="Surface Temperature", Min=DisplayTemperatureMin,Max=DisplayTemperatureMax,ColorValuePairs=_normalizedRainbow, LegendType=LegendType.Temperature, DecimalPlaces=0 } },
			{ MeshOverlay.PotentialTemperature, new MeshOverlayColors{ Title="Potential Temperature", Min=DisplayTemperatureMin,Max=DisplayTemperatureMax,ColorValuePairs=_normalizedRainbow, LegendType=LegendType.Temperature, DecimalPlaces=0 } },
			{ MeshOverlay.Pressure, new MeshOverlayColors{ Title="Pressure", Min=DisplayAirPressureMin,Max=DisplayAirPressureMax,ColorValuePairs=_normalizedRainbow, LegendType=LegendType.Pressure, DecimalPlaces=0 } },
			{ MeshOverlay.WaterTemperature, new MeshOverlayColors{ Title="Water Temperature", Min=WorldData.FreezingTemperature,Max=DisplayTemperatureMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.Temperature, DecimalPlaces=0 } },
			{ MeshOverlay.WaterCarbonDioxide, new MeshOverlayColors{Title="Water CO2",  Min=0,Max=DisplayWaterCarbonMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.PPM, DecimalPlaces=0 } },
			{ MeshOverlay.Salinity, new MeshOverlayColors{ Title="Salinity", Min=DisplaySalinityMin,Max=DisplaySalinityMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.PPM, DecimalPlaces=0 } },
			{ MeshOverlay.GroundTemperature, new MeshOverlayColors{ Title="Ground Temperature", Min=DisplayTemperatureMin,Max=DisplayTemperatureMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.Temperature, DecimalPlaces=0 } },
			{ MeshOverlay.CarbonDioxide, new MeshOverlayColors{ Title="Carbon Dioxide", Min=0,Max=DisplayCarbonDioxideMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.PPM, DecimalPlaces=0  } },
			{ MeshOverlay.HeatAbsorbed, new MeshOverlayColors{ Title="Heat Absorbed", Min=0,Max=DisplayHeatAbsorbedMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.Watts, DecimalPlaces=0 } },
			{ MeshOverlay.Rainfall, new MeshOverlayColors{ Title="Rainfall", Min=0,Max=DisplayRainfallMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.Volume, DecimalPlaces=2 } },
			{ MeshOverlay.Evaporation, new MeshOverlayColors{ Title="Evaporation", Min=0,Max=DisplayEvaporationMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.Volume, DecimalPlaces=2 } },
			{ MeshOverlay.GroundWater, new MeshOverlayColors{ Title="Ground Water", Min=0,Max=Sim.WorldData.GroundWaterMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.Mass, DecimalPlaces=0 } },
			{ MeshOverlay.GroundWaterTemperature, new MeshOverlayColors{ Title="Ground Water Temperature", Min=DisplayTemperatureMin,Max=DisplayTemperatureMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.Temperature, DecimalPlaces=0 } },
			{ MeshOverlay.CrustDepth, new MeshOverlayColors{ Title="Crust Depth", Min=DisplayCrustDepthMax,Max=0,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.None, DecimalPlaces=0 } },
			{ MeshOverlay.MagmaMass, new MeshOverlayColors{ Title="Magma", Min=0,Max=DisplayMagmaMassMax,ColorValuePairs=_normalizedRainbow,LegendType=LegendType.Mass, DecimalPlaces=0 } },
			{ MeshOverlay.DivergenceAir, new MeshOverlayColors{ Title="Divergence Air", Min=-DisplayDivergenceMax,Max=DisplayDivergenceMax,ColorValuePairs=_normalizedBlueBlackRed,LegendType=LegendType.None, DecimalPlaces=2 } },
			{ MeshOverlay.DivergenceWater, new MeshOverlayColors{ Title="Divergence Water", Min=-DisplayDivergenceMax,Max=DisplayDivergenceMax,ColorValuePairs=_normalizedBlueBlackRed,LegendType=LegendType.None, DecimalPlaces=2 } },
			{ MeshOverlay.PlateTectonics, new MeshOverlayColors{ Title="Plate Tectonics", Min=0,Max=11,ColorValuePairs=_plateColors,LegendType=LegendType.None, DecimalPlaces=0 } },
		};
		_overlayMesh = new Mesh();
		_terrainMesh = new Mesh();
		_cloudMeshFront = new Mesh();
		_cloudMeshBack = new Mesh();
		_waterMeshFront = new Mesh();
		_waterMeshBack = new Mesh();

		var cloudMeshesBack = new GameObject("Cloud Meshes Back");
		cloudMeshesBack.transform.SetParent(Planet.transform, false);

		var cloudMeshesFront = new GameObject("Cloud Meshes Front");
		cloudMeshesFront.transform.SetParent(Planet.transform, false);

		_terrainObject = new GameObject("Terrain Mesh");
		_terrainObject.transform.SetParent(Planet.transform, false);
		var terrainFilter = _terrainObject.AddComponent<MeshFilter>();
		var terrainSurfaceRenderer = _terrainObject.AddComponent<MeshRenderer>();
		var terrainCollider = _terrainObject.AddComponent<MeshCollider>();
		terrainSurfaceRenderer.material = TerrainMaterial;
		terrainFilter.mesh = terrainCollider.sharedMesh = _terrainMesh;
		_terrainMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
		_cloudMeshBack.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
		_cloudMeshFront.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;

		for (int i = 0; i < CloudLayers; i++)
		{
			var c = new GameObject("Cloud Mesh");
			c.transform.SetParent(cloudMeshesBack.transform, false);
			var cloudFilter = c.AddComponent<MeshFilter>();
			var cloudSurfaceRenderer = c.AddComponent<MeshRenderer>();
			cloudSurfaceRenderer.material = CloudMaterialBack;
			cloudFilter.sharedMesh = _cloudMeshBack;
			cloudSurfaceRenderer.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
			_cloudObjectBack.Add(c);
		}

		{
			_waterObjectFront = new GameObject("Water Mesh Front");
			_waterObjectFront.transform.SetParent(Planet.transform, false);
			var waterFilter = _waterObjectFront.AddComponent<MeshFilter>();
			var waterSurfaceRenderer = _waterObjectFront.AddComponent<MeshRenderer>();
			var waterCollider = _waterObjectFront.AddComponent<MeshCollider>();
			waterSurfaceRenderer.material = WaterMaterialFront;
			waterFilter.mesh = waterCollider.sharedMesh = _waterMeshFront;
			waterSurfaceRenderer.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
			_waterMeshFront.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
		}
		{
			_waterObjectBack = new GameObject("Water Mesh Back");
			_waterObjectBack.transform.SetParent(Planet.transform, false);
			var waterFilter = _waterObjectBack.AddComponent<MeshFilter>();
			var waterSurfaceRenderer = _waterObjectBack.AddComponent<MeshRenderer>();
			var waterCollider = _waterObjectBack.AddComponent<MeshCollider>();
			waterSurfaceRenderer.material = WaterMaterialBack;
			waterFilter.mesh = waterCollider.sharedMesh = _waterMeshBack;
			waterSurfaceRenderer.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
			_waterMeshBack.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
		}

		_overlayObject = new GameObject("Overlay Mesh");
		_overlayObject.transform.SetParent(Planet.transform, false);
		var overlayFilter = _overlayObject.AddComponent<MeshFilter>();
		var overlaySurfaceRenderer = _overlayObject.AddComponent<MeshRenderer>();
		overlaySurfaceRenderer.material = OverlayMaterial;
		overlayFilter.mesh = _overlayMesh;
		_overlayMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
		_overlayObject.SetActive(false);



		for (int i = 0; i < CloudLayers; i++)
		{
			var c = new GameObject("Cloud Mesh");
			c.transform.SetParent(cloudMeshesFront.transform, false);
			var cloudFilter = c.AddComponent<MeshFilter>();
			var cloudSurfaceRenderer = c.AddComponent<MeshRenderer>();
			cloudSurfaceRenderer.material = CloudMaterialFront;
			cloudFilter.sharedMesh = _cloudMeshFront;
			cloudSurfaceRenderer.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
			_cloudObjectFront.Add(c);
		}


		_selectionCircle = GameObject.Instantiate(SelectionCirclePrefab, Planet.transform);
		_selectionCircle.transform.localScale *= 0.02f;

		_selectionCircleNeighbors = new List<GameObject>();
		for (int i = 0; i < 6; i++)
		{
			var s = GameObject.Instantiate(SelectionCirclePrefab, Planet.transform);
			s.transform.localScale *= 0.01f;
			_selectionCircleNeighbors.Add(s);
		}

		_windArrows = new GameObject[Sim.CellCount];
		for (int i = 0; i < Sim.CellCount; i++)
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
		_renderTerrainHelper = new JobHelper(Sim.CellCount * VertsPerCell);
		_renderCloudHelper = new JobHelper(Sim.CellCount * VertsPerCloud);
		_overlayVerticesArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_overlayUVArray = new NativeArray<Vector4>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_terrainVerticesArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_terrainColorsArray1 = new NativeArray<Vector4>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_terrainColorsArray2 = new NativeArray<Vector4>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_waterVerticesArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_waterNormalsArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_waterCurrentArray = new NativeArray<Vector4>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_waterColorsArray = new NativeArray<Vector4>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_overlayColorsArray = new NativeArray<Color32>(Sim.CellCount * VertsPerCell, Allocator.Persistent);

		_cloudVerticesArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCloud, Allocator.Persistent);
		_cloudColorsArray = new NativeArray<Color32>(Sim.CellCount * VertsPerCloud, Allocator.Persistent);
		_cloudNormalsArray = new NativeArray<Vector3>(Sim.CellCount * VertsPerCloud, Allocator.Persistent);

		DisplayState = new DisplayState();
		DisplayState.Init(Sim.CellCount, ref Sim.WorldData);

		Foliage.Init(Sim.CellCount, ref Sim.StaticState);

		OnSimTick();

		//		BuildRenderState(ref Sim.ActiveSimState, ref _tempState, ref DisplayState, ref _renderStates[0], ref Sim.WorldData, ref Sim.StaticState);

		LerpMesh(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], ref _renderStates[_curRenderState]);

		const int boundsSize = 1000;
		_terrainMesh.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));
		_waterMeshFront.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));
		_waterMeshBack.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));
		_cloudMeshFront.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize)*1.1f);
		_cloudMeshBack.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize) * 0.9f);



	}
	public void OnDestroy()
	{
		Sim.OnTick -= OnSimTick;

		_normalizedRainbow.Dispose();
		_normalizedBlueBlackRed.Dispose();
		_plateColors.Dispose();

		for (int i=0;i<_renderStateCount;i++)
		{
			_renderStates[i].Dispose();
		}

		_terrainVerts.Dispose();
		_cloudVerts.Dispose();
		_overlayVerticesArray.Dispose();
		_overlayUVArray.Dispose();
		_terrainVerticesArray.Dispose();
		_terrainColorsArray1.Dispose();
		_terrainColorsArray2.Dispose();
		_waterVerticesArray.Dispose();
		_waterNormalsArray.Dispose();
		_waterColorsArray.Dispose();
		_waterCurrentArray.Dispose();
		_cloudVerticesArray.Dispose();
		_cloudColorsArray.Dispose();
		_cloudNormalsArray.Dispose();
		_overlayColorsArray.Dispose();

		Foliage.Dispose();
		DisplayState.Dispose();
	}

	public void Update()
	{
		if (Sim.TimeScale == 0)
		{
			_tickLerpTime -= Time.deltaTime * 0.5f;
		}
		else
		{
			_tickLerpTime -= Time.deltaTime * 0.5f * Sim.TimeScale;
		}

		LerpMesh(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], ref _renderStates[_curRenderState]);
		Foliage.Update(ref _renderStates[_curRenderState]);

		OverlayPostProcessingVolume.weight = (ActiveMeshOverlay != MeshOverlay.None) ? 1 : 0;
		StandardPostProcessingVolume.weight = (ActiveMeshOverlay == MeshOverlay.None) ? 1 : 0;

		if (ActiveWindOverlay != WindOverlay.None)
		{
			for (int i = 0; i < _windArrows.Length; i++)
			{
				var windHorizontal = _renderStates[_curRenderState].VelocityHorizontal[i];
				float windVertical = _renderStates[_curRenderState].VelocityVertical[i];
				float windSpeed = math.length(windHorizontal);
				var pos = _renderStates[_curRenderState].SurfacePosition[i];
				bool visible = windSpeed > 0 || windVertical != 0;
				_windArrows[i].SetActive(visible);
				if (visible)
				{
					var circle = _windArrows[i].transform.GetChild(0);
					if (windVertical > 0)
					{
						var color = Color.Lerp(Color.white, Color.red, math.min(1, windVertical / DisplayVerticalWindSpeedMin));
						var m = circle.GetChild(0).GetComponent<MeshRenderer>().material.color = color;	
					} else
					{
						var color = Color.Lerp(Color.white, Color.blue, math.min(1, -windVertical / DisplayVerticalWindSpeedMin));
						var m = circle.GetChild(0).GetComponent<MeshRenderer>().material.color = color;
					}
					circle.localScale = Vector3.one * (1 + math.abs(windVertical));
					_windArrows[i].transform.localPosition = pos;

					_windArrows[i].transform.GetChild(1).gameObject.SetActive(windSpeed > 0);
					if (windSpeed > 0)
					{
						_windArrows[i].transform.localRotation = Quaternion.LookRotation(windHorizontal / windSpeed, pos);
						_windArrows[i].transform.GetChild(1).localScale = Vector3.one * math.min(math.pow(windSpeed, 0.5f), 1.0f);
					} else
					{
						_windArrows[i].transform.localRotation = Quaternion.LookRotation(math.cross(pos, Vector3.up), pos);
					}
				}
			}
		}


		Planet.transform.SetPositionAndRotation(_renderStates[_curRenderState].Position, Quaternion.Euler(_renderStates[_curRenderState].Rotation));
		UpdateSkybox();

		SunLight.transform.rotation = Quaternion.LookRotation(Planet.transform.position - SunLight.transform.position);

		CloudMaterialBack.SetFloat("Vector1_E122B9B2", Sim.TimeScale);// sim time scale
		CloudMaterialFront.SetFloat("Vector1_E122B9B2", Sim.TimeScale);// sim time scale
		WaterMaterialFront.SetFloat("Vector1_2C57E502", _renderStates[_curRenderState].Ticks); // sim time

		for (int i = 0; i < _cloudObjectFront.Count; i++)
		{
			float h = (float)i / _cloudObjectFront.Count;
			//float opacityMultiplier = h;
			//if (opacityMultiplier < CloudCurveHeight)
			//{
			//	opacityMultiplier *= math.pow(opacityMultiplier / CloudCurveHeight, CloudCurvePower) * CloudCurveMax + (1.0f - CloudCurveMax);
			//}
			//else
			//{
			//	opacityMultiplier *= math.pow((opacityMultiplier - CloudCurveHeight) / (1.0f - CloudCurveHeight), CloudCurvePower) * CloudCurveMax + (1.0f - CloudCurveMax);
			//}
			_cloudObjectFront[i].GetComponent<MeshRenderer>().material.SetFloat("Vector1_E31C0ED2", h * CloudVerticalScale);// height multiplier
			_cloudObjectFront[i].GetComponent<MeshRenderer>().material.SetFloat("Vector1_97917029", math.pow(h, CloudLayerOpacityPower));// opacity min
		}
		for (int i = 0; i < _cloudObjectBack.Count; i++)
		{
			float h = (float)i / _cloudObjectBack.Count;
			_cloudObjectBack[i].GetComponent<MeshRenderer>().material.SetFloat("Vector1_E31C0ED2", h * CloudVerticalScale);// height multiplier
			_cloudObjectBack[i].GetComponent<MeshRenderer>().material.SetFloat("Vector1_97917029", math.pow(h, CloudLayerOpacityPower));// opacity min
		}

	}

	public void StartLerp(float lerpTime)
	{
		_tickLerpTime = lerpTime;
		_tickLerpTimeTotal = lerpTime;

		_lastRenderState = _curRenderState;
		_nextRenderState = (_curRenderState + 1) % _renderStateCount;
		_curRenderState = (_nextRenderState + 1) % _renderStateCount;

		var displayJob = default(JobHandle);
		var lastDisplay = DisplayState;
	//	if (Sim.SimSettings.CollectOverlay)
		{
			DisplayState = new DisplayState();
			DisplayState.Init(Sim.StaticState.Count, ref Sim.WorldData);
			displayJob = DisplayState.Update(ref DisplayState, ref lastDisplay, ref Sim.WorldData, ref Sim.LastTempState, ref Sim.LastSimState, ref Sim.StaticState, ref Sim.SimSettings);
		}


		var buildRenderStateJob = BuildRenderState(RunSynchronously, ref Sim.LastSimState, ref Sim.LastTempState, ref DisplayState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState, displayJob);
		var foliageJob = Foliage.Tick(ref Sim.LastTempState, displayJob);

		displayJob = JobHandle.CombineDependencies(displayJob, buildRenderStateJob, foliageJob);
		displayJob.Complete();
	//	if (Sim.SimSettings.CollectOverlay)
		{
			lastDisplay.Dispose();
		}

	}

	private void OnSimTick()
	{

		float lerpTime;
		if (Sim.TimeScale == 0)
		{
			lerpTime = 0.1f;
		}
		else
		{
			lerpTime = 0.1f + Sim.TimeTillTick / Sim.TimeScale;
		}
		StartLerp(lerpTime);


	}

	private void OnTimeScaleChanged(float ts)
	{
		StartLerp(0.05f);
	}

	public JobHandle BuildRenderState(bool sychronous, ref SimState from, ref TempState tempState, ref DisplayState display, ref RenderState to, ref WorldData worldData, ref StaticState staticState, JobHandle dependency)
	{
		to.Ticks = from.PlanetState.Ticks;
		to.Position = from.PlanetState.Position;
		to.Rotation = math.degrees(from.PlanetState.Rotation);

		MeshOverlayData meshOverlay;
		bool useMeshOverlay = GetMeshOverlayData(ActiveMeshOverlay, ref from, ref tempState, ref display, ref staticState, out meshOverlay);

		WindOverlayData windOverlayData;
		bool useWindOverlay = GetWindOverlayData(ActiveWindOverlay, ref from, ref tempState, ref display, ref staticState, ref worldData, out windOverlayData);

		var buildRenderStateJobHandle = _renderJobHelper.Schedule(
			sychronous, 64,
			new BuildRenderStateCellJob()
			{
				TerrainColor1 = to.TerrainColor1,
				TerrainColor2 = to.TerrainColor2,
				TerrainElevation = to.TerrainElevation,
				WaterColor = to.WaterColor,
				WaterElevation = to.WaterElevation,
				RWaterCurrent = to.WaterCurrent,
				CloudColor = to.CloudColor,
				CloudHeight = to.CloudHeight,
				CloudElevation = to.CloudElevation,
				VelocityHorizontal = to.VelocityHorizontal,
				VelocityVertical = to.VelocityVertical,
				SurfacePosition = to.SurfacePosition,
				OverlayColor = to.OverlayColor,

				TerrainScale = TerrainScale,
				AtmosphereScale = AtmosphereScale,
				PlanetRadius = staticState.PlanetRadius,
				CloudDropletSizeMin = worldData.rainDropMinSize,
				InverseCloudDropletSizeRange = 1.0f / (worldData.rainDropMaxSize - worldData.rainDropMinSize),
				MeshOverlayActive = useMeshOverlay,
				WindOverlayActive = useWindOverlay,
				WindVelocityMax = windOverlayData.MaxVelocity,
				WindVerticalMax = DisplayVerticalWindSpeedMax,
				WindMaskedByLand = windOverlayData.MaskLand,
				MeshOverlayMin = meshOverlay.Colors.Min,
				MeshOverlayInverseRange = meshOverlay.InverseRange,
				CloudElevationSim = tempState.CloudElevation,
				Icosphere = Sim.Icosphere.Vertices,
				Roughness = from.Roughness,
				Elevation = from.Elevation,
				CloudMass = from.CloudMass,
				CloudDropletMass = from.CloudDropletMass,
				CloudAbsorption = tempState.CloudAbsorptivity,
				IceCoverage = tempState.IceCoverage,
				FloraCoverage = tempState.FloraCoverage,
				WaterCoverage = staticState.GetSliceLayer(tempState.WaterCoverage,Sim.WorldData.SurfaceWaterLayer),
				WaterTemperature = staticState.GetSliceLayer(from.WaterTemperature,Sim.WorldData.SurfaceWaterLayer),
				WaterCurrent = staticState.GetSliceLayer(from.WaterVelocity, Sim.WorldData.SurfaceWaterLayer),
				WaterDepth = staticState.GetSliceLayer(tempState.WaterLayerDepth, Sim.WorldData.BottomWaterLayer),
				LavaMass = from.LavaMass,
				LavaTemperature = from.LavaTemperature,
				SurfaceElevation = tempState.SurfaceElevation,
				GroundWater = from.GroundWater,
				MeshOverlayData = meshOverlay.Values,
				MeshOverlayColors = meshOverlay.Colors.ColorValuePairs,
				WindOverlayData = windOverlayData.Values,
				GroundWaterMax = worldData.GroundWaterMax,
				SoilFertilityMax = DisplaySoilFertilityMax,
				LavaCrystalizationTemperature = worldData.LavaCrystalizationTemperature,
				LavaTemperatureRangeInverse = 1.0f / DisplayLavaTemperatureMax,
				DustCoverage = staticState.GetSliceAir(display.DustMass),
				DustMaxInverse = 1.0f / DisplayDustMax,
				LavaToRockMassAdjustment = worldData.LavaToRockMassAdjustment,
				WindColors = _normalizedBlueBlackRed,
				Positions = staticState.SphericalPosition,
				DisplayFloraWeight = DisplayFloraWeight,
				DisplaySandWeight = DisplaySandWeight,
				DisplaySoilWeight = DisplaySoilWeight,
				DisplayCloudHeight = CloudHeight,
				DisplayCloudMin = DisplayCloudMin,
				DisplayCloudPower = DisplayCloudPower,
				PlanktonMax = DisplayPlanktonMax,
				PlanktonPower = DisplayPlanktonPower,
				PlanktonLevels = DisplayPlanktonLevels,
				IceLevels = DisplayIceLevels,
				IceCoveragePower = DisplayIcePower,
				WaterTemperatureMax = DisplayWaterTemperatureMax,
				WaterTemperatureLevels = DisplayWaterTemperatureLevels,
				AirLayers = worldData.AirLayers - 2,
				Count = staticState.Count
			}, dependency);

		return buildRenderStateJobHandle;
	}

	public void LerpMesh(ref RenderState lastState, ref RenderState nextState, ref RenderState state)
	{

		float t = Mathf.Clamp01(1.0f - _tickLerpTime / _tickLerpTimeTotal);
		if (!LerpStates)
		{
			t = 1;
		}
		state.Ticks = (nextState.Ticks - lastState.Ticks) * t + lastState.Ticks;
		state.Rotation = new Vector3(Mathf.LerpAngle(lastState.Rotation.x, nextState.Rotation.x, t), Mathf.LerpAngle(lastState.Rotation.y, nextState.Rotation.y, t), Mathf.LerpAngle(lastState.Rotation.z, nextState.Rotation.z, t));
		state.Position = math.lerp(lastState.Position, nextState.Position, t);


		NativeList<JobHandle> dependencies = new NativeList<JobHandle>(Allocator.TempJob);
		dependencies.Add((new LerpJobfloat { Progress = t, Out = state.TerrainElevation, Start = lastState.TerrainElevation, End = nextState.TerrainElevation }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat4 { Progress = t, Out = state.TerrainColor1, Start = lastState.TerrainColor1, End = nextState.TerrainColor1 }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat4 { Progress = t, Out = state.TerrainColor2, Start = lastState.TerrainColor2, End = nextState.TerrainColor2 }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat { Progress = t, Out = state.WaterElevation, Start = lastState.WaterElevation, End = nextState.WaterElevation }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat4 { Progress = t, Out = state.WaterColor, Start = lastState.WaterColor, End = nextState.WaterColor }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat3 { Progress = t, Out = state.WaterCurrent, Start = lastState.WaterCurrent, End = nextState.WaterCurrent }).Schedule(Sim.CellCount, _batchCount));
		dependencies.Add((new LerpJobfloat3 { Progress = t, Out = state.SurfacePosition, Start = lastState.SurfacePosition, End = nextState.SurfacePosition }).Schedule(Sim.CellCount, _batchCount));

		if (true /* cloudsVisible*/)
		{
			dependencies.Add((new LerpJobfloat { Progress = t, Out = state.CloudElevation, Start = lastState.CloudElevation, End = nextState.CloudElevation }).Schedule(Sim.CellCount, _batchCount));
			dependencies.Add((new LerpJobfloat { Progress = t, Out = state.CloudHeight, Start = lastState.CloudHeight, End = nextState.CloudHeight }).Schedule(Sim.CellCount, _batchCount));
			dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.CloudColor, Start = lastState.CloudColor, End = nextState.CloudColor }).Schedule(Sim.CellCount, _batchCount));
		}
	//	if (ActiveWindOverlay != WindOverlay.None)
		{
			dependencies.Add((new LerpJobfloat3 { Progress = t, Out = state.VelocityHorizontal, Start = lastState.VelocityHorizontal, End = nextState.VelocityHorizontal }).Schedule(Sim.CellCount, _batchCount));
			dependencies.Add((new LerpJobfloat { Progress = t, Out = state.VelocityVertical, Start = lastState.VelocityVertical, End = nextState.VelocityVertical }).Schedule(Sim.CellCount, _batchCount));
		}
		if (ActiveMeshOverlay != MeshOverlay.None)
		{
			dependencies.Add((new LerpJobColor32 { Progress = t, Out = state.OverlayColor, Start = lastState.OverlayColor, End = nextState.OverlayColor }).Schedule(Sim.CellCount, _batchCount));
		}

		var getVertsHandle = _renderTerrainHelper.Schedule(
			JobType.Schedule, 64,
			new BuildTerrainVertsJob()
			{
				VOverlayPosition = _overlayVerticesArray,
				VOverlayUVs = _overlayUVArray,
				VTerrainPosition = _terrainVerticesArray,
				VTerrainColor1 = _terrainColorsArray1,
				VTerrainColor2 = _terrainColorsArray2,
				VWaterPosition = _waterVerticesArray,
				VWaterNormal = _waterNormalsArray,
				VWaterColor = _waterColorsArray,
				VWaterCurrent = _waterCurrentArray,
				VOverlayColor = _overlayColorsArray,

				TerrainElevation = state.TerrainElevation,
				TerrainColor1 = state.TerrainColor1,
				TerrainColor2 = state.TerrainColor2,
				WaterElevation = state.WaterElevation,
				WaterColor = state.WaterColor,
				WaterCurrent = state.WaterCurrent,
				OverlayColor = state.OverlayColor,
				StandardVerts = _terrainVerts,
			}, JobHandle.CombineDependencies(dependencies));

		getVertsHandle = JobHandle.CombineDependencies(getVertsHandle, _renderCloudHelper.Schedule(
			JobType.Schedule, 64,
			new BuildCloudVertsJob()
			{
				VCloudPosition = _cloudVerticesArray,
				VCloudColor = _cloudColorsArray,
				VCloudNormal = _cloudNormalsArray,

				CloudElevation = state.CloudElevation,
				CloudHeight = state.CloudHeight,
				CloudColor = state.CloudColor,
				StandardVerts = _cloudVerts,
			}, JobHandle.CombineDependencies(dependencies)));

		getVertsHandle.Complete();
		dependencies.Dispose();

		_terrainMesh.SetVertices(_terrainVerticesArray);
		_overlayMesh.SetVertices(_overlayVerticesArray);
		_waterMeshFront.SetVertices(_waterVerticesArray);
		_waterMeshFront.SetNormals(_waterNormalsArray);
		_waterMeshBack.SetVertices(_waterVerticesArray);
		_waterMeshBack.SetNormals(_waterNormalsArray);

		_cloudMeshFront.SetVertices(_cloudVerticesArray);
		_cloudMeshFront.SetColors(_cloudColorsArray);
		_cloudMeshFront.SetNormals(_cloudNormalsArray);

		_cloudMeshBack.SetVertices(_cloudVerticesArray);
		_cloudMeshBack.SetColors(_cloudColorsArray);
		_cloudMeshBack.SetNormals(_cloudNormalsArray);

		_terrainMesh.SetUVs(1, _terrainColorsArray1);
		_terrainMesh.SetUVs(2, _terrainColorsArray2);
		_waterMeshFront.SetUVs(1, _waterColorsArray);
		_waterMeshFront.SetUVs(2, _waterCurrentArray);
		_waterMeshBack.SetUVs(1, _waterColorsArray);
		_waterMeshBack.SetUVs(2, _waterCurrentArray);

		if (ActiveMeshOverlay != MeshOverlay.None)
		{
			_overlayMesh.SetColors(_overlayColorsArray);
			_overlayMesh.SetUVs(1, _overlayUVArray);
		}
		if (!_indicesInitialized)
		{
			_terrainMesh.SetTriangles(_indicesTerrain, 0);
			_overlayMesh.SetTriangles(_indicesTerrain, 0);
			_waterMeshFront.SetTriangles(_indicesTerrain, 0);
			_waterMeshBack.SetTriangles(_indicesTerrainBack, 0);
			_cloudMeshFront.SetTriangles(_indicesCloudFront, 0);
			_cloudMeshBack.SetTriangles(_indicesCloudBack, 0);
			_indicesInitialized = true;
		}

		//_terrainMesh.RecalculateBounds();
		//_waterMesh.RecalculateBounds();
		//_cloudMesh.RecalculateBounds();

		_terrainMesh.RecalculateNormals();
		//	_waterMesh.RecalculateNormals();
		//	_cloudMesh.RecalculateNormals();
		//_terrainMesh.RecalculateTangents();
		//_waterMesh.RecalculateTangents();

	}

	public void OnWaterDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		_waterObjectFront.SetActive(toggle.isOn);
		_waterObjectBack.SetActive(toggle.isOn);
	}
	public void OnCloudDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		foreach (var c in _cloudObjectFront) { c.SetActive(toggle.isOn); }
		foreach (var c in _cloudObjectBack) { c.SetActive(toggle.isOn); }
	}
	public void OnFoliageDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		Foliage.FoliageParent.SetActive(toggle.isOn);
	}
	public void OnStarsDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		_skyboxActive = toggle.isOn;
	}


	public int GetClosestVert(int triangleIndex, int vIndex)
	{
		return _indicesTerrain[triangleIndex * 3 + vIndex] / VertsPerCell;
	}

	public enum LegendType {
		None,
		Temperature,
		PPM,
		Percent,
		Pressure,
		Volume,
		Mass,
		Watts,
	}

	[Serializable]
	public struct MeshOverlayColors {
		public string Title;
		public float Min;
		public float Max;
		public LegendType LegendType;
		public int DecimalPlaces;
		public NativeArray<CVP> ColorValuePairs;
	}

	public struct MeshOverlayData {
		public MeshOverlayData(MeshOverlayColors colors, NativeSlice<float> values)
		{
			Values = values;
			Colors = colors;
			InverseRange = 1.0f / (colors.Max - colors.Min);
		}
		public MeshOverlayColors Colors;
		public NativeSlice<float> Values { get; private set; }
		public float InverseRange { get; private set; }
	}

	public struct WindOverlayData {
		public WindOverlayData(float maxVelocity, bool maskLand, NativeSlice<float3> values)
		{
			MaxVelocity = maxVelocity;
			MaskLand = maskLand;
			Values = values;
		}
		public float MaxVelocity { get; private set; }
		public bool MaskLand { get; private set; }
		public NativeSlice<float3> Values { get; private set; }
	}



	#region private functions

	private void InitVerts(Icosphere icosphere)
	{
		_terrainVertices = new Vector3[Sim.CellCount * VertsPerCell];
		_waterVertices = new Vector3[Sim.CellCount * VertsPerCell];
		_waterNormals = new Vector3[Sim.CellCount * VertsPerCell];
		_waterColors = new Color32[Sim.CellCount * VertsPerCell];

		_cloudVertices = new Vector3[Sim.CellCount * VertsPerCloud];
		_cloudColors = new Color32[Sim.CellCount * VertsPerCloud];
		_cloudNormals = new Vector3[Sim.CellCount * VertsPerCloud];

		Unity.Mathematics.Random random = new Unity.Mathematics.Random(1);

		_terrainVerts = new NativeArray<float3>(Sim.CellCount * VertsPerCell, Allocator.Persistent);
		_cloudVerts = new NativeArray<float3>(Sim.CellCount * VertsPerCloud, Allocator.Persistent);
		List<int> indices = new List<int>();
		List<int> indicesCloudFront = new List<int>();
		List<int> indicesCloudBack = new List<int>();
		List<int> indicesTerrain = new List<int>();
		List<int> indicesTerrainBack = new List<int>();
		const float maxSlope = 3;
		for (int i=0;i< icosphere.Vertices.Length;i++)
		{
			float3 pos = icosphere.Vertices[i];
			_terrainVerts[i * VertsPerCell] = pos;
			_cloudVerts[i * VertsPerCloud] = pos;
			int neighborCount = (icosphere.Neighbors[(i + 1) * MaxNeighbors - 1] >= 0) ? MaxNeighbors : (MaxNeighbors - 1);
			for (int j=0;j< neighborCount; j++)
			{
				int neighborIndex1 = icosphere.Neighbors[i * MaxNeighbors + j];
				int neighborIndex2 = icosphere.Neighbors[i * MaxNeighbors + (j + 1) % neighborCount];

				{
					float slope = random.NextFloat() * SlopeAmountTerrain;
					float3 slopePoint = (icosphere.Vertices[neighborIndex1] + icosphere.Vertices[neighborIndex2] + pos * (1 + slope)) / (maxSlope + slope);
					float slopePointLength = math.length(slopePoint);
					float3 extendedSlopePoint = slopePoint / (slopePointLength * slopePointLength);

					_terrainVerts[i * VertsPerCell + 1 + j] = extendedSlopePoint; // surface
					_terrainVerts[i * VertsPerCell + 1 + j + MaxNeighbors] = extendedSlopePoint; // wall
					_terrainVerts[i * VertsPerCell + 1 + j + MaxNeighbors * 2] = extendedSlopePoint; // wall
					_terrainVerts[i * VertsPerCell + 1 + j + MaxNeighbors * 3] = extendedSlopePoint; // corner

				}

				{
					float slope = random.NextFloat() * SlopeAmountCloud;
					float3 slopePoint = (icosphere.Vertices[neighborIndex1] + icosphere.Vertices[neighborIndex2] + pos * (1 + slope)) / (3 + slope);
					float slopePointLength = math.length(slopePoint);
					float3 extendedSlopePoint = slopePoint / (slopePointLength * slopePointLength);
					_cloudVerts[i * VertsPerCell + 1 + j] = extendedSlopePoint; // surface
					_cloudVerts[i * VertsPerCell + 1 + j + MaxNeighbors] = extendedSlopePoint; // wall
					_cloudVerts[i * VertsPerCell + 1 + j + MaxNeighbors * 2] = extendedSlopePoint; // wall
					_cloudVerts[i * VertsPerCell + 1 + j + MaxNeighbors * 3] = extendedSlopePoint; // corner
				}

				indicesTerrain.Add(i * VertsPerCell);
				indicesTerrain.Add(i * VertsPerCell + 1 + ((j + 1) % neighborCount));
				indicesTerrain.Add(i * VertsPerCell + 1 + j);

				indicesTerrainBack.Add(i * VertsPerCell);
				indicesTerrainBack.Add(i * VertsPerCell + 1 + j);
				indicesTerrainBack.Add(i * VertsPerCell + 1 + ((j + 1) % neighborCount));

				indicesCloudFront.Add(i * VertsPerCell);
				indicesCloudFront.Add(i * VertsPerCell + 1 + ((j + 1) % neighborCount));
				indicesCloudFront.Add(i * VertsPerCell + 1 + j);

				indicesCloudBack.Add(i * VertsPerCell);
				indicesCloudBack.Add(i * VertsPerCell + 1 + j);
				indicesCloudBack.Add(i * VertsPerCell + 1 + ((j + 1) % neighborCount));

				indices.Add(i * VertsPerCloud);
				indices.Add(i * VertsPerCloud + 1 + ((j + 1) % neighborCount));
				indices.Add(i * VertsPerCloud + 1 + j);

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

							indicesTerrainBack.Add(i * VertsPerCell + 1 + 2 * MaxNeighbors + ((j - 1 + neighborCount) % neighborCount));
							indicesTerrainBack.Add(neighborIndex1 * VertsPerCell + 1 + MaxNeighbors + k);
							indicesTerrainBack.Add(i * VertsPerCell + 1 + MaxNeighbors + j);

							indicesCloudFront.Add(i * VertsPerCell + 1 + 2 * MaxNeighbors + ((j - 1 + neighborCount) % neighborCount));
							indicesCloudFront.Add(i * VertsPerCell + 1 + MaxNeighbors + j);
							indicesCloudFront.Add(neighborIndex1 * VertsPerCell + 1 + MaxNeighbors + k);

							indicesCloudBack.Add(i * VertsPerCell + 1 + 2 * MaxNeighbors + ((j - 1 + neighborCount) % neighborCount));
							indicesCloudBack.Add(neighborIndex1 * VertsPerCell + 1 + MaxNeighbors + k);
							indicesCloudBack.Add(i * VertsPerCell + 1 + MaxNeighbors + j);

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

							indicesTerrainBack.Add(i * VertsPerCell + 1 + 3 * MaxNeighbors + j);
							indicesTerrainBack.Add(neighborIndex1 * VertsPerCell + 1 + 3 * MaxNeighbors + neighbor1);
							indicesTerrainBack.Add(neighborIndex2 * VertsPerCell + 1 + 3 * MaxNeighbors + k);

							indicesCloudFront.Add(i * VertsPerCell + 1 + 3 * MaxNeighbors + j);
							indicesCloudFront.Add(neighborIndex2 * VertsPerCell + 1 + 3 * MaxNeighbors + k);
							indicesCloudFront.Add(neighborIndex1 * VertsPerCell + 1 + 3 * MaxNeighbors + neighbor1);

							indicesCloudBack.Add(i * VertsPerCell + 1 + 3 * MaxNeighbors + j);
							indicesCloudBack.Add(neighborIndex1 * VertsPerCell + 1 + 3 * MaxNeighbors + neighbor1);
							indicesCloudBack.Add(neighborIndex2 * VertsPerCell + 1 + 3 * MaxNeighbors + k);

							break;
						}
					}
				}

			}
		}
		_indices = indices.ToArray();
		_indicesTerrain = indicesTerrain.ToArray();
		_indicesTerrainBack = indicesTerrainBack.ToArray();
		_indicesCloudFront = indicesCloudFront.ToArray();
		_indicesCloudBack = indicesCloudBack.ToArray();
	}

	private bool GetMeshOverlayData(MeshOverlay activeOverlay, ref SimState simState, ref TempState dependentState, ref DisplayState display, ref StaticState staticState, out MeshOverlayData overlay)
	{
		float ticksPerYear = Sim.WorldData.TicksPerSecond * 60 * 60 * 24 * 365;
		MeshOverlayColors colors;
		if (OverlayColors.TryGetValue(activeOverlay, out colors))
		{
			switch (activeOverlay)
			{
				case MeshOverlay.AbsoluteHumidity:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(dependentState.AirHumidityAbsolute, ActiveMeshLayerAir));
					return true;
				case MeshOverlay.RelativeHumidity:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(dependentState.AirHumidityRelative, ActiveMeshLayerAir));
					return true;
				case MeshOverlay.TemperatureSurface:
					overlay = new MeshOverlayData(colors, dependentState.SurfaceAirTemperatureAbsolute);
					return true;
				case MeshOverlay.PotentialTemperature:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(simState.AirTemperaturePotential, ActiveMeshLayerAir));
					return true;
				case MeshOverlay.Pressure:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(display.Pressure, ActiveMeshLayerAir));
					return true;
				case MeshOverlay.WaterTemperature:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(simState.WaterTemperature, ActiveMeshLayerWater));
					return true;
				case MeshOverlay.WaterCarbonDioxide:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(display.WaterCarbonDioxidePercent, ActiveMeshLayerWater));
					return true;
				case MeshOverlay.Salinity:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(display.Salinity, ActiveMeshLayerWater));
					return true;
				case MeshOverlay.GroundTemperature:
					overlay = new MeshOverlayData(colors, simState.GroundTemperature);
					return true;
				case MeshOverlay.CarbonDioxide:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(display.CarbonDioxidePercent, ActiveMeshLayerAir));
					return true;
				case MeshOverlay.HeatAbsorbed:
					overlay = new MeshOverlayData(colors, display.SolarRadiationAbsorbedSurface);
					return true;
				case MeshOverlay.Rainfall:
					overlay = new MeshOverlayData(colors, display.Rainfall);
					return true;
				case MeshOverlay.Evaporation:
					overlay = new MeshOverlayData(colors, display.Evaporation);
					return true;
				case MeshOverlay.GroundWater:
					overlay = new MeshOverlayData(colors, simState.GroundWater);
					return true;
				case MeshOverlay.GroundWaterTemperature:
					overlay = new MeshOverlayData(colors, simState.GroundWaterTemperature);
					return true;
				case MeshOverlay.CrustDepth:
					overlay = new MeshOverlayData(colors, simState.CrustDepth);
					return true;
				case MeshOverlay.MagmaMass:
					overlay = new MeshOverlayData(colors, simState.MagmaMass);
					return true;
				case MeshOverlay.DivergenceAir:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(display.DivergenceAir, ActiveMeshLayerAir));
					return true;
				case MeshOverlay.DivergenceWater:
					overlay = new MeshOverlayData(colors, staticState.GetSliceLayer(display.DivergenceWater, ActiveMeshLayerWater));
					return true;
				case MeshOverlay.PlateTectonics:
					overlay = new MeshOverlayData(colors, display.Plate);
					return true;
			}
		}
		overlay = new MeshOverlayData(OverlayColors[MeshOverlay.Evaporation], display.Evaporation);
		return false;

	}
	private bool GetWindOverlayData(WindOverlay activeOverlay, ref SimState simState, ref TempState dependentState, ref DisplayState displayState, ref StaticState staticState, ref WorldData worldData, out WindOverlayData overlay)
	{
		switch (activeOverlay)
		{
			case WindOverlay.Wind:
				overlay = new WindOverlayData(DisplayWindSpeedLowerAirMax, false, staticState.GetSliceLayer(simState.AirVelocity,ActiveMeshLayerAir));
				return true;
			case WindOverlay.WindCloud:
				overlay = new WindOverlayData(DisplayWindSpeedUpperAirMax, false, dependentState.CloudVelocity);
				return true;
			case WindOverlay.PGF:
				overlay = new WindOverlayData(DisplayPressureGradientForceMax, false, staticState.GetSliceLayer(displayState.PressureGradientForce,ActiveMeshLayerAir));
				return true;
			case WindOverlay.Current:
				overlay = new WindOverlayData(DisplayWindSpeedSurfaceWaterMax, true, staticState.GetSliceLayer(simState.WaterVelocity,ActiveMeshLayerWater));
				return true;
		}
		overlay = new WindOverlayData(DisplayWindSpeedDeepWaterMax, false, staticState.GetSliceLayer(simState.WaterVelocity, worldData.BottomWaterLayer));
		return false;
	}


	private void UpdateSkybox()
	{
		if (Sim.TimeScale == 0 || _skyboxActive)
		{
			if (_skyboxExposureDest != 1)
			{
				_skyboxExposureStart = _skyboxExposure;
				_skyboxExposureTime = 1.0f;
			}
			_skyboxExposureDest = 1.0f;

		}
		else
		{
			if (_skyboxExposureDest != 0)
			{
				_skyboxExposureStart = _skyboxExposure;
				_skyboxExposureTime = 1.0f;
			}
			_skyboxExposureDest = 0;
		}
		_skyboxExposureTime = math.max(0, _skyboxExposureTime - Time.deltaTime);
		_skyboxExposure = math.lerp(_skyboxExposureStart, _skyboxExposureDest, 1.0f - math.sin(_skyboxExposureTime * math.PI / 2));
		//RenderSettings.skybox.Set
		RenderSettings.skybox.SetFloat("_Exposure", _skyboxExposure);
	}

	private void OnSetActiveCell(int index)
	{
		_selectionCircle.SetActive(index >= 0);
		if (index >= 0)
		{
			//			var p = Sim.Icosphere.Vertices[index];
			var pos = _renderStates[_curRenderState].SurfacePosition[index];
			_selectionCircle.transform.localPosition = pos;
			_selectionCircle.transform.localRotation = Quaternion.LookRotation(-pos);
		}

		for (int i = 0; i < _selectionCircleNeighbors.Count; i++)
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

	private void OnSetActiveTool(GameTool tool)
	{
		SetActiveMeshOverlay(tool != null ? tool.Overlay : MeshOverlay.None);
	}

	public void SetActiveMeshOverlay(MeshOverlay o)
	{
		ActiveMeshOverlay = o;

		Sim.SimSettings.CollectOverlay = ActiveMeshOverlay != MeshOverlay.None;
		_overlayObject.SetActive(o != MeshOverlay.None);

		StartLerp(0.1f);
	}
	public void SetActiveWindOverlay(WindOverlay o)
	{
		ActiveWindOverlay = o;

		for (int i = 0; i < _windArrows.Length; i++)
		{
			var wind = _renderStates[_curRenderState].VelocityHorizontal[i];
			bool visible = ActiveWindOverlay != WindOverlay.None && !(wind.x == 0 && wind.y == 0);
			_windArrows[i].SetActive(visible);
		}
		StartLerp(0.1f);
	}
	public void SetActiveLayerAir(int l)
	{
		if (l != ActiveMeshLayerAir)
		{
			ActiveMeshLayerAir = l;

			StartLerp(0.1f);
		}
	}
	public void SetActiveLayerWater(int l)
	{
		if (l != ActiveMeshLayerWater)
		{
			ActiveMeshLayerWater = l;

			StartLerp(0.1f);
		}
	}

	#endregion
}
