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
using UnityEngine.Profiling;
using UnityEngine.VFX;

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
		FloraWater,
		CrustDepth,
		MagmaMass,
		DivergenceAir,
		DivergenceWater
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


	[Header("Display")]
	public float SlopeAmountTerrain = 3;
	public float SlopeAmountCloud = 10;
	public float TerrainScale = 100f;
	public float CloudHeight = 0.1f;
	public float DisplayCloudMin = 0.1f;
	public float AtmosphereScale = 1000f;
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
	public GameObject Planet;
	public GameObject Sun;
	public GameObject Moon;
	public GameObject SunLight;
	public VisualEffect WindEffect;
	public Material TerrainMaterial;
	public Material WaterMaterial;
	public Material CloudMaterialFront, CloudMaterialBack;
	public Material OverlayMaterial;
	public GameObject SelectionCirclePrefab;
	public GameObject WindArrowPrefab;
	public FoliageManager Foliage;

	[Header("Debug")]
	public bool RunSynchronously;

	[HideInInspector] public MeshOverlay ActiveMeshOverlay { get; private set; }
	[HideInInspector] public WindOverlay ActiveWindOverlay;
	[HideInInspector] public int ActiveMeshLayerWater = 1;
	[HideInInspector] public int ActiveMeshLayerAir = 1;
	[HideInInspector] public DisplayState DisplayState;


	private JobHelper _renderJobHelper;
	private JobHelper _renderTerrainHelper;
	private JobHelper _renderCloudHelper;

	private RenderState[] _renderStates;
	private int _curRenderState;
	private int _lastRenderState;
	private int _nextRenderState;
	private const int _renderStateCount = 3;

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
	private int[] _indices;
	private int[] _indicesTerrain;
	private int[] _indicesCloud;

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

	public void Start()
	{
		Sim.OnTick += OnSimTick;
		GameplayManager.OnSetActiveCell += OnSetActiveCell;

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
		_normalizedBlueBlackRed = new NativeArray<CVP>(new CVP[] {
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
		_terrainMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;

		_waterObject = new GameObject("Water Mesh");
		_waterObject.transform.SetParent(Planet.transform, false);
		var waterFilter = _waterObject.AddComponent<MeshFilter>();
		var waterSurfaceRenderer = _waterObject.AddComponent<MeshRenderer>();
		var waterCollider = _waterObject.AddComponent<MeshCollider>();
		waterSurfaceRenderer.material = WaterMaterial;
		waterFilter.mesh = waterCollider.sharedMesh = _waterMesh;
		waterSurfaceRenderer.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
		_waterMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;

		_cloudObject = new GameObject("Cloud Mesh");
		_cloudObject.transform.SetParent(Planet.transform, false);
		var cloudFilter = _cloudObject.AddComponent<MeshFilter>();
		var cloudSurfaceRenderer = _cloudObject.AddComponent<MeshRenderer>();
	//	cloudSurfaceRenderer.materials = new Material[] { CloudMaterialBack, CloudMaterialFront };
		cloudSurfaceRenderer.material = CloudMaterialFront;
		cloudFilter.mesh = _cloudMesh;
		cloudSurfaceRenderer.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
		_cloudMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;

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
		_waterMesh.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));
		_cloudMesh.bounds = new Bounds(Planet.transform.position, new Vector3(boundsSize, boundsSize, boundsSize));

		WindEffect.SetMesh("TerrainMesh", _terrainMesh);
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

		_terrainVerts.Dispose();
		_cloudVerts.Dispose();
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
		WaterMaterial.SetFloat("Vector1_2C57E502", _renderStates[_curRenderState].Ticks); // sim time

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
		if (Sim.SimSettings.CollectOverlay)
		{
			DisplayState = new DisplayState();
			DisplayState.Init(Sim.StaticState.Count, ref Sim.WorldData);
			displayJob = DisplayState.Update(ref DisplayState, ref lastDisplay, ref Sim.WorldData, ref Sim.LastTempState, ref Sim.LastSimState, ref Sim.StaticState, ref Sim.SimSettings);
		}


		var buildRenderStateJob = BuildRenderState(RunSynchronously, ref Sim.LastSimState, ref Sim.LastTempState, ref DisplayState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState, displayJob);
		var foliageJob = Foliage.Tick(ref Sim.LastTempState, displayJob);

		displayJob = JobHandle.CombineDependencies(displayJob, buildRenderStateJob, foliageJob);
		displayJob.Complete();
		if (Sim.SimSettings.CollectOverlay)
		{
			lastDisplay.Dispose();
		}

	}

	private void OnSimTick()
	{

		float lerpTime;
		if (Sim.TimeScale == 0)
		{
			lerpTime = Sim.TimeTillTick;
		}
		else
		{
			lerpTime = 0.1f + Sim.TimeTillTick / Sim.TimeScale;
		}
		StartLerp(lerpTime);


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
				MeshOverlayMin = meshOverlay.Min,
				MeshOverlayInverseRange = meshOverlay.InverseRange,
				CloudElevationSim = tempState.CloudElevation,
				Icosphere = Sim.Icosphere.Vertices,
				SoilFertility = from.GroundCarbon,
				Roughness = from.Roughness,
				Elevation = from.Elevation,
				CloudMass = from.CloudMass,
				CloudDropletMass = from.CloudDropletMass,
				CloudAbsorption = tempState.CloudAbsorptivity,
				IceCoverage = tempState.IceCoverage,
				FloraCoverage = tempState.FloraCoverage,
				WaterCoverage = staticState.GetSliceLayer(tempState.WaterCoverage,Sim.WorldData.SurfaceWaterLayer),
				WaterTemperature = staticState.GetSliceLayer(from.WaterTemperature,Sim.WorldData.SurfaceWaterLayer),
				PlanktonMass = staticState.GetSliceLayer(from.PlanktonMass,Sim.WorldData.SurfaceWaterLayer),
				WaterCurrent = staticState.GetSliceLayer(from.WaterVelocity, Sim.WorldData.SurfaceWaterLayer),
				WaterDepth = staticState.GetSliceLayer(tempState.WaterLayerDepth, Sim.WorldData.BottomWaterLayer),
				LavaMass = from.LavaMass,
				LavaTemperature = from.LavaTemperature,
				SurfaceElevation = tempState.SurfaceElevation,
				GroundWater = from.GroundWater,
				MeshOverlayData = meshOverlay.Values,
				MeshOverlayColors = meshOverlay.ColorValuePairs,
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
				DisplayCloudRangeInverse = 1.0f / (1.0f - DisplayCloudMin),
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
		if (ActiveWindOverlay != WindOverlay.None)
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
		_waterMesh.SetVertices(_waterVerticesArray);
		_waterMesh.SetNormals(_waterNormalsArray);

		_cloudMesh.SetVertices(_cloudVerticesArray);
		_cloudMesh.SetColors(_cloudColorsArray);
		_cloudMesh.SetNormals(_cloudNormalsArray);

		if (ActiveMeshOverlay == MeshOverlay.None)
		{
			_terrainMesh.SetUVs(1, _terrainColorsArray1);
			_terrainMesh.SetUVs(2, _terrainColorsArray2);
			_waterMesh.SetUVs(1, _waterColorsArray);
			_waterMesh.SetUVs(2, _waterCurrentArray);
		}
		else
		{
			_terrainMesh.SetColors(_overlayColorsArray);
			_waterMesh.SetColors(_overlayColorsArray);
		}


		if (!_indicesInitialized)
		{
			_terrainMesh.SetTriangles(_indicesTerrain, 0);
			_waterMesh.SetTriangles(_indicesTerrain, 0);
			_cloudMesh.SetTriangles(_indicesCloud, 0);
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
		_waterObject.SetActive(toggle.isOn);
	}
	public void OnCloudDisplayToggled(UnityEngine.UI.Toggle toggle)
	{
		_cloudObject.SetActive(toggle.isOn);
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

	public struct MeshOverlayData {
		public MeshOverlayData(float min, float max, NativeArray<CVP> colors, NativeSlice<float> values)
		{
			Values = values;
			Min = min;
			Max = max;
			ColorValuePairs = colors;
			InverseRange = 1.0f / (Max - Min);
		}
		public float Min { get; private set; }
		public float Max { get; private set; }
		public NativeSlice<float> Values { get; private set; }
		public NativeArray<CVP> ColorValuePairs { get; private set; }
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
		List<int> indicesCloud = new List<int>();
		List<int> indicesTerrain = new List<int>();
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
					float3 slopePoint = (icosphere.Vertices[neighborIndex1] + icosphere.Vertices[neighborIndex2] + pos * (1 + slope)) / (3 + slope);
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

				indicesCloud.Add(i * VertsPerCell);
				indicesCloud.Add(i * VertsPerCell + 1 + ((j + 1) % neighborCount));
				indicesCloud.Add(i * VertsPerCell + 1 + j);

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

							indicesCloud.Add(i * VertsPerCell + 1 + 2 * MaxNeighbors + ((j - 1 + neighborCount) % neighborCount));
							indicesCloud.Add(i * VertsPerCell + 1 + MaxNeighbors + j);
							indicesCloud.Add(neighborIndex1 * VertsPerCell + 1 + MaxNeighbors + k);

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

							indicesCloud.Add(i * VertsPerCell + 1 + 3 * MaxNeighbors + j);
							indicesCloud.Add(neighborIndex2 * VertsPerCell + 1 + 3 * MaxNeighbors + k);
							indicesCloud.Add(neighborIndex1 * VertsPerCell + 1 + 3 * MaxNeighbors + neighbor1);

							break;
						}
					}
				}

			}
		}
		_indices = indices.ToArray();
		_indicesTerrain = indicesTerrain.ToArray();
		_indicesCloud = indicesCloud.ToArray();
	}

	private bool GetMeshOverlayData(MeshOverlay activeOverlay, ref SimState simState, ref TempState dependentState, ref DisplayState display, ref StaticState staticState, out MeshOverlayData overlay)
	{
		float ticksPerYear = Sim.WorldData.TicksPerSecond * 60 * 60 * 24 * 365;
		switch (activeOverlay)
		{
			case MeshOverlay.AbsoluteHumidity:
				overlay = new MeshOverlayData(0, DisplayAbsoluteHumidityMax, _normalizedRainbow, staticState.GetSliceLayer(dependentState.AirHumidityAbsolute,ActiveMeshLayerAir));
				return true;
			case MeshOverlay.RelativeHumidity:
				overlay = new MeshOverlayData(0, 1.0f, _normalizedRainbow, staticState.GetSliceLayer(dependentState.AirHumidityRelative,ActiveMeshLayerAir));
				return true;
			case MeshOverlay.TemperatureSurface:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, dependentState.SurfaceAirTemperatureAbsolute);
				return true;
			case MeshOverlay.PotentialTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, staticState.GetSliceLayer(simState.AirTemperaturePotential,ActiveMeshLayerAir));
				return true;
			case MeshOverlay.Pressure:
				overlay = new MeshOverlayData(DisplayAirPressureMin, DisplayAirPressureMax, _normalizedRainbow, staticState.GetSliceLayer(display.Pressure,ActiveMeshLayerAir));
				return true;
			case MeshOverlay.WaterTemperature:
				overlay = new MeshOverlayData(WorldData.FreezingTemperature, DisplayTemperatureMax, _normalizedRainbow, staticState.GetSliceLayer(simState.WaterTemperature,ActiveMeshLayerWater));
				return true;
			case MeshOverlay.WaterCarbonDioxide:
				overlay = new MeshOverlayData(0, DisplayWaterCarbonMax, _normalizedRainbow, staticState.GetSliceLayer(display.WaterCarbonDioxidePercent,ActiveMeshLayerWater));
				return true;
			case MeshOverlay.Salinity:
				overlay = new MeshOverlayData(DisplaySalinityMin, DisplaySalinityMax, _normalizedRainbow, staticState.GetSliceLayer(display.Salinity,ActiveMeshLayerWater));
				return true;
			case MeshOverlay.GroundTemperature:
				overlay = new MeshOverlayData(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.GroundTemperature);
				return true;
			case MeshOverlay.CarbonDioxide:
				overlay = new MeshOverlayData(0, DisplayCarbonDioxideMax, _normalizedRainbow, staticState.GetSliceLayer(display.CarbonDioxidePercent,ActiveMeshLayerAir));
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
			case MeshOverlay.FloraWater:
				overlay = new MeshOverlayData(0, Sim.WorldData.FullCoverageFlora, _normalizedRainbow, simState.FloraWater);
				return true;
			case MeshOverlay.CrustDepth:
				overlay = new MeshOverlayData(DisplayCrustDepthMax, 0, _normalizedRainbow, simState.CrustDepth);
				return true;
			case MeshOverlay.MagmaMass:
				overlay = new MeshOverlayData(0, DisplayMagmaMassMax, _normalizedRainbow, simState.MagmaMass);
				return true;
			case MeshOverlay.DivergenceAir:
				overlay = new MeshOverlayData(-DisplayDivergenceMax, DisplayDivergenceMax, _normalizedBlueBlackRed, staticState.GetSliceLayer(display.DivergenceAir,ActiveMeshLayerAir));
				return true;
			case MeshOverlay.DivergenceWater:
				overlay = new MeshOverlayData(-DisplayDivergenceMax, DisplayDivergenceMax, _normalizedBlueBlackRed, staticState.GetSliceLayer(display.DivergenceWater,ActiveMeshLayerWater));
				return true;
		}
		overlay = new MeshOverlayData(0, DisplayEvaporationMax, _normalizedRainbow, display.Evaporation);
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

	public void SetActiveMeshOverlay(MeshOverlay o)
	{
		ActiveMeshOverlay = o;

		_terrainObject.GetComponent<MeshRenderer>().material = (ActiveMeshOverlay == MeshOverlay.None) ? TerrainMaterial : OverlayMaterial;
		_waterObject.GetComponent<MeshRenderer>().material = (ActiveMeshOverlay == MeshOverlay.None) ? WaterMaterial : OverlayMaterial;
		Sim.SimSettings.CollectOverlay = ActiveMeshOverlay != MeshOverlay.None;

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
