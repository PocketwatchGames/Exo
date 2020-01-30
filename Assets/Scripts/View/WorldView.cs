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
		GroundWater,
		DeepWaterTemperature,
		DeepWaterSalinity,
		ShallowWaterTemperature,
		ShallowWaterSalinity,
		LowerAirTemperature,
		LowerAirPressure,
		UpperAirTemperature,
		UpperAirPressure,
		VerticalWind,
		AbsoluteHumidity,
		RelativeHumidity,
		CloudMass,
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
	public float DisplayAbsoluteHumidityMax = 400;
	public float DisplayGroundWaterMax = 100000;
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

	static Color32 green = new Color32(0, 220, 30, 255);
	static Color32 grey = new Color32(50, 50, 80, 255);
	static Color32 brown = new Color32(100, 60, 20, 255);
	static Color32 blue = new Color32(0, 0, 255, 255);
	static Color32 white = new Color32(255, 255, 255, 255);
	static Color32 black = new Color32(0, 0, 0, 255);


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
		_curRenderState = 0;
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

		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[0], ref Sim.WorldData, ref Sim.StaticState);
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
		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState);
	}

	public void BuildRenderState(ref SimState from, ref RenderState to, ref WorldData worldData, ref StaticState staticState)
	{
		to.Ticks = from.PlanetState.Ticks;
		to.Position = from.PlanetState.Position;
		to.Rotation = math.degrees(from.PlanetState.Rotation);

		Overlay overlay;
		bool useOverlay = GetOverlay(ActiveMeshOverlay, ref from, out overlay);

		NativeArray<float2> velocity = from.AirVelocity[0];

		var buildRenderStateJob = new BuildRenderStateJob()
		{
			SurfacePosition = to.SurfacePosition,
			VelocityArrow = to.VelocityArrow,
			TerrainColor = to.TerrainColor,
			TerrainNormal = to.TerrainNormal,
			TerrainPosition = to.TerrainPosition,
			CloudColor = to.CloudColor,
			CloudNormal = to.CloudNormal,
			CloudPosition = to.CloudPosition,
			WaterColor =to.WaterColor,
			WaterNormal = to.WaterNormal,
			WaterPosition =to.WaterPosition,

			UseOverlay = useOverlay,
			TerrainScale = TerrainScale,
			PlanetRadius = staticState.PlanetRadius,
			Icosphere = Sim.Icosphere.Vertices,
			Overlay = overlay,
			CloudElevation = from.CloudElevation,
			Velocity = velocity,
			Terrain = from.CellTerrains,
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

		int batchCount = 100;
		JobHandle jobHandle = default(JobHandle);
		JobHandle.CombineDependencies(jobHandle, (new Lerpfloat3Job() { T = t, Start = lastState.SurfacePosition, End = nextState.SurfacePosition, Out = state.SurfacePosition }).Schedule(Sim.CellCount, batchCount));
		JobHandle.CombineDependencies(jobHandle, (new Lerpfloat3Job() { T = t, Start = lastState.TerrainPosition, End = nextState.TerrainPosition, Out = state.TerrainPosition }).Schedule(Sim.CellCount, batchCount));
		JobHandle.CombineDependencies(jobHandle, (new Lerpfloat3Job() { T = t, Start = lastState.TerrainNormal, End = nextState.TerrainNormal, Out = state.TerrainNormal }).Schedule(Sim.CellCount, batchCount));
		JobHandle.CombineDependencies(jobHandle, (new Lerpfloat3Job() { T = t, Start = lastState.WaterPosition, End = nextState.TerrainPosition, Out = state.TerrainPosition }).Schedule(Sim.CellCount, batchCount));
		JobHandle.CombineDependencies(jobHandle, (new Lerpfloat3Job() { T = t, Start = lastState.WaterNormal, End = nextState.WaterNormal, Out = state.WaterNormal }).Schedule(Sim.CellCount, batchCount));
		JobHandle.CombineDependencies(jobHandle, (new LerpColor32Job() { T = t, Start = lastState.TerrainColor, End = nextState.TerrainColor, Out = state.TerrainColor }).Schedule(Sim.CellCount, batchCount));
		JobHandle.CombineDependencies(jobHandle, (new LerpColor32Job() { T = t, Start = lastState.WaterColor, End = nextState.WaterColor, Out = state.WaterColor }).Schedule(Sim.CellCount, batchCount));

		if (true /* cloudsVisible*/)
		{
			JobHandle.CombineDependencies(jobHandle, (new Lerpfloat3Job() { T = t, Start = lastState.CloudPosition, End = nextState.CloudPosition, Out = state.CloudPosition }).Schedule(Sim.CellCount, batchCount));
			JobHandle.CombineDependencies(jobHandle, (new Lerpfloat3Job() { T = t, Start = lastState.CloudNormal, End = nextState.CloudNormal, Out = state.CloudNormal }).Schedule(Sim.CellCount, batchCount));
			JobHandle.CombineDependencies(jobHandle, (new LerpColor32Job() { T = t, Start = lastState.CloudColor, End = nextState.CloudColor, Out = state.CloudColor }).Schedule(Sim.CellCount, batchCount));
		}
		if (ActiveWindOverlay != WindOverlay.None)
		{
			JobHandle.CombineDependencies(jobHandle, (new Lerpfloat2Job() { T = t, Start = lastState.VelocityArrow, End = nextState.VelocityArrow, Out = state.VelocityArrow }).Schedule(Sim.CellCount, batchCount));
		}

		jobHandle.Complete();

		for (int i = 0; i < Sim.CellCount; i++)
		{
			_terrainVertices[i] = state.TerrainPosition[i];
			_terrainNormals[i] = state.TerrainNormal[i];
			_terrainColors[i] = state.TerrainColor[i];
		}

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
		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState);
	}
	public void OnHUDWindChanged(UnityEngine.UI.Dropdown dropdown)
	{
		ActiveWindOverlay = (WorldView.WindOverlay)dropdown.value;
		_tickLerpTime = 0;
		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState);
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
				return GetCellInfoGlobal(ref Sim.ActiveSimState, ref Sim.ActiveSimDependent);
			case CellInfoType.Energy:
				return GetCellInfoEnergy(ref Sim.ActiveSimState, ref Sim.ActiveSimDependent);
			case CellInfoType.Cell:
				return GetCellInfoCell(ref Sim.ActiveSimState, ref Sim.ActiveSimDependent);
			case CellInfoType.Atmosphere:
				return GetCellInfoAtmosphere(ref Sim.ActiveSimState, ref Sim.ActiveSimDependent);
			case CellInfoType.Ground:
				return GetCellInfoGround(ref Sim.ActiveSimState, ref Sim.ActiveSimDependent);
			case CellInfoType.Water:
				return GetCellInfoWater(ref Sim.ActiveSimState, ref Sim.ActiveSimDependent);
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

	public static Color32 GetTerrainColor(float roughness, float soilFertility, float waterDepth, float iceCoverage, float vegetationCoverage)
	{
		var groundColor = Color32.Lerp(grey, brown, soilFertility);
		var waterColor = Color32.Lerp(groundColor, blue, math.saturate(math.pow(waterDepth / roughness, 2)));
		var iceColor = Color32.Lerp(waterColor, white, iceCoverage);
		var vegetationColor = Color32.Lerp(groundColor, green, vegetationCoverage);
		return vegetationColor;
	}

	public static Color32 GetWaterColor(float iceCoverage)
	{
		var waterColor = blue;
		var iceColor = Color32.Lerp(waterColor, white, iceCoverage);
		return iceColor;
	}

	public static Color32 GetCloudColor(float relativeHumidity, float cloudCoverage)
	{
		var c = Color32.Lerp(white, black, relativeHumidity);
		float opacity = -math.cos(cloudCoverage * math.PI) / 2 + 0.5f;
		c.a = (byte)(255 * opacity);
		return c;
	}

	public struct Overlay {
		public Overlay(float min, float max, NativeArray<CVP> colors, NativeArray<float> values)
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



	#region private functions

	private bool GetOverlay(MeshOverlay activeOverlay, ref SimState simState, out Overlay overlay)
	{
		float ticksPerYear = Sim.WorldData.TicksPerSecond * 60 * 60 * 24 * 365;
		float maxEvap = DisplayEvaporationMax * WorldData.MassWater / ticksPerYear;
		float maxRainfall = DisplayRainfallMax * WorldData.MassWater / ticksPerYear;
		switch (activeOverlay)
		{
			case MeshOverlay.AbsoluteHumidity:
				overlay = new Overlay(0, DisplayAbsoluteHumidityMax, _normalizedRainbow, simState.AirHumidity[0]);
				return true;
			//case MeshOverlay.RelativeHumidity:
			//	overlay = new Overlay(0, 1.0f, _normalizedRainbow, simState.AirHumidity[0]);
			//	break;
			//case MeshOverlay.GroundWater:
			//	overlay = new Overlay(0, DisplayGroundWaterMax, _normalizedRainbow, simState.GroundWater[0]);
			//	break;
			//case MeshOverlay.LowerAirPressure:
			//	overlay = new Overlay(DisplayAirPressureMin, DisplayAirPressureMax, _normalizedRainbow, simState.AirPressure[0]);
			//	return true;
			case MeshOverlay.LowerAirTemperature:
				overlay = new Overlay(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.AirTemperature[0]);
				return true;
			case MeshOverlay.ShallowWaterTemperature:
				overlay = new Overlay(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.WaterTemperature[Sim.StaticState.WaterLayers-1]);
				return true;
			case MeshOverlay.GroundTemperature:
				overlay = new Overlay(DisplayTemperatureMin, DisplayTemperatureMax, _normalizedRainbow, simState.TerrainTemperature);
				return true;
			//case MeshOverlay.HeatAbsorbed:
			//	overlay = new Overlay(0, DisplayHeatAbsorbedMax, _normalizedRainbow, simState.HeatAbsorbed);
			//	return true;
			//case MeshOverlay.Rainfall:
			//	overlay = new Overlay(0, maxRainfall, _normalizedRainbow);
			//	return true;
			//case MeshOverlay.Evaporation:
			//	overlay = new Overlay(0, maxEvap, _normalizedRainbow);
			//	return true;
			//case MeshOverlay.VerticalWind:
			//	overlay = new Overlay(-DisplayVerticalWindSpeedMax, DisplayVerticalWindSpeedMax, _normalizedRainbow);
			//	return true;
		}
		overlay = new Overlay();
		return false;

	}


	private string GetCellInfoGlobal(ref SimState state, ref SimDependent dependent)
	{
		string s = "";
		s += "CO2: " + state.PlanetState.CarbonDioxide;
		s += "\nCloud Coverage: " + (dependent.DisplayPlanet.CloudCoverage * 100 * Sim.InverseCellCount).ToString("0.0") + "%";
		s += "\nGlobal Sea Level: " + (dependent.DisplayPlanet.SeaLevel * Sim.InverseCellCount).ToString("0.00");
		s += "\nOcean Coverage: " + (dependent.DisplayPlanet.OceanCoverage * 100 * Sim.InverseCellCount).ToString("0.0") + "%";
		s += "\nOcean Volume: " + (dependent.DisplayPlanet.OceanVolume / 1000000000 * Sim.InverseCellCount).ToString("0.00") + " B";
		s += "\nTemperature: " + GetTemperatureString(dependent.DisplayPlanet.Temperature * Sim.InverseCellCount, ActiveTemperatureUnits, 2);
		s += "\nCloud Mass: " + (dependent.DisplayPlanet.CloudMass).ToString("0.00");
		s += "\nWater Vapor: " + (dependent.DisplayPlanet.WaterVapor).ToString("0.00");
		s += "\nRainfall: " + (dependent.DisplayPlanet.Rainfall * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater).ToString("0.00");
		s += "\nEvaporation: " + (dependent.DisplayPlanet.Evaporation * Sim.WorldData.TicksPerYear * Sim.InverseCellCount / WorldData.MassWater).ToString("0.00");

		return s;
	}
	private string GetCellInfoEnergy(ref SimState state, ref SimDependent dependent)
	{
		string s = "";

		var totalReflected = dependent.DisplayPlanet.EnergySolarReflectedCloud + dependent.DisplayPlanet.EnergySolarReflectedAtmosphere + dependent.DisplayPlanet.EnergySolarReflectedSurface;
		var totalOutgoing = dependent.DisplayPlanet.EnergyThermalSurfaceOutAtmosphericWindow + dependent.DisplayPlanet.EnergyThermalOutAtmosphere;
		s += "Delta: " + ConvertTileEnergyToWatts((dependent.DisplayPlanet.EnergyIncoming - totalReflected - totalOutgoing) * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Incoming: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergyIncoming * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Reflected: " + ConvertTileEnergyToWatts((totalReflected) * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Reflected Cloud: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergySolarReflectedCloud * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Reflected Atmos: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergySolarReflectedAtmosphere * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Reflected Surf: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergySolarReflectedSurface * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Abs Atm Total: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergySolarAbsorbedAtmosphere * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Abs Clouds: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergySolarAbsorbedCloud * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Abs Surface Total: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergySolarAbsorbedSurface * Sim.InverseCellCount).ToString("0.0");
		s += "\nS Abs Ocean: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergySolarAbsorbedOcean * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Outgoing: " + ConvertTileEnergyToWatts(totalOutgoing * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Out Atm Window: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergyThermalSurfaceOutAtmosphericWindow * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Out Radiation: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergyThermalOutAtmosphere * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Surface Radiation: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergyThermalSurfaceRadiation * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Atm Absorbed: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergyThermalAbsorbedAtmosphere * Sim.InverseCellCount).ToString("0.0");
		s += "\nT Back Radiation: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergyThermalBackRadiation * Sim.InverseCellCount).ToString("0.0");
		s += "\nEvapotranspiration: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergyEvapotranspiration * Sim.InverseCellCount).ToString("0.0");
		s += "\nSurface Conduction: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergySurfaceConduction * Sim.InverseCellCount).ToString("0.0");
		s += "\nOcean Radiation: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergyThermalOceanRadiation * Sim.InverseCellCount).ToString("0.0");
		s += "\nOcean Conduction: " + ConvertTileEnergyToWatts(dependent.DisplayPlanet.EnergyOceanConduction * Sim.InverseCellCount).ToString("0.0");



		return s;
	}
	private string GetCellInfoCell(ref SimState state, ref SimDependent dependent)
	{
		if (ActiveCellIndex < 0)
			return "";

		var terrain = state.CellTerrains[ActiveCellIndex];
		var coord = Sim.StaticState.Coordinate[ActiveCellIndex];
		string s = "";
		s += "Index: " + ActiveCellIndex + "\n";
		s += "Coord: (" + math.degrees(coord.x).ToString("0.0") + ", " + math.degrees(coord.y).ToString("0.0") + ")\n";
		s += "Surface: " + (terrain.Elevation + dependent.WaterAndIceDepth).ToString("0.000") + " m\n";
		s += "Elevation: " + terrain.Elevation.ToString("0.000") + " m\n";
		s += "H2O Depth: " + dependent.WaterDepth.ToString("0.000") + " m\n";
		s += "Ice Depth: " + (state.IceMass[ActiveCellIndex] / WorldData.MassIce).ToString("0.000") + " m\n";
		s += "Rainfall: " + (dependent.CellDisplays[ActiveCellIndex].Rainfall / WorldData.MassWater * 1000).ToString("0.000") + " mm\n";
		s += "Evaporation: " + (dependent.CellDisplays[ActiveCellIndex].Evaporation / WorldData.MassWater * 1000).ToString("0.000") + " mm\n";
		//		s += "Condensation: " + (display.Condensation / WorldData.MassWater * 1000000).ToString("0.000") + " nm3\n";
		return s;
	}
	private string GetCellInfoAtmosphere(ref SimState state, ref SimDependent dependent)
	{
		if (ActiveCellIndex < 0)
			return "";

		string s = "";
		s += "Surface Temp: " + GetTemperatureString(state.AirTemperature[0][ActiveCellIndex], ActiveTemperatureUnits, 0) + "\n";
		s += "Surface Pressure: " + dependent.AirPressure[0][ActiveCellIndex].ToString("0") + " Pa\n";
		s += "Surface Wind Horz: (" + state.AirVelocity[0][ActiveCellIndex].x.ToString("0.0") + ", " + state.AirVelocity[0][ActiveCellIndex].y.ToString("0.0") + ") m/s\n";
		s += "Surface Wind Vert: " + wind.WindVertical.ToString("0.0") + " m/s\n";
		s += "Surface Humidity: " + (dependent.RelativeHumidity[0][ActiveCellIndex] * 100).ToString("0.0") + "%\n";
		s += "Clouds: " + (state.CloudMass[ActiveCellIndex] / WorldData.MassWater).ToString("0.000") + " m3\n";
		s += "Droplets: " + (state.CloudDropletMass[ActiveCellIndex] * 1000000 / (WorldData.MassWater * state.CloudMass[ActiveCellIndex])).ToString("0.000") + " nm3\n";
		return s;
	}
	private string GetCellInfoGround(ref SimState state, ref SimDependent dependent)
	{
		if (ActiveCellIndex < 0)
			return "";

		var terrain = state.CellTerrains[ActiveCellIndex];
		string s = "";
		s += "Fertility: " + terrain.SoilFertility + "\n";
		s += "Temp: " + GetTemperatureString(state.TerrainTemperature[ActiveCellIndex], ActiveTemperatureUnits, 0);
		s += "H2O Volume: " + (terrain.GroundWater / cell.WaterMass).ToString("0.0") + " m3\n";
		s += "H2O Depth: " + terrain.GroundWaterDepth.ToString("0") + " m\n";
		s += "Roughness: " + terrain.Roughness.ToString("0") + " m\n";
		return s;
	}
	private string GetCellInfoWater(ref SimState state, ref SimDependent dependent)
	{
		if (ActiveCellIndex < 0)
			return "";

		string s = "";
		float saltMass = state.WaterSaltMass[Sim.StaticState.WaterLayers - 1][ActiveCellIndex];
		s += "Surface Temp: " + GetTemperatureString(state.WaterTemperature[Sim.StaticState.WaterLayers-1][ActiveCellIndex], ActiveTemperatureUnits, 0) + "\n";
		s += "Surface Salinity: " + (1000000 * saltMass / (saltMass + cell.WaterMass)).ToString("0.0") + " ppm\n";
		return s;
	}

	#endregion
}
