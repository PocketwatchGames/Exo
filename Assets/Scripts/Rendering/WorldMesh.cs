using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;
using Unity.Burst;

public class WorldMesh : MonoBehaviour {

	public enum MeshOverlay {
		None,
		TerrainTemperature,
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
		CloudCoalescence,
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

	[Header("Display")]
	public float ElevationScale = 0.002f;
	public float MinZoom;
	public float MaxZoom;
	public float Zoom { get { return MinZoom + MaxZoom * (float)Mathf.Pow(ZoomLevel, 3); } }
	public float ZoomLevel = 0.5f;
	public float CameraMoveSpeed = 2;
	public float CameraZoomSpeed = 2;
//	public TemperatureDisplayType TemperatureDisplay;
	public float minPressure = 300;
	public float maxPressure = 600;
	public float MinElevation = -11000;
	public float MaxElevation = 10000;
	public float MaxDepth = 11000;
	public float maxHumidity = 50;
	public float maxRainfall = 5.0f;
	public float maxCloudColor = 300.0f;
	public float MaxEnergyAbsorbed = 300;
	public float MinSalinity = 0;
	public float MaxSalinity = 50;
	public float waterDepthThreshold = 10;
	public float DisplayMaxWindSpeedLowerAtm = 50;
	public float DisplayMaxWindSpeedUpperAtm = 250;
	public float DisplayMaxWindSpeedSurfaceWater = 5;
	public float DisplayMaxWindSpeedDeepWater = 0.5f;
	public float DisplayMaxVerticalWindSpeed = 1.0f;
	public float MaxEvap = 5.0f;
	public float DisplayMaxCanopy = 1000;
	public float DisplayMinTemperature = 223;
	public float DisplayMaxTemperature = 323;


	public bool LerpStates = true;
	public float DistanceToSun = 100;
	public float TerrainScale = 0.00005f;
	public MeshOverlay ActiveMeshOverlay;
	public WindOverlay ActiveWindOverlay;

	public Material TerrainMaterial;
	public Material WaterMaterial;
	public Material CloudMaterial;

	static Color32 green = new Color32(0, 220, 30, 255);
	static Color32 grey = new Color32(50, 50, 80, 255);
	static Color32 brown = new Color32(100, 60, 20, 255);
	static Color32 blue = new Color32(0, 0, 255, 255);
	static Color32 white = new Color32(255, 255, 255, 255);
	static Color32 black = new Color32(0, 0, 0, 255);


	public WorldSim Sim;
	private RenderState[] _renderStates;
	private int _curRenderState;
	private int _lastRenderState;
	private int _nextRenderState;
	private const int _renderStateCount = 3;

	private Mesh _terrainMesh;
	private Mesh _waterMesh;
	private Mesh _cloudMesh;

	private GameObject _terrainObject;
	private GameObject _waterObject;
	private GameObject _cloudObject;

	private float _tickLerpTime = 0;
	private float _tickLerpTimeTotal = 1;


	bool _indicesInitialized;
	int[] indices;

	public void Start()
	{
		Sim.OnTick += OnSimTick;

		if (_terrainMesh)
		{
			GameObject.Destroy(_terrainMesh);
		}
		_terrainMesh = new Mesh();
		_cloudMesh = new Mesh();
		_waterMesh = new Mesh();

		_terrainObject = new GameObject("Terrain Mesh");
		_terrainObject.transform.SetParent(gameObject.transform, false);
		var terrainFilter = _terrainObject.AddComponent<MeshFilter>();
		var terrainSurfaceRenderer = _terrainObject.AddComponent<MeshRenderer>();
		terrainSurfaceRenderer.material = TerrainMaterial;
		terrainFilter.mesh = _terrainMesh;

		_waterObject = new GameObject("Water Mesh");
		_waterObject.transform.SetParent(gameObject.transform, false);
		var waterFilter = _waterObject.AddComponent<MeshFilter>();
		var waterSurfaceRenderer = _waterObject.AddComponent<MeshRenderer>();
		waterSurfaceRenderer.material = WaterMaterial;
		waterFilter.mesh = _waterMesh;

		_cloudObject = new GameObject("Cloud Mesh");
		_cloudObject.transform.SetParent(gameObject.transform, false);
		var cloudFilter = _cloudObject.AddComponent<MeshFilter>();
		var cloudSurfaceRenderer = _cloudObject.AddComponent<MeshRenderer>();
		cloudSurfaceRenderer.material = CloudMaterial;
		cloudFilter.mesh = _cloudMesh;



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

		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[0], Sim.StaticState);
		UpdateMesh(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], ref _renderStates[_curRenderState]);
	}
	public void OnDestroy()
	{
		Sim.OnTick -= OnSimTick;
	}

	public void Update()
	{
		_tickLerpTime -= Time.deltaTime * Sim.TimeScale;
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
		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[_nextRenderState], Sim.StaticState);
	}

	public void BuildRenderState(ref SimState from, ref RenderState to, StaticState staticState)
	{
		to.Ticks = from.Ticks;
		to.OrbitSpeed = from.OrbitSpeed;
		to.SpinAngle = from.SpinAngle;
		to.SpinSpeed = from.SpinSpeed;
		to.TiltAngle = from.TiltAngle;
		for (int i = 0; i < from.Count; i++)
		{
			ref var fromCell = ref from.Cells[i];
			if (ActiveMeshOverlay == MeshOverlay.None)
			{
				to.TerrainColor[i] = GetTerrainColor(staticState, fromCell);
				to.WaterColor[i] = GetWaterColor(staticState, fromCell);
			} else
			{
				var overlayColor = GetOverlayColor(fromCell);
				to.TerrainColor[i] = overlayColor;
				to.WaterColor[i] = overlayColor;
			}
			to.CloudColor[i] = GetCloudColor(staticState, fromCell);
			to.TerrainNormal[i] = Sim.Icosphere.Vertices[i];
			to.WaterNormal[i] = Sim.Icosphere.Vertices[i];
			to.CloudNormal[i] = Sim.Icosphere.Vertices[i];
			to.TerrainPosition[i] = Sim.Icosphere.Vertices[i] * ((fromCell.Elevation + fromCell.Roughness) * TerrainScale + staticState.Radius) / staticState.Radius;
			to.WaterPosition[i] = Sim.Icosphere.Vertices[i] * ((fromCell.Elevation + fromCell.WaterDepth) * TerrainScale + staticState.Radius) / staticState.Radius * math.clamp(fromCell.WaterDepth / fromCell.Roughness, 0, 1);
			to.CloudPosition[i] = Sim.Icosphere.Vertices[i] * (fromCell.CloudElevation * TerrainScale + staticState.Radius) / staticState.Radius;
		}

	}

	public void UpdateMesh(ref RenderState lastState, ref RenderState nextState, ref RenderState state)
	{
		float t = Mathf.Clamp01(1.0f - _tickLerpTime / _tickLerpTimeTotal);
		if (!LerpStates)
		{
			t = 1;
		}
		state.Ticks = (nextState.Ticks - lastState.Ticks) * t + lastState.Ticks;
		state.TiltAngle = Mathf.LerpAngle(lastState.TiltAngle, nextState.TiltAngle, t);
		for (int i = 0; i < state.CloudColor.Length; i++)
		{
			state.TerrainColor[i] = Color32.Lerp(lastState.TerrainColor[i], nextState.TerrainColor[i], t);
			state.TerrainPosition[i] = Vector3.Lerp(lastState.TerrainPosition[i], nextState.TerrainPosition[i], t);
			state.TerrainNormal[i] = Vector3.Lerp(lastState.TerrainNormal[i], nextState.TerrainNormal[i], t);
			state.WaterColor[i] = Color32.Lerp(lastState.WaterColor[i], nextState.WaterColor[i], t);
			state.WaterPosition[i] = Vector3.Lerp(lastState.WaterPosition[i], nextState.WaterPosition[i], t);
			state.WaterNormal[i] = Vector3.Lerp(lastState.WaterNormal[i], nextState.WaterNormal[i], t);
			state.CloudColor[i] = Color32.Lerp(lastState.CloudColor[i], nextState.CloudColor[i], t);
			state.CloudPosition[i] = Vector3.Lerp(lastState.CloudPosition[i], nextState.CloudPosition[i], t);
			state.CloudNormal[i] = Vector3.Lerp(lastState.CloudNormal[i], nextState.CloudNormal[i], t);
		}

		_terrainMesh.vertices = state.TerrainPosition;
		_terrainMesh.normals = state.TerrainNormal;
		_terrainMesh.colors32 = state.TerrainColor;

		_waterMesh.vertices = state.WaterPosition;
		_waterMesh.normals = state.WaterNormal;
		_waterMesh.colors32 = state.WaterColor;

		_cloudMesh.vertices = state.CloudPosition;
		_cloudMesh.normals = state.CloudNormal;
		_cloudMesh.colors32 = state.CloudColor;

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

		var time = WorldTime.GetTime(state.Ticks, state.SpinSpeed);
		var spin = Quaternion.Euler(state.TiltAngle, time * 360, 0);
		float orbitAngle = state.Ticks * state.OrbitSpeed * 360;
		transform.SetPositionAndRotation(new Vector3(math.cos(orbitAngle), 0, math.sin(orbitAngle)) * DistanceToSun, spin);
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
		ActiveMeshOverlay = (WorldMesh.MeshOverlay)dropdown.value;
		_tickLerpTime = 0;
		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[_nextRenderState], Sim.StaticState);
	}
	public void OnHUDWindChanged(UnityEngine.UI.Dropdown dropdown)
	{
		ActiveWindOverlay = (WorldMesh.WindOverlay)dropdown.value;
		_tickLerpTime = 0;
		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[_nextRenderState], Sim.StaticState);
	}


	#region private functions

	private Color32 GetTerrainColor(StaticState staticState, SimStateCell cell)
	{
		var groundColor = Color32.Lerp(grey, brown, cell.SoilFertility);
		var waterColor = Color32.Lerp(groundColor, blue, math.clamp(math.pow(cell.WaterDepth / cell.Roughness, 2), 0, 1));
		var iceColor = Color32.Lerp(waterColor, white, cell.IceMass);
		var vegetationColor = Color32.Lerp(groundColor, green, cell.Vegetation);
		return vegetationColor;
	}

	private Color32 GetWaterColor(StaticState staticState, SimStateCell cell)
	{
		return Color32.Lerp(blue, white, cell.IceMass);
	}

	private Color32 GetCloudColor(StaticState staticState, SimStateCell cell)
	{
		var c = Color32.Lerp(white, black, cell.RelativeHumidity);
		float opacity = -math.cos(cell.CloudCoverage * math.PI)/ 2 + 0.5f;
		c.a = (byte)(255 * opacity);
		return c;
	}

	private Color32 GetOverlayColor(SimStateCell cell)
	{
		switch (ActiveMeshOverlay)
		{
			case MeshOverlay.LowerAirTemperature:
				return Lerp(NormalizedRainbow, cell.AirTemperature, DisplayMinTemperature, DisplayMaxTemperature);
			default:
				return black;
		}
	}

	static List<CVP> NormalizedRainbow = new List<CVP> {
											new CVP(Color.black, 0),
											new CVP(Color.white, 0.1667f),
											new CVP(Color.blue, 0.3333f),
											new CVP(Color.green, 0.5f),
											new CVP(Color.yellow, 0.6667f),
											new CVP(Color.red, 0.8333f),
											new CVP(Color.magenta, 1) };
	struct CVP {
		public Color Color;
		public float Value;
		public CVP(Color c, float v) { Color = c; Value = v; }
	};

	Color Lerp(List<CVP> colors, float value, float min, float max)
	{
		return Lerp(colors, (value - min) / (max - min));
	}
	Color Lerp(List<CVP> colors, float value)
	{
		for (int i = 0; i < colors.Count - 1; i++)
		{
			if (value < colors[i + 1].Value)
			{
				return Color.Lerp(colors[i].Color, colors[i + 1].Color, (value - colors[i].Value) / (colors[i + 1].Value - colors[i].Value));
			}
		}
		return colors[colors.Count - 1].Color;
	}


	#endregion
}
