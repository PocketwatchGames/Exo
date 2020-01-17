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

	[Header("Display")]
	public bool LerpStates = true;
	public float TerrainScale = 100f;
//	public TemperatureDisplayType TemperatureDisplay;
	public float MinElevation = -11000;
	public float MaxElevation = 10000;
	public float MaxDepth = 11000;
	public float maxCloudColor = 300.0f;
	public float WaterDepthThreshold = 10;

	public float RRDisplayEnergyAborsobedMax = 300;
	public float DisplayRainfallMax = 5.0f;
	public float DisplayMinSalinity = 0;
	public float DisplayMaxSalinity = 50;
	public float DisplayMaxWindSpeedLowerAtm = 50;
	public float DisplayMaxWindSpeedUpperAtm = 250;
	public float DisplayMaxWindSpeedSurfaceWater = 5;
	public float DisplayMaxWindSpeedDeepWater = 0.5f;
	public float DisplayMaxVerticalWindSpeed = 1.0f;
	public float DisplayEvaporationMax = 5.0f;
	public float DisplayMaxCanopy = 1000;
	public float DisplayTemperatureMin = 223;
	public float DisplayTemperatureMax = 323;
	public float DisplayAbsoluteHumidityMax = 400;
	public float DisplayGroundWaterMax = 100000;
	public float DisplayAirPressureMin = 97000;
	public float DisplayAirPressureMax = 110000;



	[Header("References")]
	public WorldSim Sim;
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

	private bool _indicesInitialized;
	private int[] indices;

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

		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[0], ref Sim.WorldData, ref Sim.StaticState);
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
		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState);
	}

	public void BuildRenderState(ref SimState from, ref RenderState to, ref WorldData worldData, ref StaticState staticState)
	{
		to.Ticks = from.PlanetState.Ticks;
		to.Position = from.PlanetState.Position;
		to.Rotation = from.PlanetState.Rotation;
		for (int i = 0; i < from.Cells.Length; i++)
		{
			ref var fromCell = ref from.Cells[i];
			if (ActiveMeshOverlay == MeshOverlay.None)
			{
				to.TerrainColor[i] = GetTerrainColor(ref worldData, ref staticState, ref fromCell);
				to.WaterColor[i] = GetWaterColor(ref worldData, ref staticState, ref fromCell);
			} else
			{
				var overlayColor = GetOverlayColor(ref fromCell);
				to.TerrainColor[i] = overlayColor;
				to.WaterColor[i] = overlayColor;
			}
			to.CloudColor[i] = GetCloudColor(ref staticState, ref fromCell);
			to.TerrainNormal[i] = Sim.Icosphere.Vertices[i];
			to.WaterNormal[i] = Sim.Icosphere.Vertices[i];
			to.CloudNormal[i] = Sim.Icosphere.Vertices[i];
			to.TerrainPosition[i] = Sim.Icosphere.Vertices[i] * ((fromCell.Elevation + fromCell.Roughness) * TerrainScale + staticState.PlanetRadius) / staticState.PlanetRadius;
			to.WaterPosition[i] = Sim.Icosphere.Vertices[i] * ((fromCell.Elevation + fromCell.WaterDepth) * TerrainScale + staticState.PlanetRadius) / staticState.PlanetRadius * math.clamp(fromCell.WaterDepth / fromCell.Roughness, 0, 1);
			to.CloudPosition[i] = Sim.Icosphere.Vertices[i] * (fromCell.CloudElevation * TerrainScale + staticState.PlanetRadius) / staticState.PlanetRadius;
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
		state.Rotation = Quaternion.Lerp(lastState.Rotation, nextState.Rotation, t);
		state.Position = math.lerp(lastState.Position, nextState.Position, t);
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

		transform.SetPositionAndRotation(state.Position, state.Rotation);
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
		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState);
	}
	public void OnHUDWindChanged(UnityEngine.UI.Dropdown dropdown)
	{
		ActiveWindOverlay = (WorldMesh.WindOverlay)dropdown.value;
		_tickLerpTime = 0;
		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[_nextRenderState], ref Sim.WorldData, ref Sim.StaticState);
	}


	#region private functions

	private Color32 GetTerrainColor(ref WorldData worldData, ref StaticState staticState, ref SimStateCell cell)
	{
		var groundColor = Color32.Lerp(grey, brown, cell.SoilFertility);
		var waterColor = Color32.Lerp(groundColor, blue, math.clamp(math.pow(cell.WaterDepth / cell.Roughness, 2), 0, 1));
		var iceColor = Color32.Lerp(waterColor, white, math.clamp(cell.IceMass / (WorldData.MassIce * worldData.FullIceCoverage), 0, 1));
		var vegetationColor = Color32.Lerp(groundColor, green, math.clamp(cell.Vegetation / worldData.FullCanopyCoverage, 0, 1));
		return vegetationColor;
	}

	private Color32 GetWaterColor(ref WorldData worldData, ref StaticState staticState, ref SimStateCell cell)
	{
		var waterColor = blue;
		var iceColor = Color32.Lerp(waterColor, white, math.clamp(cell.IceMass / (WorldData.MassIce * worldData.FullIceCoverage), 0, 1));
		return iceColor;
	}

	private Color32 GetCloudColor(ref StaticState staticState, ref SimStateCell cell)
	{
		var c = Color32.Lerp(white, black, cell.RelativeHumidity);
		float opacity = -math.cos(cell.CloudCoverage * math.PI)/ 2 + 0.5f;
		c.a = (byte)(255 * opacity);
		return c;
	}

	private Color32 GetOverlayColor(ref SimStateCell cell)
	{
		switch (ActiveMeshOverlay)
		{
			case MeshOverlay.AbsoluteHumidity:
				return Lerp(NormalizedRainbow, cell.AirWaterMass, 0, DisplayAbsoluteHumidityMax);
			case MeshOverlay.RelativeHumidity:
				return Lerp(NormalizedRainbow, cell.RelativeHumidity, 0, 1.0f);
			case MeshOverlay.GroundWater:
				return Lerp(NormalizedRainbow, cell.GroundWater, 0, DisplayGroundWaterMax);
			case MeshOverlay.LowerAirPressure:
				return Lerp(NormalizedRainbow, cell.AirPressure, DisplayAirPressureMin, DisplayAirPressureMax);
			case MeshOverlay.LowerAirTemperature:
				return Lerp(NormalizedRainbow, cell.AirTemperature, DisplayTemperatureMin, DisplayTemperatureMax);
			case MeshOverlay.ShallowWaterTemperature:
				return Lerp(NormalizedRainbow, cell.WaterTemperature, DisplayTemperatureMin, DisplayTemperatureMax);
			case MeshOverlay.GroundTemperature:
				return Lerp(NormalizedRainbow, Atmosphere.GetLandTemperature(ref Sim.WorldData, cell.GroundEnergy, cell.GroundWater, cell.SoilFertility, cell.Vegetation), DisplayTemperatureMin, DisplayTemperatureMax);
			case MeshOverlay.VerticalWind:
				return Lerp(NormalizedBlueBlackRed, cell.WindVertical, -DisplayMaxVerticalWindSpeed, DisplayMaxVerticalWindSpeed);
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
	static List<CVP> NormalizedBlueBlackRed = new List<CVP> {
											new CVP(Color.blue, 0),
											new CVP(Color.black, 0.5f),
											new CVP(Color.red, 1) };
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
