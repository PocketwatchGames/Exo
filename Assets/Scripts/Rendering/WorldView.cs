using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;
using Unity.Burst;

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
	public float DisplayHeatAbsorbedMax = 1000;

	float inverseMaxEvapMass;
	float inverseMaxRainfall;


	[Header("References")]
	public WorldSim Sim;
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

	private Mesh _terrainMesh;
	private Mesh _waterMesh;
	private Mesh _cloudMesh;

	private GameObject _terrainObject;
	private GameObject _waterObject;
	private GameObject _cloudObject;
	private GameObject _selectionCircle;

	private GameObject[] _windArrows;

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
		_selectionCircle.transform.localScale *= 0.1f;

		_windArrows = new GameObject[Sim.Icosphere.Vertices.Count];
		for (int i = 0; i < Sim.Icosphere.Vertices.Count;i++)
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

		BuildRenderState(ref Sim.ActiveSimState, ref _renderStates[0], ref Sim.WorldData, ref Sim.StaticState);
		UpdateMesh(ref _renderStates[_lastRenderState], ref _renderStates[_nextRenderState], ref _renderStates[_curRenderState]);
	}
	public void OnDestroy()
	{
		Sim.OnTick -= OnSimTick;
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

		UpdateWind();
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
		inverseMaxEvapMass = worldData.TicksPerSecond * 60 * 60 * 24 * 365 / (DisplayEvaporationMax * WorldData.MassWater);
		inverseMaxRainfall = worldData.TicksPerSecond * 60 * 60 * 24 * 365 / (DisplayRainfallMax * WorldData.MassWater);



		to.Ticks = from.PlanetState.Ticks;
		to.Position = from.PlanetState.Position;
		to.Rotation = math.degrees(from.PlanetState.Rotation);
		for (int i = 0; i < from.Cells.Length; i++)
		{
			ref var fromCell = ref from.Cells[i];
			ref var fromWind = ref from.Wind[i];
			if (ActiveMeshOverlay == MeshOverlay.None)
			{
				to.TerrainColor[i] = GetTerrainColor(ref worldData, ref staticState, ref fromCell);
				to.WaterColor[i] = GetWaterColor(ref worldData, ref staticState, ref fromCell);
			} else
			{
				var overlayColor = GetOverlayColor(ref fromCell, ref fromWind, from.DisplayCells[i]);
				to.TerrainColor[i] = overlayColor;
				to.WaterColor[i] = overlayColor;
			}
			to.CloudColor[i] = GetCloudColor(ref staticState, ref fromCell);
			to.TerrainNormal[i] = Sim.Icosphere.Vertices[i];
			to.WaterNormal[i] = Sim.Icosphere.Vertices[i];
			to.CloudNormal[i] = Sim.Icosphere.Vertices[i];
			to.SurfacePosition[i] = Sim.Icosphere.Vertices[i] * ((fromCell.Elevation + math.max(fromCell.Roughness, fromCell.WaterAndIceDepth)) * TerrainScale + staticState.PlanetRadius) / staticState.PlanetRadius;
			to.TerrainPosition[i] = Sim.Icosphere.Vertices[i] * ((fromCell.Elevation + fromCell.Roughness) * TerrainScale + staticState.PlanetRadius) / staticState.PlanetRadius;
			to.WaterPosition[i] = Sim.Icosphere.Vertices[i] * ((fromCell.Elevation + fromCell.WaterDepth) * TerrainScale + staticState.PlanetRadius) / staticState.PlanetRadius * math.saturate(fromCell.WaterDepth / fromCell.Roughness);
			to.CloudPosition[i] = Sim.Icosphere.Vertices[i] * (fromCell.CloudElevation * TerrainScale + staticState.PlanetRadius) / staticState.PlanetRadius;
			to.WaterDepth[i] = fromCell.WaterDepth;
			to.Wind[i] = fromWind.WindSurface;
			to.Current[i] = fromWind.CurrentSurface;
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
		state.Rotation = new Vector3(Mathf.LerpAngle(lastState.Rotation.x, nextState.Rotation.x, t), Mathf.LerpAngle(lastState.Rotation.y, nextState.Rotation.y, t), Mathf.LerpAngle(lastState.Rotation.z, nextState.Rotation.z, t));
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
			state.SurfacePosition[i] = Vector3.Lerp(lastState.SurfacePosition[i], nextState.SurfacePosition[i], t);
			state.Wind[i] = Vector2.Lerp(lastState.Wind[i], nextState.Wind[i], t);
			state.Current[i] = Vector2.Lerp(lastState.Current[i], nextState.Current[i], t);
			state.WaterDepth[i] = math.lerp(lastState.WaterDepth[i], nextState.WaterDepth[i], t);
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

		Planet.transform.SetPositionAndRotation(state.Position, Quaternion.Euler(state.Rotation));
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
				return GetCellInfoGlobal(ref Sim.ActiveSimState);
			case CellInfoType.Energy:
				return GetCellInfoEnergy(ref Sim.ActiveSimState);
			case CellInfoType.Cell:
				return GetCellInfoCell(ref Sim.ActiveSimState);
			case CellInfoType.Atmosphere:
				return GetCellInfoAtmosphere(ref Sim.ActiveSimState);
			case CellInfoType.Ground:
				return GetCellInfoGround(ref Sim.ActiveSimState);
			case CellInfoType.Water:
				return GetCellInfoWater(ref Sim.ActiveSimState);
		}
		return "";
	}

	public int GetClosestVert(int triangleIndex)
	{
		return indices[triangleIndex * 3];
	}

	public void SetActiveCell(int index, bool locked)
	{
		_selectionCircle.SetActive(index >= 0);
		ActiveCellIndex = index;
		ActiveCellLocked = locked;
		if (index >= 0)
		{
			//			var p = Sim.Icosphere.Vertices[index];
			var pos = _renderStates[_curRenderState].SurfacePosition[index];
			_selectionCircle.transform.localPosition = pos;
			_selectionCircle.transform.localRotation = Quaternion.LookRotation(-pos);
		}
	}

	float ConvertTileEnergyToWatts(float energy)
	{
		return energy * 1000 / Sim.WorldData.SecondsPerTick;
	}



	#region private functions


	private Color32 GetTerrainColor(ref WorldData worldData, ref StaticState staticState, ref SimCell cell)
	{
		var groundColor = Color32.Lerp(grey, brown, cell.SoilFertility);
		var waterColor = Color32.Lerp(groundColor, blue, math.saturate(math.pow(cell.WaterDepth / cell.Roughness, 2)));
		var iceColor = Color32.Lerp(waterColor, white, math.saturate(cell.IceMass / (WorldData.MassIce * worldData.FullIceCoverage)));
		var vegetationColor = Color32.Lerp(groundColor, green, math.saturate(cell.Vegetation / worldData.FullCanopyCoverage));
		return vegetationColor;
	}

	private Color32 GetWaterColor(ref WorldData worldData, ref StaticState staticState, ref SimCell cell)
	{
		var waterColor = blue;
		var iceColor = Color32.Lerp(waterColor, white, math.saturate(cell.IceMass / (WorldData.MassIce * worldData.FullIceCoverage)));
		return iceColor;
	}

	private Color32 GetCloudColor(ref StaticState staticState, ref SimCell cell)
	{
		var c = Color32.Lerp(white, black, cell.RelativeHumidity);
		float opacity = -math.cos(cell.CloudCoverage * math.PI)/ 2 + 0.5f;
		c.a = (byte)(255 * opacity);
		return c;
	}

	private Color32 GetOverlayColor(ref SimCell cell, ref SimWind wind, DisplayCell displayCell)
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
			case MeshOverlay.HeatAbsorbed:
				return Lerp(NormalizedRainbow, displayCell.Heat, 0, DisplayHeatAbsorbedMax);
			case MeshOverlay.Rainfall:
				return Lerp(NormalizedRainbow, displayCell.Rainfall * inverseMaxRainfall, 0, DisplayRainfallMax);
			case MeshOverlay.Evaporation:
				return Lerp(NormalizedRainbow, displayCell.Evaporation * inverseMaxEvapMass, 0, DisplayEvaporationMax);
			case MeshOverlay.VerticalWind:
				return Lerp(NormalizedBlueBlackRed, wind.WindVertical, -DisplayMaxVerticalWindSpeed, DisplayMaxVerticalWindSpeed);
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


	private string GetCellInfoGlobal(ref SimState state)
	{
		float inverseCellCount = 1.0f / state.Cells.Length;
		string s = "";
		s += "CO2: " + state.PlanetState.CarbonDioxide;
		s += "\nCloud Coverage: " + (state.DisplayPlanet.CloudCoverage * 100 * inverseCellCount).ToString("0.0") + "%";
		s += "\nGlobal Sea Level: " + (state.DisplayPlanet.SeaLevel * inverseCellCount).ToString("0.00");
		s += "\nOcean Coverage: " + (state.DisplayPlanet.OceanCoverage * 100 * inverseCellCount).ToString("0.0") + "%";
		s += "\nOcean Volume: " + (state.DisplayPlanet.OceanVolume / 1000000000 * inverseCellCount).ToString("0.00") + " B";
		s += "\nTemperature: " + GetTemperatureString(state.DisplayPlanet.Temperature * inverseCellCount, ActiveTemperatureUnits, 2);
		s += "\nAtmospheric Mass: " + (state.DisplayPlanet.AtmosphericMass / 1000).ToString("0") + " K";
		s += "\nCloud Mass: " + (state.DisplayPlanet.CloudMass).ToString("0.00");
		s += "\nWater Vapor: " + (state.DisplayPlanet.WaterVapor).ToString("0.00");
		s += "\nRainfall: " + (state.DisplayPlanet.Rainfall * Sim.WorldData.TicksPerYear * inverseCellCount / WorldData.MassWater).ToString("0.00");
		s += "\nEvaporation: " + (state.DisplayPlanet.Evaporation * Sim.WorldData.TicksPerYear * inverseCellCount / WorldData.MassWater).ToString("0.00");

		return s;
	}
	private string GetCellInfoEnergy(ref SimState state)
	{
		string s = "";

		float inverseCellCount = 1.0f / state.Cells.Length;
		var totalReflected = state.DisplayPlanet.EnergySolarReflectedCloud + state.DisplayPlanet.EnergySolarReflectedAtmosphere + state.DisplayPlanet.EnergySolarReflectedSurface;
		var totalOutgoing = state.DisplayPlanet.EnergyThermalOutAtmosphericWindow + state.DisplayPlanet.EnergyThermalOutAtmosphere;
		s += "Delta: " + ConvertTileEnergyToWatts((state.DisplayPlanet.EnergyIncoming - totalReflected - totalOutgoing) * inverseCellCount).ToString("0.0");
		s += "\nS Incoming: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergyIncoming * inverseCellCount).ToString("0.0");
		s += "\nS Reflected: " + ConvertTileEnergyToWatts((totalReflected) * inverseCellCount).ToString("0.0");
		s += "\nS Reflected Cloud: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergySolarReflectedCloud * inverseCellCount).ToString("0.0");
		s += "\nS Reflected Atmos: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergySolarReflectedAtmosphere * inverseCellCount).ToString("0.0");
		s += "\nS Reflected Surf: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergySolarReflectedSurface * inverseCellCount).ToString("0.0");
		s += "\nS Abs Atm Total: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergySolarAbsorbedAtmosphere * inverseCellCount).ToString("0.0");
		s += "\nS Abs Clouds: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergySolarAbsorbedCloud * inverseCellCount).ToString("0.0");
		s += "\nS Abs Surface Total: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergySolarAbsorbedSurface * inverseCellCount).ToString("0.0");
		s += "\nS Abs Ocean: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergySolarAbsorbedOcean * inverseCellCount).ToString("0.0");
		s += "\nT Outgoing: " + ConvertTileEnergyToWatts(totalOutgoing * inverseCellCount).ToString("0.0");
		s += "\nT Out Atm Window: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergyThermalOutAtmosphericWindow * inverseCellCount).ToString("0.0");
		s += "\nT Out Radiation: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergyThermalOutAtmosphere * inverseCellCount).ToString("0.0");
		s += "\nT Surface Radiation: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergyThermalSurfaceRadiation * inverseCellCount).ToString("0.0");
		s += "\nT Atm Absorbed: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergyThermalAbsorbedAtmosphere * inverseCellCount).ToString("0.0");
		s += "\nT Back Radiation: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergyThermalBackRadiation * inverseCellCount).ToString("0.0");
		s += "\nEvapotranspiration: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergyEvapotranspiration * inverseCellCount).ToString("0.0");
		s += "\nSurface Conduction: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergySurfaceConduction * inverseCellCount).ToString("0.0");
		s += "\nOcean Radiation: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergyThermalOceanRadiation * inverseCellCount).ToString("0.0");
		s += "\nOcean Conduction: " + ConvertTileEnergyToWatts(state.DisplayPlanet.EnergyOceanConduction * inverseCellCount).ToString("0.0");



		return s;
	}
	private string GetCellInfoCell(ref SimState state)
	{
		if (ActiveCellIndex < 0)
			return "";

		var cell = state.Cells[ActiveCellIndex];
		var display = state.DisplayCells[ActiveCellIndex];
		var coord = Sim.StaticState.Coordinate[ActiveCellIndex];
		string s = "";
		s += "Index: " + ActiveCellIndex + "\n";
		s += "Coord: (" + math.degrees(coord.x).ToString("0.0") + ", " + math.degrees(coord.y).ToString("0.0") + ")\n";
		s += "Surface: " + (cell.Elevation + cell.WaterAndIceDepth).ToString("0.000") + " m\n";
		s += "Elevation: " + cell.Elevation.ToString("0.000") + " m\n";
		s += "H2O Depth: " + cell.WaterDepth.ToString("0.000") + " m\n";
		s += "Ice Depth: " + (cell.IceMass / WorldData.MassIce).ToString("0.000") + " m\n";
		s += "Rainfall: " + (display.Rainfall / WorldData.MassWater * 1000).ToString("0.000") + " mm\n";
		s += "Evaporation: " + (display.Evaporation / WorldData.MassWater * 1000).ToString("0.000") + " mm\n";
		//		s += "Condensation: " + (display.Condensation / WorldData.MassWater * 1000000).ToString("0.000") + " nm3\n";
		return s;
	}
	private string GetCellInfoAtmosphere(ref SimState state)
	{
		if (ActiveCellIndex < 0)
			return "";

		var cell = state.Cells[ActiveCellIndex];
		var wind = state.Wind[ActiveCellIndex];
		string s = "";
		s += "Temp: " + GetTemperatureString(cell.AirTemperature, ActiveTemperatureUnits, 0) + "\n";
		s += "Pressure: " + cell.AirPressure.ToString("0") + " Pa\n";
		s += "Wind Horz: (" + wind.WindSurface.x.ToString("0.0") + ", " + wind.WindSurface.y.ToString("0.0") + ") m/s\n";
		s += "Wind Vert: " + wind.WindVertical.ToString("0.0") + " m/s\n";
		s += "Humidity: " + (cell.RelativeHumidity * 100).ToString("0.0") + "%\n";
		s += "Clouds: " + (cell.CloudMass / WorldData.MassWater).ToString("0.000") + " m3\n";
		s += "Droplets: " + (cell.CloudDropletMass * 1000000 / (WorldData.MassWater * cell.CloudMass)).ToString("0.000") + " nm3\n";
		return s;
	}
	private string GetCellInfoGround(ref SimState state)
	{
		if (ActiveCellIndex < 0)
			return "";

		var cell = state.Cells[ActiveCellIndex];
		string s = "";
		s += "Fertility: " + cell.SoilFertility + "\n";
		s += "Temp: " + GetTemperatureString(Atmosphere.GetLandTemperature(ref Sim.WorldData, cell.GroundEnergy, cell.GroundWater, cell.SoilFertility, math.saturate(cell.Vegetation / Sim.WorldData.FullCanopyCoverage)), ActiveTemperatureUnits, 0);
		s += "H2O Volume: " + (cell.GroundWater / cell.WaterMass).ToString("0.0") + " m3\n";
		s += "H2O Depth: " + cell.GroundWaterDepth.ToString("0") + " m\n";
		s += "Roughness: " + cell.Roughness.ToString("0") + " m\n";
		return s;
	}
	private string GetCellInfoWater(ref SimState state)
	{
		if (ActiveCellIndex < 0)
			return "";

		var cell = state.Cells[ActiveCellIndex];
		string s = "";
		s += "Temp: " + GetTemperatureString(cell.WaterTemperature, ActiveTemperatureUnits, 0) + "\n";
		s += "Salinity: " + (1000000 * cell.SaltMass / (cell.SaltMass + cell.WaterMass)).ToString("0.0") + " ppm\n";
		return s;
	}
	private void UpdateWind()
	{
		switch (ActiveWindOverlay)
		{
			case WindOverlay.LowerAirWind:
				UpdateArrowsWind();
				break;
			case WindOverlay.ShallowWaterCurrent:
				UpdateArrowsCurrent();
				break;
		}
	}
	private void UpdateArrowsWind()
	{
		for (int i=0;i<_windArrows.Length;i++)
		{
			var wind = _renderStates[_curRenderState].Wind[i];

			float speed = wind.magnitude / DisplayWindMax;
			bool visible = speed > 0.01f;
			_windArrows[i].SetActive(visible);
			if (visible)
			{
				var pos = _renderStates[_curRenderState].SurfacePosition[i];
				_windArrows[i].transform.localPosition = pos;
				_windArrows[i].transform.localRotation = Quaternion.LookRotation(-pos, Vector3.Cross(new Vector3(wind.x, wind.y, 0), pos));
				_windArrows[i].transform.GetChild(1).localScale = Vector3.one * math.min(1, speed);
			}
		}
	}
	private void UpdateArrowsCurrent()
	{
		for (int i = 0; i < _windArrows.Length; i++)
		{
			var wind = _renderStates[_curRenderState].Current[i];

			float speed = wind.magnitude / DisplayCurrentMax;
			bool visible = speed > 0.01f && _renderStates[_curRenderState].WaterDepth[i] >= DisplayCurrentMinDepth;
			_windArrows[i].SetActive(visible);
			if (visible)
			{
				var pos = _renderStates[_curRenderState].SurfacePosition[i];
				_windArrows[i].transform.localPosition = pos;
				var forward = new Vector3(wind.x, wind.y, 0);
				_windArrows[i].transform.localRotation = Quaternion.LookRotation(-pos, Vector3.Cross(new Vector3(wind.x, wind.y, 0), pos));
				_windArrows[i].transform.GetChild(1).localScale = Vector3.one * math.min(1, speed);
			}

		}
	}

	#endregion
}
