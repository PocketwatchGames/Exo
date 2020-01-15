using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;
using Unity.Mathematics;
using Unity.Entities;
using Unity;
using Unity.Burst;
using Unity.Jobs;
using Unity.Collections;

public class RenderState {
	public Vector3[] terrainVertices;
	public Vector3[] terrainNormals;
	public Color32[] terrainColors;
	public Vector3[] waterVertices;
	public Vector3[] waterNormals;
	public Color32[] waterColors;
	public Vector3[] cloudVertices;
	public Vector3[] cloudNormals;
	public Color32[] cloudColors;
	public float SpinAngle;
	public float TiltAngle;
	public float Scale;

	public RenderState(int verts)
	{
		terrainVertices = new Vector3[verts];
		terrainNormals = new Vector3[verts];
		terrainColors = new Color32[verts];
		waterVertices = new Vector3[verts];
		waterNormals = new Vector3[verts];
		waterColors = new Color32[verts];
		cloudVertices = new Vector3[verts];
		cloudNormals = new Vector3[verts];
		cloudColors = new Color32[verts];
	}
}

public class SimToRenderCellSystem : ComponentSystem {

	public List<Vector3> SphereVertices;
	public StaticState StaticState;
	public RenderState RenderState;
	public float Scale;

	static Color32 green = new Color32(20, 255, 30, 255);
	static Color32 brown = new Color32(220, 150, 70, 255);
	static Color32 blue = new Color32(0, 0, 255, 255);
	static Color32 white = new Color32(255, 255, 255, 255);
	static Color32 black = new Color32(0, 0, 0, 255);

	private Color32 GetTerrainColor(StaticState staticState, CellComponent cell)
	{
		return Color32.Lerp(brown, green, cell.Vegetation);
	}

	private Color32 GetWaterColor(StaticState staticState, CellComponent cell)
	{
		return Color32.Lerp(blue, white, cell.Ice);
	}

	private Color32 GetCloudColor(StaticState staticState, CellComponent cell)
	{
		var humidity = cell.RelativeHumidity;
		var c = Color32.Lerp(white, black, humidity);
		float opacity = humidity > 0.5f ? math.pow(humidity, 0.25f) * 0.75f : math.pow(humidity, 4);
		c.a = (byte)(255 * opacity);
		return c;
	}


	[BurstCompile]
	protected override void OnUpdate()
	{
		int i = 0;
		Entities.ForEach((ref CellComponent cell) => {
			RenderState.terrainVertices[i] = SphereVertices[i] * (cell.Elevation + StaticState.Radius) * Scale;
			RenderState.terrainColors[i] = GetTerrainColor(StaticState, cell);
			RenderState.terrainNormals[i] = SphereVertices[i];

			RenderState.waterVertices[i] = SphereVertices[i] * (cell.WaterElevation + StaticState.Radius) * Scale;
			RenderState.waterColors[i] = GetWaterColor(StaticState, cell);
			RenderState.waterNormals[i] = SphereVertices[i];

			RenderState.cloudVertices[i] = SphereVertices[i] * (cell.CloudElevation + StaticState.Radius) * Scale;
			RenderState.cloudColors[i] = GetCloudColor(StaticState, cell);
			RenderState.cloudNormals[i] = SphereVertices[i];

			i++;
		});
	}
}


public struct WorldRenderer {

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

	public MeshOverlay ActiveMeshOverlay;
	public WindOverlay ActiveWindOverlay;
	public Icosphere Icosphere;


	public void Init(int subdivisions)
	{
		Icosphere = new Icosphere(subdivisions);
	}


	public void InitMesh(Mesh terrainMesh, Mesh waterMesh, Mesh cloudMesh)
	{

		int vertexCount = Icosphere.Polygons.Count * 3;

		int[] indices = new int[vertexCount];

		for (int i = 0; i < Icosphere.Polygons.Count; i++)
		{
			var poly = Icosphere.Polygons[i];

			indices[i * 3 + 0] = poly.m_Vertices[0];
			indices[i * 3 + 1] = poly.m_Vertices[1];
			indices[i * 3 + 2] = poly.m_Vertices[2];
		}

		terrainMesh.SetTriangles(indices, 0);
		waterMesh.SetTriangles(indices, 0);
		cloudMesh.SetTriangles(indices, 0);

	}


	public void LerpRenderState(RenderState lastState, RenderState nextState, float t, RenderState state)
	{
		for (int i = 0; i < state.cloudColors.Length; i++)
		{
			state.cloudColors[i] = Color32.Lerp(lastState.cloudColors[i], nextState.cloudColors[i], t);
			state.terrainColors[i] = Color32.Lerp(lastState.terrainColors[i], nextState.terrainColors[i], t);
			state.waterColors[i] = Color32.Lerp(lastState.waterColors[i], nextState.waterColors[i], t);
			state.cloudNormals[i] = Vector3.Lerp(lastState.cloudNormals[i], nextState.cloudNormals[i], t);
			state.terrainNormals[i] = Vector3.Lerp(lastState.terrainNormals[i], nextState.terrainNormals[i], t);
			state.waterNormals[i] = Vector3.Lerp(lastState.waterNormals[i], nextState.waterNormals[i], t);
			state.cloudVertices[i] = Vector3.Lerp(lastState.cloudVertices[i], nextState.cloudVertices[i], t);
			state.terrainVertices[i] = Vector3.Lerp(lastState.terrainVertices[i], nextState.terrainVertices[i], t);
			state.waterVertices[i] = Vector3.Lerp(lastState.waterVertices[i], nextState.waterVertices[i], t);
		}
	}

	public void BuildRenderState(StaticState staticState, SimState simState, RenderState renderState, float scale)
	{
		renderState.SpinAngle = simState.SpinAngle;
		renderState.TiltAngle = simState.TiltAngle;

		var s = simState.World.GetExistingSystem<SimToRenderCellSystem>();
		s.RenderState = renderState;
		s.StaticState = staticState;
		s.Scale = scale;
		s.SphereVertices = Icosphere.Vertices;
		s.Update();
	}

	public void UpdateMesh(Mesh terrainMesh, Mesh waterMesh, Mesh cloudMesh, StaticState staticState, RenderState renderstate, float scale)
	{
		terrainMesh.vertices = renderstate.terrainVertices;
		terrainMesh.normals = renderstate.terrainNormals;
		terrainMesh.colors32 = renderstate.terrainColors;

		waterMesh.vertices = renderstate.waterVertices;
		waterMesh.normals = renderstate.waterNormals;
		waterMesh.colors32 = renderstate.waterColors;

		cloudMesh.vertices = renderstate.cloudVertices;
		cloudMesh.normals = renderstate.cloudNormals;
		cloudMesh.colors32 = renderstate.cloudColors;

		terrainMesh.RecalculateBounds();
		waterMesh.RecalculateBounds();
		cloudMesh.RecalculateBounds();

	}

}
