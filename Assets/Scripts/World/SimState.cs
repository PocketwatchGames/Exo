using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Entities;
using Unity.Collections;
using UnityEngine;

public class SimState {

	public SimState(int count, string name)
	{
		World = new World(name);
		var cells = new NativeArray<Entity>(count, Allocator.Temp);
		var cellArchetype = World.EntityManager.CreateArchetype(typeof(CellComponent), typeof(RenderCellComponent));
		World.EntityManager.CreateEntity(cellArchetype, cells);
		var entities = World.EntityManager.GetAllEntities();

		Cells = new NativeArray<CellComponent>(count, Allocator.Persistent);
		for (int i=0;i<count; i++)
		{
			Cells[i] = World.EntityManager.GetComponentData<CellComponent>(entities[i]);
		}
	}
	public NativeArray<CellComponent> Cells;
	public World World;
	public int Ticks;
	public float Gravity;
	public float SpinAngle;
	public float SpinSpeed;
	public float OrbitSpeed;
	public float TiltAngle;


	public void Update(SimState lastState)
	{
		Ticks = lastState.Ticks+1;
		Gravity = lastState.Gravity;
		SpinAngle = lastState.SpinAngle;
		SpinSpeed = lastState.SpinSpeed;
		OrbitSpeed = lastState.OrbitSpeed;
		TiltAngle = lastState.TiltAngle;

		for (int i=0;i<lastState.Cells.Length;i++)
		{
			Cells[i] = lastState.Cells[i];
		}
		World.Update();
		World.GetExistingSystem<WorldSimSystem>()?.Update();
	}
}


public struct CellComponent : IComponentData {
	public float Elevation;
	public float WaterElevation;
	public float Vegetation;
	public float Ice;
	public float CloudCoverage;
	public float CloudElevation;
	public float RelativeHumidity;
}

public struct RenderCellComponent : IComponentData {
	public Color32 CloudColor;
	public Color32 WaterColor;
	public Color32 TerrainColor;
	public float WaterElevation;
	public float TerrainElevation;
	public float CloudElevation;
}