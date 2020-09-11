using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

[Serializable]
public class Species
{
	public string Name;
	public RenderTexture Portrait;
	public Mesh Mesh;
}

public struct SpeciesState {
	public float Population;
}