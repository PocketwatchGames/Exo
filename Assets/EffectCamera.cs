using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using UnityEngine.VFX;

public class EffectCamera : MonoBehaviour
{
	public WorldView View;
	public VisualEffect RainEffect;
	public VisualEffect WindEffect;
	public Texture2D PositionRenderTexture;
	public Texture2D _rainRenderTexture;

	[HideInInspector] public Texture2D PixelTexture;

	private bool _positionTextureInitialized;

	// Start is called before the first frame update
	void Start()
    {
		int rtSize = (int)math.pow(2, math.ceil(math.log(math.sqrt(View.Sim.StaticState.Count)) / math.log(2)));
		_rainRenderTexture = new Texture2D(rtSize, rtSize, TextureFormat.ARGB32, 0, true );
		_rainRenderTexture.filterMode = FilterMode.Point;

		RainEffect.SetTexture("RainTexture", _rainRenderTexture);
		RainEffect.SetTexture("PositionTexture", PositionRenderTexture);
		WindEffect.SetTexture("PositionTexture", PositionRenderTexture);
		WindEffect.SetTexture("WindTex", _rainRenderTexture);
		//WindEffect.SetMesh("TerrainMesh", _terrainMesh);


		PixelTexture = new Texture2D(1, 1);
		PixelTexture.SetPixel(0, 0, Color.white);

	}

	// Update is called once per frame
	void Update()
    {
		GL.PushMatrix();
		GL.LoadOrtho();
		GL.sRGBWrite = false;

		int effectTextureWidth = _rainRenderTexture.width;
		if (!_positionTextureInitialized)
		{
			_positionTextureInitialized = true;

//			RenderTexture.active = PositionRenderTexture;
//			float2 screenToRT = new float2(1.0f / 64, 1.0f / 64);
			Color32[] colors = new Color32[PositionRenderTexture.width * PositionRenderTexture.height];
			for (int i = 0; i < PositionRenderTexture.width; i++)
			{
				for (int j = 0; j < PositionRenderTexture.height; j++)
				{
					float2 latlon = new float2((float)i / PositionRenderTexture.width * math.PI * 2, (float)j / PositionRenderTexture.height * math.PI);

					float sinY = math.sin(latlon.y);
					Ray ray = new Ray(Vector3.zero, new Vector3(math.sin(latlon.x) * sinY, math.cos(latlon.y), math.cos(latlon.x) * sinY));
					RaycastHit hitInfo;
					bool hit = View.Sim.Icosphere.MeshCollider.Raycast(ray, out hitInfo, 10);
					Debug.Assert(hit);
					if (hit)
					{
						int tIndex = 2;
						if (hitInfo.barycentricCoordinate.x > hitInfo.barycentricCoordinate.y && hitInfo.barycentricCoordinate.x > hitInfo.barycentricCoordinate.z)
						{
							tIndex = 0;
						}
						else if (hitInfo.barycentricCoordinate.y > hitInfo.barycentricCoordinate.x && hitInfo.barycentricCoordinate.y > hitInfo.barycentricCoordinate.z)
						{
							tIndex = 1;
						}
						int pIndex = View.Sim.Icosphere.Polygons[hitInfo.triangleIndex].m_Vertices[tIndex];

						int row = pIndex / effectTextureWidth;
						float2 pos = new float2(pIndex - row * effectTextureWidth, row);
//						Rect r = new Rect(i * screenToRT.x, j * screenToRT.y, 1 * screenToRT.x, 1 * screenToRT.y);
						colors[i+j*PixelTexture.width] = new Color32((byte)pos.x, (byte)pos.y, 1, 1);
//						Graphics.DrawTexture(r, PixelTexture, new Rect(0, 0, 1, 1), 0, 0, 0, 0, c);
					}
				}
			}
			PixelTexture.SetPixels32(colors);

		}

		{
			//			RenderTexture.active = _rainRenderTexture;
			//			GL.Clear(false, true, new Color(0, 0, 0, 0));
			//			float2 screenToRT = new float2(1.0f / _rainRenderTexture.width, 1.0f / _rainRenderTexture.height);
			Color32[] colors = new Color32[(int)math.ceil((float)View.Sim.StaticState.Count / _rainRenderTexture.width) * _rainRenderTexture.width];
			for (int i = 0; i < View.Sim.StaticState.Count; i++)
			{
				float3 w = View.CurRenderState.VelocityHorizontal[i] / View.DisplayWindSpeedLowerAirMax;
				float speed = math.length(w);
				if (speed > 1)
				{
					w /= speed;
				}
				w = w * 0.5f + 0.5f;
				int row = i / effectTextureWidth;
//				Rect r = new Rect((i - row * effectTextureWidth) * screenToRT.x, row * screenToRT.y, 1 * screenToRT.x, 1 * screenToRT.y);
				colors[i] = new Color(w.x, w.y, w.z, 1 /* math.saturate(DisplayState.Rainfall[i] / DisplayRainfallMax)*/);
//				Graphics.DrawTexture(r, PixelTexture, new Rect(0, 0, 1, 1), 0, 0, 0, 0, c);
			}
			_rainRenderTexture.SetPixels32(0, 0, _rainRenderTexture.width, (int)math.ceil((float)colors.Length / _rainRenderTexture.width), colors, 0);
		}

		GL.PopMatrix();
		RenderTexture.active = null;
		GL.sRGBWrite = true;
	}
}
