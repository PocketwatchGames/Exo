Shader "Custom/CloudShaderFront"
{
    Properties
    {
		_MainTex("Albedo (RGB)", 2D) = "white" {}
		_NormalMap("Normal Map", 2D) = "bump" {}
		_AlphaTest("Alpha Test", Range(0,1)) = 0.1
		_HardEdge("Hard Edges", Range(0,100)) = 10
		_Alpha("Alpha", Range(0,1)) = 0.5
		_TextureScale("Texture Scale", Range(0,10)) = 0.1
		_TimeScale("Time Scale", Range(0,0.1)) = 0.001
		_ViewAlphaMin("XRay Alpha Min", Range(0, 1)) = 0.075
		_ViewDistMin("XRay Alpha Dist", float) = 2
		_ViewDistScale("XRay Alpha Dist Scale", float) = 0.2
	}
	SubShader
	{

		Tags { "RenderType" = "Opaque" }
		LOD 200
		Blend SrcAlpha OneMinusSrcAlpha

		Pass{
			ZWrite On
			ColorMask 0
		}


		CGPROGRAM
		// Physically based Standard lighting model, and enable shadows on all light types
		#pragma surface surf StandardSpecular alpha fullforwardshadows addshadow 

		// Use shader model 3.0 target, to get nicer looking lighting
		#pragma target 3.0

		#include "CloudShader.cginc"


		float _ViewAlphaMin;
		float _ViewDistMin;
		float _ViewDistScale;

		void surf(Input IN, inout SurfaceOutputStandardSpecular o)
		{
			standardPass(IN, o);
			float3 forward = mul(unity_CameraToWorld, float3(0, 0, 1));

			float dist = distance(IN.worldPos, _WorldSpaceCameraPos);
			o.Alpha *= max(_ViewAlphaMin,saturate(1 + dot(WorldNormalVector(IN, o.Normal), forward) / (max(0, dist + _ViewDistMin) * _ViewDistScale)));

		}
		ENDCG

    }
    FallBack "Diffuse"
}
