Shader "Custom/CloudShaderBack"
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
	}
    SubShader
    {

        Tags { "RenderType"="Opaque" }
        LOD 200
		Blend SrcAlpha OneMinusSrcAlpha


		Cull Front
		ZWrite Off

		CGPROGRAM

		// Physically based Standard lighting model, and enable shadows on all light types
		#pragma surface surf StandardSpecular alpha fullforwardshadows addshadow 
		#pragma vertex vert
		// Use shader model 3.0 target, to get nicer looking lighting
		#pragma target 3.0

		#include "CloudShader.cginc"


		void vert(inout appdata_full v) {
			v.normal *= -1;
		}

		void surf(Input IN, inout SurfaceOutputStandardSpecular o)
		{
			standardPass(IN, o);
		}
		ENDCG

    }
    FallBack "Diffuse"
}
