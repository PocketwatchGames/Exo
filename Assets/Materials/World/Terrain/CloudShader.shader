Shader "Custom/CloudShader"
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

		CGINCLUDE


		sampler2D _MainTex;
		sampler2D _NormalMap;
		float _AlphaTest;
		float _HardEdge;
		float _Alpha;
		float _TextureScale;
		float _TimeScale;
		uniform float _GameTime;

		struct Input
		{
			float2 uv_MainTex;
			float4 color : Color;
			float3 worldPos;
			float3 worldNormal;
		};


		struct TriplanarUV {
			float2 x, y, z;
		};

		TriplanarUV GetTriplanarUV(float3 p) {
			TriplanarUV triUV;
			triUV.x = p.zy;
			triUV.y = p.xz;
			triUV.z = p.xy;
			return triUV;
		}

		float3 GetTriplanarWeights(float3 normal) {
			float3 triW = abs(normal);
			return triW / (triW.x + triW.y + triW.z);
		}

		void standardPass(Input IN, inout SurfaceOutputStandardSpecular o)
		{
			//TriplanarUV triUV1 = GetTriplanarUV(localPos * _TextureScale + _GameTime * _TimeScale);
			//float4 texX1 = tex2D(_MainTex, triUV1.x);
			//float4 texY1 = tex2D(_MainTex, triUV1.y);
			//float4 texZ1 = tex2D(_MainTex, triUV1.z);
			//float3 triW = GetTriplanarWeights(IN.worldNormal);
			//float4 tex = triW.x * texX + triW.y * texY + triW.z * texZ;
			//float4 tex1 = (texX1 + texY1 + texZ1) / 3;

			float3 localPos = mul(unity_WorldToObject, float4(IN.worldPos, 1)).xyz;

			TriplanarUV triUV1 = GetTriplanarUV(localPos * _TextureScale + _GameTime * _TimeScale);
			float4 texX1 = tex2D(_MainTex, triUV1.x);
			float4 texY1 = tex2D(_MainTex, triUV1.y);
			float4 texZ1 = tex2D(_MainTex, triUV1.z);

			TriplanarUV triUV2 = GetTriplanarUV(localPos * _TextureScale - _GameTime * _TimeScale + 0.4);
			float4 texX2 = tex2D(_MainTex, triUV2.x);
			float4 texY2 = tex2D(_MainTex, triUV2.y);
			float4 texZ2 = tex2D(_MainTex, triUV2.z);
			float4 tex = (texX1 + texY1 + texZ1 + texX2 + texY2 + texZ2) / 6;

			float4 normalX1 = tex2D(_NormalMap, triUV1.x);
			float4 normalY1 = tex2D(_NormalMap, triUV1.y);
			float4 normalZ1 = tex2D(_NormalMap, triUV1.z);
			float4 normalX2 = tex2D(_NormalMap, triUV2.x);
			float4 normalY2 = tex2D(_NormalMap, triUV2.y);
			float4 normalZ2 = tex2D(_NormalMap, triUV2.z);
			float4 normal = (normalX1 + normalY1 + normalZ1 + normalX2 + normalY2 + normalZ2) / 6;

			float dist = distance(IN.worldPos, _WorldSpaceCameraPos);
			float alpha = tex.x * IN.color.a;
			alpha = saturate((alpha - _AlphaTest) * _HardEdge) * _Alpha;
			clip(alpha);

			o.Albedo = 1;
			o.Alpha = saturate(alpha);
			o.Specular *= _Alpha;
			o.Normal = UnpackNormal(normal);
		}

		ENDCG


        Tags { "RenderType"="Transparent" "Queue"="Transparent" }
        LOD 200
		Blend SrcAlpha OneMinusSrcAlpha

		Cull Front

        CGPROGRAM
        // Physically based Standard lighting model, and enable shadows on all light types
		#pragma surface surf StandardSpecular alpha
		#pragma vertex vert
        // Use shader model 3.0 target, to get nicer looking lighting
        #pragma target 3.0

		void vert(inout appdata_full v) {
			v.normal *= -1;
		}

        void surf (Input IN, inout SurfaceOutputStandardSpecular o)
        {
			standardPass(IN, o);
		}
        ENDCG

		Cull Back
		Pass{
			ZWrite On
			ColorMask 0
		}


		CGPROGRAM
		// Physically based Standard lighting model, and enable shadows on all light types
		#pragma surface surf StandardSpecular alpha

		// Use shader model 3.0 target, to get nicer looking lighting
		#pragma target 3.0

		void surf(Input IN, inout SurfaceOutputStandardSpecular o)
		{
			standardPass(IN, o);
		}
		ENDCG

    }
    FallBack "Diffuse"
}
