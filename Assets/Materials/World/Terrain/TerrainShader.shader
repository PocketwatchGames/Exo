Shader "Custom/TerrainShader"
 {
    Properties{
		_UVScale("UV Scale", Range(0, 100)) = 1
		_HeightBlend("Height Blend Threshold", Range(0, 1)) = 0.025
		_ColorBase("Color Base", Color) = (1,1,1,1)
		_Color1("Color 1", Color) = (1,1,1,1)
		_Color2("Color 2", Color) = (1,1,1,1)
		_Color3("Color 3", Color) = (1,1,1,1)
		_Color4("Color 4", Color) = (1,1,1,1)
		_Color5("Color 5", Color) = (1,1,1,1)
		_Color6("Color 6", Color) = (1,1,1,1)
		_GlossBase("Gloss Base", Range(0,1)) = 0.5
		_Gloss1("Gloss 1", Range(0,1)) = 0.5
		_Gloss2("Gloss 2", Range(0,1)) = 0.5
		_Gloss3("Gloss 3", Range(0,1)) = 0.5
		_Gloss4("Gloss 4", Range(0,1)) = 0.5
		_Gloss5("Gloss 5", Range(0,1)) = 0.5
		_Gloss6("Gloss 6", Range(0,1)) = 0.5
		_HeightMap("Height Map", 2D) = "white" {}
	}
    SubShader{
        Tags{ "RenderType" = "Opaque"}
        LOD 200

        CGPROGRAM
        // Physically based Standard lighting model, and enable shadows on all light types
#pragma surface surf Standard fullforwardshadows addshadow

        // Use shader model 3.0 target, to get nicer looking lighting
#pragma target 3.0

        sampler2D _HeightMap;

        struct Input {
            float2 uv_HeightMap;
			float4 height1 : TEXCOORD1;
			float4 height2 : TEXCOORD2;
			float3 worldNormal;
			float3 worldPos;
		};

		fixed4 _ColorBase;
		fixed4 _Color1;
		fixed4 _Color2;
		fixed4 _Color3;
		fixed4 _Color4;
		fixed4 _Color5;
		fixed4 _Color6;

		half _GlossBase;
		half _Gloss1;
		half _Gloss2;
		half _Gloss3;
		half _Gloss4;
		half _Gloss5;
		half _Gloss6;

		float _UVScale;
		float _HeightBlend;

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

		void surf(Input IN, inout SurfaceOutputStandard o)
        {                        
			float4 c;
			float gloss;
			float2 uv = IN.uv_HeightMap;

			float3 localPos = mul(unity_WorldToObject, float4(IN.worldPos,1)).xyz;
			TriplanarUV triUV = GetTriplanarUV(localPos*_UVScale);
			float4 heightX = tex2D(_HeightMap, triUV.x);
			float4 heightY = tex2D(_HeightMap, triUV.y);
			float4 heightZ = tex2D(_HeightMap, triUV.z);


			//float4 height = (heightX + heightY + heightZ) / 3;
			float3 triW = GetTriplanarWeights(IN.worldNormal);
			float4 height = triW.x * heightX + triW.y * heightY + triW.z * heightZ;

			float4 w1 = IN.height1 * height;
			float4 w2 = IN.height2 * height;
			float maxVal = max(w1.x, max(w1.y, max(w1.z, max(w1.w, max(w2.x, max(w2.y, max(w2.z, w2.w)))))));
			float baseVal = 1 - maxVal;
			maxVal = max(maxVal, baseVal);
			float minVal = max(0, maxVal - _HeightBlend);
			w1 = max(0, w1 - minVal);
			w2 = max(0, w2 - minVal);
			baseVal = max(0, baseVal - minVal);
			float totalHeight = w1.x + w1.y + w1.z + w1.w + w2.x + w2.y + w2.z + w2.w + baseVal;
			w1 /= totalHeight;
			w2 /= totalHeight;
			baseVal /= totalHeight;

			c.rgb = _Color1 * w1.x + _Color2 * w1.y + _Color3 * w1.z + _Color4 * w1.w + _Color5 * w2.x + _Color6 * w2.y + _ColorBase * baseVal;
			gloss = _Gloss1 * w1.x + _Gloss2 * w1.y + _Gloss3 * w1.z + _Gloss4 * w1.w + _Gloss5 * w2.x + _Gloss6 * w2.y + _GlossBase * baseVal;

			o.Albedo = c.rgb;
            o.Smoothness = gloss;
		}
        ENDCG
        }
        FallBack "Diffuse"
}