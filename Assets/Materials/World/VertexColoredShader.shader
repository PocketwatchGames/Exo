Shader "Custom/VertexColoredShader"
 {
    Properties{
		_Color1("Color Rock", Color) = (1,1,1,1)
		_Color2("Color Sand", Color) = (1,1,1,1)
		_Color3("Color Dirt", Color) = (1,1,1,1)
		_Color4("Color Grass", Color) = (1,1,1,1)
		_Color5("Color Water", Color) = (1,1,1,1)
		_Gloss1("Gloss Rock", Range(0,1)) = 0.5
		_Gloss2("Gloss Sand", Range(0,1)) = 0.5
		_Gloss3("Gloss Dirt", Range(0,1)) = 0.5
		_Gloss4("Gloss Grass", Range(0,1)) = 0.5
		_Gloss5("Gloss Water", Range(0,1)) = 0.5
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
            float4 color : Color;
			float3 worldPos : SV_POSITION;
        };

		fixed4 _Color1;
		fixed4 _Color2;
		fixed4 _Color3;
		fixed4 _Color4;
		fixed4 _Color5;

		half _Gloss1;
		half _Gloss2;
		half _Gloss3;
		half _Gloss4;
		half _Gloss5;

		void surf(Input IN, inout SurfaceOutputStandard o)
        {                        
			float4 c;
			float gloss;
			float2 uv = IN.uv_HeightMap;
			float4 height = tex2D(_HeightMap, uv);
			float height5 = 1 - (IN.color.x + IN.color.y + IN.color.z + IN.color.w);

			float4 w = IN.color * height;
			float maxVal = max(w.x, max(w.y, max(w.z, max(w.w, height5))));
			float minVal = max(0, maxVal - 0.1);
			w = max(0, w- minVal);
			height5 = max(0, height5 - minVal);
			float totalHeight = w.x + w.y + w.z + w.w + height5;
			w /= totalHeight;
			height5 /= totalHeight;

			c.rgb = _Color1 * w.x + _Color2 * w.y + _Color3 * w.z + _Color4 * w.w + _Color5 * height5;
			gloss = _Gloss1 * w.x + _Gloss2 * w.y + _Gloss3 * w.z + _Gloss4 * w.w + _Gloss5 * height5;

			o.Albedo = c.rgb;
            o.Smoothness = gloss;
        }
        ENDCG
        }
        FallBack "Diffuse"
}