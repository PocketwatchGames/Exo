Shader "Custom/OverlayShader"
 {
    Properties{
    }
    SubShader{
        Tags{ "RenderType" = "Opaque"}
        LOD 200
		Lighting Off

        CGPROGRAM
        // Physically based Standard lighting model, and enable shadows on all light types
#pragma surface surf Unlit

        // Use shader model 3.0 target, to get nicer looking lighting
#pragma target 3.0

		half4 LightingUnlit(SurfaceOutput s, fixed3 lightDir, fixed atten) {
		   return fixed4(0,0,0,0);
		 }

		struct Input {
			float4 color : Color;
		};

		void surf(Input IN, inout SurfaceOutput o)
        {                        
            o.Albedo = IN.color.rgb;
        }
        ENDCG
        }
        FallBack "Diffuse"
}