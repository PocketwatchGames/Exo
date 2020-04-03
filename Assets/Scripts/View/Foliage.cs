using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using System.Xml.Serialization;

public class Foliage : MonoBehaviour
{
	[XmlIgnore]
	public float3 Position;
	private float3 _baseScale;
	public int CellIndex;
	private float _curScale;
	private float _growthSpeed;
	private float _growthDelay;
	const float treeGrowthSpeed = 4;
	const float treeGrowthSpeedRange = 4;
	const float treeGrowthDelayRange = 1.5f;

	public void InitPosition(GameObject parent, int cellIndex, int treeIndex, float3 pos, float elevation, float3 baseScale, float scale, float perturbMax, Unity.Mathematics.Random random)
	{
		CellIndex = cellIndex;
		_baseScale = baseScale;
		_curScale = scale;
		float3 forward = math.cross(pos, new float3(0, 1, 0));
		float3 right = math.cross(forward, pos);
		float2 perturb = perturbMax * (new float2(random.NextFloat(), random.NextFloat()) * 2 - 1);		

		Position = pos + perturb.x * right + perturb.y * forward;
		var rot = Quaternion.FromToRotation(Vector3.up, pos) * Quaternion.AngleAxis(random.NextFloat(360), Vector3.up);
		transform.localPosition = Position * elevation;
		transform.localScale = baseScale * scale;
		transform.SetParent(parent.transform, false);
		transform.localRotation = rot;
		_growthSpeed = treeGrowthSpeed + random.NextFloat() * treeGrowthSpeedRange;
		_growthDelay = random.NextFloat() * treeGrowthDelayRange;
	}

	public bool UpdatePosition(float elevation, float scale)
	{
		if (_growthDelay > 0)
		{
			_growthDelay -= Time.deltaTime;
			return false;
		}
		else
		{
			_curScale += (scale - _curScale) * math.min(1, Time.deltaTime * _growthSpeed);
			transform.localPosition = Position * elevation;
			transform.localScale = _baseScale * _curScale;
			return (scale < _curScale && _curScale < 0.01f);
		}
	}

}
