using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;
using Unity.Jobs;

public static class Utils {
	public static float Sqr(float x) { return x * x; }
	public static int Sqr(int x) { return x * x; }

	static public float RepeatExclusive(float x, float y)
	{
		while (x < 0)
		{
			x += y;
		}
		while (x >= y)
		{
			x -= y;
		}
		return x;
	}

	static public float WrapAngle(float x)
	{
		while (x < math.PI)
		{
			x += math.PI * 2;
		}
		while (x >= math.PI)
		{
			x -= math.PI * 2;
		}
		return x;
	}

}

public struct Optional<T> {
	public bool hasValue { get; private set; }
	public T value { get; private set; }

	public Optional(T t) {
		value = t;
		hasValue = true;
	}
	public void Set(T t)
	{
		value = t;
		hasValue = true;
	}
	public void Clear()
	{
		hasValue = false;
	}
}

[Serializable]
public class JobHelper {

	public bool Enabled = true;
	[NonSerialized]
	public bool Async = true;
	private int _cellCount;

	public JobHelper(int cellCount)
	{
		_cellCount = cellCount;
	}

	public static int DefaultBatchCount = 100;
	public JobHandle Run<T>(T job, JobHandle dependences = default(JobHandle)) where T : struct, IJobParallelFor
	{
		return Run(job, Enabled, Async, _cellCount, DefaultBatchCount, dependences);
	}

	public static JobHandle Run<T>(T job, bool enabled, bool async, int count, int batchCount, JobHandle dependences) where T : struct, IJobParallelFor
	{
		if (enabled)
		{
			if (async)
			{
				return job.Schedule(count, batchCount, dependences);
			}
			dependences.Complete();
			job.Run(count);
			return default(JobHandle);
		}
		return dependences;
	}

}