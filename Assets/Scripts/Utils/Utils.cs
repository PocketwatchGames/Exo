using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

public static class Utils {
	public static float Sqr(float x) { return x * x; }
	public static int Sqr(int x) { return x * x; }
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