using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Unity.Mathematics;

static class WorldTime {
	public static float GetTime(int ticks, float spinSpeed)
	{
		return math.modf((float)ticks * spinSpeed, out var i);
	}
	public static float GetDays(int ticks, float spinSpeed)
	{
		return (int)(ticks * spinSpeed);
	}
}
