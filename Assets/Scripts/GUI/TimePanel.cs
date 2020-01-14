using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TimePanel : MonoBehaviour
{
	public WorldComponent World;
	public UnityEngine.UI.Text Time;
	public UnityEngine.UI.Text Date;
	public UnityEngine.UI.Text Ticks;
	public UnityEngine.UI.Text TimeScale;

	// Start is called before the first frame update
	void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
		TimeScale.text = "x" + World.TimeScale;

		var rs = World.RenderState;
		Ticks.text = "[" + rs.Ticks + "]";
		Time.text = ((int)(WorldTime.GetTime(rs.Ticks, rs.SpinSpeed) / rs.SpinSpeed)).ToString() + "\"";
		float days = WorldTime.GetDays(rs.Ticks, rs.SpinSpeed);
		Date.text = ((int)(days / rs.SpinSpeed)).ToString() + "\\" + ((int)(days / rs.OrbitSpeed)).ToString("X4");
	}

	public void OnTimescaleChanged(float t)
	{
		World.TimeScale = t;
	}

	public void OnTimeStep()
	{
		World.StepTime();
	}

}
