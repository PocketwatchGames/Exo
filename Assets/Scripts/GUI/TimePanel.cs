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

		var state = World.ActiveSimState;
		Ticks.text = "[" + state.Ticks + "]";
		Time.text = ((int)(WorldTime.GetTime(state.Ticks, state.SpinSpeed) / state.SpinSpeed)).ToString() + "\"";
		float days = WorldTime.GetDays(state.Ticks, state.SpinSpeed);
		Date.text = ((int)days).ToString() + "\\" + ((int)(state.Ticks * state.OrbitSpeed)).ToString("X4");
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
