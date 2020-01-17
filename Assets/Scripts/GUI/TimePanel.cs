using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TimePanel : MonoBehaviour
{
	public WorldSim World;
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
		Ticks.text = "[" + state.PlanetState.Ticks + "]";
		Time.text = ((int)(WorldTime.GetTime(state.PlanetState.Ticks, state.PlanetState.SpinSpeed) / state.PlanetState.SpinSpeed)).ToString() + "\"";
		float days = WorldTime.GetDays(state.PlanetState.Ticks, state.PlanetState.SpinSpeed);
		Date.text = ((int)days).ToString() + "\\" + ((int)(state.PlanetState.Ticks * state.PlanetState.OrbitSpeed)).ToString("X4");
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
