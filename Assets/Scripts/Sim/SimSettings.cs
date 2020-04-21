using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

[Serializable]
public struct SimSettings {
	public bool MakeAirIncompressible;
	public bool CheckForDegeneracy;
	public bool CollectGlobals;
	public bool CollectOverlay;
	public bool LogState;
	public int LogStateIndex;
}
