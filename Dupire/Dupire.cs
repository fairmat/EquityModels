using System;
using System.Collections.Generic;
using DVPLI;

namespace Dupire
{
	/// <summary>
	/// Dupire context: created after Parsing
	/// </summary>
	internal class DupireContext
	{
		internal double s0;
		internal IFunction r;
		internal IFunction q;
		internal IFunction localVol;
	}

	/// <summary>
	/// Implements Dupire local volatiltiy model simulation
	/// </summary>
	[Serializable]
    public class Dupire:  IExtensibleProcess, IParsable,IMarkovSimulator
    {
		[SettingDescription("S0")]
		public IModelParameter s0;//scalar
		[SettingDescription("Forwarding Risk Free Rate (Zr)")]
		public IModelParameter r;//1d function
		public IModelParameter q;//1d function
		[SettingDescription("Local Volatility")]
		public IModelParameter localVol;//2d function

		[NonSerialized] private DupireContext context;



        public Dupire ()
        {
        }

		  #region IExtensibleProcess implementation
        /// <summary>
        /// Gets a value indicating whether FullSimulation is implemented, in this case it does
        /// so it always returns true.
        /// </summary>
        public bool ImplementsFullSimulation
        {
            get
            {
                return false;
            }
        }

        /// <summary>
        /// Gets a value indicating whether a Markov based simulation is implemented, in this
        /// case it doesn't so it always returns false.
        /// </summary>
        public bool ImplementsMarkovBasedSimulation
        {
            get
            {
                return true;
            }
        }

        /// <summary>
        /// Gets the ProcessInfo for this plug-in, in this case Historical Simulator.
        /// </summary>
        public ProcessInfo ProcessInfo
        {
            get
            {
                return new ProcessInfo("Dupire Local Volatility Model");
            }
        }

        /// <summary>
        /// Called by Simulator after parse.
        /// Initializes here time-dependant but not state dependent variables.
        /// </summary>
        /// <param name="simulationDates">
        /// The dates at which the process realizations will be requested.
        /// </param>
        public void Setup(double[] simulationDates)
        {
			//todo: creates context class 
        }

        /// <summary>
        /// Gets the information required in order to allow the simulation to run.
        /// </summary>
        public SimulationInfo SimulationInfo
        {
            get
            {
                SimulationInfo simulationInfo = new SimulationInfo();
                simulationInfo.NoiseSize = 1;
                simulationInfo.LatentSize = 0;
                simulationInfo.StateSize = 1;
                return simulationInfo;
            }
        }

        /// <summary>
        /// Creates a list of all the sub-objects that can be edited.
        /// </summary>
        /// <param name="recursive">The parameter is not used.</param>
        /// <returns>
        /// The created list with all the sub objects that can be edited (empty in this case).
        /// </returns>
        public List<IExportable> ExportObjects(bool recursive)
        {
            return new List<IExportable>();
        }
        #endregion // IExtensibleProcess implementation

        #region IParsable implementation
        /// <summary>
        /// Parses the process (in this case nothing has to be done).
        /// </summary>
        /// <param name="context">The project representing the context of the parsing.</param>
        /// <returns>true if the the parsing caused errors; otherwise false.</returns>
        public bool Parse(IProject context)
        {
            return false;
        }
        #endregion // IParsable implementation

		#region IMarkovSimulator implementation
		public unsafe void a (int i, double* x, double* a)
		{
			throw new System.NotImplementedException ();
		}

		public unsafe void b (int i, double* x, double* b)
		{
			throw new System.NotImplementedException ();
		}

		public void isLog (ref bool[] isLog)
		{
			isLog[0]=true;
		}

		public DynamicInfo DynamicInfo {
			get {
				return  new DynamicInfo(true,false,true,true);
			}
		}

		public double[] x0 {
			get {
				return  new double[]{context.s0};
			}
		}
		#endregion
    }
}

