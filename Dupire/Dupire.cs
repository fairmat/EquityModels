/* Copyright (C) 2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Enrico Degiuli (enrico.degiuli@fairmat.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

using System;
using System.Collections.Generic;
using DVPLDOM;
using DVPLI;

namespace Dupire
{
    /// <summary>
    /// Dupire context: created after Parsing.
    /// </summary>
    internal class DupireContext
    {
        internal double s0;
        internal IFunction r;
        internal IFunction q;
        internal IFunction localVol;
    }

    /// <summary>
    /// Implements Dupire local volatility model simulation.
    /// </summary>
    [Serializable]
    public class DupireProcess : IExtensibleProcess, IParsable, IMarkovSimulator, IEstimationResultPopulable
    {
        [SettingDescription("S0")]
        public IModelParameter s0;//scalar
        [SettingDescription("Time Dependent Risk Free Rate (Zero Rate)")]
        public IModelParameter r;//1d function
        [SettingDescription("Time Dependent Continuous Dividend Yield")]
        public IModelParameter q;//1d function
        [SettingDescription("Local Volatility")]
        public IModelParameter localVol;//2d function

        [NonSerialized]
        private DupireContext context;
        [NonSerialized]
        private double[] mu;
        [NonSerialized]
        private double[] simDates;

        public DupireProcess()
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
                return new ProcessInfo("Dupire Local Volatility");
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
            mu = new double[simulationDates.Length];
            double[] zr = new double[simulationDates.Length];
            for (int i = 0; i < simulationDates.Length; i++)
                zr[i] = ZR(simulationDates[i]);
            for (int i = 0; i < simulationDates.Length - 1; i++)
                mu[i] = (zr[i + 1] * simulationDates[i + 1] - zr[i] * simulationDates[i])
                    / (simulationDates[i + 1] - simulationDates[i]);
            mu[simulationDates.Length - 1] = mu[simulationDates.Length - 2];
            this.simDates = simulationDates;
            // la superficie della local vol va precalcolata?
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
            bool errors = false;
            errors = s0.Parse(context);
            errors = BoolHelper.AddBool(errors, q.Parse(context));
            errors = BoolHelper.AddBool(errors, r.Parse(context));
            errors = BoolHelper.AddBool(errors, localVol.Parse(context));

            this.context = new DupireContext();
            this.context.s0 = s0.fV();
            this.context.q = q.fVRef() as IFunction;
            this.context.r = r.fVRef() as IFunction;
            this.context.localVol = localVol.fVRef() as IFunction;

            return errors;
        }
        #endregion // IParsable implementation

        #region IMarkovSimulator implementation
        public unsafe void a(int i, double* x, double* a)
        {
            a[0] = mu[i] - 0.5 * System.Math.Pow(context.localVol.Evaluate(simDates[i], x[0]), 2.0);
        }

        public unsafe void b(int i, double* x, double* b)
        {
            b[0] = context.localVol.Evaluate(simDates[i], x[0]);
        }

        /// <summary>
        /// Sets the passed array with a Boolean stating if the process
        /// must be simulated as a log-normal process. In this case it's
        /// a log-normal process so the first component is set to true.
        /// </summary>
        /// <param name="isLog">
        /// A reference to the array to be set with the required information.
        /// </param>
        public void isLog(ref bool[] isLog)
        {
            isLog[0] = true;
        }

        /// <summary>
        /// Gets details about the structure of the functions A and B of the Markov
        /// process.
        /// In this case drift is time dependant and not state dependant and
        /// volatility is state and time dependent.
        /// </summary>
        public DynamicInfo DynamicInfo
        {
            get
            {
                return new DynamicInfo(true, false, true, true);
            }
        }

        /// <summary>
        /// Gets the starting point for the process.
        /// </summary>
        public double[] x0
        {
            get
            {
                return new double[] { context.s0 };
            }
        }
        #endregion

        /// <summary>
        /// Helper function to make functions easier to read.
        /// Just returns the value of the zero rate at position t.
        /// </summary>
        /// <param name="t">The position where to get the value of the zero rate from.</param>
        /// <returns>The value of the zero rate at position t.</returns>
        private double ZR(double t)
        {
            return (context.r.Evaluate(t) - context.q.Evaluate(t));
        }

        #region IEstimationResultPopulable implementation
        void IEstimationResultPopulable.Populate(IStochasticProcess Container, EstimationResult Estimate)
        {
            bool found;
            s0 = new ModelParameter(PopulateHelper.GetValue("S0", Estimate.Names, Estimate.Values, out found), s0.Description);
            PFunction rFunc = Estimate.Objects[0] as PFunction;
            PFunction rFuncDest = r.fVRef() as PFunction;
            rFuncDest.Expr = rFunc.Expr;

            PFunction qFunc = Estimate.Objects[0] as PFunction;
            PFunction qFuncDest = q.fVRef() as PFunction;
            qFuncDest.Expr = qFunc.Expr;

            PFunction2D.PFunction2D localVolSrc = Estimate.Objects[2] as PFunction2D.PFunction2D;
            PFunction2D.PFunction2D localVolDest = localVol.fVRef() as PFunction2D.PFunction2D;
            localVolDest.Expr = localVolSrc.Expr;
            localVolDest.Interpolation = localVolSrc.Interpolation;
        }
        #endregion
    }
}