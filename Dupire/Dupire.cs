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
    public class DupireProcess : IExtensibleProcess, IParsable,
                                 IMarkovSimulator, IEstimationResultPopulable
    {
        #region Serialized Parameters

        [SettingDescription("S0")]
        public IModelParameter s0; // scalar

        /// <summary>
        /// The time Dependent Risk Free Rate (Zero Rate) (1D Function).
        /// </summary>
        [SettingDescription("Time Dependent Risk Free Rate (Zero Rate)")] 
        [ExternalSymbolReference("ZR", typeof(PFunction))]
        public IModelParameter r;

        /// <summary>
        /// The time Dependent Continuous Dividend Yield (1D Function).
        /// </summary>
        [SettingDescription("Time Dependent Continuous Dividend Yield")]
        [ExternalSymbolReference("DividendYield", typeof(PFunction))]
        public IModelParameter q;

        /// <summary>
        /// The Local volatility (2D Function).
        /// </summary>
        [SettingDescription("Local Volatility")] 
        [ExternalSymbolReference("LocalVolatility", typeof(PFunction2D.PFunction2D))]
        public IModelParameter localVol;

        #endregion Serialized Parameters

        [NonSerialized]
        private DupireContext context;
        /// <summary>
        /// Temporary calculated Function fitting the zero rate values.
        /// </summary>
        [NonSerialized]
        private Function zrCurve;

        /// <summary>
        /// Temporary calculated Function fitting the dividend yield values.
        /// </summary>
        [NonSerialized]
        private Function dyCurve;
        [NonSerialized]
        private double[] mu;
        [NonSerialized]
        private double[] simDates;

        public DupireProcess()
        {
            this.s0 = new ModelParameter(string.Empty, "S0");
            this.r = new ModelParameter(string.Empty, "Time Dependent Risk Free Rate (Zero Rate)");
            this.q = new ModelParameter(string.Empty, "Time Dependent Continuous Dividend Yield");
            this.localVol = new ModelParameter(string.Empty, "Local Volatility");
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
            // todo: creates context class
            this.mu = new double[simulationDates.Length];
            Vector time = new Vector(1);
            for (int i = 0; i < simulationDates.Length; i++)
            {
                time[0] = simulationDates[i];
                this.mu[i] =  this.zrCurve.Partial(time, 0) * time[0] + this.zrCurve.Evaluate(time) - this.dyCurve.Evaluate(time);
            }
            this.simDates = simulationDates;
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
                simulationInfo.DefaultComponent = -1;
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
            errors = this.s0.Parse(context);
            errors = BoolHelper.AddBool(errors, this.q.Parse(context));
            errors = BoolHelper.AddBool(errors, this.r.Parse(context));
            errors = BoolHelper.AddBool(errors, this.localVol.Parse(context));

            this.context = new DupireContext();
            this.context.s0 = this.s0.fV();
            this.context.q = this.q.fVRef() as IFunction;
            this.context.r = this.r.fVRef() as IFunction;
            this.context.localVol = this.localVol.fVRef() as IFunction;

            return RetrieveCurve(context, errors);
        }

        /// <summary>
        /// Retrieves zr and dy curve from the model.
        /// </summary>
        /// <param name="p_Context">The parameter is not used.</param>
        /// <param name="errors">
        /// The current status of errors, if errors happened previously to this call.
        /// </param>
        /// <returns>True if there have been errors during this call or before.</returns>
        private bool RetrieveCurve(IProject p_Context, bool errors)
        {
            // check the presence of a zero rate curve
            object zrref = Engine.Parser.EvaluateAsReference(this.r.Expression);
            if (!Engine.Parser.GetParserError())
            {
                this.zrCurve = zrref as Function;
                if (this.zrCurve == null)
                {
                    errors = true;
                    p_Context.AddError("Cannot find the zero rate curve " +
                                       this.r.Expression);
                }
            }
            else
                errors = true;

            // check the presence of a dividend yield curve
            object dyref = Engine.Parser.EvaluateAsReference(this.q.Expression);
            if (!Engine.Parser.GetParserError())
            {
                this.dyCurve = dyref as Function;
                if (this.dyCurve == null)
                {
                    errors = true;
                    p_Context.AddError("Cannot find the dividend yield curve " +
                                       this.q.Expression);
                }
            }
            else
                errors = true;

            return errors;
        }

        #endregion // IParsable implementation

        #region IMarkovSimulator implementation
        public unsafe void a(int i, double* x, double* a)
        {
            a[0] = this.mu[i] - 0.5 * System.Math.Pow(this.context.localVol.Evaluate(this.simDates[i], x[0]), 2.0);
        }

        public unsafe void b(int i, double* x, double* b)
        {
            b[0] = this.context.localVol.Evaluate(this.simDates[i], x[0]);
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
                return new double[] { this.context.s0 };
            }
        }
        #endregion

        #region IEstimationResultPopulable implementation

        public void Populate(IStochasticProcess container, EstimationResult estimate)
        {
            bool found;
            this.s0 = new ModelParameter(PopulateHelper.GetValue("S0", estimate.Names, estimate.Values, out found), this.s0.Description);

            bool errors = RetrieveCurve(container.Context, false);
            if (!errors)
            {
                PFunction rFunc = estimate.Objects[0] as PFunction;
                PFunction rFuncDest = this.r.fVRef() as PFunction;
                rFuncDest.Expr = rFunc.Expr;

                PFunction qFunc = estimate.Objects[1] as PFunction;
                PFunction qFuncDest = this.q.fVRef() as PFunction;
                qFuncDest.Expr = qFunc.Expr;
                //Calibrator assumes dividend yield is a step constant function, the simulation model must be coherent with that assumption. 
                qFuncDest.m_Function.iType = DVPLUtils.EInterpolationType.ZERO_ORDER_LEFT;


                PFunction2D.PFunction2D localVolSrc = estimate.Objects[2] as PFunction2D.PFunction2D;
                PFunction2D.PFunction2D localVolDest = this.localVol.fVRef() as PFunction2D.PFunction2D;
                localVolDest.Expr = localVolSrc.Expr;
                localVolDest.Interpolation = localVolSrc.Interpolation;
            }
        }

        #endregion
    }
}
