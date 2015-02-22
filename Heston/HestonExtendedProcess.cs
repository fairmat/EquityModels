/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Enrico Degiuli (enrico.degiuli@fairmat.com)
 *            Matteo Tesser (matteo.tesser@fairmat.com)
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
using System.Runtime.Serialization;
using DVPLDOM;
using DVPLI;

namespace HestonExtended
{
    /// <summary>
    /// Implementation of the Heston Stochastic Process with time varying drift.
    /// </summary>
    [Serializable]
    public unsafe class HestonExtendedProcess : IExtensibleProcess, IMarkovSimulator, IParsable,
                                                IEstimationResultPopulable, IGreeksDerivativesInfo,
                                                IExportableContainer
    {
        #region Serialized Parameters
        /// <summary>
        /// Previously referenced the zr,
        /// now it's just kept for compatibility.
        /// </summary>
        public IModelParameter MuReference;

        /// <summary>
        /// Reference to the zero rate curve.
        /// </summary>
        [OptionalField(VersionAdded = 2)]
        [ExternalSymbolReference("RateCurve", typeof(PFunction))]
        public IModelParameter zrReference;

        /// <summary>
        /// Reference to the dividend yield curve.
        /// </summary>
        [OptionalField(VersionAdded = 2)]
        [ExternalSymbolReference("DividendCurve", typeof(PFunction))]
        public IModelParameter dyReference;

        /// <summary>
        /// The speed of mean reversion of the process.
        /// </summary>
        public IModelParameter k;

        /// <summary>
        /// The mean reversion level of the process.
        /// The volatility process will reach this value as time goes to infinity.
        /// </summary>
        public IModelParameter theta;

        /// <summary>
        /// The "volatility of volatility" regulates the variance of the volatility process.
        /// </summary>
        public IModelParameter sigma;

        /// <summary>
        /// Starting value for the stock process.
        /// </summary>
        public IModelParameter S0;

        /// <summary>
        /// Starting value for the volatility process.
        /// </summary>
        public IModelParameter V0;

        #endregion Serialized Parameters

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

        /// <summary>
        /// Temporary drift values, used to optimize the simulation.
        /// </summary>
        [NonSerialized]
        private double[] mu;

        /// <summary>
        /// A string representing this process to the user.
        /// </summary>
        internal static string extendedHestonDescription = "Heston (with time dependent drift)";

        /// <summary>
        /// Initializes the zr and dy reference to empty values, in case they are not initialized.
        /// </summary>
        private void InitOptionalFields()
        {
            if (this.zrReference == null)
                this.zrReference = new ModelParameter(string.Empty, "zero rate curve");

            if (this.dyReference == null)
                this.dyReference = new ModelParameter(string.Empty, "dividend yield curve");
        }

        /// <summary>
        /// Initializes a new instance of the HestonExtendedProcess class.
        /// This is the default constructor setting all parameters to zero/empty.
        /// </summary>
        public HestonExtendedProcess()
        {
            InitOptionalFields();
            this.k = new ModelParameter(0.0, "k");
            this.theta = new ModelParameter(0.0, "theta");
            this.sigma = new ModelParameter(0.0, "sigma");
            this.S0 = new ModelParameter(0.0, "S0");
            this.V0 = new ModelParameter(0.0, "V0");
        }

        /// <summary>
        /// Initializes optional fields after deserialization.
        /// </summary>
        /// <param name='context'>
        /// The parameter is not used.
        /// </param>
        [OnDeserialized]
        public void OnDeserialized(StreamingContext context)
        {
            // Setup non initialized fields, if needed.
            InitOptionalFields();
        }

        #region IGreeksDerivativesInfo

        /// <summary>
        /// Gets the factors for Delta Greek derivative.
        /// </summary>
        /// <returns>
        /// An array of model parameters containing the starting
        /// value for the stock process.
        /// </returns>
        public IModelParameter[] GetDeltaFactors()
        {
            return new IModelParameter[] { this.S0 };
        }

        /// <summary>
        /// Gets the factors for Vega Greek derivative.
        /// </summary>
        /// <returns>
        /// An array of model parameters containing the starting
        /// value for the volatility process and the mean reversion level.
        /// </returns>
        public IModelParameter[] GetVegaFactors()
        {
            return new IModelParameter[] { this.V0, this.theta };
        }

        #endregion IGreeksDerivativesInfo

        /// <summary>
        /// Sets some default values for the parameters.
        /// * k = 2.5
        /// * theta = 0.4
        /// * sigma = 0.2
        /// * S0 = 1
        /// * V0 = 0.3.
        /// </summary>
        public void DefaultInstance()
        {
            this.k = new ModelParameter(2.5, "k");
            this.theta = new ModelParameter(0.4, "theta");
            this.sigma = new ModelParameter(0.2, "sigma");
            this.S0 = new ModelParameter(1, "S0");
            this.V0 = new ModelParameter(0.3, "V0");
        }

        #region IExtensibleProcess Members

        /// <summary>
        /// Gets a value indicating whether FullSimulation is implemented, in this
        /// case it doesn't so it always returns false.
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
        /// case it does so it always returns true.
        /// </summary>
        public bool ImplementsMarkovBasedSimulation
        {
            get
            {
                return true;
            }
        }

        /// <summary>
        /// Gets the ProcessInfo for this plugin, in this case the text contained in
        /// <see cref="extendedHestonDescription"/>.
        /// </summary>
        public ProcessInfo ProcessInfo
        {
            get
            {
                return new ProcessInfo(extendedHestonDescription);
            }
        }

        /// <summary>
        /// Called by Simulator after parse.
        /// Initializes here time-dependant but not state dependent variables.
        /// </summary>
        /// <param name='simulationDates'>
        /// The dates at which the process realizations will be requested.
        /// </param>
        public void Setup(double[] simulationDates)
        {
            // Precalculates drift values
            this.mu = new double[simulationDates.Length];
            Vector time = new Vector(1);
            for (int i = 0; i < simulationDates.Length; i++)
            {
                time[0] = simulationDates[i];
                //Istantaneous growth rate is d/dt[z(t)*t] hence z(t)+d/dt[z(t)]*t
                this.mu[i] = this.zrCurve.Partial(time, 0) * time[0] + this.zrCurve.Evaluate(time) - this.dyCurve.Evaluate(time);
            }
        }

        /// <summary>
        /// Gets the information required for the simulator in order to control
        /// the plug-in through the <see cref="IMarkovSimulator"/> interface.
        /// The Heston Extended Process has:
        /// * 1 latent components in the state components.
        /// * 2 components of noise.
        /// * 2 state components.
        /// </summary>
        public SimulationInfo SimulationInfo
        {
            get
            {
                SimulationInfo si = new SimulationInfo();
                si.LatentSize = 1;
                si.NoiseSize = 2;
                si.StateSize = 2;
                si.DefaultComponent = 0;
                return si;
            }
        }

        #endregion

        #region IExportableContainer Members

        /// <summary>
        /// Creates a list of all the sub-objects that can be edited.
        /// </summary>
        /// <param name="recursive">
        /// The parameter is not used.
        /// </param>
        /// <returns>
        /// The created list with all the sub objects that can be edited.
        /// </returns>
        public List<IExportable> ExportObjects(bool recursive)
        {
            List<IExportable> parameters = new List<IExportable>();
            parameters.Add(this.S0);
            parameters.Add(this.V0);
            parameters.Add(this.k);
            parameters.Add(this.theta);
            parameters.Add(this.sigma);
            parameters.Add(this.zrReference);
            parameters.Add(this.dyReference);
            return parameters;
        }

        #endregion

        /// <summary>
        /// Associate the process to a zero rate and a dividend yield defined in the Fairmat model
        /// (e.g. @zr1).
        /// </summary>
        /// <param name='zrstring'>
        /// The zero rate reference.
        /// </param>
        /// <param name='dystring'>
        /// The dividend yield reference.
        /// </param>
        public void SetCurveReference(string zrstring, string dystring)
        {
            this.zrReference = new ModelParameter(zrstring);
            this.dyReference = new ModelParameter(dystring);
        }

        #region IMarkovSimulator Members

        /// <summary>
        /// Gets details about the structure of the functions A and B of the Markov
        /// process.
        /// In this case drift and volatility are only state dependant and not time
        /// dependant.
        /// </summary>
        public DynamicInfo DynamicInfo
        {
            get
            {
                return new DynamicInfo(false, true, false, true);
            }
        }

        /// <summary>
        /// Gets the starting point for the process.
        /// In this case it's in the first dimension the value
        /// of S0 and in the second dimension the value of V0.
        /// </summary>
        public double[] x0
        {
            get
            {
                double[] r0 = new double[2];
                r0[0] = this.S0.fV();
                r0[1] = this.V0.fV();
                return r0;
            }
        }

        /// <summary>
        /// This function defines the drift in the Heston process.
        /// </summary>
        /// <remarks>
        /// Heston operates on two dimensions.
        /// </remarks>
        /// <param name="i">The time step of the simulation.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="a">The output of the function.</param>
        public void a(int i, double* x, double* a)
        {
            double m = Math.Max(0, x[1]);
            a[0] = this.mu[i] - 0.5 * m;
            a[1] = this.k.fV() * (this.theta.fV() - m);
        }

        /// <summary>
        /// This function defines the volatility in the Heston process.
        /// </summary>
        /// <remarks>Heston operates on two dimensions.</remarks>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="b">The output of the function.</param>
        public void b(int i, double* x, double* b)
        {
            double m = Math.Sqrt(Math.Max(0, x[1]));
            b[0] = m;
            b[1] = this.sigma.fV() * m;
        }

        /// <summary>
        /// Sets the passed array with a Boolean stating if the process
        /// must be simulated as a log-normal process. Here the first component
        /// must be simulated as a log-normal process the second not.
        /// </summary>
        /// <param name="isLog">
        /// A reference to the array to be set with the required information.
        /// </param>
        public void isLog(ref bool[] isLog)
        {
            isLog[0] = true;
            isLog[1] = false;
        }

        #endregion IMarkovSimulator Members

        #region IParsable Members

        /// <summary>
        /// Ensure the parameters are correct.
        /// </summary>
        /// <param name='p_Context'>
        /// The underlying project.
        /// </param>
        /// <returns>
        /// False if the parameter were correct.
        /// </returns>
        public bool Parse(IProject p_Context)
        {
            bool errors = false;
            List<IExportable> list = this.ExportObjects(false);
            foreach (IExportable parameter in list)
            {
                if (parameter is IParsable)
                    BoolHelper.AddBool(errors, (parameter as IParsable).Parse(p_Context));
            }

            if (this.zrReference.Expression.IndexOf("@") == -1)
            {
                p_Context.AddError(this.zrReference.Expression +
                                   " is not a reference to a zero rate curve");
            }

            if (this.dyReference.Expression.IndexOf("@") == -1)
            {
                p_Context.AddError(this.dyReference.Expression +
                                   " is not a reference to a dividend yield curve");
            }

            return RetrieveCurve(p_Context, errors);
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
            object zrref = Engine.Parser.EvaluateAsReference(this.zrReference.Expression);
            if (!Engine.Parser.GetParserError())
            {
                this.zrCurve = zrref as Function;
                if (this.zrCurve == null)
                {
                    errors = true;
                    p_Context.AddError("Cannot find the zero rate curve " +
                                       this.zrReference.Expression);
                }
            }
            else
                errors = true;

            // check the presence of a dividend yield curve
            object dyref = Engine.Parser.EvaluateAsReference(this.dyReference.Expression);
            if (!Engine.Parser.GetParserError())
            {
                this.dyCurve = dyref as Function;
                if (this.dyCurve == null)
                {
                    errors = true;
                    p_Context.AddError("Cannot find the dividend yield curve " +
                                       this.dyReference.Expression);
                }
            }
            else
                errors = true;

            return errors;
        }

        #endregion

        /// <summary>
        /// Populate editable fields from name and value vectors
        /// specific to the Heston extended process.
        /// </summary>
        /// <param name="stocProcess">
        /// The stochastic process which is being referenced to.
        /// </param>
        /// <param name="estimate">
        /// The estimation result which contains values and names of parameters.
        /// It will be searched for S0, kappa, theta, sigma, V0 and rho.
        /// </param>
        public void Populate(IStochasticProcess stocProcess, EstimationResult estimate)
        {
            bool found;
            this.S0 = new ModelParameter(PopulateHelper.GetValue("S0", estimate.Names, estimate.Values, out found), this.S0.Description);
            this.k = new ModelParameter(PopulateHelper.GetValue("kappa", estimate.Names, estimate.Values, out found), this.k.Description);
            this.theta = new ModelParameter(PopulateHelper.GetValue("theta", estimate.Names, estimate.Values, out found), this.theta.Description);
            this.sigma = new ModelParameter(PopulateHelper.GetValue("sigma", estimate.Names, estimate.Values, out found), this.sigma.Description);
            this.V0 = new ModelParameter(PopulateHelper.GetValue("V0", estimate.Names, estimate.Values, out found), this.V0.Description);

            int index = stocProcess.NoiseIndex;
            ProjectProcess prj = stocProcess.Context as ProjectProcess;

            // Update the correlation matrix.
            prj.Processes.r.Set(index, index + 1, (RightValue)PopulateHelper.GetValue("rho", estimate.Names, estimate.Values, out found));
            bool errors = RetrieveCurve(stocProcess.Context, false);
            if (!errors)
            {
                this.zrCurve.Expr = (estimate.Objects[0] as Matrix).ToArray();
                this.dyCurve.Expr = (estimate.Objects[1] as Matrix).ToArray();
                //Calibrator assumes dividend yield is a step constant function, the simulation model must be coherent with that assumption. 
                (this.dyCurve as PFunction).m_Function.iType = DVPLUtils.EInterpolationType.ZERO_ORDER_LEFT;
            }
        }
    }
}
