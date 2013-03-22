/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Safe Khampol (safe.khampol@gmail.com)
 *            Matteo Tesser (matteo.tesser@fairmat.com)
 *            Enrico Degiuli (enrico.degiuli@fairmat.com)
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

namespace VarianceGamma
{
    /// <summary>
    /// Implementation of the Variance Gamma model simulation.
    /// </summary>
    [Serializable]
    public class VarianceGamma : IExtensibleProcess, IFullSimulator, IParsable,
                                 IPopulable, IGreeksDerivativesInfo, ISerializable
    {
        #region SerializedFields

        /// <summary>
        /// Version of Variance Gamma, used to mark versions of serialized objects.
        /// </summary>
        private static int varianceGammaVersion = 1;

        /// <summary>
        /// VG starting value.
        /// </summary>
        private IModelParameter s0;

        /// <summary>
        /// VG theta parameter.
        /// </summary>
        private IModelParameter theta;

        /// <summary>
        /// VG sigma parameter.
        /// </summary>
        private IModelParameter sigma;

        /// <summary>
        /// VG nu parameter (adjusting gamma distribution parameter).
        /// </summary>
        private IModelParameter nu;

        /// <summary>
        /// Risk free rate.
        /// </summary>
        private IModelParameter rate;

        /// <summary>
        /// Continuous dividend yield.
        /// </summary>
        private IModelParameter dividend;

        #endregion

        /// <summary>
        /// Gamma distribution generator.
        /// </summary>
        [NonSerialized]
        private Fairmat.Statistics.Gamma gamma;

        /// <summary>
        /// Drift of the exponent composed by risk free rate, dividend yield and omega.
        /// </summary>
        [NonSerialized]
        private double drift;

        /// <summary>
        /// Defines the description of the S0 model variable.
        /// </summary>
        private const string s0Description = "S0";

        /// <summary>
        /// Defines the description of the theta model variable.
        /// </summary>
        private const string thetaDescription = "Theta";

        /// <summary>
        /// Defines the description of the sigma model variable.
        /// </summary>
        private const string sigmaDescription = "Sigma";

        /// <summary>
        /// Defines the description of the nu model variable.
        /// </summary>
        private const string nuDescription = "Nu";

        /// <summary>
        /// Defines the description of the risk free rate parameter.
        /// </summary>
        private const string rateDescription = "r";

        /// <summary>
        /// Defines the description of the dividend yield parameter.
        /// </summary>
        private const string dividendDescription = "q";

        /// <summary>
        /// Initializes a new instance of the VarianceGamma class.
        /// This is the default constructor.
        /// </summary>
        public VarianceGamma()
            : this(100, 0.1, 0.1, 0.1, 0.02, 0.01)
        {
        }

        /// <summary>
        /// Initializes a new instance of the VarianceGamma class.
        /// Defines a new process with given parameter values.
        /// </summary>
        /// <param name="s0">Starting value.</param>
        /// <param name="theta">VG theta parameter.</param>
        /// <param name="sigma">VG sigma parameter.</param>
        /// <param name="nu">VG nu parameter.</param>
        /// <param name="r">Risk free rate.</param>
        /// <param name="q">Continuous dividend yield.</param>
        public VarianceGamma(double s0, double theta, double sigma, double nu, double r, double q)
        {
            this.s0 = new ModelParameter(s0, s0Description);
            this.theta = new ModelParameter(theta, thetaDescription);
            this.sigma = new ModelParameter(sigma, sigmaDescription);
            this.nu = new ModelParameter(nu, nuDescription);
            this.rate = new ModelParameter(r, rateDescription);
            this.dividend = new ModelParameter(q, dividendDescription);
        }

        #region IParsable Members

        /// <summary>
        /// Ensure the parameters are correct.
        /// </summary>
        /// <param name='context'>
        /// The underlying project.
        /// </param>
        /// <returns>
        /// False if there were no parse errors.
        /// </returns>
        public bool Parse(IProject context)
        {
            bool errors = false;
            BoolHelper.AddBool(errors, this.s0.Parse(context));
            BoolHelper.AddBool(errors, this.theta.Parse(context));
            BoolHelper.AddBool(errors, this.sigma.Parse(context));
            BoolHelper.AddBool(errors, this.nu.Parse(context));
            BoolHelper.AddBool(errors, this.rate.Parse(context));
            BoolHelper.AddBool(errors, this.dividend.Parse(context));
            return errors;
        }

        #endregion

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
        /// Gets the ProcessInfo for this plugin, in this case Variance Gamma.
        /// </summary>
        public ProcessInfo ProcessInfo
        {
            get
            {
                return new ProcessInfo("Variance Gamma");
            }
        }

        /// <summary>
        /// Gets the information required in order to allow the simulation to run.
        /// </summary>
        public SimulationInfo SimulationInfo
        {
            get
            {
                SimulationInfo s = new SimulationInfo();
                s.LatentSize = 0;
                s.NoiseSize = 1;
                s.StateDescription = new string[] { "Index value" };
                s.StateSize = 1;
                return s;
            }
        }

        /// <summary>
        /// Called by Simulator after parse.
        /// Initializes here time-dependant but not state dependent variables.
        /// </summary>
        /// <param name='dates'>
        /// The dates at which the process realizations will be requested.
        /// </param>
        public void Setup(double[] dates)
        {
            // Note that gamma distribution in Fairmat has the parametrization with shape
            // alpha and rate beta (see "Gamma distribution" in wikipedia).
            // In some paper on variance gamma this formulas are different
            // due to different gamma parametrization.
            double omega = Math.Log(1.0 - this.nu.fV() * this.theta.fV() - 0.5 * this.sigma.fV() * this.sigma.fV() * this.nu.fV()) / this.nu.fV();
            this.drift = this.rate.fV() - this.dividend.fV() + omega;
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
            parameters.Add(this.s0);
            parameters.Add(this.theta);
            parameters.Add(this.sigma);
            parameters.Add(this.nu);
            parameters.Add(this.rate);
            parameters.Add(this.dividend);
            return parameters;
        }

        #endregion

        #region IPopulable Members

        /// <summary>
        /// Populate editable fields from name and value vectors specific to VG.
        /// </summary>
        /// <param name="names">
        /// An array with the names of the variables.
        /// </param>
        /// <param name="values">The values associated to the parameters in names.</param>
        public void Populate(string[] names, double[] values)
        {
            bool found = false;
            this.theta = new ModelParameter(PopulateHelper.GetValue("S0", names, values, out found), thetaDescription);
            this.sigma = new ModelParameter(PopulateHelper.GetValue("theta", names, values, out found), sigmaDescription);
            this.theta = new ModelParameter(PopulateHelper.GetValue("sigma", names, values, out found), thetaDescription);
            this.theta = new ModelParameter(PopulateHelper.GetValue("nu", names, values, out found), thetaDescription);
            this.theta = new ModelParameter(PopulateHelper.GetValue("rate", names, values, out found), thetaDescription);
            this.theta = new ModelParameter(PopulateHelper.GetValue("dividend", names, values, out found), thetaDescription);
        }

        #endregion

        #region IGreeksDerivativesInfo implementation

        /// <summary>
        /// Gets the factors for Delta Greek derivative.
        /// </summary>
        /// <returns>
        /// Null as the functionality is not implemented.
        /// </returns>
        public IModelParameter[] GetDeltaFactors()
        {
            return null;
        }

        /// <summary>
        /// Gets the factors for Vega Greek derivative.
        /// </summary>
        /// <returns>
        /// A model parameter containing the sigma value of VG.
        /// </returns>
        public IModelParameter[] GetVegaFactors()
        {
            return new IModelParameter[] { this.sigma };
        }

        #endregion

        #region Serialization

        #region ISerializable Members

        /// <summary>
        /// Initializes the serialization information of the current object.
        /// </summary>
        /// <param name="info">The SerializationInfo to populate with data.</param>
        /// <param name="context">The StreamingContext that contains contextual information
        /// about the destination.</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("varianceGammaVersion", varianceGammaVersion);
            info.AddValue("s0", this.s0);
            info.AddValue("theta", this.theta);
            info.AddValue("sigma", this.sigma);
            info.AddValue("nu", this.nu);
            info.AddValue("rate", this.rate);
            info.AddValue("dividend", this.dividend);
        }

        #endregion

        /// <summary>
        /// Initializes a new instance of the VarianceGamma class.
        /// Initializes the object based on the serialized data.
        /// </summary>
        /// <param name="info">The SerializationInfo that holds the serialized object data.</param>
        /// <param name="context">The StreamingContext that contains contextual
        /// information about the source.</param>
        protected VarianceGamma(SerializationInfo info, StreamingContext context)
        {
            int readVersion = info.GetInt32("varianceGammaVersion");
            this.s0 = (IModelParameter)ObjectSerialization.GetValue2(info, "s0", typeof(IModelParameter));
            this.theta = (IModelParameter)ObjectSerialization.GetValue2(info, "theta", typeof(IModelParameter));
            this.sigma = (IModelParameter)ObjectSerialization.GetValue2(info, "sigma", typeof(IModelParameter));
            this.nu = (IModelParameter)ObjectSerialization.GetValue2(info, "nu", typeof(IModelParameter));
            this.rate = (IModelParameter)ObjectSerialization.GetValue2(info, "rate", typeof(IModelParameter));
            this.dividend = (IModelParameter)ObjectSerialization.GetValue2(info, "dividend", typeof(IModelParameter));
        }

        #endregion

        #region IFullSimulator
        /// <summary>
        /// Manage the simulation of the variance gamma process
        /// </summary>
        /// <param name="Dates">Simulation dates</param>
        /// <param name="Noise">Gaussian noise for a single path</param>
        /// <param name="OutDynamic">Single path process realization</param>
        public void Simulate(double[] Dates, IReadOnlyMatrixSlice Noise, IMatrixSlice OutDynamic)
        { 
            int steps = OutDynamic.R;
            double GammaNoise, dt;
            OutDynamic[0, 0] = this.s0.fV();
            for (int i = 1; i < steps; i++)
            {
                dt = Dates[i] - Dates[i - 1];
                this.gamma = new Fairmat.Statistics.Gamma(dt / this.nu.fV(), 1.0 / this.nu.fV());
                GammaNoise = this.gamma.Draw(Engine.Generator);
                OutDynamic[i, 0] = OutDynamic[i - 1, 0] * Math.Exp(this.drift * dt +
                    this.theta.fV() * GammaNoise + this.sigma.fV() * Math.Sqrt(GammaNoise) * Noise[i - 1, 0]);
            }
        }

        #endregion
    }
}
