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
using DVPLDOM;
using HestonEstimator;
using DVPLI;


namespace Heston
{
    /// <summary>
    /// Implementation of the Heston Stochastic Process with a constant drift.
    /// </summary>
    [Serializable]
    public unsafe class HestonProcess : IExtensibleProcess, IMarkovSimulator, IParsable,
                                        IGreeksDerivativesInfo, IEstimationResultPopulable,
                                        IOpenCLCode, IExportableContainer, IPlainVanillaPricing, IForwardStartingPricing
    {
        #region Serialized Parameters

        /// <summary>
        /// The risk free rate.
        /// </summary>
        public IModelParameter r;

        /// <summary>
        /// The dividend yield rate of the stock.
        /// </summary>
        public IModelParameter q;

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

        /// <summary>
        /// Correlation parameter
        /// </summary>
        public IModelParameter rho;

        #endregion Serialized Parameters


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
        /// value for the volatility process.
        /// </returns>
        public IModelParameter[] GetVegaFactors()
        {
            return new IModelParameter[] { this.V0 };
        }

        /// <summary>
        /// Initializes a new instance of the HestonProcess class.
        /// This is the default constructor and sets all parameters to 0.0.
        /// </summary>
        public HestonProcess()
        {
            this.r = new ModelParameter(0.0, "r");
            this.q = new ModelParameter(0.0, "q");
            this.k = new ModelParameter(0.0, "k");
            this.theta = new ModelParameter(0.0, "theta");
            this.sigma = new ModelParameter(0.0, "sigma");
            this.S0 = new ModelParameter(0.0, "S0");
            this.V0 = new ModelParameter(0.0, "V0");
            this.rho = new ModelParameter(0.0, "rho");
        }

        /// <summary>
        /// Sets some default values for the parameters.
        /// * r = 0.05
        /// * q = 0.03
        /// * k = 2.5
        /// theta = 0.4
        /// sigma = 0.2
        /// S0 = 100
        /// V0 = 0.3.
        /// </summary>
        public void DefaultInstance()
        {
            this.r = new ModelParameter(0.05, "r");
            this.q = new ModelParameter(0.03, "q");
            this.k = new ModelParameter(2.5, "k");
            this.theta = new ModelParameter(0.4, "theta");
            this.sigma = new ModelParameter(0.2, "sigma");
            this.S0 = new ModelParameter(100, "S0");
            this.V0 = new ModelParameter(0.3, "V0");
            this.rho = new ModelParameter(0.0, "rho");
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
        /// Gets the ProcessInfo for this plugin, in this case Heston.
        /// </summary>
        public ProcessInfo ProcessInfo
        {
            get
            {
                return new ProcessInfo("Heston");
            }
        }

        /// <summary>
        /// Called by Simulator after parse.
        /// The method implements only an interface requisite but does nothing.
        /// </summary>
        /// <param name='simulationDates'>The parameter is not used.</param>
        public void Setup(double[] simulationDates)
        {
        }

        /// <summary>
        /// Gets the information required for the simulator in order to control
        /// the plug-in through the <see cref="IMarkovSimulator"/> interface.
        /// The Heston Process has:
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
                si.DefaultComponent = 0; //the equity is defined at component zero
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
            return new List<IExportable>
            {
                this.S0,
                this.V0,
                this.r,
                this.q,
                this.k,
                this.theta,
                this.sigma,
                this.rho
            };
        }

        #endregion

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
                return new DynamicInfo(p_a_td: false, p_a_sd: true, p_b_td: false, p_b_sd: true);
            }
        }

        /// <summary>
        /// This function defines the drift in the Heston Markov process.
        /// </summary>
        /// <remarks>
        /// Heston operates on two dimensions.
        /// </remarks>
        /// <param name="i">The time step of the simulation.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="a">The output of the function.</param>
        [Obsolete("This method is obsolete, use the ab method instead.")]
        public void a(int i, double* x, double* a)
        {
            double m = Math.Max(0, x[1]);
            a[0] = this.r.fV() - this.q.fV() - 0.5 * m;
            a[1] = this.k.fV() * (this.theta.fV() - m);
        }

        /// <summary>
        /// This function defines the volatility in the Heston Markov process.
        /// </summary>
        /// <remarks>Heston operates on two dimensions.</remarks>
        /// <param name="i">The parameter is not used.</param>
        /// <param name="x">The state vector at the previous state.</param>
        /// <param name="b">The output of the function.</param>
        [Obsolete("This method is obsolete, use the ab method instead.")]
        public void b(int i, double* x, double* b)
        {
            double m = Math.Sqrt(Math.Max(0, x[1]));
            b[0] = m;
            b[1] = this.sigma.fV() * m;
        }

        public void ab(int i, double* x, double* a, double* b)
        {
            double V = Math.Max(0, x[1]);
            double radqV = Math.Sqrt(V);
            a[0] = this.r.fV() - this.q.fV() - 0.5 * V;
            a[1] = this.k.fV() * (this.theta.fV() - V);
            b[0] = radqV;
            b[1] = this.sigma.fV() * radqV;
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

        #endregion

        #region IOpenCLCode implementation

        /// <summary>
        /// Gets the arguments needed for an OpenCL simulation.
        /// q, k, r, theta, sigma are used in this context.
        /// </summary>
        public List<Tuple<string, object>> Arguments
        {
            get
            {
                List<Tuple<string, object>> args = new List<Tuple<string, object>>();
                args.Add(new Tuple<string, object>("q", this.q));
                args.Add(new Tuple<string, object>("k", this.k));
                args.Add(new Tuple<string, object>("r", this.r));
                args.Add(new Tuple<string, object>("theta", this.theta));
                args.Add(new Tuple<string, object>("sigma", this.sigma));
                return args;
            }
        }

        /// <summary>
        /// Gets the OpenCL code used to calculate A and B.
        /// </summary>
        public Dictionary<string, string> Code
        {
            get
            {
                Dictionary<string, string> sources = new Dictionary<string, string>();
                sources.Add("B", "double m = sqrt(max(0.0, x[1])); b[0] = m; b[1] = sigma*m;");
                sources.Add("A", "double m = max(0.0, x[1]); a[0] = r-q-0.5*m; a[1] = k*(theta-m);");
                return sources;
            }
        }

        /// <summary>
        /// Gets a value indicating whether the plugin OpenCL implementation is usable.
        /// This plugin can always run through the OpenCL simulator so it always returns true.
        /// </summary>
        public bool IsOpenCLUsable
        {
            get
            {
                return true;
            }
        }

        #endregion IOpenCLCode implementation

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

            return errors;
        }

        #endregion

        #region IPlainVanillaPricing Members

        /// <summary>
        /// Calculates the discount factor for a given time period.
        /// </summary>
        /// <param name="t1">The start of the time period.</param>
        /// <param name="t2">The end of the time period.</param>
        /// <returns>The discount factor, calculated as e^(-r*(t2-t1)), where r is the risk-free interest rate.</returns>
        public double GetDiscountFactor(int component, double t1, double t2)
        {
            // component in the context of Heston process is not used
            return Math.Exp(-this.r.fV() * (t2 - t1));
        }


        /// <summary>
        /// Calculate the price of a european call option.
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strike">
        /// The strike of the call option
        /// <param name="timeToMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding the call option
        /// </param>
        /// <returns>The call price</returns>

        public GreeksDerivatives Call(int component, double strike, double timeToMaturity, Dictionary<string, object> additionalInformation)
        {
            GreeksDerivatives result = new GreeksDerivatives();

            var kappaParam = this.k.fV();
            var thetaParam = this.theta.fV();
            var sigmaParam = this.sigma.fV();
            var rhoParam = this.rho.fV();
            var v0Param = this.V0.fV();
            var s0Param = this.S0.fV();
            var rParam = this.r.fV();
            var qParam = this.q.fV();

            AnalyticalPricingFunctionsValuationMode requestedResult =
                AttributesUtility.RetrieveAttributeOrDefaultValue(
                    additionalInformation,
                    AnalyticalPricingFunctions.GreekNameKey,
                    AnalyticalPricingFunctionsValuationMode.Price
                );


            switch (requestedResult)
            {
                default:
                case AnalyticalPricingFunctionsValuationMode.Price:
                    result.MarkToMarket = HestonCall.HestonCallPrice(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Delta:
                    var delta = HestonDelta.DeltaCall(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Deltas = (Vector)(delta);
                    break;

                case AnalyticalPricingFunctionsValuationMode.Gamma:

                    var gamma = HestonGamma.GammaCall(
                        kappa: kappaParam,
                        theta: thetaParam,
                        sigma: sigmaParam,
                        rho: rhoParam,
                        v0: v0Param,
                        s0: s0Param,
                        T: timeToMaturity,
                        K: strike,
                        r: rParam,
                        q: qParam);

                    result.Gammas = (Vector)(gamma);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Rho:

                    var rho = HestonRho.RhoCall(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Rho = rho;
                    break;


                case AnalyticalPricingFunctionsValuationMode.Theta:

                    var theta = HestonNumericalGreeks.ThetaCall(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Theta = theta;
                    break;

                case AnalyticalPricingFunctionsValuationMode.Vega:
                    var vega = HestonVega.VegaCall(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);
                    result.Vegas = (Vector)vega;
                    break;
                
                case AnalyticalPricingFunctionsValuationMode.All:
                    Func< Dictionary<string, object>, GreeksDerivatives > myFunc = 
                        (Dictionary<string, object>  info) => Call(component, strike, timeToMaturity, info);
                    return AnalyticalPricingFunctions.ComputeAllGreeks(myFunc); 
            }

            
            return result;


        }


        /// <summary>
        /// Calculate the price of a european put option.
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strike">
        /// The strike of the put option
        /// <param name="timeToMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding the put option
        /// </param>
        /// <returns>The put price</returns>
        public GreeksDerivatives Put(int component, double strike, double timeToMaturity, Dictionary<string, object> additionalInformation)
        {
            GreeksDerivatives result = new GreeksDerivatives();

            var kappaParam = this.k.fV();
            var thetaParam = this.theta.fV();
            var sigmaParam = this.sigma.fV();
            var rhoParam = this.rho.fV();
            var v0Param = this.V0.fV();
            var s0Param = this.S0.fV();
            var rParam = this.r.fV();
            var qParam = this.q.fV();



            AnalyticalPricingFunctionsValuationMode requestedResult =
                AttributesUtility.RetrieveAttributeOrDefaultValue(
                    additionalInformation,
                    AnalyticalPricingFunctions.GreekNameKey,
                    AnalyticalPricingFunctionsValuationMode.Price
                );


            switch (requestedResult)
            {
                default:
                case AnalyticalPricingFunctionsValuationMode.Price:
                    result.MarkToMarket = HestonCall.HestonPutPrice(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Delta:
                    var delta = HestonDelta.DeltaPut(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Deltas = (Vector)(delta);
                    break;

                case AnalyticalPricingFunctionsValuationMode.Gamma:

                    var gamma = HestonGamma.GammaPut(
                        kappa: kappaParam,
                        theta: thetaParam,
                        sigma: sigmaParam,
                        rho: rhoParam,
                        v0: v0Param,
                        s0: s0Param,
                        T: timeToMaturity,
                        K: strike,
                        r: rParam,
                        q: qParam);

                    result.Gammas = (Vector)(gamma);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Rho:

                    var rho = HestonRho.RhoPut(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Rho = rho;
                    break;


                case AnalyticalPricingFunctionsValuationMode.Theta:

                    var theta = HestonNumericalGreeks.ThetaPut(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Theta = theta;
                    break;

                case AnalyticalPricingFunctionsValuationMode.Vega:
                    var vega = HestonVega.VegaPut(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);
                    result.Vegas = (Vector)vega;
                    break;

                case AnalyticalPricingFunctionsValuationMode.All:
                    Func<Dictionary<string, object>, GreeksDerivatives> myFunc =
                        (Dictionary<string, object> info) => Put(component, strike, timeToMaturity, info);
                    return AnalyticalPricingFunctions.ComputeAllGreeks(myFunc);
            }


            return result;

        }



        /// <summary>
        /// Calculate the price of a european digital call option.
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strike">
        /// The strike of the digital call option
        /// <param name="timeToMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding the digital call option
        /// </param>
        /// <returns>The digital call price</returns>

        public GreeksDerivatives DigitalCall(int component, double strike, double timeToMaturity, Dictionary<string, object> additionalInformation)
        {
            GreeksDerivatives result = new GreeksDerivatives();

            var kappaParam = this.k.fV();
            var thetaParam = this.theta.fV();
            var sigmaParam = this.sigma.fV();
            var rhoParam = this.rho.fV();
            var v0Param = this.V0.fV();
            var s0Param = this.S0.fV();
            var rParam = this.r.fV();
            var qParam = this.q.fV();



            AnalyticalPricingFunctionsValuationMode requestedResult =
                AttributesUtility.RetrieveAttributeOrDefaultValue(
                    additionalInformation,
                    AnalyticalPricingFunctions.GreekNameKey,
                    AnalyticalPricingFunctionsValuationMode.Price
                );


            switch (requestedResult)
            {
                default:
                case AnalyticalPricingFunctionsValuationMode.Price:
                    result.MarkToMarket = HestonDigital.HestonDigitalCallPrice(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Delta:
                    var delta = HestonDelta.DeltaDigitalCall(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Deltas = (Vector)(delta);
                    break;

                case AnalyticalPricingFunctionsValuationMode.Gamma:

                    var gamma = HestonGamma.GammaDigitalCall(
                        kappa: kappaParam,
                        theta: thetaParam,
                        sigma: sigmaParam,
                        rho: rhoParam,
                        v0: v0Param,
                        s0: s0Param,
                        T: timeToMaturity,
                        K: strike,
                        r: rParam,
                        q: qParam);

                    result.Gammas = (Vector)(gamma);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Rho:

                    var rho = HestonRho.RhoDigitalCall(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Rho = rho;
                    break;


                case AnalyticalPricingFunctionsValuationMode.Theta:

                    var theta = HestonNumericalGreeks.ThetaDCall(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Theta = theta;
                    break;

                case AnalyticalPricingFunctionsValuationMode.Vega:
                    var vega = HestonVega.VegaDigitalCall(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);
                    result.Vegas = (Vector)vega;
                    break;

                case AnalyticalPricingFunctionsValuationMode.All:
                    Func<Dictionary<string, object>, GreeksDerivatives> myFunc =
                        (Dictionary<string, object> info) => DigitalCall(component, strike, timeToMaturity, info);
                    return AnalyticalPricingFunctions.ComputeAllGreeks(myFunc);
            }


            return result;



        }

        /// <summary>
        /// Calculate the price of a european digital put option.
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strike">
        /// The strike of the digital put option
        /// <param name="timeToMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding the digital put option
        /// </param>
        /// <returns>The digital put price</returns>
        public GreeksDerivatives DigitalPut(int component, double strike, double timeToMaturity, Dictionary<string, object> additionalInformation)
        {
            GreeksDerivatives result = new GreeksDerivatives();

            var kappaParam = this.k.fV();
            var thetaParam = this.theta.fV();
            var sigmaParam = this.sigma.fV();
            var rhoParam = this.rho.fV();
            var v0Param = this.V0.fV();
            var s0Param = this.S0.fV();
            var rParam = this.r.fV();
            var qParam = this.q.fV();



            AnalyticalPricingFunctionsValuationMode requestedResult =
                AttributesUtility.RetrieveAttributeOrDefaultValue(
                    additionalInformation,
                    AnalyticalPricingFunctions.GreekNameKey,
                    AnalyticalPricingFunctionsValuationMode.Price
                );


            switch (requestedResult)
            {
                default:
                case AnalyticalPricingFunctionsValuationMode.Price:
                    result.MarkToMarket = HestonDigital.HestonDigitalPutPrice(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Delta:
                    var delta = HestonDelta.DeltaDigitalPut(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Deltas = (Vector)(delta);
                    break;

                case AnalyticalPricingFunctionsValuationMode.Gamma:

                    var gamma = HestonGamma.GammaDigitalPut(
                        kappa: kappaParam,
                        theta: thetaParam,
                        sigma: sigmaParam,
                        rho: rhoParam,
                        v0: v0Param,
                        s0: s0Param,
                        T: timeToMaturity,
                        K: strike,
                        r: rParam,
                        q: qParam);

                    result.Gammas = (Vector)(gamma);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Rho:

                    var rho = HestonRho.RhoDigitalPut(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Rho = rho;
                    break;


                case AnalyticalPricingFunctionsValuationMode.Theta:

                    var theta = HestonNumericalGreeks.ThetaDPut(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);

                    result.Theta = theta;
                    break;

                case AnalyticalPricingFunctionsValuationMode.Vega:
                    var vega = HestonVega.VegaDigitalPut(
                       kappa: kappaParam,
                       theta: thetaParam,
                       sigma: sigmaParam,
                       rho: rhoParam,
                       v0: v0Param,
                       s0: s0Param,
                       T: timeToMaturity,
                       K: strike,
                       r: rParam,
                       q: qParam);
                    result.Vegas = (Vector)vega;
                    break;

                case AnalyticalPricingFunctionsValuationMode.All:
                    Func<Dictionary<string, object>, GreeksDerivatives> myFunc =
                        (Dictionary<string, object> info) => DigitalPut(component, strike, timeToMaturity, info);
                    return AnalyticalPricingFunctions.ComputeAllGreeks(myFunc);
            }


            return result;
   
        }

        /// <summary>
        /// Calculate the price of a swap option
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strike">
        /// The strike of the swap option
        /// <param name="swapMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding the swap option
        /// </param>
        /// <returns>The swap price</returns>
        public GreeksDerivatives Swap(int component, double strike, double swapMaturity, Dictionary<string, object> additionalInformation = null)
        {
            throw new NotImplementedException();
        }

        #endregion

        #region IForwardStartingPricing Members

        /// <summary>
        /// Calculate the price of a forward starting call option.
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strikeFraction">
        /// The strike percentage of the option
        /// <param name="fsTime">
        /// The strike fixing date of the option 
        /// <param name="timeToMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding option
        /// </param>
        /// <returns>The option price</returns>
        public GreeksDerivatives FSCall(int component, double strikeFraction, double fsTime, double timeToMaturity, Dictionary<string, object> additionalInformation = null)
        {
            var result = new GreeksDerivatives();

            AnalyticalPricingFunctionsValuationMode requestedResult =
                AttributesUtility.RetrieveAttributeOrDefaultValue(
                    additionalInformation,
                    AnalyticalPricingFunctions.GreekNameKey,
                    AnalyticalPricingFunctionsValuationMode.Price
                );


            switch (requestedResult)
            {
                default:
                case AnalyticalPricingFunctionsValuationMode.Price:
                    result.MarkToMarket = HestonForwardApproximated.HestonForwardCallPrice(
                    kappa: k.fV(),
                    theta: this.theta.fV(),
                    sigma: this.sigma.fV(),
                    rho: this.rho.fV(),
                    v0: this.V0.fV(),
                    s0: this.S0.fV(),
                    T: timeToMaturity,
                    K: strikeFraction,
                    r: this.r.fV(),
                    q: this.q.fV(),
                    T0: fsTime
                    );

                    break;

                case AnalyticalPricingFunctionsValuationMode.Delta:
                    var delta = HestonForwardApproximated.FSCallCalculateDelta(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Deltas = (Vector)(delta);
                    break;

                case AnalyticalPricingFunctionsValuationMode.Gamma:

                    var gamma = HestonForwardApproximated.FSCallCalculateGamma(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Gammas = (Vector)(gamma);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Rho:
                    var rho = HestonForwardApproximated.FSCallCalculateRho(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Rho = rho;
                    break;


                case AnalyticalPricingFunctionsValuationMode.Theta:

                    var theta = HestonForwardApproximated.FSCallCalculateTheta(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );
                    result.Theta = theta;
                    break;

                case AnalyticalPricingFunctionsValuationMode.Vega:
                    var vega = HestonForwardApproximated.FSCallCalculateVega(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Vegas = (Vector)vega;
                    break;

                case AnalyticalPricingFunctionsValuationMode.All:
                    return HestonForwardApproximated.HestonForwardDigitalCallWithGreeks(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );
            }

            return result;


        }

        /// <summary>
        /// Calculate the price of a forward starting put option.
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strikeFraction">
        /// The strike percentage of the option
        /// <param name="fsTime">
        /// The strike fixing date of the option 
        /// <param name="timeToMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding option
        /// </param>
        /// <returns>The option price</returns>
        public GreeksDerivatives FSPut(int component, double strikeFraction, double fsTime, double timeToMaturity, Dictionary<string, object> additionalInformation = null)
        {
            var result = new GreeksDerivatives();

            AnalyticalPricingFunctionsValuationMode requestedResult =
                AttributesUtility.RetrieveAttributeOrDefaultValue(
                    additionalInformation,
                    AnalyticalPricingFunctions.GreekNameKey,
                    AnalyticalPricingFunctionsValuationMode.Price
                );


            switch (requestedResult)
            {
                default:
                case AnalyticalPricingFunctionsValuationMode.Price:
                    result.MarkToMarket = HestonForwardApproximated.HestonForwardPutPrice(
                    kappa: k.fV(),
                    theta: this.theta.fV(),
                    sigma: this.sigma.fV(),
                    rho: this.rho.fV(),
                    v0: this.V0.fV(),
                    s0: this.S0.fV(),
                    T: timeToMaturity,
                    K: strikeFraction,
                    r: this.r.fV(),
                    q: this.q.fV(),
                    T0: fsTime
                    );

                    break;

                case AnalyticalPricingFunctionsValuationMode.Delta:
                    var delta = HestonForwardApproximated.FSPutCalculateDelta(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Deltas = (Vector)(delta);
                    break;

                case AnalyticalPricingFunctionsValuationMode.Gamma:

                    var gamma = HestonForwardApproximated.FSPutCalculateGamma(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Gammas = (Vector)(gamma);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Rho:
                    var rho = HestonForwardApproximated.FSPutCalculateRho(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Rho = rho;
                    break;


                case AnalyticalPricingFunctionsValuationMode.Theta:

                    var theta = HestonForwardApproximated.FSPutCalculateTheta(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );
                    result.Theta = theta;
                    break;

                case AnalyticalPricingFunctionsValuationMode.Vega:
                    var vega = HestonForwardApproximated.FSPutCalculateVega(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Vegas = (Vector)vega;
                    break;

                case AnalyticalPricingFunctionsValuationMode.All:
                    return HestonForwardApproximated.HestonForwardPutWithGreeks(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );
            }

            return result;


        }

        /// <summary>
        /// Calculate the price of a forward starting digital call option.
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strikeFraction">
        /// The strike percentage of the option
        /// <param name="fsTime">
        /// The strike fixing date of the option 
        /// <param name="timeToMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding option
        /// </param>
        /// <returns>The option price</returns>
        public GreeksDerivatives FSDigitalCall(int component, double strikeFraction, double fsTime, double timeToMaturity, Dictionary<string, object> additionalInformation = null)
        {


            var result = new GreeksDerivatives();

            AnalyticalPricingFunctionsValuationMode requestedResult =
                AttributesUtility.RetrieveAttributeOrDefaultValue(
                    additionalInformation,
                    AnalyticalPricingFunctions.GreekNameKey,
                    AnalyticalPricingFunctionsValuationMode.Price
                );


            switch (requestedResult)
            {
                default:
                case AnalyticalPricingFunctionsValuationMode.Price:
                    result.MarkToMarket = HestonForwardApproximated.HestonForwardDigitalCallPrice(
                    kappa: k.fV(),
                    theta: this.theta.fV(),
                    sigma: this.sigma.fV(),
                    rho: this.rho.fV(),
                    v0: this.V0.fV(),
                    s0: this.S0.fV(),
                    T: timeToMaturity,
                    K: strikeFraction,
                    r: this.r.fV(),
                    q: this.q.fV(),
                    T0: fsTime
                    );

                    break;

                case AnalyticalPricingFunctionsValuationMode.Delta:
                    var delta = HestonForwardApproximated.FSDCallCalculateDelta(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Deltas = (Vector)(delta);
                    break;

                case AnalyticalPricingFunctionsValuationMode.Gamma:

                    var gamma = HestonForwardApproximated.FSDCallCalculateGamma(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Gammas = (Vector)(gamma);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Rho:
                    var rho = HestonForwardApproximated.FSDCallCalculateRho(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Rho = rho;
                    break;


                case AnalyticalPricingFunctionsValuationMode.Theta:

                    var theta = HestonForwardApproximated.FSDCallCalculateTheta(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );
                    result.Theta = theta;
                    break;

                case AnalyticalPricingFunctionsValuationMode.Vega:
                    var vega = HestonForwardApproximated.FSDCallCalculateVega(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Vegas = (Vector)vega;
                    break;

                case AnalyticalPricingFunctionsValuationMode.All:
                    return HestonForwardApproximated.HestonForwardDigitalCallWithGreeks(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );
            }

            return result;

        }



        /// <summary>
        /// Calculate the price of a forward starting digital put option.
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strikeFraction">
        /// The strike percentage of the option
        /// <param name="fsTime">
        /// The strike fixing date of the option 
        /// <param name="timeToMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding option
        /// </param>
        /// <returns>The option price</returns>
        public GreeksDerivatives FSDigitalPut(int component, double strikeFraction, double fsTime, double timeToMaturity, Dictionary<string, object> additionalInformation = null)
        {
            var result = new GreeksDerivatives();

            AnalyticalPricingFunctionsValuationMode requestedResult =
                AttributesUtility.RetrieveAttributeOrDefaultValue(
                    additionalInformation,
                    AnalyticalPricingFunctions.GreekNameKey,
                    AnalyticalPricingFunctionsValuationMode.Price
                );


            switch (requestedResult)
            {
                default:
                case AnalyticalPricingFunctionsValuationMode.Price:
                    result.MarkToMarket = HestonForwardApproximated.HestonForwardDigitalPutPrice(
                    kappa: k.fV(),
                    theta: this.theta.fV(),
                    sigma: this.sigma.fV(),
                    rho: this.rho.fV(),
                    v0: this.V0.fV(),
                    s0: this.S0.fV(),
                    T: timeToMaturity,
                    K: strikeFraction,
                    r: this.r.fV(),
                    q: this.q.fV(),
                    T0: fsTime
                    );

                    break;

                case AnalyticalPricingFunctionsValuationMode.Delta:
                    var delta = HestonForwardApproximated.FSDPutCalculateDelta(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Deltas = (Vector)(delta);
                    break;

                case AnalyticalPricingFunctionsValuationMode.Gamma:

                    var gamma = HestonForwardApproximated.FSDPutCalculateGamma(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Gammas = (Vector)(gamma);
                    break;


                case AnalyticalPricingFunctionsValuationMode.Rho:
                    var rho = HestonForwardApproximated.FSDPutCalculateRho(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Rho = rho;
                    break;


                case AnalyticalPricingFunctionsValuationMode.Theta:

                    var theta = HestonForwardApproximated.FSDPutCalculateTheta(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );
                    result.Theta = theta;
                    break;

                case AnalyticalPricingFunctionsValuationMode.Vega:
                    var vega = HestonForwardApproximated.FSDPutCalculateVega(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );

                    result.Vegas = (Vector)vega;
                    break;

                case AnalyticalPricingFunctionsValuationMode.All:
                    return HestonForwardApproximated.HestonForwardDigitalPutWithGreeks(
                        kappa: k.fV(),
                        theta: this.theta.fV(),
                        sigma: this.sigma.fV(),
                        rho: this.rho.fV(),
                        v0: this.V0.fV(),
                        s0: this.S0.fV(),
                        T: timeToMaturity,
                        K: strikeFraction,
                        r: this.r.fV(),
                        q: this.q.fV(),
                        T0: fsTime
                        );
            }

            return result;
        }

        /// <summary>
        /// Calculate the price of a forward starting swap option.
        /// </summary>
        /// <param name="component">
        /// The component of the process
        /// </param>
        /// <param name="strikeFraction">
        /// The strike percentage of the option
        /// <param name="fsTime">
        /// The strike fixing date of the option 
        /// <param name="timeToMaturity">
        /// The time to maturity of the option: T-t
        /// <param name="additionalInformation">"
        /// Additional information regarding option
        /// </param>
        /// <returns>The option price</returns>
        public GreeksDerivatives FSSwap(int component, double strikeFraction, double fsTime, double swapMaturity, Dictionary<string, object> additionalInformation = null)
        {
            throw new NotImplementedException();
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
        /// It will be searched for S0, kappa, theta, sigma, V0, r, q and rho.
        /// </param>
        public void Populate(IStochasticProcess stocProcess, EstimationResult estimate)
        {
            bool found;
            this.S0 = new ModelParameter(PopulateHelper.GetValue("S0", estimate.Names, estimate.Values, out found), this.S0.Description);
            this.k = new ModelParameter(PopulateHelper.GetValue("kappa", estimate.Names, estimate.Values, out found), this.k.Description);
            this.theta = new ModelParameter(PopulateHelper.GetValue("theta", estimate.Names, estimate.Values, out found), this.theta.Description);
            this.sigma = new ModelParameter(PopulateHelper.GetValue("sigma", estimate.Names, estimate.Values, out found), this.sigma.Description);
            this.V0 = new ModelParameter(PopulateHelper.GetValue("V0", estimate.Names, estimate.Values, out found), this.V0.Description);
            this.r = new ModelParameter(PopulateHelper.GetValue("r", estimate.Names, estimate.Values, out found), this.r.Description);
            this.q = new ModelParameter(PopulateHelper.GetValue("q", estimate.Names, estimate.Values, out found), this.q.Description);

            var rhoEstimate = PopulateHelper.GetValue("rho", estimate.Names, estimate.Values, out found);
            this.rho = new ModelParameter(rhoEstimate, this.rho.Description);

            int index = stocProcess.NoiseIndex;
            ProjectProcess prj = stocProcess.Context as ProjectProcess;

            // Updates the correlation.
            prj.Processes.r.Set(index, index + 1, (RightValue)rhoEstimate);

        }

    }
}
