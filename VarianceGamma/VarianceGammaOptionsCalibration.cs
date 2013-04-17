/* Copyright (C) 2013 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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
using DVPLI;
using Fairmat.MarketData;
using Fairmat.Optimization;
using Mono.Addins;
using Fairmat.Math;

namespace VarianceGamma
{
    /// <summary>
    /// Calibrates the Variance Gamma model against observed option prices.
    /// </summary>
    [Extension("/Fairmat/Estimator")]
    public class VarianceGammaOptionsCalibration : IEstimator, IDescription
    {
        /// <summary>
        /// A vector of Maturities.
        /// </summary>
        private Vector m;

        /// <summary>
        /// A vector of Strikes.
        /// </summary>
        private Vector k;

        /// <summary>
        /// A matrix with the Call prices.
        /// </summary>
        private Matrix cp;

        /// <summary>
        /// Dividend yield
        /// </summary>
        private double q;

        /// <summary>
        /// Process starting value
        /// </summary>
        private double s0;

        /// <summary>
        /// Risk free rate
        /// </summary>
        private double r;

        /// <summary>
        /// Gets the value requested by the interface ProvidesTo,
        /// returning VarianceGamma as the type.
        /// </summary>
        public Type ProvidesTo
        {
            get
            {
                return typeof(VarianceGamma);
            }
        }

        /// <summary>
        /// Gets the description of the implemented calibration function.
        /// </summary>
        public string Description
        {
            get
            {
                return "Variance Gamma Options calibrator";
            }
        }

        /// <summary>
        /// Gets the types required by the estimator in order to work:
        /// EquitySpotMarketData and CallPriceMarketData are
        /// the required type for this estimator.
        /// </summary>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="multivariateRequest">The parameter is not used.</param>
        /// <returns>
        /// An array containing the type EquitySpotMarketData and CallPriceMarketData.
        /// </returns>
        public EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { new EstimateRequirement(typeof(EquitySpotMarketData)), 
                                               new EstimateRequirement(typeof(CallPriceMarketData)), 
                                               new EstimateRequirement(typeof(DiscountingCurveMarketData)) };
        }

        /// <summary>
        /// Attempts to solve the Variance Gamma Optimization problem using
        /// <see cref="Heston.VarianceGammaOptimizationProblem"/>.
        /// </summary>
        /// <param name="data">
        /// The data to be used in order to perform the optimization.
        /// </param>
        /// <param name="settings">The parameter is not used.</param>
        /// <returns>The results of the optimization.</returns>
        public EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            EquitySpotMarketData espmd = data[0] as EquitySpotMarketData;
            CallPriceMarketData cpmd = data[1] as CallPriceMarketData;
            DiscountingCurveMarketData dcmd = data[2] as DiscountingCurveMarketData;
            EquityCalibrationData ecd = new EquityCalibrationData(cpmd, dcmd);
            this.s0 = espmd.Price;
            this.r = espmd.RiskFreeRate;
            this.q = espmd.DividendYield;

            this.k = ecd.Hdata.Strike;
            this.m = ecd.Hdata.Maturity;
            this.cp = ecd.Hdata.CallPrice;

            Vector x0 = (Vector)new double[] { 0.1, 0.1, 0.1 };
            IOptimizationAlgorithm algorithm = new QADE();
            OptimizationSettings optimizationSettings = new DESettings();

            // Maximum number of iteration allowed.
            optimizationSettings.MaxIter = 50;

            // Positive integer values print debug info.
            optimizationSettings.Verbosity = 1;

            // Tolerance.
            optimizationSettings.epsilon = 10e-5;
            var solution = algorithm.Minimize(new VarianceGammaOptimizationProblem(this.q, this.s0, this.k,
                                                                          this.r, this.cp, this.m),
                                                                          optimizationSettings, x0);

            if (solution.errors)
                return new EstimationResult(solution.message);

            var er = new EstimationResult();
            er.Names= new string[]{"S0","theta","sigma","nu","rate","dividend"};
            er.Values = new double[] {this.s0,solution.x[0],solution.x[1],solution.x[2],this.r,this.q};
            return er;
        }

        /// <summary>
        /// Calculate the price of a call option.
        /// Notation and formulas follow Whitley A. "Pricing of European, Bermudan and American Options
        /// under the Exponential Variance Gamma Process" (2009) apart from a correction in the last line
        /// where we substituted q (dividend yield) in place of r (risk free rate).
        /// </summary>
        /// <param name="theta">VG theta parameter</param>
        /// <param name="sigma">VG sigma parameter</param>
        /// <param name="nu">VG nu parameter</param>
        /// <param name="t">Call maturity</param>
        /// <param name="k">Call strike</param>
        /// <param name="q">Dividend yield</param>
        /// <param name="s0">Process starting value</param>
        /// <param name="r">Risk free rate</param>
        /// <returns>Call price value</returns>
        public static double VGCall(double theta, double sigma, double nu, double t, double k, double q, double s0, double r)
        {
            double omega = Math.Log(1 - 0.5 * sigma * sigma * nu - theta * nu) / nu;
            double zeta = (Math.Log(s0 / k) + omega * t) / sigma;
            double v = 1.0 - nu * (theta + 0.5 * sigma * sigma);
            double a1 = zeta * Math.Sqrt(v / nu);
            double b1 = (theta + sigma * sigma) * Math.Sqrt(nu / v) / sigma;
            double a2 = zeta / Math.Sqrt(nu);
            double b2 = theta * Math.Sqrt(nu) / sigma;
            double c = t / nu;
            double p1 = Psi(a1, b1, c);
            double p2 = Psi(a2, b2, c);
            return s0 * Math.Exp(-q * t) * p1 - k * Math.Exp(-r * t) * p2;
        }

        /// <summary>
        /// Psi function that enters in the call price formula
        /// Notation and formulas follow Whitley A. "Pricing of European, Bermudan and American Options
        /// under the Exponential Variance Gamma Process" (2009)
        /// </summary>
        /// <param name="a">First variable</param>
        /// <param name="b">Second variable</param>
        /// <param name="c">Third variable</param>
        /// <returns>PsiBH value</returns>
        private static double Psi(double a, double b, double c)
        {
            if (c <= 0 && (c - Math.Floor(c)) == 0)
                throw new Exception("Psi function is not defined on negative integers");
            double u = b / Math.Sqrt(2.0 + b * b);
            double d = Math.Abs(a) * Math.Sqrt(2.0 + b * b);
            double fi = Math.Pow(d, c + 0.5) * Math.Exp(Math.Sign(a) * d) * Math.Pow(1.0 + u, c) * BesselK(c + 0.5, d) / (Math.Sqrt(2.0 * Math.PI) * Gamma(c) * c);
            double se = -Math.Sign(a) * Math.Pow(d, c + 0.5) * Math.Exp(Math.Sign(a) * d) * Math.Pow(1.0 + u, 1.0 + c) * BesselK(c - 0.5, d) / (Math.Sqrt(2.0 * Math.PI) * Gamma(c) * (1.0 + c));
            double th = Math.Sign(a) * Math.Pow(d, c + 0.5) * Math.Exp(Math.Sign(a) * d) * Math.Pow(1 + u, c) * BesselK(c - 0.5, d) / (Math.Sqrt(2.0 * Math.PI) * Gamma(c) * c);
            double first = CH(c, 1.0 - c, 1.0 + c, 0.5 * (1.0 + u), -Math.Sign(a) * d * (1.0 + u));
            double second = CH(1.0 + c, 1.0 - c, 2.0 + c, 0.5 * (1.0 + u), -Math.Sign(a) * d * (1.0 + u));
            double third = CH(c, 1.0 - c, 1.0 + c, 0.5 * (1.0 + u), -Math.Sign(a) * d * (1.0 + u));
            return fi * first + se * second + th * third;
        }

        /// <summary>
        /// Simple wrapper for the Fairmat provided factorial function.
        /// </summary>
        /// <param name="n">The value to use to evaluate the factorial function.</param>
        /// <returns>The value of the factorial function.</returns>
        private static double Factorial(int n)
        {
            return Fairmat.Math.SpecialFunctions.Fact(n);
        }

        /// <summary>
        /// Simple wrapper for the Fairmat provided Gamma function.
        /// It handles some special cases required by Variance Gamma.
        /// </summary>
        /// <param name="z">The value to use to evaluate the Gamma function.</param>
        /// <returns>The value of the Gamma function.</returns>
        private static double Gamma(double z)
        {
            if (z >= 0)
                return Fairmat.Math.SpecialFunctions.Gamma(z);
            else
            {
                if (Math.Floor(z) == z)
                    throw new Exception("Gamma function is not defined on negative integers");
                else
                    return Math.PI / (Fairmat.Math.SpecialFunctions.Gamma(1.0 - z) * Math.Sin(Math.PI * z));
            }
        }

        /// <summary>
        /// Calculates the modified Bessel function of the first kind
        /// </summary>
        /// <param name="nu">Parameter</param>
        /// <param name="z">Variable</param>
        /// <returns>Bessel function of the first kind value</returns>
        private static double BesselI(double nu, double z)
        {
            double s = 0;
            for (int k = 0; k <= 20; k++)
            {
                s = s + Math.Pow(z / 2.0, 2.0 * k + nu) / (Factorial(k) * Gamma(k + nu + 1.0));
            }
            return s;
        }

        /// <summary>
        /// Calculates the modified Bessel function of the second kind
        /// </summary>
        /// <param name="nu">Parameter</param>
        /// <param name="z">Variable</param>
        /// <returns>Bessel function of the second kind value</returns>
        private static double BesselK(double nu, double z)
        {
            if ((nu - Math.Floor(nu)) == 0.0)
                throw new Exception("BesselK function is not defined on integers");
            return (Math.PI / 2.0) * (BesselI(-nu, z) - BesselI(nu, z)) / Math.Sin(nu * Math.PI);
        }

        /// <summary>
        /// Calculates the confluent hypergeometric function of two variables
        /// </summary>
        /// <param name="a">First parameter</param>
        /// <param name="b">Second parameter</param>
        /// <param name="c">Third parameter</param>
        /// <param name="x">First variable</param>
        /// <param name="y">Second variable</param>
        /// <returns>Confluent hypergeometric function of two variables value</returns>
        private static double CH(double a, double b, double c, double x, double y)
        {
            // Throw exception if z is not greater than zero.
            if (Math.Abs(x) > 1)
            {
                throw new ArgumentException("Value |x| must be smaller than 1");
            }
            else if (((a - Math.Floor(a)) == 0.0 && a <= 0) || ((b - Math.Floor(b)) == 0.0 && b <= 0) || ((c - Math.Floor(c)) == 0.0 && c <= 0))
            {
                throw new ArgumentException("a, b, c cannot be negative integers");
            }
            else
            {
                int numMax = 50;
                double s = 0;
                double err = 1;
                while (err > Math.Pow(10, -6))
                {
                    for (int m = 0; m <= numMax; m++)
                    {
                        for (int n = 0; n <= numMax; n++)
                        {
                            double term = ((Gamma(a + m + n) / Gamma(a)) * (Gamma(b + m) / Gamma(b))) / ((Gamma(c + m + n) / Gamma(c)) * Factorial(m) * Factorial(n)) * Math.Pow(x, m) * Math.Pow(y, n);
                            double snew = s + term;
                            err = Math.Abs(snew - s);
                            s = snew;
                        }
                    }
                    numMax = numMax + 1;
                }
                return s;
            }
        }
    }
}
