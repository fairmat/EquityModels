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
using DVPLI;
using Fairmat.Math;

namespace HestonEstimator
{
    /// <summary>
    /// Handles the calculation of the call prices on a given Heston model instance.
    /// </summary>
    public class HestonCall : IIntegrable
    {
        #region Model Parameters
        /// <summary>
        /// Heston volatility mean reversion speed parameter.
        /// </summary>
        private double kappa;

        /// <summary>
        /// Heston volatility mean reversion level parameter.
        /// </summary>
        private double theta;

        /// <summary>
        /// Heston volatility of volatility parameter.
        /// </summary>
        private double sigma;

        /// <summary>
        /// Correlation between the two Wiener processes in the Heston dynamics.
        /// </summary>
        private double rho;

        /// <summary>
        /// Starting value for the stock process.
        /// </summary>
        private double v0;

        /// <summary>
        /// Starting value for the volatility process.
        /// </summary>
        private double s0;

        #endregion Model Parameters

        /// <summary>
        /// Risk free rate.
        /// </summary>
        internal double rate;

        /// <summary>
        /// Dividend yield.
        /// </summary>
        internal double dividend;

        /// <summary>
        /// Call option maturity.
        /// </summary>
        internal double T;

        /// <summary>
        /// Call option strike value.
        /// </summary>
        internal double K;

        /// <summary>
        /// Defines the row of the hestonCallPrice matrix on which
        /// perform a partial calculation of the objective function.
        /// </summary>
        internal int row;

        /// <summary>
        /// Keeps track of partial summation of objective function during Heston calibration.
        /// </summary>
        internal double sum = 0;

        /// <summary>
        /// Heston model call price matrix.
        /// </summary>
        internal Matrix hestonCallPrice;

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonCall"/> class.
        /// </summary>
        public HestonCall()
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonCall"/> class.
        /// </summary>
        /// <param name='problem'>
        /// HestonCallOptimizationProblem at which call price calculations are to be linked.
        /// </param>
        internal HestonCall(HestonCallOptimizationProblem problem)
        {
            hestonCallPrice = new Matrix(problem.callMarketPrice.R, problem.callMarketPrice.C);
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonCall"/> class.
        /// </summary>
        /// <param name='problem'>
        /// HestonCallOptimizationProblem at which call price calculations are to be linked.
        /// </param>
        /// <param name='x'>
        /// Vector of Heston model parameters.
        /// </param>
        /// <param name='s0'>
        /// Starting value for the stock process.
        /// </param>
        internal HestonCall(HestonCallOptimizationProblem problem, Vector x, double s0)
            : this(problem)
        {
            this.kappa = x[0];
            this.theta = x[1];
            this.sigma = x[2];
            this.rho = x[3];
            this.v0 = x[4];
            this.s0 = s0;
        }

        /// <summary>
        /// Calculates the Heston model call price by using local variables.
        /// </summary>
        /// <remarks>
        /// The indefinite integral present in the formula is performed with extremes [1E-8, 1000].
        /// This choice seems to be sufficient to correctly estimate the indefinite integral.
        /// </remarks>
        /// <returns>The price of a call option within the Heston.</returns>
        internal double HestonCallPrice()
        {
            double F = this.s0 * Math.Exp((this.rate - this.dividend) * this.T);
            double firstTerm = 0.5 * (F - this.K);
            double a = 1E-8;
            double b = 1000.0;
            Integrate integrate = new Integrate(this);

            // The second term of this expressions approximates the integral in the interval [0,a].
            double integral = integrate.AdaptLobatto(a, b) + a * IntegrandFunc(a / 2.0);
            double call = Math.Exp(-this.rate * this.T) * (firstTerm + integral / Math.PI);
            return call;
        }

        /// <summary>
        /// Calculates the Heston model call price.
        /// </summary>
        /// <remarks>
        /// Note that this function is a wrapper for the other one which sets
        /// all the variables locally before calling it.</remarks>
        /// <param name="x">Vector of Heston model parameters.</param>
        /// <param name="s0">Starting value for the stock process.</param>
        /// <param name="T">Maturity of the call option which price is to be calculated.</param>
        /// <param name="K">Strike of the call option which price is to be calculated.</param>
        /// <param name="r">Risk free rate.</param>
        /// <param name="q">Dividend yield.</param>
        /// <returns>The price of a call option within the Heston model.</returns>
        public double HestonCallPrice(Vector x, double s0, double T, double K, double r, double q)
        {
            this.kappa = x[0];
            this.theta = x[1];
            this.sigma = x[2];
            this.rho = x[3];
            this.v0 = x[4];

            this.s0 = s0;
            this.T = T;
            this.K = K;

            this.rate = r;
            this.dividend = q;

            return HestonCallPrice();
        }

        /// <summary>
        /// Calculates value of the integrand function that
        /// appears in the Heston model call price formula.
        /// </summary>
        /// <param name="u"> Value at which the integrand function is to be calculated.</param>
        /// <returns>
        /// The value of the integrand function.
        /// </returns>
        public double IntegrandFunc(double u)
        {
            Complex I = Complex.I;
            Complex Iu = Complex.I * u;
            Complex A = Complex.Exp(-Iu * Math.Log(this.K));
            Complex complexVal1 = A * Phi(u - I, this.kappa, this.theta, this.sigma, this.rho, this.s0, this.v0, this.rate - this.dividend, this.T) / Iu;
            Complex complexVal2 = A * Phi(u, this.kappa, this.theta, this.sigma, this.rho, this.s0, this.v0, this.rate - this.dividend, this.T) / Iu;
            return complexVal1.Re - this.K * complexVal2.Re;
        }

        /// <summary>
        /// Calculates value of the integrand function that appears
        /// in the Heston model call price formula.
        /// </summary>
        /// <param name="u"> Value at which the integrand function is to be calculated.</param>
        /// <param name="kappa">Heston volatility mean reversion speed parameter.</param>
        /// <param name="theta">Heston volatility mean reversion level parameter.</param>
        /// <param name="sigma">Heston volatility of volatility parameter.</param>
        /// <param name="rho">
        /// Correlation between the two Wiener processes in the Heston dynamics.
        /// </param>
        /// <param name="s0">Starting value for the stock process.</param>
        /// <param name="v0">Starting value for the volatility process.</param>
        /// <param name="rate">Risk free rate.</param>
        /// <param name="dividend">Dividend yield.</param>
        /// <param name="T">Maturity of the call option which price is to be calculated.</param>
        /// <param name="K">Strike of the call option which price is to be calculated.</param>
        /// <returns>
        /// The value of the integrand function.
        /// </returns>
        public double IntegrandFunc(double u, double kappa, double theta, double sigma, double rho, double s0, double v0, double rate, double dividend, double T, double K)
        {
            this.kappa = kappa;
            this.theta = theta;
            this.sigma = sigma;
            this.rho = rho;
            this.s0 = s0;
            this.v0 = v0;
            this.rate = rate;
            this.dividend = dividend;
            this.T = T;
            this.K = K;

            return IntegrandFunc(u);
        }

        /// <summary>
        /// Calculates the Heston characteristic function.
        /// </summary>
        /// <remarks>
        /// The calculation does not follow the original paper by Heston
        /// but uses the form found in the paper
        /// <br />
        /// Schoutens, W., Simons E. and Tistaert, J. (2004)
        /// A perfect calibration! Now what?, Wilmott Magazine, March 2005, 66–78.
        /// <br />
        /// As stated in:
        /// Albrecher, H., Mayer, Ph., Schoutens, W and Tistaert, J. (2007)
        /// The Little Heston Trap. Wilmott Magazine, January Issue, 83-92.
        /// <br />
        /// This form does not have the discontinuity issue that affects the form
        /// found in the original Heston paper.
        /// </remarks>
        /// <param name="u">
        /// Complex value at which the characteristic function is to be calculated.
        /// </param>
        /// <param name="kappa">Heston volatility mean reversion speed parameter.</param>
        /// <param name="theta">Heston volatility mean reversion level parameter.</param>
        /// <param name="sigma">Heston volatility of volatility parameter.</param>
        /// <param name="rho">
        /// Correlation between the two Wiener processes in the Heston dynamics.
        /// </param>
        /// <param name="s0">Starting value for the stock process.</param>
        /// <param name="v0">Starting value for the volatility process.</param>
        /// <param name="r">Risk free rate.</param>
        /// <param name="T">Time at which the characteristic function is to be calculated.</param>
        /// <returns>The value of the characteristic function.</returns>
        public Complex Phi(Complex u, double kappa, double theta, double sigma, double rho, double s0, double v0, double r, double T)
        {
            Complex d, par, g, A, B, val;
            Complex I = Complex.I;
            double ss = sigma * sigma;
            Complex tmp1 = I * rho * sigma * u;
            d = Complex.Sqrt(Complex.Pow(tmp1 - kappa, 2.0) + ss * (I * u + u * u));
            Complex tmp2 = kappa - tmp1;
            par = tmp2 - d;
            g = par / (tmp2 + d);

            Complex edT = Complex.Exp(-d * T);
            Complex numArg = 1.0 - g * edT;
            A = (theta * kappa) * (par * T - 2.0 * Complex.Log(numArg / (1.0 - g))) / (sigma * sigma);
            B = v0 * (par * (1.0 - edT) / numArg) / ss;

            val = Complex.Exp(I * u * (Math.Log(s0) + r * T) + A + B);
            return val;
        }

        /// <summary>
        /// Calculates Heston characteristic function with input u real.
        /// </summary>
        /// <param name="u"> Value at which the characteristic function is to be calculated.</param>
        /// <param name="kappa">Heston volatility mean reversion speed parameter.</param>
        /// <param name="theta">Heston volatility mean reversion level parameter.</param>
        /// <param name="sigma">Heston volatility of volatility parameter.</param>
        /// <param name="rho">
        /// Correlation between the two Wiener processes in the Heston dynamics.
        /// </param>
        /// <param name="v0">Starting value for the volatility process.</param>
        /// <param name="s0">Starting value for the stock process.</param>
        /// <param name="r">Risk free rate.</param>
        /// <param name="T">Time at which the characteristic function is to be calculated.</param>
        /// <returns>The value of the characteristic function.</returns>
        internal Complex Phi(double u, double kappa, double theta, double sigma, double rho, double v0, double s0, double r, double T)
        {
            Complex Cu = new Complex(u);
            return Phi(Cu, kappa, theta, sigma, rho, v0, s0, r, T);
        }
    }
}
