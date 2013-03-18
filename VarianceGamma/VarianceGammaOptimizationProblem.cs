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
using DVPLI;
using Fairmat.Optimization;

namespace VarianceGamma
{
    /// <summary>
    /// Describes the Variance Gamma Optimization problem.
    /// </summary>
    public class VarianceGammaOptimizationProblem : IOptimizationProblem
    {
        /// <summary>
        /// The dividend yield.
        /// </summary>
        private double q;

        /// <summary>
        /// The actual value of stock.
        /// </summary>
        private double s0;

        /// <summary>
        /// The risk free rate.
        /// </summary>
        private double r;

        /// <summary>
        /// The strike Prices (unrolled).
        /// </summary>
        private Vector k;

        /// <summary>
        /// The maturities (unrolled).
        /// </summary>
        private Vector m;

        /// <summary>
        /// The observed call prices (unrolled).
        /// </summary>
        private Matrix cp;

        /// <summary>
        /// Initializes a new instance of the VarianceGammaOptimizationProblem class.
        /// </summary>
        /// <param name="q">The dividend yield.</param>
        /// <param name="s0">The actual value of stock.</param>
        /// <param name="k">The strike Prices (unrolled).</param>
        /// <param name="r">The risk free rate.</param>
        /// <param name="cp">The observed call prices (unrolled).</param>
        /// <param name="m">The Maturities (unrolled).</param>
        public VarianceGammaOptimizationProblem(double q, double s0, Vector k, double r, Matrix cp, Vector m)
        {
            this.q = q;
            this.s0 = s0;
            this.k = k;
            this.r = r;
            this.cp = cp;
            this.m = m;
        }

        /// <summary>
        /// Calibration objective function.
        /// </summary>
        /// <param name='x'>
        /// The vector of parameters.
        /// </param>
        /// <returns>
        /// Objective function value.
        /// </returns>
        public double Obj(Vector x)
        {
            return VGDiff(x, this.q, this.s0, this.k, this.r, this.cp, this.m);
        }

        /// <summary>
        /// This method is unused but part of the interface.
        /// </summary>
        /// <param name="x">The parameter is not used.</param>
        /// <exception cref="NotImplementedException">
        /// The exception is always thrown upon calling this method.
        /// </exception>
        /// <returns>Nothing. The function always throws a NotImplementedException.</returns>
        public Vector Grad(Vector x)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Gets the bounds for the optimization.
        /// </summary>
        public Bounds Bounds
        {
            get
            {
                Bounds b = new Bounds();
                b.Lb = new Vector() { -0.2, 0.01, 0.01 };
                b.Ub = new Vector() { 0.2, 0.2, 1 };
                return b;
            }
        }

        /// <summary>
        /// Gets a value indicating whether there are non linear constraints in this
        /// optimization problem. In this case there are not.
        /// </summary>
        public bool HasNonLinearConstraints
        {
            get
            {
                return false;
            }
        }

        /// <summary>
        /// Gets The linear constrains for Variance Gamma.
        /// Being A: [0 0 0] and b: [0].
        /// </summary>
        public LinearConstraints LinearIneqConstraints
        {
            get
            {
                LinearConstraints linearInq = new LinearConstraints();
                linearInq.A = new Matrix(1, 3);
                linearInq.b = new Vector(1);
                linearInq.A[0, 0] = 0;
                linearInq.A[0, 1] = 0;
                linearInq.A[0, 2] = 0;
                linearInq.b[0] = 0;
                return linearInq;
            }
        }

        /// <summary>
        /// This method is unused but part of the interface.
        /// </summary>
        /// <param name="x">The parameter is not used.</param>
        /// <exception cref="NotImplementedException">
        /// The exception is always thrown upon calling this method.
        /// </exception>
        /// <returns>Nothing. The function always throws a NotImplementedException.</returns>
        public Vector G(Vector x)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Returns the root mean square error between market and model call prices
        /// </summary>
        /// <param name="x">Parameter vector</param>
        /// <param name="q">Dividend yield</param>
        /// <param name="s0">Index starting value</param>
        /// <param name="k">Call strike vector</param>
        /// <param name="r">Short rate</param>
        /// <param name="cp">Call price matrix</param>
        /// <param name="m">Call maturity vector</param>
        /// <returns>Call price root mean square error</returns>
        public static double VGDiff(Vector x, double q, double s0, Vector k, double r, Matrix cp, Vector m)
        {
            double y = 0;
            for (int i = 0; i < m.Length; i++)
            {
                for (int j = 0; j < k.Length; j++)
                {
                    double residual = Math.Pow(cp[i, j] - VarianceGammaOptionsCalibration.VGCall(x[0], x[1], x[2], m[i], k[j], q, s0, r), 2);
                    if (residual > Math.Pow(10, 10))
                    {
                        y = y + 1000 * Math.Sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
                    }
                    else
                    {
                        y = y + residual;
                    }
                }
            }

            return Math.Sqrt(y / (m.Length * k.Length));
        }
    }
}
