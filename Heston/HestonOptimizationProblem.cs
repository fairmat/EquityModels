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
using System.Threading.Tasks;
using DVPLI;
using DVPLI.TaskScheduler;
using Fairmat.MarketData;
using Fairmat.Optimization;

namespace HestonEstimator
{
    /// <summary>
    /// Describes the Heston optimization problem.
    /// </summary>
    public class HestonCallOptimizationProblem : IOptimizationProblem
    {
        /// <summary>
        /// The call market price matrix.
        /// </summary>
        internal Matrix callMarketPrice;

        /// <summary>
        /// The maturity vector relative to the callMarketPrice matrix.
        /// </summary>
        private Vector maturity;

        /// <summary>
        /// The strike vector relative to callMarketPrice matrix.
        /// </summary>
        private Vector strike;

        /// <summary>
        /// The rate vector.
        /// </summary>
        private Vector rate;

        /// <summary>
        /// The dividend yield vector.
        /// </summary>
        private Vector dividendYield;

        /// <summary>
        /// The drift vector.
        /// </summary>
        private Vector drift;

        /// <summary>
        /// Establish whether to use the boundary penalty function or not.
        /// </summary>
        public bool useBoundPenalty = false;

        /// <summary>
        /// Establish whether to use the Feller penalty function or not.
        /// </summary>
        public bool useFellerPenalty = false;

        /// <summary>
        /// Small value used in the boundary penalty function.
        /// </summary>
        private double smallValue = 1e-4;

        /// <summary>
        /// Value that weights the boundary penalty function.
        /// </summary>
        private double k1 = 1e10;

        /// <summary>
        /// Value that weights the Feller inequality penalty function.
        /// </summary>
        private double k2 = 1e6;

        /// <summary>
        /// The number of call option on which calibration is performed.
        /// </summary>
        internal int numCall;

        /// <summary>
        /// Process starting value.
        /// </summary>
        private double s0;

        /// <summary>
        /// Initializes a new instance of the HestonCallOptimizationProblem class using the
        /// EquityCalibrationData data structure.
        /// </summary>
        /// <param name="equityCalData">
        /// An EquityCalibrationData object containing market data for calibration.
        /// </param>
        /// <param name="matBound">
        /// A vector containing the minimum and maximum values
        /// for maturities to be used in calibration.
        /// </param>
        /// <param name="strikeBound">
        /// A vector containing the minimum and maximum values
        /// for strikes to be used in calibration.
        /// </param>
        public HestonCallOptimizationProblem(EquityCalibrationData equityCalData, Vector matBound, Vector strikeBound)
        {
            SetVariables(equityCalData.Hdata.CallPrice, equityCalData.Hdata.Maturity,
                         equityCalData.Hdata.Strike, equityCalData.CallMatrixRiskFreeRate,
                         equityCalData.CallMatrixDividendYield, equityCalData.Hdata.S0,
                         matBound, strikeBound);
        }

        /// <summary>
        /// Initializes a new instance of the HestonCallOptimizationProblem class.
        /// </summary>
        /// <param name="callMarketPrice">A matrix containing call option market prices.</param>
        /// <param name="maturity">
        /// Vector of call option maturities relative to callMarketPrice matrix.
        /// </param>
        /// <param name="strike">
        /// Vector of call option strikes relative to callMarketPrice matrix.
        /// </param>
        /// <param name="rate">
        /// Vector of zero coupon bond rates calculated relative to maturity vector maturities.
        /// </param>
        /// <param name="dividendYield">
        /// Vector of dividend yield rates calculated relative to maturity vector maturities.
        /// </param>
        /// <param name="s0">Index/Equity value at the time of calibration.</param>
        /// <param name="matBound">
        /// A vector containing the minimum and maximum values
        /// for maturities to be used in calibration.
        /// </param>
        /// <param name="strikeBound">
        /// A vector containing the minimum and maximum values
        /// for strikes to be used in calibration.
        /// </param>
        public HestonCallOptimizationProblem(Matrix callMarketPrice, Vector maturity, Vector strike, Vector rate, Vector dividendYield, double s0, Vector matBound, Vector strikeBound)
        {
            SetVariables(callMarketPrice, maturity, strike, rate,
                         dividendYield, s0, matBound, strikeBound);
        }

        /// <summary>
        /// Sets several variables used to solve the optimization problem.
        /// </summary>
        /// <param name="callMarketPrice">A matrix containing call option market prices.</param>
        /// <param name="maturity">
        /// Vector of call option maturities relative to callMarketPrice matrix.
        /// </param>
        /// <param name="strike">
        /// Vector of call option strikes relative to callMarketPrice matrix.
        /// </param>
        /// <param name="rate">
        /// Vector of zero coupon bond rates calculated relative to maturity vector maturities.
        /// </param>
        /// <param name="dividendYield">
        /// Vector of dividend yield rates calculated relative to maturity vector maturities.
        /// </param>
        /// <param name="s0">Index/Equity value at the time of calibration.</param>
        /// <param name="matBound">
        /// A vector containing the minimum and maximum values
        /// for maturities to be used in calibration.
        /// </param>
        /// <param name="strikeBound">
        /// A vector containing the minimum and maximum values
        /// for strikes to be used in calibration.</param>
        private void SetVariables(Matrix callMarketPrice, Vector maturity, Vector strike, Vector rate, Vector dividendYield, double s0, Vector matBound, Vector strikeBound)
        {
            Console.WriteLine("SetVariables");
            this.s0 = s0;
            this.rate = rate;
            this.dividendYield = dividendYield;
            Vector drift = this.rate - this.dividendYield;
            int[] matI = new int[2];
            int[] strikeI = new int[2];
            Console.WriteLine("maturity = " + maturity);
            Console.WriteLine("matBound = " + matBound);
            Console.WriteLine("strike = " + strike);
            Console.WriteLine("S0_ = " + s0);
            Console.WriteLine("strikeBound = " + strikeBound);
            matI = this.FindExtremes(maturity, matBound);
            strikeI = this.FindExtremes(strike, s0 * strikeBound);

            int numMat = matI[1] - matI[0] + 1;
            int numStrike = strikeI[1] - strikeI[0] + 1;
            this.maturity = new Vector(numMat);
            this.drift = new Vector(numMat);
            this.strike = new Vector(numStrike);
            this.callMarketPrice = new Matrix(numMat, numStrike) - 1;
            this.numCall = 0;
            Console.WriteLine("Call = " + callMarketPrice);
            Console.WriteLine("Mat = " + maturity);
            Console.WriteLine("drift = " + drift);
            for (int i = 0; (i < numMat) && ((matI[0] + i) < callMarketPrice.R); i++)
            {
                this.maturity[i] = maturity[matI[0] + i];
                this.drift[i] = drift[matI[0] + i];
                for (int j = 0; (j < numStrike) && ((strikeI[0] + j) < callMarketPrice.C); j++)
                {
                    this.callMarketPrice[i, j] = callMarketPrice[matI[0] + i, strikeI[0] + j];
                    if (this.callMarketPrice[i, j] != -1)
                        this.numCall++;
                }
            }

            for (int j = 0; j < numStrike; j++)
                this.strike[j] = strike[strikeI[0] + j];
        }

        /// <summary>
        /// Finds the couple of integer {i,j} so that the values
        /// vector[i], vector[i+1],..., vector[j-1], vector[j]
        /// are all included in the interval (bound[0], bound[1])
        /// while vector[i-1] and vector[j+1] are not.
        /// </summary>
        /// <param name='vector'>
        /// Vector to be filtered.
        /// </param>
        /// <param name='bounds'>
        /// Bounds to be used for the filtering.
        /// </param>
        /// <returns>
        /// Integer array representing the index couple.
        /// </returns>
        private int[] FindExtremes(Vector vector, Vector bounds)
        {
            int[] result = new int[2];
            int i = 0;

            while (vector[i] < bounds[0])
                i++;
            result[0] = i;
            i = vector.Length - 1;
            while (vector[i] > bounds[1])
                i--;
            result[1] = i;
            return result;
        }

        #region IOptimizationProblem Members

        /// <summary>
        /// Gets the bounds for the optimization.
        /// </summary>
        public Bounds Bounds
        {
            get
            {
                Bounds bounds = new Bounds();

                // The order of parameters is k, theta, sigma, rho, S0, V0.
                bounds.Lb = (Vector)new double[] { 0, 0, 0.01, -1, 0 };
                bounds.Ub = (Vector)new double[] { 10, 1, 5, 1, 1 };
                return bounds;
            }
        }

        /// <summary>
        /// This method is unused but part of the interface.
        /// </summary>
        /// <param name="x">The parameter is not used.</param>
        /// <returns>Null as it's unused.</returns>
        public DVPLI.Vector G(DVPLI.Vector x)
        {
            return null;
        }

        /// <summary>
        /// This method is unused but part of the interface.
        /// </summary>
        /// <param name="x">The parameter is not used.</param>
        /// <returns>Nothing the function always throws a NotImplementedException.</returns>
        public DVPLI.Vector Grad(DVPLI.Vector x)
        {
            throw new NotImplementedException();
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
        /// Gets null as we have no linear constrains defined.
        /// </summary>
        public virtual LinearConstraints LinearIneqConstraints
        {
            get
            {
                return null;
            }
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
        public double Obj(DVPLI.Vector x)
        {
            double sum = 0;

            HestonCall hc = new HestonCall(this, x, this.s0);
            if (Engine.MultiThread)
            {
                // Instantiate parallel computation if enabled.
                List<Task> tl = new List<Task>();

                // Context contains both input parameters and outputs.
                List<HestonCall> context = new List<HestonCall>();
                for (int r = 0; r < this.callMarketPrice.R; r++)
                {
                    context.Add(hc);
                    hc.T = this.maturity[r];
                    hc.rate = this.rate[r];
                    hc.dividend = this.dividendYield[r];
                    hc.row = r;
                    tl.Add(Task.Factory.StartNew(this.CalculateSingleRow, hc));
                }

                tsScheduler.WaitTaskList(tl);
                for (int r = 0; r < this.callMarketPrice.R; r++)
                    sum += context[r].sum;
            }
            else
            {
                // Sequential version of the code, used when parallel computation is disabled.
                for (int r = 0; r < this.callMarketPrice.R; r++)
                {
                    hc.T = this.maturity[r];
                    hc.rate = this.rate[r];
                    hc.dividend = this.dividendYield[r];
                    hc.row = r;
                    this.CalculateSingleRow(hc);
                    sum += hc.sum;
                }
            }

            sum = sum / ((double)this.numCall);
            if (this.useBoundPenalty)
                sum = sum + this.BoundPenalty(x);

            if (this.useFellerPenalty)
                sum = sum + this.FellerPenalty(x);
            return sum;
        }

        /// <summary>
        /// Calculates a single row of quadratic difference matrix.
        /// </summary>
        /// <param name='context'>
        /// An object of type <see cref="HestonCall"/> containing the context.
        /// </param>
        private void CalculateSingleRow(object context)
        {
            HestonCall hc = context as HestonCall;
            int r = hc.row;
            for (int c = 0; c < this.callMarketPrice.C; c++)
            {
                if (this.callMarketPrice[r, c] != -1)
                {
                    hc.K = this.strike[c];
                    hc.hestonCallPrice[r, c] = hc.HestonCallPrice();
                    hc.sum += Math.Pow(hc.hestonCallPrice[r, c] - this.callMarketPrice[r, c], 2);
                }
            }

            return;
        }

        #endregion

        /// <summary>
        /// Penalty function relative to bounds.
        /// </summary>
        /// <param name='x'>
        /// Vector of parameters.
        /// </param>
        /// <returns>
        /// Penalty value.
        /// </returns>
        private double BoundPenalty(Vector x)
        {
            Vector t1, t2, t3;
            t1 = Bounds.Lb + this.smallValue - x;
            t2 = -Bounds.Ub + this.smallValue + x;

            t3 = t1 * t1 * (t1 > 0) + t2 * t2 * (t2 > 0);
            return this.k1 * t3.Sum();
        }

        /// <summary>
        /// Penalty function relative to Feller condition.
        /// </summary>
        /// <param name='x'>
        /// Vector of parameters.
        /// </param>
        /// <returns>
        /// Penalty value.
        /// </returns>
        private double FellerPenalty(Vector x)
        {
            double result;
            result = Math.Max(0, x[2] * x[2] - 2 * x[0] * x[1]);
            return this.k2 * result * result;
        }
    }
}
