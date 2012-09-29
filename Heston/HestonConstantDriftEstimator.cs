/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Matteo Tesser (matteo.tesser@fairmat.com)
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
using DVPLDOM;
using DVPLI;
using Fairmat.MarketData;
using Fairmat.Optimization;

namespace HestonEstimator
{
    /// <summary>
    /// Implements and resolves the Heston optimization problem with a constant drift.
    /// </summary>
    [Mono.Addins.Extension("/Fairmat/Estimator")]
    public class HestonConstantDriftEstimator : IEstimatorEx, IMenuItemDescription
    {
        #region IMenuItemDescription Members

        /// <summary>
        /// Gets the description of the implemented calibration function.
        /// </summary>
        public string Description
        {
            get
            {
                return "Heston constant drift calibration";
            }
        }

        /// <summary>
        /// Gets the tooltip for the implemented calibration function.
        /// </summary>
        public string ToolTipText
        {
            get
            {
                return "Calibrates Heston fixing risk free rate and " +
                       "dividend yield given by the specified maturity";
            }
        }

        #endregion IMenuItemDescription Members

        #region IEstimatorEx

        /// <summary>
        /// Gets the default settings to be used by the estimator.
        /// </summary>
        public IEstimationSettings DefaultSettings
        {
            get
            {
                return UserSettings.GetSettings(typeof(HestonEstimationSettings)) as HestonEstimationSettings;
            }
        }

        #endregion IEstimatorEx

        #region IEstimator Members

        /// <summary>
        /// Gets the value requested by the interface ProvidesTo,
        /// returning HestonProcess as the type.
        /// </summary>
        public Type ProvidesTo
        {
            get
            {
                return typeof(Heston.HestonProcess);
            }
        }

        /// <summary>
        /// Gets the types required by the estimator in order to work:
        /// InterestRateMarketData and CallPriceMarketData are
        /// the required type for this estimator.
        /// </summary>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="multivariateRequest">The parameter is not used.</param>
        /// <returns>
        /// An array containing the type InterestRateMarketData and CallPriceMarketData.
        /// </returns>
        public Type[] GetRequirements(IEstimationSettings settings, bool multivariateRequest)
        {
            return new Type[] { typeof(InterestRateMarketData), typeof(CallPriceMarketData) };
        }

        /// <summary>
        /// Attempts to solve the Heston optimization problem using
        /// <see cref="Heston.HestonOptimizationProblem"/> with a constant drift.
        /// </summary>
        /// <param name="marketData">
        /// The data to be used in order to perform the optimization.
        /// </param>
        /// <param name="settings">The settings to be used for the optimization.</param>
        /// <returns>The results of the optimization.</returns>
        public unsafe EstimationResult Estimate(List<object> marketData, IEstimationSettings settings)
        {
            InterestRateMarketData interestDataSet = (InterestRateMarketData)marketData[0];
            CallPriceMarketData callDataSet = (CallPriceMarketData)marketData[1];
            EquityCalibrationData equityCalData = new EquityCalibrationData(callDataSet, interestDataSet);
            HestonEstimationSettings localSettings = (HestonEstimationSettings)settings;

            equityCalData.SetToSpecificMaturity(localSettings.Maturity);

            // Creates the context.
            Document doc = new Document();
            ProjectROV prj = new ProjectROV(doc);
            doc.Part.Add(prj);

            // Optimization problem instance.
            Vector matBound = new Vector(2);
            Vector strikeBound = new Vector(2);
            matBound[0] = 0.0;
            matBound[1] = 2.0;
            strikeBound[0] = 0.79;
            strikeBound[1] = 1.21;

            HestonCallOptimizationProblem problem = new HestonCallOptimizationProblem(equityCalData, matBound, strikeBound);
            Console.WriteLine("Optimization based on " + problem.numCall + " call options");

            IOptimizationAlgorithm solver = new QADE();
            IOptimizationAlgorithm solver2 = new SteepestDescent();

            DESettings o = new DESettings();

            o.NP = 40;
            o.MaxIter = 40;
            o.Verbosity = 1;

            // If true the optimization algorithm will operate in parallel.
            o.Parallel = true;

            // If true the objective function will be calculated in parallel.
            Engine.MultiThread = true;
            SolutionInfo solution = null;

            Vector x0 = new Vector(new double[] { 5.0, 0.1, 0.8, -0.7, 0.05 });

            EstimationResult result = new EstimationResult();

            // GA
            solution = solver.Minimize(problem, o, x0);

            o.epsilon = 10e-8;
            o.options = "qn";
            o.h = 10e-8;

            o.MaxIter = 2000;
            if (solution != null)
                solution = solver2.Minimize(problem, o, solution.x);
            else
            {
                solution = solver2.Minimize(problem, o, x0);
            }

            string[] names = new string[] { "S0", "kappa", "theta", "sigma",
                                            "rho", "V0", "r", "q" };
            Vector param = new Vector(8);
            param[0] = callDataSet.S0;
            param[Range.New(1, 5)] = solution.x;
            param[6] = equityCalData.Rate[0];
            param[7] = equityCalData.DividendYield[0];
            result = new EstimationResult(names, param);
            return result;
        }

        #endregion
    }
}
