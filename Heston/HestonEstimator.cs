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
    /// Implements and resolves the Heston optimization problem.
    /// </summary>
    [Mono.Addins.Extension("/Fairmat/Estimator")]
    public class CallEstimator : IEstimator, IDescription
    {
        #region IEstimator Members

        /// <summary>
        /// Gets the value requested by the interface ProvidesTo,
        /// returning HestonExtendedProcess as the type.
        /// </summary>
        public Type ProvidesTo
        {
            get
            {
                return typeof(HestonExtended.HestonExtendedProcess);
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
        public EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { new EstimateRequirement(typeof(DiscountingCurveMarketData)), 
                                               new EstimateRequirement(typeof(CallPriceMarketData)) };
        }

        /// <summary>
        /// Attempts to solve the Heston optimization problem using
        /// <see cref="Heston.HestonOptimizationProblem"/>.
        /// </summary>
        /// <param name="marketData">Data to be used in order to perform the optimization.</param>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="controller">IController.</param>
        /// <returns>The results of the optimization.</returns>
        public unsafe EstimationResult Estimate(List<object> marketData, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            DiscountingCurveMarketData interestDataSet = (DiscountingCurveMarketData)marketData[0];
            CallPriceMarketData callDataSet = (CallPriceMarketData)marketData[1];
            EquityCalibrationData equityCalData = new EquityCalibrationData(callDataSet, interestDataSet);

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
            o.controller = controller;
            o.NP = 40;
            o.MaxIter = 40;
            o.Verbosity = 1;

            // If true the optimization algorithm will operate in parallel.
            o.Parallel = Engine.MultiThread;

            // If true the objective function will be calculated in parallel.
            Engine.MultiThread = true;
            SolutionInfo solution = null;

            Vector x0 = new Vector(new double[] { 5.0, 0.1, 0.8, -0.7, 0.05 });

            // GA
            solution = solver.Minimize(problem, o, x0);
            if (solution.errors)
                return null;
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
            if (solution.errors)
                return null;

            string[] names = new string[] { "S0", "kappa", "theta", "sigma", "rho", "V0" };
            Vector param = new Vector(6);
            param[0] = callDataSet.S0;
            param[Range.New(1, 5)] = solution.x;
            EstimationResult result = new EstimationResult(names, param);

            // In the following the two function describing the ZR and dividend yields are created
            Matrix zerorate = new Matrix(interestDataSet.Durations.Length, 2);
            zerorate[Range.All, 0] = interestDataSet.Durations;
            zerorate[Range.All, 1] = interestDataSet.Values;

            Matrix dividendYield = new Matrix(equityCalData.MaturityDY.Length, 2);
            dividendYield[Range.All, 0] = equityCalData.MaturityDY;
            dividendYield[Range.All, 1] = equityCalData.DividendYield;

            result.Objects = new object[2];
            result.Objects[0] = zerorate;
            result.Objects[1] = dividendYield;
            return result;
        }

        #endregion

        public string Description
        {
            get { return "Calibrate against call options"; }
        }
    }
}
