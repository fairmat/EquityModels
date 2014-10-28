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
        public virtual Type ProvidesTo
        {
            get
            {
                return typeof(HestonExtended.HestonExtendedProcess);
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
            return new EstimateRequirement[] { new EstimateRequirement(typeof(DiscountingCurveMarketData),MarketRequirement.TickerMarket), 
                                               new EstimateRequirement(typeof(CallPriceMarketData)) };
        }

        protected virtual void Setup(EquityCalibrationData equityCalData, IEstimationSettings settings)
        {
           
        }

        protected virtual HestonCallOptimizationProblem NewOptimizationProblem(EquityCalibrationData equityCalData, Vector matBound, Vector strikeBound)
        {
            return new HestonCallOptimizationProblem(equityCalData, matBound, strikeBound);
        }


        /// <summary>
        /// Attempts to solve the Heston optimization problem using
        /// <see cref="Heston.HestonOptimizationProblem"/>.
        /// </summary>
        /// <param name="marketData">Data to be used in order to perform the optimization.</param>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="controller">IController.</param>
        /// <returns>The results of the optimization.</returns>
        public EstimationResult Estimate(List<object> marketData, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            DateTime t0 = DateTime.Now;
            var interestDataSet = (CurveMarketData)marketData[0];
            CallPriceMarketData callDataSet = (CallPriceMarketData)marketData[1];
            EquityCalibrationData equityCalData = new EquityCalibrationData(callDataSet, interestDataSet);
            

            Setup(equityCalData, settings);

            // Creates the context.
            Document doc = new Document();
            ProjectROV prj = new ProjectROV(doc);
            doc.Part.Add(prj);

            // Optimization problem instance.
            Vector matBound = new Vector(2);
            Vector strikeBound = new Vector(2);
            matBound[0] = 0.0;// .25;
            matBound[1] = 6;// 10; //Up to 6Y maturities
            strikeBound[0] = 0.7;// 0.5;
            strikeBound[1] = 1.3;//1.5;

            HestonCallOptimizationProblem problem = NewOptimizationProblem(equityCalData, matBound, strikeBound);
            int totalOpts = problem.numCall + problem.numPut;
            Console.WriteLine("Calibration based on "+totalOpts+ " options. (" + problem.numCall + " call options and "+problem.numPut+" put options).");

            //IOptimizationAlgorithm solver = new  QADE();
            IOptimizationAlgorithm solver = new MultiLevelSingleLinkage();
            IOptimizationAlgorithm solver2 = new SteepestDescent();

            DESettings o = new DESettings();
            o.controller = controller;
            
            // If true the optimization algorithm will operate in parallel.
            o.Parallel = Engine.MultiThread;
            o.h = 10e-7;
            o.epsilon = 10e-7;
          
            SolutionInfo solution = null;

            double minObj=double.MaxValue;
            Vector minX= null;
            int Z = 1;
            //if (problem.GetType() == typeof(Heston.HestonCallSimulationOptimizationProblem))
            //    Z = 2;

            for(int z=0;z<Z;z++)
            {
                o.NP = 50;
                o.MaxIter =  25;
                o.Verbosity = 1;

            Vector x0 = null;// new Vector(new double[] { 0.5, 0.5, 0.8, -0.5, 0.05 });

            // GA
            solution = solver.Minimize(problem, o, x0);
            if (solution.errors)
                return null;
            
            o.options = "qn";
            o.MaxIter = 500;// 1000;

            if (solution != null)
                solution = solver2.Minimize(problem, o, solution.x);
            else
            {
                solution = solver2.Minimize(problem, o, x0);
            }
            if (solution.errors)
                return null;

            if (solution.obj < minObj)
            {
                minObj = solution.obj;
                minX = solution.x.Clone();
            }
            }


            
            solution.obj = minObj;
            solution.x = minX;

            //Displays pricing error structure
            HestonCallOptimizationProblem.displayPricingError = true;
            problem.Obj(solution.x);
            HestonCallOptimizationProblem.displayPricingError = false;
            Console.WriteLine("Calibration Time (s)\t" + (DateTime.Now - t0).TotalSeconds);

            return BuildEstimate(interestDataSet, callDataSet, equityCalData, solution);
        }

        protected virtual EstimationResult BuildEstimate(CurveMarketData interestDataSet, CallPriceMarketData callDataSet, EquityCalibrationData equityCalData, SolutionInfo solution)
        {
            string[] names = new string[] { "S0", "kappa", "theta", "sigma", "rho", "V0" };
            Vector param = new Vector(6);
            param[0] = callDataSet.S0;
            param[Range.New(1, 5)] = solution.x;
            var result = new EstimationResult(names, param);

            // In the following the two function describing the ZR and dividend yields are created
            //Matrix zerorate = new Matrix(interestDataSet.Durations.Length, 2);
            //zerorate[Range.All, 0] = interestDataSet.Durations;
            //zerorate[Range.All, 1] = interestDataSet.Values;

            //Matrix dividendYield = new Matrix(equityCalData.MaturityDY.Length, 2);
            //dividendYield[Range.All, 0] = equityCalData.MaturityDY;
            //dividendYield[Range.All, 1] = equityCalData.DividendYield;

            Matrix zerorate = new Matrix((equityCalData.zrFunc as PFunction).Expr);
            Matrix dividendYield = new Matrix((equityCalData.dyFunc as PFunction).Expr);
            result.Objects = new object[2];
            result.Objects[0] = zerorate;
            result.Objects[1] = dividendYield;
            result.Fit = solution.obj;
            Console.WriteLine(result);
            return result;
        }

        #endregion

        public virtual string Description
        {
            get { return "Calibrate against options (closed form)"; }
        }



    }
}
