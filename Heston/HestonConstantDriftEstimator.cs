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
    public class HestonConstantDriftEstimator : CallEstimator, IEstimatorEx, IMenuItemDescription
    {
        /// <summary>
        /// If set to True, dividends are calibrated during the optimization process.
        /// </summary>
        static internal bool impliedDividends = false;

        #region IMenuItemDescription Members

        /// <summary>
        /// Gets the description of the implemented calibration function.
        /// </summary>
        public override string Description
        {
            get
            {
                return "Heston constant drift calibration";
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
        public override Type ProvidesTo
        {
            get
            {
                return typeof(Heston.HestonProcess);
            }
        }



        static double Mean(IFunction f, double a, double b)
        {
            double sum = 0;
            int N = 20;
            double h = b - a;
            double dt = h / N;
            for (int z = 0; z < N; z++)
                sum += f.Evaluate(a + (dt * z)) * dt;
            return sum / h;
        }



        protected override void Setup(EquityCalibrationData equityCalData, IEstimationSettings  settings)
        {
            //No need to used constant values. In the optimization procedure we can use the term structure
            //as it is.


            /*
             //equityCalData.dyFunc= LeastSquaresDividendCalibration(equityCalData.Hdata, equityCalData.zrFunc);

            if (settings != null) //uses settings if provided.
            {
                var localSettings = (HestonEstimationSettings)settings;

                equityCalData.SetToSpecificMaturity(localSettings.Maturity);
            }
            else
                equityCalData.SetToSpecificMaturity(1);
            */
            return;
        }

        double DY(EquityCalibrationData equityCalData)
        {
            double dy= 0.5*(equityCalData.dyFunc.Evaluate(1) + equityCalData.dyFunc.Evaluate(2));
            Console.WriteLine("Call/Put Parity Dividend\t" + dy);
            return dy;
        }

        protected override EstimationResult BuildEstimate(DiscountingCurveMarketData interestDataSet,CallPriceMarketData callDataSet, EquityCalibrationData equityCalData, SolutionInfo solution)
        {
            string[] names = new string[] { "S0", "kappa", "theta", "sigma",
                                            "rho", "V0", "r", "q" };
            Vector param = new Vector(8);
            param[0] = callDataSet.S0;
            param[Range.New(1, 5)] = solution.x[Range.New(0, 4)];
            param[6] = equityCalData.zrFunc.Evaluate(TheoreticalModelsSettings.ConstantDYRFMaturity);
            if (impliedDividends)
                param[7] = solution.x[Range.End];// equityCalData.dyFunc.Evaluate(TheoreticalModelsSettings.ConstantDYRFMaturity); 
            else
                param[7] = DY(equityCalData);
            
            var result = new EstimationResult(names, param);
            result.Fit = HestonCallOptimizationProblem.avgPricingError;
            Console.WriteLine(result);
            return result;
        }

        #endregion
    }
}
