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


        protected override void Setup(EquityCalibrationData equityCalData, IEstimationSettings  settings)
        {
            var localSettings = (HestonEstimationSettings)settings;

            equityCalData.SetToSpecificMaturity(localSettings.Maturity);
        }


        protected override EstimationResult BuildEstimate(DiscountingCurveMarketData interestDataSet,CallPriceMarketData callDataSet, EquityCalibrationData equityCalData, SolutionInfo solution)
        {
            EstimationResult result = new EstimationResult();
            string[] names = new string[] { "S0", "kappa", "theta", "sigma",
                                            "rho", "V0", "r", "q" };
            Vector param = new Vector(8);
            param[0] = callDataSet.S0;
            param[Range.New(1, 5)] = solution.x;
            param[6] = equityCalData.ShortRate;
            param[7] = equityCalData.DividendYield[0];
            result = new EstimationResult(names, param);
            return result;
        }

        #endregion
    }
}
