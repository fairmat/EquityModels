/* Copyright (C) 2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Matteo Tesser (matteo.tesser@fairmat.com)
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
using System.IO;
using System.Linq;
using DVPLDOM;
using DVPLI;

namespace HistoricalSimulator
{
    //[Mono.Addins.Extension("/Fairmat/Estimator")]
    public class HistoricalSimulatorCalibrator: IEstimatorEx
    {
        public IEstimationSettings DefaultSettings
        {
            get { throw new NotImplementedException(); }
        }

        public EstimationResult Estimate(List<object> data, IEstimationSettings settings = null, IController controller = null, Dictionary<string, object> properties = null)
        {
            EstimationResult r = new EstimationResult();
            r.Objects = new object[] { };
            return r;
        }

        public EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { new EstimateRequirement(typeof(DVPLI.MarketDataTypes.Scalar[])) };
        }

        public Type ProvidesTo
        {
            get { return typeof(HistoricalSimulator); }
        }
    }
}
