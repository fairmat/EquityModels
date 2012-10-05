/* Copyright (C) 2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Enrico Degiuli (enrico.degiuli@fairmat.com)
 *            Matteo Tesser  (matteo.tesser@fairmat.com) 
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

namespace Dupire
{
    public enum LocalVolatilityCalculation
    {
        /// <summary>
        /// todo:find better description
        /// </summary>
        Method1,
        /// <summary>
        /// Uses quantlib metho
        /// </summary>
        QuantLib,
    }

    /// <summary>
    /// Calass containing calibration options for the Dupire process 
    /// </summary>
    [SettingsContainer("Calibration Settings for Heston")]
    [Mono.Addins.Extension("/Fairmat/UserSettings")]
    [Serializable]
    public class DupireCalibrationSettings:IEstimationSettings
    {
        [SettingDescription("Local volatility calculation method")]
        public LocalVolatilityCalculation  LocalVolatilityCalculation= LocalVolatilityCalculation.Method1;
    }
}

