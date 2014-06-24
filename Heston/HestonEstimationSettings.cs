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
using DVPLI;

namespace HestonEstimator
{
    /// <summary>
    /// Maintains some settings to be used by the Heston calibrator.
    /// These settings will be available from the plugin settings tab in Fairmat.
    /// </summary>
    [SettingsContainer("Calibration Settings for Heston",SettingType= SettingType.Calibration)]
    [Mono.Addins.Extension("/Fairmat/UserSettings")]
    [Serializable]
    public class HestonEstimationSettings : IEstimationSettings
    {
        /// <summary>
        /// Maintains the value of the maturity to fix risk free rate and dividend yield,
        /// this is used by <see cref="HestonConstantDriftEstimator"/>.
        /// </summary>
        [RangeSettingDescription("Maturity to fix risk free rate and dividend yield", 0.0, 10)]
        public double Maturity=1;
    }
}
