/* Copyright (C) 2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Francesco Biondi (francesco.biondi@fairmat.com)
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

namespace HistoricalSimulator
{
    /// <summary>
    /// Operating mode describes how the plug-in must use historical data.
    /// </summary>
    public enum OperatingMode
    {
        /// <summary>
        /// Translates historical realizations forward starting from the start date.
        /// </summary>
        TranslateHistoricalRealizationsForward = 0,
		/// <summary>
		/// Simulate using historical increments
		/// </summary>
		Bootstrap=1,
    }
}
