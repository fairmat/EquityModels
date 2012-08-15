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
using DVPLI;

namespace HistoricalSimulator
{
    /// <summary>
    /// Implements simulation by bootstrapping historical data.
    /// </summary>
    internal class Bootstrap
    {
        private Vector[] returns;

        internal Bootstrap(List<Tuple<DateTime, Vector>> samples)
        {
            CreateReturnSamples(samples);
        }

        // This could be done in setup.
        internal void CreateReturnSamples(List<Tuple<DateTime, Vector>> samples)
        {
            this.returns = new Vector[samples.Count - 1];
            for (int z = 1; z < samples.Count; z++)
            {
                // Finds normalization coefficient.
                TimeSpan delta = samples[z].Item1 - samples[z - 1].Item1;
                double dt = delta.TotalDays / 365;
                this.returns[z - 1] = Vector.Log(samples[z].Item2 / samples[z - 1].Item2) / dt;
            }
        }

        /// <summary>
        /// Simulate a realization of the stochastic process by bootstrapping historical data.
        /// </summary>
        /// <remarks>
        /// This function is called once for realization.
        /// </remarks>
        /// <param name="dates">
        /// The dates (in years fractions) at which the process must be simulated.
        /// </param>
        /// <param name="outDynamic">Where the dynamic should be written.</param>
        internal void Simulate(double[] dates, IMatrixSlice outDynamic)
        {
            // Number of components.
            int components = this.returns[0].Length;
            for (int c = 0; c < components; c++)
                outDynamic[0, c] = 1;

            // Assumes equispaced dates, otherwise we can normalize returns.
            for (int i = 1; i < dates.Length; i++)
            {
                double dt = dates[i] - dates[i - 1];

                // Selects a vector of returns.
                double u = Engine.Generator.Uniform();
                int z = (int)(u * this.returns.Length);
                for (int c = 0; c < components; c++)
                    outDynamic[i, c] = outDynamic[i - 1, c] * Math.Exp(this.returns[z][c] * dt);
            }
        }
    }
}
