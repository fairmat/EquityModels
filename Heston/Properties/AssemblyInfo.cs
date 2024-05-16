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

using System.Reflection;
using Mono.Addins;

[assembly: Addin("Heston model", "1.1.2", Category = "Stochastic Process")]
[assembly: AddinDependency("Fairmat", "1.0")]
[assembly: AddinAuthor("Fairmat SRL")]
[assembly: AddinDescription("The Heston model simulates equity or index prices taking into " +
                            "account stochastic volatility effects. The main feature of the " +
                            "model is that the price process follows a geometric brownian " +
                            "motion with a stochastic volatility while the volatility follows " +
                            "a square root mean reverting process. Usually the correlation is " +
                            "negative, so that a lowering in the stock price is correlated with " +
                            "an increasing in the volatility. Once installed the plug-in offers " +
                            "the possibility of using two new processes, the Heston process and " +
                            "the Heston time dependent drift process and to calibrate them to a " +
                            "series of call prices.")]
[assembly: AssemblyTrademark("Fairmat")]
[assembly: AssemblyCulture("")]
