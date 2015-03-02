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
using System.Runtime.InteropServices;
using Mono.Addins;

// The following lines tell that the assembly is an addin
[assembly: AssemblyVersion("1.1.1.0")]
[assembly: AssemblyFileVersion("1.1.1.0")]
[assembly: Mono.Addins.Addin("Heston model", "1.1.1", Category = "Stochastic Process")]
[assembly: Mono.Addins.AddinDependency("Fairmat", "1.0")]
[assembly: Mono.Addins.AddinAuthor("Fairmat SRL")]
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

// General Information about an assembly is controlled through the following
// set of attributes. Change these attribute values to modify the information
// associated with an assembly.
[assembly: AssemblyTitle("Heston")]
[assembly: AssemblyDescription("The Heston model simulates equity or index prices taking into " +
                               "account stochastic volatility effects. The main feature of the " +
                               "model is that the price process follows a geometric brownian " +
                               "motion with a stochastic volatility while the volatility follows " +
                               "a square root mean reverting process. Usually the correlation " +
                               "is negative, so that a lowering in the stock price is " +
                               "correlated with an increasing in the volatility. Once installed " +
                               "the plug-in offers the possibility of using two new processes, " +
                               "the Heston process and the Heston time dependent drift process " +
                               "and to calibrate them to a series of call prices.")]
[assembly: AssemblyConfiguration("")]
[assembly: AssemblyCompany("Fairmat SRL")]
[assembly: AssemblyCopyright("Copyright © Fairmat SRL 2009-2015")]
[assembly: AssemblyProduct("Heston")]
[assembly: AssemblyTrademark("Fairmat")]
[assembly: AssemblyCulture("")]

// Setting ComVisible to false makes the types in this assembly not visible
// to COM components.  If you need to access a type in this assembly from
// COM, set the ComVisible attribute to true on that type.
[assembly: ComVisible(false)]

// The following GUID is for the ID of the typelib if this project is exposed to COM
[assembly: Guid("0b44bce0-4c4b-4bea-a5c9-25df658972b9")]

