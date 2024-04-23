/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Enrico Degiuli (enrico.degiuli@fairmat.com)
 *            Matteo Tesser (matteo.tesser@fairmat.com)
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
using HestonEstimator;
using NUnit.Framework;

namespace Heston
{
    /// <summary>
    /// This test compares the price of a call calculated in c# with a benchmark value.
    /// </summary>
    [TestFixture]
    public class TestHestonCallPrice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void Test()
        {
            double k = 0.9;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 1.0;
            double v0 = 0.3;
            double rho = -0.8;

            // Calculates the theoretical value of the call.
            Vector param = new Vector(5);
            param[0] = kappa;
            param[1] = theta;
            param[2] = sigma;
            param[3] = rho;
            param[4] = v0;
            double fairmatPrice = HestonCall.HestonCallPrice(param, s0, tau, k, rate, dy);
            double tol = 1e-3;
            double benchmarkPrice = 0.339537359104676;

            Console.WriteLine("Theoretical Benchmark  Price = " + benchmarkPrice.ToString());
            Console.WriteLine("Theoretical Fairmat    Price = " + fairmatPrice);

            Assert.Less(Math.Abs(fairmatPrice - benchmarkPrice), tol);
        }
    }
}
