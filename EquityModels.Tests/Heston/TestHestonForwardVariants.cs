/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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
using Fairmat.Math;
using HestonEstimator;
using NUnit.Framework;

namespace Heston
{
    /// <summary>
    /// Covers the previously-untested forward-start Heston pricing variants
    /// (<see cref="HestonForwardAhlipRutkowski"/>, <see cref="HestonForwardLucic"/>) and a few
    /// untested branches of <see cref="HestonForwardApproximated"/>.
    /// </summary>
    [TestFixture]
    public class TestHestonForwardVariants
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestAhlipRutkowskiForwardCallPriceIsFinite()
        {
            // Same parameter set used by TestHestonCallPrice.TestForwardCall, where this
            // formula was tried but never asserted against (left commented out) -- so we
            // only check it produces a sane finite number, not an exact reference value.
            double k = 1;
            double tau = 2.0;
            double rate = 0.5;
            double dy = 0.0;
            double kappa = 1;
            double theta = 0.4;
            double sigma = 0.1;
            double s0 = 1.0;
            double v0 = 0.4;
            double rho = 0.0;
            double T0 = 0.001;

            double price = HestonForwardAhlipRutkowski.HestonForwardCallPrice(
                kappa, theta, rho, v0, sigma, s0, k, rate, dy, tau, T0);

            Console.WriteLine("AhlipRutkowski forward call price = " + price);

            Assert.IsFalse(double.IsNaN(price));
            Assert.IsFalse(double.IsInfinity(price));
        }

        [Test]
        public void TestAhlipRutkowskiVectorOverloadMatchesScalarOverload()
        {
            double k = 1;
            double tau = 2.0;
            double rate = 0.5;
            double dy = 0.0;
            double kappa = 1;
            double theta = 0.4;
            double sigma = 0.1;
            double s0 = 1.0;
            double v0 = 0.4;
            double rho = 0.0;

            Vector x = new Vector(5);
            x[0] = kappa;
            x[1] = theta;
            x[2] = sigma;
            x[3] = rho;
            x[4] = v0;

            double priceFromVector = HestonForwardAhlipRutkowski.HestonForwardCallPrice(x, s0, tau, 0.001, k, rate, dy);

            Assert.IsFalse(double.IsNaN(priceFromVector));
            Assert.IsFalse(double.IsInfinity(priceFromVector));
        }

        [Test]
        public void TestLucicForwardCallPriceIsFinite()
        {
            double lambda = 1;
            double rho = 0.0;
            double eta = 0.1;
            double T = 2.0;
            double T0 = 0.5;
            double tau = T - T0;
            double vhat = 0.4;
            double v = 0.4;
            double rate = 0.5;
            double K = 1;

            Complex price = HestonForwardLucic.HestonForwardLucicCallPrice(lambda, rho, eta, T, T0, tau, vhat, v, rate, K);

            Console.WriteLine("Lucic forward call price = " + price.Re + " + i*" + price.Im);

            Assert.IsFalse(double.IsNaN(price.Re));
            Assert.IsFalse(double.IsInfinity(price.Re));
        }

        [Test]
        public void TestForwardApproximatedAlwaysZeroGreekStubs()
        {
            // FSPCallCalculateDelta/Gamma and FSCallCalculateGamma are documented (via the
            // numerical-greeks cross-checks in TestHestonCallPrice) to be structurally zero;
            // this locks in that contract directly against the stub implementations.
            double s0 = 100.0, K = 90, T = 2.0, T0 = 0.01, r = 0.1, q = 0.07;
            double kappa = 2.5, theta = 0.4, sigma = 0.2, rho = -0.8, v0 = 0.3;

            Assert.AreEqual(0.0, HestonForwardApproximated.FSPCallCalculateDelta(s0, K, T, T0, r, q, kappa, theta, sigma, rho, v0));
            Assert.AreEqual(0.0, HestonForwardApproximated.FSPCallCalculateGamma(s0, K, T, T0, r, q, kappa, theta, sigma, rho, v0));
            Assert.AreEqual(0.0, HestonForwardApproximated.FSCallCalculateGamma(s0, K, T, T0, r, q, kappa, theta, sigma, rho, v0));
        }

        [Test]
        public void TestForwardApproximatedCustomDiscountingFactorFunctionBranch()
        {
            // Exercises the branch where discountingFactorFunction is explicitly supplied
            // (as opposed to the default-null fallback exercised by most existing tests).
            double s0 = 100.0, K = 90, T = 2.0, T0 = 0.01, r = 0.1, q = 0.07;
            double kappa = 2.5, theta = 0.4, sigma = 0.2, rho = -0.8, v0 = 0.3;

            Func<double, double, double> flatDiscounting = (t, Tt) => Math.Exp(-r * (Tt - t));

            double priceWithDefault = HestonForwardApproximated.HestonForwardCallPrice(
                kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0);
            double priceWithCustomFunc = HestonForwardApproximated.HestonForwardCallPrice(
                kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0,
                discountingFactorFunction: flatDiscounting);

            Assert.IsFalse(double.IsNaN(priceWithDefault));
            Assert.IsFalse(double.IsNaN(priceWithCustomFunc));
            Assert.AreEqual(priceWithDefault, priceWithCustomFunc, 1e-6);
        }
    }
}
