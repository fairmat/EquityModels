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
using HestonEstimator;
using NUnit.Framework;

namespace Heston
{
    /// <summary>
    /// Directly unit tests the internals of <see cref="HestonNumericalGreeks"/>:
    /// the shared finite-difference core (<see cref="HestonNumericalGreeks.GreeksBumper"/>)
    /// and a couple of branch splits not exercised by the higher-level tests in
    /// TestHestonCallPrice.cs.
    /// </summary>
    [TestFixture]
    public class TestHestonNumericalGreeksInternals
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestGreeksBumperCentralDifferenceMatchesAnalyticDerivative()
        {
            // f(x) = x^2, f'(x) = 2x. At x=5 the analytic derivative is 10.
            Func<double, double> f = x => x * x;
            double result = HestonNumericalGreeks.GreeksBumper(bumpPercentage: 0.0001, init: 5.0, fun: f);

            Assert.AreEqual(10.0, result, 1e-3);
        }

        [Test]
        public void TestGreeksBumperOneSidedDifferenceWithUnbumpedPrice()
        {
            // Same function, but forcing the one-sided (forward) difference branch by
            // supplying unbumpedPrice = f(init).
            Func<double, double> f = x => x * x;
            double init = 5.0;
            double result = HestonNumericalGreeks.GreeksBumper(bumpPercentage: 0.0001, init: init, fun: f, unbumpedPrice: f(init));

            Assert.AreEqual(10.0, result, 1e-2);
        }

        [Test]
        public void TestThetaFSPCallTpNullVsSuppliedBranches()
        {
            double kappa = 2.5, theta = 0.4, sigma = 0.2, rho = -0.8, v0 = 0.3;
            double s0 = 100.0, K = 90, T = 2.0, T0 = 0.01, r = 0.1, q = 0.07;

            Engine.Verbose = 0;

            double thetaWithoutTp = HestonNumericalGreeks.ThetaFSPCall(
                kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma,
                s0: s0, T: T, K: K, r: r, q: q, T0: T0);

            double thetaWithTp = HestonNumericalGreeks.ThetaFSPCall(
                kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma,
                s0: s0, T: T, K: K, r: r, q: q, T0: T0, Tp: T);

            Console.WriteLine("Theta (Tp=null) = " + thetaWithoutTp);
            Console.WriteLine("Theta (Tp=T)    = " + thetaWithTp);

            Assert.IsFalse(double.IsNaN(thetaWithoutTp));
            Assert.IsFalse(double.IsNaN(thetaWithTp));
        }

        [Test]
        public void TestRhoCallCustomDiscountingFactorFunctionBranch()
        {
            double kappa = 2.5, theta = 0.4, sigma = 0.2, rho = -0.8, v0 = 0.3;
            double s0 = 100.0, K = 90, T = 2.0, r = 0.1, q = 0.07;

            // A non-flat curve that nevertheless matches the flat rate r at T (so the price
            // itself is identical for both branches; only the code path differs).
            double[] tenors = new double[] { 0.1, 0.2, 0.5, 1, 2, 5 };
            double[] rates = new double[] { 0.01, 0.02, 0.05, r, r, 0.2 };
            var curve = new DVPLDOM.PFunction((Vector)tenors, (Vector)rates);
            curve.m_Function.iType = DVPLUtils.EInterpolationType.ZERO_ORDER;
            Func<double, double, double> nonFlatDiscounting = (t, Tt) => Math.Exp(-curve.Evaluate(Tt) * Tt);

            double rhoDefault = HestonNumericalGreeks.RhoCall(
                kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma,
                s0: s0, T: T, K: K, r: r, q: q);

            double rhoCustom = HestonNumericalGreeks.RhoCall(
                kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma,
                s0: s0, T: T, K: K, r: r, q: q, discountingFactorFunction: nonFlatDiscounting);

            Console.WriteLine("Rho (default flat curve)   = " + rhoDefault);
            Console.WriteLine("Rho (custom supplied curve) = " + rhoCustom);

            // RhoCall's DiscountingFunctionAsFunctionOfIR closure re-derives a spread from
            // the supplied curve and holds it fixed while bumping r, so Rho ends up
            // independent of the curve's shape -- this pins down that (documented-by-code,
            // not obvious) invariance while still exercising the custom-function branch.
            Assert.IsFalse(double.IsNaN(rhoDefault));
            Assert.IsFalse(double.IsNaN(rhoCustom));
            Assert.AreEqual(rhoDefault, rhoCustom, 1e-6);
        }
    }
}
