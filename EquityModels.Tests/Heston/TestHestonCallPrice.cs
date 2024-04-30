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
using Accord.Math;
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
            DVPLI.Vector param = new DVPLI.Vector(5);
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


        [Test]
        public void TestComparisonCarrMadan()
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
            
            double fairmatPrice = HestonCall.HestonCallPrice(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double carrMadanPrice = HestonCall.HestonCallPriceCarrMadan(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double tol = 1e-3;
            

            Console.WriteLine("Theoretical CarrMadan  Price = " + carrMadanPrice);
            Console.WriteLine("Theoretical Fairmat    Price = " + fairmatPrice);

            Assert.Less(Math.Abs(fairmatPrice - carrMadanPrice), tol);
        }



        [Test]
        public void TestQuantLib()
        {
            double k = 90;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;

            // Calculates the theoretical value of the call.

            double fairmatPrice = HestonCall.HestonCallPrice(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double quantlibPrice = 34.08981109547019;

            var relativeDifference = Math.Abs(fairmatPrice - quantlibPrice) / quantlibPrice;
            Console.WriteLine("Theoretical QuantLib  Price = " + quantlibPrice);
            Console.WriteLine("Theoretical Fairmat    Price = " + fairmatPrice);
            Assert.Less(relativeDifference, 1e-2);
        }


        [Test]
        public void TestGreeksCall()
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

            // Calculates the greeks.
            Engine.Verbose = 0;

            var analyticalDelta = HestonDelta.DeltaCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            var analyticalGamma = HestonGamma.GammaCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            var analyticalRho = HestonRho.RhoCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            var analyticalVega = HestonVega.VegaCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);


            (double delta, double gamma) = HestonNumericalGreeks.DeltaGammaCall(bumpPercentage : 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedDelta = 0.6459985446;
            Assert.AreEqual(expectedDelta, delta, 1e-3);
            Assert.AreEqual(analyticalDelta, delta, 1e-3);



            double expectedGamma = 0.3264832;
            Assert.AreEqual(expectedGamma, gamma, 1e-3);
            Assert.AreEqual(analyticalGamma, gamma, 1e-3);


            double rhoGreek = HestonNumericalGreeks.RhoCall(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedRho = 0.61272535;
            Assert.AreEqual(expectedRho, rhoGreek, 1e-3);
            Assert.AreEqual(analyticalRho, rhoGreek, 1e-3);


            double vega = HestonNumericalGreeks.VegaCall(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedVega = 0.06372181708;
            Assert.AreEqual(expectedVega, vega, 1e-3);
            // Assert.AreEqual(analyticalVega, vega, 1e-3);


            double thetaGreek = HestonNumericalGreeks.ThetaCall(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);  
            double expectedTheta = 0.0494223;
            Assert.AreEqual(expectedTheta, thetaGreek, 1e-3);


        }

        [Test]
        public void TestGreeksPut()
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

            // Calculates the greeks.
            Engine.Verbose = 0;

            var analyticalDelta = HestonDelta.DeltaPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            var analyticalGamma = HestonGamma.GammaPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            var analyticalRho = HestonRho.RhoPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            var analyticalVega = HestonVega.VegaPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);


            (double delta, double gamma) = HestonNumericalGreeks.DeltaGammaPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedDelta = -0.22334975741;
            Assert.AreEqual(expectedDelta, delta, 1e-3);
            Assert.AreEqual(analyticalDelta, delta, 1e-3);

            double expectedGamma = 0.3264832;
            Assert.AreEqual(expectedGamma, gamma, 1e-3);
            Assert.AreEqual(analyticalGamma, gamma, 1e-3);


            double rhoGreek = HestonNumericalGreeks.RhoPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedRho = -0.86098999518;
            Assert.AreEqual(expectedRho, rhoGreek, 1e-3);
            Assert.AreEqual(analyticalRho, rhoGreek, 1e-3);


            double vega = HestonNumericalGreeks.VegaPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedVega = 0.06372181708;
            Assert.AreEqual(expectedVega, vega, 1e-3);
            // Assert.AreEqual(analyticalVega, vega, 1e-3);


            double thetaGreek = HestonNumericalGreeks.ThetaPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedTheta = 0.03659261685613;
            Assert.AreEqual(expectedTheta, thetaGreek, 1e-3);


        }


    }
}
