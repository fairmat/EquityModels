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
using System.Numerics;
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
        public void TestForwardCall()
        {
            /*
            double k = 0.9;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 1.0;
            double v0 = 0.3;
            double rho = 0.8;
            double T0 = 0.00001;
            */
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

            // Calculates the theoretical value of the call.
            DVPLI.Vector param = new DVPLI.Vector(5);
            param[0] = kappa;
            param[1] = theta;
            param[2] = sigma;
            param[3] = rho;
            param[4] = v0;
            double fairmatPrice = HestonCall.HestonCallPrice(param, s0, tau, k, rate, dy);

            double lambda = kappa;
            double vhat = theta;
            double eta = sigma;

            double tol = 1e-3;
            double benchmarkPrice = 0.339537359104676;

            Console.WriteLine("Theoretical Benchmark  Price = " + benchmarkPrice.ToString());
            Console.WriteLine("Theoretical Fairmat    Price = " + fairmatPrice);

            var fairmatPriceApprox = HestonForwardApproximated.HestonForwardCallPrice(x: param, s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var fairmatPriceApprox2 = HestonForwardApproximated.HestonForwardPercentageCallPrice(x: param, s0, T: tau, T0: T0, K: k, r: rate, q: dy);

            //Assert.Less(Math.Abs(fairmatPrice - benchmarkPrice), tol);
            Assert.Less(Math.Abs(fairmatPrice - fairmatPriceApprox), tol);
            Assert.Less(Math.Abs(fairmatPrice - fairmatPriceApprox2), tol);
        }

        [Test]
        public void TestForwardPut()
        {
            /*
            double k = 0.9;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 1.0;
            double v0 = 0.3;
            double rho = 0.8;
            double T0 = 0.00001;
            */
            double k = 1;
            double tau = 2.0;

            double rate = 0.5;
            double dy = 0.0;
            double kappa = 1;
            double theta = 0.4;
            double sigma = 0.1;
            double s0 = 1.0;
            double v0 = 0.4;
            double rho = 0.2;
            double T0 = 0.001;

            // Calculates the theoretical value of the call.
            DVPLI.Vector param = new DVPLI.Vector(5);
            param[0] = kappa;
            param[1] = theta;
            param[2] = sigma;
            param[3] = rho;
            param[4] = v0;
            double fairmatPrice = HestonCall.HestonPutPrice(param, s0, tau, k, rate, dy);
            //double fairmatPriceF = HestonForwardAhlipRutkowski.HestonForwardCallPrice(kappa, theta, rho, v0, sigma, s0, k, rate, dy, tau, T0);
            // in lucic paper we have 
            // lambda * (vhat - v )
            // in our notation we have
            //kappa(theta-v)

            double lambda = kappa;
            double vhat = theta;
            double eta = sigma;

            double tol = 1e-3;
            double benchmarkPrice = 0.339537359104676;

            Console.WriteLine("Theoretical Benchmark  Price = " + benchmarkPrice.ToString());
            Console.WriteLine("Theoretical Fairmat    Price = " + fairmatPrice);
            //Console.WriteLine("Theoretical Fairmat forward  Price = " + fairmatPriceF);
            var fairmatPriceApprox = HestonForwardApproximated.HestonForwardPutPrice(x: param, s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var fairmatPriceApprox2 = HestonForwardApproximated.HestonForwardPercentagePutPrice(x: param, s0, T: tau, T0: T0, K: k, r: rate, q: dy);

            Assert.Less(Math.Abs(fairmatPriceApprox - fairmatPriceApprox2), tol);

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
            double s0 = 100.0;

            double k = 0.9*s0;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double v0 = 0.3;
            double rho = -0.8;
            var Tp = tau + 0.5;

            double[] tenors = [0.1, 0.2, 0.5, 1, 2];
            double[] rates = [0.1, 0.1, 0.1, 0.1, 0.1];
            var df = new DVPLDOM.PFunction((DVPLI.Vector)tenors, (DVPLI.Vector)rates);
            df.m_Function.iType = DVPLUtils.EInterpolationType.ZERO_ORDER;

            Func<double, double, double> discountingFunction = (t, T) => Math.Exp(-df.Evaluate(T)*T) ;

            // Calculates the greeks.
            Engine.Verbose = 0;

            var analyticalDelta = HestonDelta.DeltaCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);
            var analyticalGamma = HestonGamma.GammaCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);
            var analyticalRho = HestonRho.RhoCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);
            var analyticalVega = HestonVega.VegaCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);


            (double delta, double gamma) = HestonNumericalGreeks.DeltaGammaCall(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction:discountingFunction, Tp: Tp);
            double expectedDelta = 0.61450217824972952;
            Assert.AreEqual(expectedDelta, delta, 1e-3);
            Assert.AreEqual(analyticalDelta, delta, 1e-3);


            double expectedGamma = 0.00310550;
            Assert.AreEqual(expectedGamma, gamma, 1e-3);
            Assert.AreEqual(analyticalGamma, gamma, 1e-3);


            double rhoGreek = HestonNumericalGreeks.RhoCall(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, Tp:Tp);
            double expectedRho = 42.1301853;
            Assert.AreEqual(expectedRho, rhoGreek, 1e-3);
            Assert.AreEqual(analyticalRho, rhoGreek, 1e-3);


            double vega = HestonNumericalGreeks.VegaCall(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, timeToPaymentDate:Tp);
            double expectedVega = 6.061406719;
            Assert.AreEqual(expectedVega, vega, 1e-3);
            Assert.AreEqual(analyticalVega, vega, 1e-3);


            double thetaGreek = HestonNumericalGreeks.ThetaCall(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime:Tp);
            double expectedTheta = -4.701290525161;
            Assert.AreEqual(expectedTheta, thetaGreek, 1e-3);


        }

        /*
        [Test]
        public void TestGreeksPut()
        {
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            var k = 0.9 * s0;
            
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

            double expectedGamma = 0.003264832;
            Assert.AreEqual(expectedGamma, gamma, 1e-3);
            Assert.AreEqual(analyticalGamma, gamma, 1e-3);


            double rhoGreek = HestonNumericalGreeks.RhoPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedRho = -86.098999518;
            Assert.AreEqual(expectedRho, rhoGreek, 1e-3);
            Assert.AreEqual(analyticalRho, rhoGreek, 1e-3);


            double vega = HestonNumericalGreeks.VegaPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedVega = 6.372181708;
            Assert.AreEqual(expectedVega, vega, 1e-3);
            Assert.AreEqual(analyticalVega, vega, 1e-3);


            double thetaGreek = HestonNumericalGreeks.ThetaPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy);
            double expectedTheta =-3.659261685613;
            Assert.AreEqual(expectedTheta, thetaGreek, 1e-3);


        }
        */

        [Test]
        public void TestGreeksPut()
        {
            double s0 = 100.0;

            double k = 0.9 * s0;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double v0 = 0.3;
            double rho = -0.8;
            var Tp = tau + 0.5;

            double[] tenors = [0.1, 0.2, 0.5, 1, 2];
            double[] rates = [0.1, 0.1, 0.1, 0.1, 0.1];
            var df = new DVPLDOM.PFunction((DVPLI.Vector)tenors, (DVPLI.Vector)rates);
            df.m_Function.iType = DVPLUtils.EInterpolationType.ZERO_ORDER;

            Func<double, double, double> discountingFunction = (t, T) => Math.Exp(-df.Evaluate(T) * T);

            // Calculates the greeks.
            Engine.Verbose = 0;

            var analyticalDelta = HestonDelta.DeltaPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);
            var analyticalGamma = HestonGamma.GammaPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);
            var analyticalRho = HestonRho.RhoPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);
            var analyticalVega = HestonVega.VegaPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);


            (double delta, double gamma) = HestonNumericalGreeks.DeltaGammaPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);
            double expectedDelta = -0.21245695569364287;
            Assert.AreEqual(expectedDelta, delta, 1e-3);
            Assert.AreEqual(analyticalDelta, delta, 1e-3);


            double expectedGamma = 0.0032647285578946135;
            Assert.AreEqual(expectedGamma, gamma, 1e-3);
            Assert.AreEqual(analyticalGamma, gamma, 1e-3);


            double rhoGreek = HestonNumericalGreeks.RhoPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);
            double expectedRho = -91.752036270751347;
            Assert.AreEqual(expectedRho, rhoGreek, 1e-3);
            Assert.AreEqual(analyticalRho, rhoGreek, 1e-3);


            double vega = HestonNumericalGreeks.VegaPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);
            double expectedVega = 6.3721817085079149;
            Assert.AreEqual(expectedVega, vega, 1e-3);
            Assert.AreEqual(analyticalVega, vega, 1e-3);


            double thetaGreek = HestonNumericalGreeks.ThetaPut(bumpPercentage: 0.001, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, timeToPaymentDate: Tp);
            double expectedTheta = -3.4807973873052234;
            Assert.AreEqual(expectedTheta, thetaGreek, 1e-3);
        }

        [Test]
        public void TestGreeksFSCall()
        {
            double k = 0.9;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;
            double T0 = 0.01;

            // Calculates the greeks.
            Engine.Verbose = 0;

            (var numericalDelta, var numericalGamma) = HestonNumericalGreeks.DeltaGammaFSCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy); 
            var numericalRho = HestonNumericalGreeks.RhoFSCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalVega = HestonNumericalGreeks.VegaFSCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalTheta = HestonNumericalGreeks.ThetaFSCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);


            var analyticalGreeks = HestonForwardApproximated.HestonForwardCallWithGreeks(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);

            Assert.AreEqual(analyticalGreeks.Deltas[0], numericalDelta, 1e-3);
            Assert.AreEqual(analyticalGreeks.Gammas[0], numericalGamma, 1e-3);
            Assert.AreEqual(analyticalGreeks.Theta, numericalTheta, 1e-3);
            Assert.AreEqual(analyticalGreeks.Vegas[0], numericalVega, 1e-3);
            Assert.AreEqual(analyticalGreeks.Rho, numericalRho, 1e-3);


        }

        [Test]
        public void TestGreeksFSPCall()
        {
            double k = 0.9;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;
            double T0 = 0.01;
            var Tp = tau + 0.1;
            
            double[] tenors = [0.1, 0.2, 0.5, 1, 2];
            double[] rates = [0.01, 0.015, 0.02, 0.05, 0.06];
            var df = new DVPLDOM.PFunction((DVPLI.Vector)tenors, (DVPLI.Vector)rates);
            Func<double, double, double> discountingFunction = (t, T) => Math.Exp(-df.Evaluate(T) * T);

            // Calculates the greeks.
            Engine.Verbose = 0;

            (var numericalDelta, var numericalGamma) = HestonNumericalGreeks.DeltaGammaFSPCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, timeToPaymentDate: Tp);
            var numericalRho = HestonNumericalGreeks.RhoFSPCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, Tp: Tp);
            var numericalVega = HestonNumericalGreeks.VegaFSPCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, Tp: Tp);
            var numericalTheta = HestonNumericalGreeks.ThetaFSPCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, Tp: Tp);


            var analyticalGreeksPercentage = HestonForwardApproximated.HestonForwardPercentageCallWithGreeks(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);

            var expectedDelta = 0;
            Assert.AreEqual(analyticalGreeksPercentage.Deltas[0], numericalDelta, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Deltas[0], expectedDelta, 1e-5);

            var expectedGamma = 0;
            Assert.AreEqual(analyticalGreeksPercentage.Gammas[0], numericalGamma, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Gammas[0], expectedGamma, 1e-5);

            var expectedVega = 0.066917570932606418;
            Assert.AreEqual(analyticalGreeksPercentage.Vegas[0], numericalVega, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Vegas[0], expectedVega, 1e-5);

            var expectedTheta = 0.013138142550550569;
            Assert.AreEqual(analyticalGreeksPercentage.Theta, numericalTheta, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Theta, expectedTheta, 1e-5);

            var expectedRho = 0.6156029386666928;
            Assert.AreEqual(analyticalGreeksPercentage.Rho, numericalRho, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Rho, expectedRho, 1e-5);
        }

        [Test]
        public void TestGreeksFSPut()
        {
            double k = 0.9;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;
            double T0 = 1.0;

            double[] tenors = [0.1, 0.2, 0.5, 1, 2];
            double[] rates = [0.01, 0.015, 0.02, 0.05, 0.06];
            var df = new DVPLDOM.PFunction((DVPLI.Vector)tenors, (DVPLI.Vector)rates);
            Func<double, double, double> discountingFunction = (t, T) => df.Evaluate(T);

            // Calculates the greeks.
            Engine.Verbose = 0;

            (var numericalDelta, var numericalGamma) = HestonNumericalGreeks.DeltaGammaFSPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalRho = HestonNumericalGreeks.RhoFSPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalVega = HestonNumericalGreeks.VegaFSPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalTheta = HestonNumericalGreeks.ThetaFSPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);


            var analyticalGreeks = HestonForwardApproximated.HestonForwardPutWithGreeks(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);


            Assert.AreEqual(analyticalGreeks.Deltas[0], numericalDelta, 1e-3);
            Assert.AreEqual(analyticalGreeks.Gammas[0], numericalGamma, 1e-3);
            Assert.AreEqual(analyticalGreeks.Theta, numericalTheta, 1e-3);
            Assert.AreEqual(analyticalGreeks.Vegas[0], numericalVega, 1e-3);
            Assert.AreEqual(analyticalGreeks.Rho, numericalRho, 1e-3);


        }

        [Test]
        public void TestGreeksFSPPut()
        {
            double k = 0.9;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;
            double T0 = 0.01;
            var Tp = tau + 0.1;

            double[] tenors = [0.1, 0.2, 0.5, 1, 2];
            double[] rates = [0.01, 0.015, 0.02, 0.05, 0.06];
            var df = new DVPLDOM.PFunction((DVPLI.Vector)tenors, (DVPLI.Vector)rates);
            Func<double, double, double> discountingFunction = (t, T) => df.Evaluate(T);

            // Calculates the greeks.
            Engine.Verbose = 0;

            (var numericalDelta, var numericalGamma) = HestonNumericalGreeks.DeltaGammaFSPPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, timeToPaymentDate: Tp);
            var numericalRho = HestonNumericalGreeks.RhoFSPPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, Tp: Tp);
            var numericalVega = HestonNumericalGreeks.VegaFSPPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, Tp: Tp);
            var numericalTheta = HestonNumericalGreeks.ThetaFSPPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, Tp: Tp);


            var analyticalGreeksPercentage = HestonForwardApproximated.HestonForwardPercentagePutWithGreeks(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy, discountingFactorFunction: discountingFunction, paymentTime: Tp);

            var expectedDelta = 0;
            Assert.AreEqual(analyticalGreeksPercentage.Deltas[0], numericalDelta, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Deltas[0], expectedDelta, 1e-5);

            var expectedGamma = 0;
            Assert.AreEqual(analyticalGreeksPercentage.Gammas[0], numericalGamma, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Gammas[0], expectedGamma, 1e-5);

            var expectedVega = 0.0046398412977064164;
            Assert.AreEqual(analyticalGreeksPercentage.Vegas[0], numericalVega, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Vegas[0], expectedVega, 1e-5);

            var expectedTheta = -0.0036849534534085654;
            Assert.AreEqual(analyticalGreeksPercentage.Theta, numericalTheta, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Theta, expectedTheta, 1e-5);

            var expectedRho = -0.0654834715356255;
            Assert.AreEqual(analyticalGreeksPercentage.Rho, numericalRho, 1e-5);
            Assert.AreEqual(analyticalGreeksPercentage.Rho, expectedRho, 1e-5);
        }

        [Test]
        public void TestGreeksFSDigitalPut()
        {
            double k = 0.9;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;
            double T0 = 1.0;

            // Calculates the greeks.
            Engine.Verbose = 0;

            (var numericalDelta, var numericalGamma) = HestonNumericalGreeks.DeltaGammaFSDPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalRho = HestonNumericalGreeks.RhoFSDPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalVega = HestonNumericalGreeks.VegaFSDPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalTheta = HestonNumericalGreeks.ThetaFSDPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var analyticalGreeks = HestonForwardApproximated.HestonForwardDigitalPutWithGreeks(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);


            Assert.AreEqual(analyticalGreeks.Deltas[0], numericalDelta, 1e-3);
            Assert.AreEqual(analyticalGreeks.Gammas[0], numericalGamma, 1e-3);
            Assert.AreEqual(analyticalGreeks.Theta, numericalTheta, 1e-3);
            Assert.AreEqual(analyticalGreeks.Vegas[0], numericalVega, 1e-3);
            Assert.AreEqual(analyticalGreeks.Rho, numericalRho, 1e-3);


        }

        [Test]
        public void TestGreeksFSDigitalCall()
        {
            double k = 0.9;
            double tau = 2.0;

            double rate = 0.1;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;
            double T0 = 1.0;

            // Calculates the greeks.
            Engine.Verbose = 0;

            (var numericalDelta, var numericalGamma) = HestonNumericalGreeks.DeltaGammaFSDCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalRho = HestonNumericalGreeks.RhoFSDCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalVega = HestonNumericalGreeks.VegaFSDCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);
            var numericalTheta = HestonNumericalGreeks.ThetaFSDCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);


            var analyticalGreeks = HestonForwardApproximated.HestonForwardDigitalCallWithGreeks(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: tau, T0: T0, K: k, r: rate, q: dy);


            Assert.AreEqual(analyticalGreeks.Deltas[0], numericalDelta, 1e-3);
            Assert.AreEqual(analyticalGreeks.Gammas[0], numericalGamma, 1e-3);
            Assert.AreEqual(analyticalGreeks.Theta, numericalTheta, 1e-3);
            Assert.AreEqual(analyticalGreeks.Vegas[0], numericalVega, 1e-3);
            Assert.AreEqual(analyticalGreeks.Rho, numericalRho, 1e-3);


        }


    }
}
