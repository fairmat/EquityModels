using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Heston;
using HestonEstimator;
using Fairmat.Math;
using System.Diagnostics;

namespace EquityModels.Tests.Heston
{
    [TestFixture]
    public class TestHestonForward
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        //[Test]
        public void TestApproximatedApproach()
        {
            double k = 1.0;
            double tau = 2.0;

            double rate = 0.01;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;
            double T0 = 0.5;


            var result = HestonForwardApproximated.HestonForwardCallPrice(
                K: k,
                T: tau,
                T0: T0,
                r: rate,
                q: dy,
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                s0: s0,
                v0: v0,
                rho: rho
                );



            var expected = 22.872215559994689d;


        }

        //[Test]
        public void TestApproximatedApproachVersusMC()
        {
            double k = 1.0;
            double tau = 1.50;

            double rate = 0.0;
            double dy = 0.0;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.9;
            double T0 = 1.0;


            var result = HestonForwardApproximated.HestonForwardCallPrice(
                K: k,
                T: tau,
                T0: T0,
                r: rate,
                q: dy,
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                s0: s0,
                v0: v0,
                rho: rho
                );

            var resultJJ = HestonForwardJJ.CalculateFwdCallPrice(
                K: k,
                T: tau,
                T0: T0,
                r: rate,
                q: dy,
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                s0: s0,
                v0: v0,
                rho: rho
                );

            resultJJ = s0 * Math.Exp(-dy * T0) * resultJJ;

            Assert.AreEqual(result, resultJJ, 1e-2);
        }

        //[Test]
        public void TestForwardAhlipRutkowski()
        {
            double k = 1.0;
            double tau = 2.0;

            double rate = 0.01;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.2;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;
            double T0 = 0.5;


            var result = HestonForwardAhlipRutkowski.HestonForwardCallPrice(
                K: k,
                T: tau,
                T0: T0,
                r: rate,
                q: dy,
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                s0: s0,
                v0: v0,
                rho: rho
                );

            var expected = 0.0;

            Assert.AreEqual(expected, result, 1e-10);

        }

        //[Test]
        public void TestCallPriceAhlipRutkowski()
        {
            double strike = 0.1;
            double tau = 0.5;

            double rate = 0.01;
            double dy = 0.00;
            double kappa = 2.5;
            double theta = 0.2;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.15;
            double rho = -0.8;


            var resultFromFwd = HestonForwardAhlipRutkowski.CallPrice(
                K: strike,
                T: tau,
                T0: 0.0,
                r: rate,
                q: dy,
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                s0: s0,
                v0: v0,
                rho: rho
                );

            // in the context of fwd starting k is a percentage of s0
            var resultFromHeston = HestonCall.HestonCallPrice(
                K: strike*s0,
                T: tau,
                r: rate,
                q: dy,
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                s0: s0,
                v0: v0,
                rho: rho);

            Assert.AreEqual(resultFromHeston, resultFromFwd, 1e-10);

        }
        
        //[Test]
        public void TestCallPriceJJ()
        {
            double strike = 1.0;
            double tau = 1;
            double T0 = 0.0;
            double rate = 0.0;
            double dy = 0.00;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.5; 


            var resultFromFwd = HestonForwardJJ.CalculateFwdCallPrice(
                K: strike,
                T: tau,
                T0: T0,
                r: rate,
                q: dy,
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                s0: s0,
                v0: v0,
                rho: rho
                );
            
            // in the context of fwd starting k is a percentage of s0
            var resultFromHeston = HestonCall.HestonCallPrice(
                K: strike * s0,
                T: tau,
                r: rate,
                q: dy,
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                s0: s0,
                v0: v0,
                rho: rho);

            Assert.AreEqual(resultFromHeston, s0 * Math.Exp(-dy * T0) * resultFromFwd, 1e-10);

        }

        //[Test]
        public void TestValidationFwdJJ()
        {
            double strike = 1.0;
            double tau = 1;
            double T0 = 0.1;
            double rate = 0.0;
            double dy = 0.00;
            double kappa = 2.1;
            double theta = 0.03;
            double sigma = 0.1;
            double s0 = 1;
            double v0 = 0.05;
            double rho = -0.2;

            var mp = new HestonEstimator.ModelParameters
            {
                kappa = kappa,
                theta = theta,
                sigma = sigma,
                rho = rho,
                v0 = v0,
                s0 = s0,
                r = rate,
                q = dy
            };

            var u = new Complex(0.5, 1); 
            var b = HestonForwardJJ.GetLittleB(mp, T0);
            var d = HestonForwardJJ.D(mp, u);
            var littleGamma = HestonForwardJJ.GetLittleGamma(mp, u);
            var cfHeston = HestonForwardJJ.CfHestonFwd(mp, u, T0, tau);
            var intPhi = HestonForwardJJ.GetIntegrandFunction(mp, 2.0, 0.5, mp.kappa, T0, tau);

            Assert.AreEqual(0.0002254949452, b, 1e-10);
            Assert.AreEqual((2.1129597440980485), d.Re, 1e-10);
            Assert.AreEqual(0.0199719848510477, d.Im, 1e-10);
            //Assert.AreEqual(-0.037404360915137914, intPhi, 1e-10); 



            var resultFromFwd = HestonForwardJJ.CalculateFwdCallPrice(
                K: strike,
                T: tau,
                T0: T0,
                r: rate,
                q: dy,
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                s0: s0,
                v0: v0,
                rho: rho
                );

            // in the context of fwd starting k is a percentage of s0
            
            Assert.AreEqual(0.07597741087461451, resultFromFwd, 1e-5);

        }

        [TestFixture]
        public class HestonForwardAhlipRutkowskiTests
        {
            [Test]
            public void TestQ2()
            {
                
                // Arrange
                var mp = new HestonForwardAhlipRutkowski.ModelParameters
                {
                    rho = -0.8,
                    sigma = 0.2, 
                    kappa = 0.5
                };
                double u = 0.5;

                // Act
                Complex result = HestonForwardAhlipRutkowski.q2(u, mp);

                // Assert
                // Here you should check if the result is what you expect.
                // This depends on what you expect the output to be for the given input.
                // For example, if you expect the result to be (0.25 + 0i), you could write:
                Assert.AreEqual(new Complex(0, 2), result);
            }

            [Test]
            public void TestQ1()
            {

                // Arrange
                var mp = new HestonForwardAhlipRutkowski.ModelParameters
                {
                    rho = -0.8,
                    sigma = 0.2,
                    kappa = 0.5
                };
                double u = 0.5;

                // Act
                Complex result = HestonForwardAhlipRutkowski.q1(u, mp);

                // Assert
                // Here you should check if the result is what you expect.
                // This depends on what you expect the output to be for the given input.
                // For example, if you expect the result to be (0.25 + 0i), you could write:
                Assert.AreEqual(new Complex(0.045, 1.25), result);
            }
        }

    }
}
