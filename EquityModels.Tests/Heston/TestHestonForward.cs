﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Heston;
using HestonEstimator;

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

        [Test]
        public void TestApproximatedApproach()
        {
            double k = 0.9;
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

            var expected = 26.1556170228032;


            Assert.AreEqual(expected, result, 1e-10);

        }

        [Test]
        public void TestForwardAhlipRutkowski()
        {
            double k = 0.9;
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
    }
}
