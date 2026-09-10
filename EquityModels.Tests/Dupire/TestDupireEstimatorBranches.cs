/* Copyright (C) 2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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
using DVPLDOM;
using DVPLI;
using Fairmat.MarketData;
using NUnit.Framework;

namespace Dupire
{
    /// <summary>
    /// Covers <see cref="DupireEstimator"/> branches not touched by the happy-path
    /// end-to-end Monte-Carlo calibration test in DupireCalibration.cs: settings/requirements
    /// metadata, the null-settings and invalid-enum failure paths, and the QuantLib branch.
    /// </summary>
    [TestFixture]
    public class TestDupireEstimatorBranches
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static List<object> LoadMarketData(out CallPriceMarketData hData)
        {
            InterestRateMarketData iData = (InterestRateMarketData)ObjectSerialization.ReadFromFile("../../../TestData/IRMD-sample.bin");
            hData = (CallPriceMarketData)ObjectSerialization.ReadFromFile("../../../TestData/CallData-sample.bin");

            var marketData = new List<object>();
            marketData.Add(iData.DiscountingCurve);
            marketData.Add(hData);
            return marketData;
        }

        [Test]
        public void TestDefaultSettingsAndRequirements()
        {
            var estimator = new DupireEstimator();

            Assert.IsInstanceOf<DupireCalibrationSettings>(estimator.DefaultSettings);

            var requirements = estimator.GetRequirements(null, null);
            Assert.AreEqual(2, requirements.Length);
            Assert.AreEqual(typeof(DiscountingCurveMarketData), requirements[0].MarketDataType);
            Assert.AreEqual(typeof(CallPriceMarketData), requirements[1].MarketDataType);

            Assert.AreEqual(typeof(DupireProcess), estimator.ProvidesTo);
        }

        [Test]
        public void TestEstimateNullSettingsThrows()
        {
            var marketData = LoadMarketData(out _);
            var estimator = new DupireEstimator();

            Assert.Throws<NullReferenceException>(() => estimator.Estimate(marketData, settings: null));
        }

        [Test]
        public void TestEstimateInvalidEnumThrowsNotImplemented()
        {
            var marketData = LoadMarketData(out _);
            var estimator = new DupireEstimator();
            var settings = new DupireCalibrationSettings
            {
                LocalVolatilityCalculation = (LocalVolatilityCalculation)99
            };

            Assert.Throws<NotImplementedException>(() => estimator.Estimate(marketData, settings));
        }

        [Test, Category("BigTest")]
        public void TestEstimateMethod1HappyPath()
        {
            var marketData = LoadMarketData(out CallPriceMarketData hData);
            var estimator = new DupireEstimator();
            var settings = new DupireCalibrationSettings
            {
                LocalVolatilityCalculation = LocalVolatilityCalculation.Method1
            };

            EstimationResult result = estimator.Estimate(marketData, settings);

            Assert.IsNotNull(result);
            Assert.AreEqual(4, result.Objects.Length);
            Assert.AreEqual("S0", result.Names[0]);
            Assert.AreEqual(hData.S0, result.Values[0], 1e-8);
            Assert.IsInstanceOf<PFunction>(result.Objects[0]);
            Assert.IsInstanceOf<PFunction>(result.Objects[1]);
            Assert.IsInstanceOf<PFunction2D.PFunction2D>(result.Objects[2]);
        }

        [Test, Category("BigTest")]
        public void TestEstimateQuantLibBranch()
        {
            // DupireEstimator.Estimate has a comment noting the QuantLib path "does not work
            // correctly"; with the shared sample fixture it fails inside FitImplVolModel
            // because too few (Maturity,Strike) points pass the 0.01 volatility filter to fit
            // the quadratic implied-vol model. This test locks in that documented behavior
            // (and exercises the QuantLib branch/FitImplVolModel code path) rather than
            // asserting a successful calibration.
            var marketData = LoadMarketData(out _);
            var estimator = new DupireEstimator();
            var settings = new DupireCalibrationSettings
            {
                LocalVolatilityCalculation = LocalVolatilityCalculation.QuantLib
            };

            Assert.Throws<Exception>(() => estimator.Estimate(marketData, settings));
        }
    }
}
