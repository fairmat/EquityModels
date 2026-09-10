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
using System.Collections.Generic;
using DVPLDOM;
using DVPLI;
using Fairmat.MarketData;
using HestonEstimator;
using NUnit.Framework;

namespace Heston
{
    /// <summary>
    /// Covers <see cref="CallEstimator"/> (HestonEstimator.cs), which previously had no direct
    /// test coverage: metadata (DefaultSettings/GetRequirements/ProvidesTo), the fast
    /// deterministic dummy-calibration path, and (as a heavier BigTest) the real
    /// QADE + SteepestDescent optimization path.
    /// </summary>
    [TestFixture]
    public class TestHestonEstimator
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
            marketData.Add(new DVPLI.MarketDataTypes.Scalar { Value = hData.S0 });
            return marketData;
        }

        [Test]
        public void TestDefaultSettingsAndRequirements()
        {
            var estimator = new CallEstimator();

            Assert.IsInstanceOf<HestonCalibrationSettings>(estimator.DefaultSettings);

            var requirements = estimator.GetRequirements(null, null);
            Assert.AreEqual(3, requirements.Length);
            Assert.AreEqual(typeof(DiscountingCurveMarketData), requirements[0].MarketDataType);
            Assert.AreEqual(typeof(CallPriceMarketData), requirements[1].MarketDataType);
            Assert.AreEqual(typeof(DVPLI.MarketDataTypes.Scalar), requirements[2].MarketDataType);

            Assert.AreEqual(typeof(HestonExtended.HestonExtendedProcess), estimator.ProvidesTo);
            Assert.IsFalse(string.IsNullOrEmpty(estimator.Description));
            Assert.IsFalse(string.IsNullOrEmpty(estimator.ToolTipText));
        }

        [Test]
        public void TestDummyCalibrationReturnsFixedSolution()
        {
            var marketData = LoadMarketData(out CallPriceMarketData hData);
            var estimator = new CallEstimator();
            var settings = new HestonCalibrationSettings { DummyCalibration = true };

            EstimationResult result = estimator.Estimate(marketData, settings);

            Assert.IsNotNull(result);
            CollectionAssert.AreEqual(new[] { "S0", "kappa", "theta", "sigma", "rho", "V0" }, result.Names);

            double[] expected = new double[] { hData.S0, 0.5, 0.5, 0.8, -0.5, 0.05 };
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], result.Values[i], 1e-8);
        }

        [Test, Category("BigTest")]
        public void TestEstimateEndToEndWithSeed()
        {
            var marketData = LoadMarketData(out _);
            var estimator = new CallEstimator();
            var settings = new HestonCalibrationSettings();

            var properties = new Dictionary<string, object>
            {
                { "RandomSeed", 42 },
                { "NP", 8 },
                { "MaxIter", 5 },
            };

            EstimationResult result = estimator.Estimate(marketData, settings, controller: null, properties: properties);

            Assert.IsNotNull(result);
            Assert.GreaterOrEqual(result.Fit, 0);
            Assert.AreEqual(2, result.Objects.Length);
        }
    }
}
