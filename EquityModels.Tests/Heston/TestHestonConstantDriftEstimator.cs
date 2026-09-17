using System.Collections.Generic;
using DVPLI;
using DVPLI.MarketDataTypes;
using Fairmat.MarketData;
using Heston;
using NUnit.Framework;

namespace HestonEstimator
{
    /// <summary>
    /// Tests <see cref="HestonConstantDriftEstimator"/>.
    /// </summary>
    [TestFixture]
    public class TestHestonConstantDriftEstimator
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestDescription()
        {
            Assert.AreEqual("Heston constant drift calibration", new HestonConstantDriftEstimator().Description);
        }

        [Test]
        public void TestProvidesTo()
        {
            Assert.AreEqual(typeof(Heston.HestonProcess), new HestonConstantDriftEstimator().ProvidesTo);
        }

        [Test]
        public void TestDefaultSettings()
        {
            Assert.IsInstanceOf<HestonEstimationSettings>(new HestonConstantDriftEstimator().DefaultSettings);
        }

        /// <summary>
        /// Exercises Setup() and BuildEstimate() through the base class' Estimate() pipeline,
        /// using DummyCalibration to skip the (slow) stochastic/deterministic optimizer and
        /// go straight to a known solution vector.
        /// </summary>
        [Test]
        public void TestBuildEstimateWithDummyCalibration()
        {
            var interestData = (InterestRateMarketData)ObjectSerialization.ReadFromFile("../../../TestData/IRMD-sample.bin");
            var callData = (CallPriceMarketData)ObjectSerialization.ReadFromFile("../../../TestData/CallData-sample.bin");
            var spot = new Scalar(callData.S0);

            var estimator = new HestonConstantDriftEstimator();
            var settings = new HestonCalibrationSettings { DummyCalibration = true };
            var marketData = new List<object> { interestData.DiscountingCurve, callData, spot };

            EstimationResult result = estimator.Estimate(marketData, settings);

            Assert.IsNotNull(result);
            CollectionAssert.AreEqual(
                new[] { "S0", "kappa", "theta", "sigma", "rho", "V0", "r", "q" },
                result.Names);

            bool found;
            Assert.AreEqual(callData.S0, PopulateHelper.GetValue("S0", result.Names, result.Values, out found));

            // Dummy solution injected by CallEstimator.Estimate's DummyCalibration path.
            Assert.AreEqual(0.5, PopulateHelper.GetValue("kappa", result.Names, result.Values, out found));
            Assert.AreEqual(0.5, PopulateHelper.GetValue("theta", result.Names, result.Values, out found));
            Assert.AreEqual(0.8, PopulateHelper.GetValue("sigma", result.Names, result.Values, out found));
            Assert.AreEqual(-0.5, PopulateHelper.GetValue("rho", result.Names, result.Values, out found));
            Assert.AreEqual(0.05, PopulateHelper.GetValue("V0", result.Names, result.Values, out found));

            // avgPricingError is never touched on the dummy path, so it keeps its default value.
            Assert.AreEqual(0.0, result.Fit);

            // q comes from DY(equityCalData) (call/put parity implied dividend) since
            // HestonConstantDriftEstimator.impliedDividends defaults to false.
            double q = PopulateHelper.GetValue("q", result.Names, result.Values, out found);
            Assert.IsFalse(double.IsNaN(q));
        }
    }
}
