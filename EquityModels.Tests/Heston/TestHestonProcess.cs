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
using HestonEstimator;
using NUnit.Framework;

namespace Heston
{
    /// <summary>
    /// Covers <see cref="HestonProcess"/>, which previously had no direct test coverage of
    /// its own lifecycle/dispatch members (only indirectly, through subclasses/static pricing
    /// helpers exercised elsewhere). Split into context-free members (no live Document/Project
    /// needed) and context-dependent members (Call/Put/.../Populate), which require
    /// Project.ActiveProject and Option.CurrentSolving to be populated -- these are normally
    /// only set by the DVPL solver mid-valuation, so the tests below set them directly on the
    /// same public static fields the solver uses, mirroring a real valuation context without
    /// running a full Monte Carlo solve.
    /// </summary>
    [TestFixture]
    public class TestHestonProcess
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        #region Context-free members

        [Test]
        public void TestDefaultInstanceSetsDocumentedDefaults()
        {
            var process = new HestonProcess();
            process.DefaultInstance();

            Assert.AreEqual(0.05, process.r.fV(), 1e-12);
            Assert.AreEqual(0.03, process.q.fV(), 1e-12);
            Assert.AreEqual(2.5, process.k.fV(), 1e-12);
            Assert.AreEqual(0.4, process.theta.fV(), 1e-12);
            Assert.AreEqual(0.2, process.sigma.fV(), 1e-12);
            Assert.AreEqual(100, process.S0.fV(), 1e-12);
            Assert.AreEqual(0.3, process.V0.fV(), 1e-12);
            Assert.AreEqual(0.0, process.rho.fV(), 1e-12);
        }

        [Test]
        public void TestGetDeltaAndVegaFactors()
        {
            var process = new HestonProcess();
            process.DefaultInstance();

            CollectionAssert.AreEqual(new[] { process.S0 }, process.GetDeltaFactors());
            CollectionAssert.AreEqual(new[] { process.V0 }, process.GetVegaFactors());
        }

        [Test]
        public void TestLifecycleMetadata()
        {
            var process = new HestonProcess();

            Assert.IsFalse(process.ImplementsFullSimulation);
            Assert.IsTrue(process.ImplementsMarkovBasedSimulation);
            Assert.AreEqual("Heston", process.ProcessInfo.ProcessType);

            Assert.DoesNotThrow(() => process.Setup(new double[] { 0.0, 1.0 }));

            var simInfo = process.SimulationInfo;
            Assert.AreEqual(1, simInfo.LatentSize);
            Assert.AreEqual(2, simInfo.NoiseSize);
            Assert.AreEqual(2, simInfo.StateSize);
        }

        [Test]
        public void TestExportObjectsReturnsAllEightParameters()
        {
            var process = new HestonProcess();
            process.DefaultInstance();

            var exported = process.ExportObjects(false);

            Assert.AreEqual(8, exported.Count);
            CollectionAssert.Contains(exported, process.S0);
            CollectionAssert.Contains(exported, process.V0);
            CollectionAssert.Contains(exported, process.r);
            CollectionAssert.Contains(exported, process.q);
            CollectionAssert.Contains(exported, process.k);
            CollectionAssert.Contains(exported, process.theta);
            CollectionAssert.Contains(exported, process.sigma);
            CollectionAssert.Contains(exported, process.rho);
        }

        [Test]
        public void TestSwapAndFSSwapThrowNotImplemented()
        {
            var process = new HestonProcess();
            process.DefaultInstance();

            Assert.Throws<NotImplementedException>(() => process.Swap(0, 100, 1.0));
            Assert.Throws<NotImplementedException>(() => process.FSSwap(0, 1.0, 0.1, 1.0));
        }

        [Test]
        public void TestParseDoesNotThrow()
        {
            var process = new HestonProcess();
            process.DefaultInstance();

            Assert.DoesNotThrow(() => process.Parse(null));
        }

        #endregion

        #region Context-dependent members

        private static (Document doc, ProjectROV rov) CreateSolvingContext(double flatRate)
        {
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);

            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.RiskFree;
            rfi.m_deterministicRF = flatRate;

            return (doc, rov);
        }

        [Test]
        public void TestCallAndPutPriceThroughLiveContext()
        {
            var (doc, rov) = CreateSolvingContext(0.1);

            var process = new HestonProcess();
            process.DefaultInstance();
            process.r = (ModelParameter)0.1;
            process.q = (ModelParameter)0.07;
            process.S0 = (ModelParameter)100.0;

            var op = new OptionTree(rov);
            rov.Map.Root = op;
            op.CalculateDiscountingModelCache();

            var previousActiveProject = Project.ActiveProject;
            Project.ActiveProject = rov;
            Option.CurrentSolving = op;
            try
            {
                var call = process.Call(0, strike: 90, timeToMaturity: 2.0, additionalInformation: null);
                var put = process.Put(0, strike: 90, timeToMaturity: 2.0, additionalInformation: null);

                double expectedCall = HestonCall.HestonCallPrice(
                    kappa: process.k.fV(), theta: process.theta.fV(), sigma: process.sigma.fV(),
                    rho: process.rho.fV(), v0: process.V0.fV(), s0: process.S0.fV(),
                    T: 2.0, K: 90, r: process.r.fV(), q: process.q.fV());

                Assert.AreEqual(expectedCall, call.MarkToMarket, 1e-6);
                Assert.IsFalse(double.IsNaN(put.MarkToMarket));
                Assert.Greater(call.MarkToMarket, 0);
                Assert.Greater(put.MarkToMarket, 0);
            }
            finally
            {
                Option.CurrentSolving = null;
                Project.ActiveProject = previousActiveProject;
            }
        }

        [Test]
        public void TestCallValuationModeBranches()
        {
            var (doc, rov) = CreateSolvingContext(0.1);

            var process = new HestonProcess();
            process.DefaultInstance();
            process.r = (ModelParameter)0.1;
            process.q = (ModelParameter)0.07;
            process.S0 = (ModelParameter)100.0;

            var op = new OptionTree(rov);
            rov.Map.Root = op;
            op.CalculateDiscountingModelCache();

            var previousActiveProject = Project.ActiveProject;
            Project.ActiveProject = rov;
            Option.CurrentSolving = op;
            try
            {
                GreeksDerivatives Price(AnalyticalPricingFunctionsValuationMode mode)
                {
                    var info = new Dictionary<string, object> { { AnalyticalPricingFunctions.GreekNameKey, mode } };
                    return process.Call(0, strike: 90, timeToMaturity: 2.0, additionalInformation: info);
                }

                var delta = Price(AnalyticalPricingFunctionsValuationMode.Delta);
                var gamma = Price(AnalyticalPricingFunctionsValuationMode.Gamma);
                var rho = Price(AnalyticalPricingFunctionsValuationMode.Rho);
                var theta = Price(AnalyticalPricingFunctionsValuationMode.Theta);
                var vega = Price(AnalyticalPricingFunctionsValuationMode.Vega);
                var all = Price(AnalyticalPricingFunctionsValuationMode.All);

                Assert.IsFalse(double.IsNaN(delta.Deltas[0]));
                Assert.IsFalse(double.IsNaN(gamma.Gammas[0]));
                Assert.IsFalse(double.IsNaN(rho.Rho));
                Assert.IsFalse(double.IsNaN(theta.Theta));
                Assert.IsFalse(double.IsNaN(vega.Vegas[0]));

                Assert.AreEqual(delta.Deltas[0], all.Deltas[0], 1e-6);
                Assert.AreEqual(gamma.Gammas[0], all.Gammas[0], 1e-6);
                Assert.AreEqual(rho.Rho, all.Rho, 1e-6);
                Assert.AreEqual(theta.Theta, all.Theta, 1e-6);
                Assert.AreEqual(vega.Vegas[0], all.Vegas[0], 1e-6);
            }
            finally
            {
                Option.CurrentSolving = null;
                Project.ActiveProject = previousActiveProject;
            }
        }

        [Test]
        public void TestPopulateUpdatesProcessFieldsAndCorrelation()
        {
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);

            var process = new HestonProcess();
            var stocProcess = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(stocProcess);

            var estimate = new EstimationResult(
                new[] { "S0", "kappa", "theta", "sigma", "rho", "V0" },
                new[] { 120.0, 3.0, 0.5, 0.25, -0.6, 0.35 });

            process.Populate(stocProcess, estimate);

            Assert.AreEqual(120.0, process.S0.fV(), 1e-8);
            Assert.AreEqual(3.0, process.k.fV(), 1e-8);
            Assert.AreEqual(0.5, process.theta.fV(), 1e-8);
            Assert.AreEqual(0.25, process.sigma.fV(), 1e-8);
            Assert.AreEqual(0.35, process.V0.fV(), 1e-8);
            Assert.AreEqual(-0.6, process.rho.fV(), 1e-8);
        }

        #endregion
    }
}
