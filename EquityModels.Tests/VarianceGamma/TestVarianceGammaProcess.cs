using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Runtime.Serialization.Formatters.Binary;
using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace VarianceGamma
{
    /// <summary>
    /// Tests the pure/deterministic members of <see cref="VarianceGamma"/> directly
    /// (not through a full Document/ProjectROV Monte Carlo simulation).
    /// </summary>
    [TestFixture]
    public class TestVarianceGammaProcess
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static double GetDrift(VarianceGamma process)
        {
            FieldInfo driftField = typeof(VarianceGamma).GetField("drift", BindingFlags.NonPublic | BindingFlags.Instance);
            return (double)driftField.GetValue(process);
        }

        [Test]
        public void TestConstructorWithDefaultValues()
        {
            VarianceGamma process = new VarianceGamma();

            List<IExportable> exported = process.ExportObjects(false);

            Assert.AreEqual(100.0, (exported[0] as ModelParameter).fV());
            Assert.AreEqual(0.1, (exported[1] as ModelParameter).fV());
            Assert.AreEqual(0.1, (exported[2] as ModelParameter).fV());
            Assert.AreEqual(0.1, (exported[3] as ModelParameter).fV());
            Assert.AreEqual(0.02, (exported[4] as ModelParameter).fV());
            Assert.AreEqual(0.01, (exported[5] as ModelParameter).fV());
        }

        [Test]
        public void TestConstructorWithExplicitValues()
        {
            VarianceGamma process = new VarianceGamma(120.0, -0.2, 0.25, 0.6, 0.05, 0.03);

            List<IExportable> exported = process.ExportObjects(false);

            Assert.AreEqual(6, exported.Count);
            Assert.AreEqual(120.0, (exported[0] as ModelParameter).fV());
            Assert.AreEqual(-0.2, (exported[1] as ModelParameter).fV());
            Assert.AreEqual(0.25, (exported[2] as ModelParameter).fV());
            Assert.AreEqual(0.6, (exported[3] as ModelParameter).fV());
            Assert.AreEqual(0.05, (exported[4] as ModelParameter).fV());
            Assert.AreEqual(0.03, (exported[5] as ModelParameter).fV());
        }

        [Test]
        public void TestProcessInfo()
        {
            VarianceGamma process = new VarianceGamma();

            Assert.AreEqual("Variance Gamma", process.ProcessInfo.ProcessType);
        }

        [Test]
        public void TestSimulationInfo()
        {
            VarianceGamma process = new VarianceGamma();

            SimulationInfo info = process.SimulationInfo;

            Assert.AreEqual(0, info.LatentSize);
            Assert.AreEqual(1, info.NoiseSize);
            Assert.AreEqual(1, info.StateSize);
            CollectionAssert.AreEqual(new string[] { "Index value" }, info.StateDescription);
        }

        [Test]
        public void TestSimpleFlags()
        {
            VarianceGamma process = new VarianceGamma();

            Assert.IsFalse(process.ImplementsFullSimulation);
            Assert.IsTrue(process.ImplementsMarkovBasedSimulation);
        }

        [Test]
        public void TestGetDeltaFactors()
        {
            VarianceGamma process = new VarianceGamma();

            Assert.IsNull(process.GetDeltaFactors());
        }

        [Test]
        public void TestGetVegaFactors()
        {
            VarianceGamma process = new VarianceGamma(120.0, -0.2, 0.25, 0.6, 0.05, 0.03);

            IModelParameter[] vegaFactors = process.GetVegaFactors();
            List<IExportable> exported = process.ExportObjects(false);

            Assert.AreEqual(1, vegaFactors.Length);
            Assert.AreSame(exported[2], vegaFactors[0]);
        }

        [Test]
        public void TestParseReturnsNoErrorsForValidParameters()
        {
            VarianceGamma process = new VarianceGamma(120.0, -0.2, 0.25, 0.6, 0.05, 0.03);

            bool errors = process.Parse(null);

            Assert.IsFalse(errors);
        }

        [Test]
        [TestCase(100.0, 0.1, 0.1, 0.1, 0.02, 0.01)]
        [TestCase(120.0, -0.2, 0.25, 0.6, 0.05, 0.03)]
        [TestCase(80.0, 0.05, 0.15, 0.4, 0.03, 0.0)]
        public void TestSetupComputesDrift(double s0, double theta, double sigma, double nu, double r, double q)
        {
            VarianceGamma process = new VarianceGamma(s0, theta, sigma, nu, r, q);

            process.Setup(new double[] { 0.0, 1.0 });

            double omega = Math.Log(1.0 - nu * theta - 0.5 * sigma * sigma * nu) / nu;
            double expectedDrift = r - q + omega;

            Assert.AreEqual(expectedDrift, GetDrift(process), 1e-12);
        }

        [Test]
        public void TestSimulateSetsFirstRowToS0()
        {
            VarianceGamma process = new VarianceGamma(120.0, -0.2, 0.25, 0.6, 0.05, 0.03);
            double[] dates = new double[] { 0.0 };
            Matrix noise = new Matrix(0, 1);
            Matrix outDynamic = new Matrix(1, 1);

            process.Simulate(dates, noise, outDynamic);

            Assert.AreEqual(120.0, outDynamic[0, 0]);
        }

        [Test]
        public void TestIsSerializable()
        {
            VarianceGamma process = new VarianceGamma(120.0, -0.2, 0.25, 0.6, 0.05, 0.03);

            BinaryFormatter formatter = new BinaryFormatter();
            VarianceGamma deserialized;
            using (MemoryStream stream = new MemoryStream())
            {
                formatter.Serialize(stream, process);
                stream.Position = 0;
                deserialized = (VarianceGamma)formatter.Deserialize(stream);
            }

            List<IExportable> original = process.ExportObjects(false);
            List<IExportable> roundTripped = deserialized.ExportObjects(false);

            Assert.AreEqual(original.Count, roundTripped.Count);
            for (int i = 0; i < original.Count; i++)
            {
                Assert.AreEqual((original[i] as ModelParameter).fV(), (roundTripped[i] as ModelParameter).fV());
            }
        }
    }
}
