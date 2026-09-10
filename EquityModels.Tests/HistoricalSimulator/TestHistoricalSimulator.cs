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
using System.IO;
using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace HistoricalSimulator
{
    /// <summary>
    /// Covers <see cref="HistoricalSimulator"/>, which previously had zero test coverage:
    /// the file-based Setup/Simulate happy path, the missing-file and malformed-file
    /// defensive branches, the unknown-date fallback, and the lifecycle/metadata members.
    /// </summary>
    [TestFixture]
    public class TestHistoricalSimulator
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static string WriteCsvFixture(DateTime firstDate)
        {
            string path = Path.Combine(Path.GetTempPath(), "HistoricalSimulatorTest_" + Guid.NewGuid().ToString("N") + ".csv");
            string[] lines = new string[]
            {
                firstDate.ToString("yyyy-MM-dd") + ";1.1;2.2;3.3",
                firstDate.AddDays(1).ToString("yyyy-MM-dd") + ";1.2;2.3;3.4",
                firstDate.AddDays(2).ToString("yyyy-MM-dd") + ";1.3;2.4;3.5",
            };
            File.WriteAllLines(path, lines);
            return path;
        }

        private static Document CreateActiveDocumentWithSimulationStart(DateTime simulationStartDate)
        {
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            doc.SimulationStartDate = simulationStartDate;
            return doc;
        }

        [Test]
        public void TestSetupAndSimulate_TranslateForward()
        {
            DateTime firstDate = new DateTime(2020, 1, 1);
            string path = WriteCsvFixture(firstDate);
            try
            {
                CreateActiveDocumentWithSimulationStart(firstDate);

                var sim = new global::HistoricalSimulator.HistoricalSimulator
                {
                    FilePath = path,
                    StartDate = firstDate,
                    OperatingMode = OperatingMode.TranslateHistoricalRealizationsForward
                };

                sim.Setup(new double[] { 0.0 });

                Assert.AreEqual(3, sim.SimulationInfo.StateSize);

                var outDynamic = new Matrix(1, 2);
                var noise = new Matrix(1, 2);
                sim.Simulate(new double[] { 0.0 }, noise, outDynamic);

                // Simulate only copies Length-1 columns of the matched row (drops the last one).
                Assert.AreEqual(1.1, outDynamic[0, 0], 1e-8);
                Assert.AreEqual(2.2, outDynamic[0, 1], 1e-8);
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void TestSimulateFallsBackToIndexZeroForUnknownDate()
        {
            DateTime firstDate = new DateTime(2020, 1, 1);
            string path = WriteCsvFixture(firstDate);
            try
            {
                CreateActiveDocumentWithSimulationStart(firstDate);

                var sim = new global::HistoricalSimulator.HistoricalSimulator
                {
                    FilePath = path,
                    StartDate = firstDate,
                    OperatingMode = OperatingMode.TranslateHistoricalRealizationsForward
                };

                sim.Setup(new double[] { 0.0 });

                var outDynamic = new Matrix(1, 2);
                var noise = new Matrix(1, 2);
                // 99.0 was never registered by Setup -> falls back to dateIndex = 0.
                sim.Simulate(new double[] { 99.0 }, noise, outDynamic);

                Assert.AreEqual(1.1, outDynamic[0, 0], 1e-8);
                Assert.AreEqual(2.2, outDynamic[0, 1], 1e-8);
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void TestSetupMissingFileLeavesStateSizeZero()
        {
            CreateActiveDocumentWithSimulationStart(DateTime.Today);

            var sim = new global::HistoricalSimulator.HistoricalSimulator
            {
                FilePath = Path.Combine(Path.GetTempPath(), "does_not_exist_" + Guid.NewGuid().ToString("N") + ".csv"),
                OperatingMode = OperatingMode.TranslateHistoricalRealizationsForward
            };

            Assert.DoesNotThrow(() => sim.Setup(new double[] { 0.0 }));
            Assert.AreEqual(0, sim.SimulationInfo.StateSize);
        }

        [Test]
        public void TestSetupMalformedFileSwallowsException()
        {
            string path = Path.Combine(Path.GetTempPath(), "HistoricalSimulatorTest_malformed_" + Guid.NewGuid().ToString("N") + ".csv");
            File.WriteAllLines(path, new[] { "not-a-date;1;2" });
            try
            {
                CreateActiveDocumentWithSimulationStart(DateTime.Today);

                var sim = new global::HistoricalSimulator.HistoricalSimulator
                {
                    FilePath = path,
                    OperatingMode = OperatingMode.TranslateHistoricalRealizationsForward
                };

                Assert.DoesNotThrow(() => sim.Setup(new double[] { 0.0 }));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void TestSetupAndSimulate_Bootstrap()
        {
            DateTime firstDate = new DateTime(2020, 1, 1);
            string path = WriteCsvFixture(firstDate);
            try
            {
                CreateActiveDocumentWithSimulationStart(firstDate);

                var sim = new global::HistoricalSimulator.HistoricalSimulator
                {
                    FilePath = path,
                    StartDate = firstDate,
                    OperatingMode = OperatingMode.Bootstrap
                };

                sim.Setup(new double[] { 0.0, 1.0 });

                var outDynamic = new Matrix(2, 3);
                var noise = new Matrix(2, 3);
                Assert.DoesNotThrow(() => sim.Simulate(new double[] { 0.0, 1.0 }, noise, outDynamic));

                for (int c = 0; c < 3; c++)
                    Assert.AreEqual(1.0, outDynamic[0, c], 1e-8);
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void TestParseAlwaysFalse()
        {
            var sim = new global::HistoricalSimulator.HistoricalSimulator();
            Assert.IsFalse(sim.Parse(null));
        }

        [Test]
        public void TestExportObjectsEmpty()
        {
            var sim = new global::HistoricalSimulator.HistoricalSimulator();
            var exported = sim.ExportObjects(true);
            Assert.IsNotNull(exported);
            Assert.AreEqual(0, exported.Count);
        }

        [Test]
        public void TestImplementsFlags()
        {
            var sim = new global::HistoricalSimulator.HistoricalSimulator();
            Assert.IsTrue(sim.ImplementsFullSimulation);
            Assert.IsFalse(sim.ImplementsMarkovBasedSimulation);
        }

        [Test]
        public void TestProcessInfo()
        {
            var sim = new global::HistoricalSimulator.HistoricalSimulator();
            Assert.AreEqual("Historical Simulator", sim.ProcessInfo.ProcessType);
        }

        [Test]
        public void TestDefaultStartDateIsToday()
        {
            var sim = new global::HistoricalSimulator.HistoricalSimulator();
            Assert.AreEqual(DateTime.Now.Date, sim.StartDate);
        }
    }
}
