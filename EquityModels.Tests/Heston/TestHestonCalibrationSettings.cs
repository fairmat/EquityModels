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

using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using DVPLI;
using Heston;
using NUnit.Framework;

namespace HestonEstimator
{
    /// <summary>
    /// Tests <see cref="Heston.HestonCalibrationSettings"/>.
    /// </summary>
    [TestFixture]
    public class TestHestonCalibrationSettings
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestDefaultValues()
        {
            var settings = new HestonCalibrationSettings();

            Assert.AreEqual(0.4, settings.MinStrike);
            Assert.AreEqual(1.6, settings.MaxStrike);
            Assert.AreEqual(1.0 / 12, settings.MinMaturity);
            Assert.AreEqual(6, settings.MaxMaturity);
            Assert.IsFalse(settings.DummyCalibration);
        }

        [Test]
        public void TestIsEstimationSettings()
        {
            var settings = new HestonCalibrationSettings();

            Assert.IsInstanceOf<IEstimationSettings>(settings);
        }

        [Test]
        public void TestPropertiesAreSettable()
        {
            var settings = new HestonCalibrationSettings
            {
                MinStrike = 0.2,
                MaxStrike = 1.8,
                MinMaturity = 0.5,
                MaxMaturity = 10,
                DummyCalibration = true
            };

            Assert.AreEqual(0.2, settings.MinStrike);
            Assert.AreEqual(1.8, settings.MaxStrike);
            Assert.AreEqual(0.5, settings.MinMaturity);
            Assert.AreEqual(10, settings.MaxMaturity);
            Assert.IsTrue(settings.DummyCalibration);
        }

        [Test]
        public void TestIsSerializable()
        {
            var settings = new HestonCalibrationSettings
            {
                MinStrike = 0.3,
                MaxStrike = 1.7,
                MinMaturity = 0.2,
                MaxMaturity = 8,
                DummyCalibration = true
            };

            var formatter = new BinaryFormatter();
            using (var stream = new MemoryStream())
            {
                formatter.Serialize(stream, settings);
                stream.Position = 0;
                var deserialized = (HestonCalibrationSettings)formatter.Deserialize(stream);

                Assert.AreEqual(settings.MinStrike, deserialized.MinStrike);
                Assert.AreEqual(settings.MaxStrike, deserialized.MaxStrike);
                Assert.AreEqual(settings.MinMaturity, deserialized.MinMaturity);
                Assert.AreEqual(settings.MaxMaturity, deserialized.MaxMaturity);
                Assert.AreEqual(settings.DummyCalibration, deserialized.DummyCalibration);
            }
        }
    }
}
