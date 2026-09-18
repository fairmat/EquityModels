using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace HistoricalSimulator
{
    /// <summary>
    /// Tests <see cref="HistoricalSimulatorSymbolChoice"/>.
    /// </summary>
    [TestFixture]
    public class TestHistoricalSimulatorSymbolChoice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestDescription()
        {
            HistoricalSimulatorSymbolChoice choice = new HistoricalSimulatorSymbolChoice();

            Assert.AreEqual("Historical Simulator", choice.Description);
        }

        [Test]
        public void TestCreateInstance()
        {
            HistoricalSimulatorSymbolChoice choice = new HistoricalSimulatorSymbolChoice();

            IEditable instance = choice.CreateInstance();

            Assert.IsNotNull(instance);
            StochasticProcessExtendible extendible = instance as StochasticProcessExtendible;
            Assert.IsNotNull(extendible);
            Assert.IsInstanceOf<HistoricalSimulator>(extendible.Plugin);
        }
    }
}
