using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace Dupire
{
    /// <summary>
    /// Tests <see cref="DupireSymbolChoice"/>.
    /// </summary>
    [TestFixture]
    public class TestDupireSymbolChoice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestDescription()
        {
            DupireSymbolChoice choice = new DupireSymbolChoice();

            Assert.AreEqual("Equity/Dupire Local Volatility Model", choice.Description);
        }

        [Test]
        public void TestCreateInstance()
        {
            DupireSymbolChoice choice = new DupireSymbolChoice();

            IEditable instance = choice.CreateInstance();

            Assert.IsNotNull(instance);
            StochasticProcessExtendible extendible = instance as StochasticProcessExtendible;
            Assert.IsNotNull(extendible);
            Assert.IsInstanceOf<DupireProcess>(extendible.Plugin);
        }
    }
}
