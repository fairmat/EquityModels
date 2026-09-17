using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace HestonExtended
{
    /// <summary>
    /// Tests <see cref="HestonExtendedChoice"/>.
    /// </summary>
    [TestFixture]
    public class TestHestonExtendedChoice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestDescription()
        {
            HestonExtendedChoice choice = new HestonExtendedChoice();

            Assert.AreEqual("Equity/Heston (with time dependent drift)", choice.Description);
        }

        [Test]
        public void TestCreateInstance()
        {
            HestonExtendedChoice choice = new HestonExtendedChoice();

            IEditable instance = choice.CreateInstance();

            Assert.IsNotNull(instance);
            StochasticProcessExtendible extendible = instance as StochasticProcessExtendible;
            Assert.IsNotNull(extendible);
            Assert.IsInstanceOf<HestonExtendedProcess>(extendible.Plugin);
        }
    }
}
