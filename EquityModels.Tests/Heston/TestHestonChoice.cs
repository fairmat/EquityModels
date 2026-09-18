using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace Heston
{
    /// <summary>
    /// Tests <see cref="HestonChoice"/>.
    /// </summary>
    [TestFixture]
    public class TestHestonChoice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestDescription()
        {
            HestonChoice choice = new HestonChoice();

            Assert.AreEqual("Equity/Heston", choice.Description);
        }

        [Test]
        public void TestCreateInstance()
        {
            HestonChoice choice = new HestonChoice();

            IEditable instance = choice.CreateInstance();

            Assert.IsNotNull(instance);
            StochasticProcessExtendible extendible = instance as StochasticProcessExtendible;
            Assert.IsNotNull(extendible);
            Assert.IsInstanceOf<HestonProcess>(extendible.Plugin);
        }
    }
}
