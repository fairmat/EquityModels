using DVPLDOM;
using DVPLI;
using NUnit.Framework;

namespace VarianceGamma
{
    /// <summary>
    /// Tests <see cref="VarianceGammaChoice"/>.
    /// </summary>
    [TestFixture]
    public class TestVarianceGammaChoice
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestDescription()
        {
            VarianceGammaChoice choice = new VarianceGammaChoice();

            Assert.AreEqual("Equity/Variance Gamma Model", choice.Description);
        }

        [Test]
        public void TestCreateInstance()
        {
            VarianceGammaChoice choice = new VarianceGammaChoice();

            IEditable instance = choice.CreateInstance();

            Assert.IsNotNull(instance);
            StochasticProcessExtendible extendible = instance as StochasticProcessExtendible;
            Assert.IsNotNull(extendible);
            Assert.IsInstanceOf<VarianceGamma>(extendible.Plugin);
        }
    }
}
