using System;
using DVPLDOM;
using DVPLI;
using DVPLSolver;
using HestonEstimator;
using NUnit.Framework;

namespace Heston
{
    [TestFixture]
    public class TestHestonDigital
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestTheoreticalPrice()
        {
            double k = 90;
            double tau = 2.0;

            double rate = 0.01;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.8;

            // Calculates the theoretical value of the call.
            Vector param = new Vector(5);
            param[0] = kappa;
            param[1] = theta;
            param[2] = sigma;
            param[3] = rho;
            param[4] = v0;
            var hestonDigital = new HestonDigital();
            double fairmatPrice = hestonDigital.HestonDigitalCallPrice(param, s0, tau, k, rate, dy);
            double tol = 1e-3;
            double benchmarkPrice = 0.331117450319339;

            Console.WriteLine("Theoretical Benchmark  Price = " + benchmarkPrice.ToString());
            Console.WriteLine("Theoretical Fairmat    Price = " + fairmatPrice);

            Assert.Less(Math.Abs(fairmatPrice - benchmarkPrice), tol);
        }


        [TestCase(1,100.0)]
        [TestCase(5, 100.0)]
        [TestCase(1, 80.0)]
        [TestCase(.5, 80.0)]
        public void TestTheoreticalPriceVsMC(double tau, double k)
        {

            double rate = 0.01;
            double dy = 0.07;
            double kappa = 2.5;
            double theta = 0.4;
            double sigma = 0.2;
            double s0 = 100.0;
            double v0 = 0.3;
            double rho = -0.3;

            // Calculates the theoretical value of the call.
            Vector param = new Vector(5);
            param[0] = kappa;
            param[1] = theta;
            param[2] = sigma;
            param[3] = rho;
            param[4] = v0;
            var hestonDigital = new HestonDigital();
            double fairmatPrice = hestonDigital.HestonDigitalCallPrice(param, s0, tau, k, rate, dy);


            // set up MC simu

            Engine.MultiThread = true;
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            doc.DefaultProject.NMethods.m_UseAntiteticPaths = true;

            int n_sim = 50000;
            int n_steps = 512;
      

            ModelParameter pStrike = new ModelParameter(k, "strike");
            pStrike.VarName = "strike";
            rov.Symbols.Add(pStrike);

            ModelParameter pRate = new ModelParameter(rate, "rfrate");
            pRate.VarName = "rfrate";
            rov.Symbols.Add(pRate);

            AFunction payoff = new AFunction(rov);
            payoff.VarName = "payoff";
            payoff.m_IndependentVariables = 1;
            payoff.m_Value = (RightValue)("iif(x1 > strike ; 1 ; 0)");
            rov.Symbols.Add(payoff);

            HestonProcess process = new HestonProcess();
            process.r = (ModelParameter)rate;
            process.q = (ModelParameter)dy;
            process.k = (ModelParameter)kappa;
            process.theta = (ModelParameter)theta;
            process.sigma = (ModelParameter)sigma;
            process.S0 = (ModelParameter)s0;
            process.V0 = (ModelParameter)v0;


            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

            // Set the discounting.
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.RiskFree;
            rfi.m_deterministicRF = rate;

            OptionTree op = new OptionTree(rov);
            op.PayoffInfo.PayoffExpression = "payoff(v1a)";
            op.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)tau;
            op.PayoffInfo.European = true;
            rov.Map.Root = op;

            rov.NMethods.Technology = ETechType.T_SIMULATION;
            rov.NMethods.PathsNumber = n_sim;
            rov.NMethods.SimulationSteps = n_steps;

            ROVSolver solver = new ROVSolver();
            solver.BindToProject(rov);
            solver.DoValuation(-1);

            if (rov.HasErrors)
            {
                Console.WriteLine(rov.m_RuntimeErrorList[0]);
            }

            Assert.IsFalse(rov.HasErrors);

            ResultItem price = rov.m_ResultList[0] as ResultItem;

            double samplePrice = price.value;
            double sampleDevSt = price.stdDev / Math.Sqrt((double)n_sim);

           
            Console.WriteLine("Theoretical Price = " + fairmatPrice.ToString());
            Console.WriteLine("Monte Carlo Price = " + samplePrice);
            Console.WriteLine("Standard Deviation = " + sampleDevSt.ToString());
            double tol = 4.0 * sampleDevSt;

            Assert.Less(Math.Abs(fairmatPrice - samplePrice), tol);

        }
    }
}
