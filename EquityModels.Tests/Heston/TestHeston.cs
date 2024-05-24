/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Enrico Degiuli (enrico.degiuli@fairmat.com)
 *            Matteo Tesser (matteo.tesser@fairmat.com)
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
using DVPLDOM;
using DVPLI;
using DVPLSolver;
using HestonEstimator;
using HestonExtended;
using NUnit.Framework;

namespace Heston
{
    /// <summary>
    /// This test compares the price of a call option optained through monte carlo simulation
    /// and the price of the same option obtained through analytic results.
    /// </summary>
    [TestFixture]
    public class TestHeston
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Category("BigTest")]
        [TestCase(100, 5, 0.1)]
        [TestCase(80,1,0.01)]
        [TestCase(150,0.5,0.1)]
        public void TestCall(double strike, double tau, double r)
        {
            Engine.MultiThread = true;
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            doc.DefaultProject.NMethods.m_UseAntiteticPaths = true;

            int n_sim = 50000;
            int n_steps = 512;
           
            double rate = r;
            double dy = 0.00;

            ModelParameter pStrike = new ModelParameter(strike, "strike");
            pStrike.VarName = "strike";
            rov.Symbols.Add(pStrike);

            ModelParameter pRate = new ModelParameter(rate, "rfrate");
            pRate.VarName = "rfrate";
            rov.Symbols.Add(pRate);

            AFunction payoff = new AFunction(rov);
            payoff.VarName = "payoff";
            payoff.m_IndependentVariables = 1;
            payoff.m_Value = (RightValue)("max(x1 - strike ; 0)");
            rov.Symbols.Add(payoff);

            HestonProcess process = new HestonProcess();
            process.r = (ModelParameter)rate;
            process.q = (ModelParameter)dy;
            process.k = (ModelParameter)2.5;
            process.theta = (ModelParameter)0.4;
            process.sigma = (ModelParameter)0.2;
            process.S0 = (ModelParameter)100.0;
            process.V0 = (ModelParameter)0.3;


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

            // Calculates the theoretical value of the call.
            Vector param = new Vector(5);
            param[0] = process.k.V();
            param[1] = process.theta.V();
            param[2] = process.sigma.V();
            param[3] = 0.0;
            param[4] = process.V0.V();
            double thPrice = HestonCall.HestonCallPrice(param, process.S0.V(),
                                                        tau, strike, rate, dy);
            Console.WriteLine("Theoretical Price = " + thPrice.ToString());
            Console.WriteLine("Monte Carlo Price = " + samplePrice);
            Console.WriteLine("Standard Deviation = " + sampleDevSt.ToString());
            double tol = 4.0 * sampleDevSt;

            Assert.Less(Math.Abs(thPrice - samplePrice), tol);
        }

        [Category("BigTest")]
        [TestCase(100, 5, 0.1)]
        [TestCase(80, 1, 0.01)]
        [TestCase(150, 0.5, 0.1)]
        public void TestPut(double strike, double tau, double r)
        {
            Engine.MultiThread = true;
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            doc.DefaultProject.NMethods.m_UseAntiteticPaths = true;

            int n_sim = 50000;
            int n_steps = 512;

            double rate = r;
            double dy = 0.05;

            ModelParameter pStrike = new ModelParameter(strike, "strike");
            pStrike.VarName = "strike";
            rov.Symbols.Add(pStrike);

            ModelParameter pRate = new ModelParameter(rate, "rfrate");
            pRate.VarName = "rfrate";
            rov.Symbols.Add(pRate);

            AFunction payoff = new AFunction(rov);
            payoff.VarName = "payoff";
            payoff.m_IndependentVariables = 1;
            payoff.m_Value = (RightValue)("max( -x1 + strike ; 0)");
            rov.Symbols.Add(payoff);

            HestonProcess process = new HestonProcess();
            process.r = (ModelParameter)rate;
            process.q = (ModelParameter)dy;
            process.k = (ModelParameter)2.5;
            process.theta = (ModelParameter)0.4;
            process.sigma = (ModelParameter)0.2;
            process.S0 = (ModelParameter)100.0;
            process.V0 = (ModelParameter)0.3;


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

            // Calculates the theoretical value of the put.
            Vector param = new Vector(5);
            param[0] = process.k.V();
            param[1] = process.theta.V();
            param[2] = process.sigma.V();
            param[3] = 0.0;
            param[4] = process.V0.V();
            double thPrice = HestonCall.HestonPutPrice(param, process.S0.V(),
                                                        tau, strike, rate, dy);

            Console.WriteLine("Theoretical Price = " + thPrice.ToString());
            Console.WriteLine("Monte Carlo Price = " + samplePrice);
            Console.WriteLine("Standard Deviation = " + sampleDevSt.ToString());
            double tol = 4.0 * sampleDevSt;

            Assert.Less(Math.Abs(thPrice - samplePrice), tol);
        }



        [TestCase(1.0, 1, 0.5, 0.05, 19.41)]
        public void TestFwdApproximation(double strike, double T, double T0,  double r, double expectedResult)
        {
            double rate = r;
            double dy = 0.00;

            var process = new HestonProcess();
            process.r = (ModelParameter)rate;
            process.q = (ModelParameter)dy;
            process.k = (ModelParameter)2.5;
            process.theta = (ModelParameter)0.4;
            process.sigma = (ModelParameter)0.2;
            process.S0 = (ModelParameter)100.0;
            process.V0 = (ModelParameter)0.3;
            process.rho = (ModelParameter)(-0.5);


           

            // Calculates the theoretical value of the call.
            Vector param = new Vector(5);
            param[0] = process.k.V();
            param[1] = process.theta.V();
            param[2] = process.sigma.V();
            param[3] = process.rho.V(); 
            param[4] = process.V0.V();
            double thPrice = HestonForwardApproximated.HestonForwardCallPrice(param, process.S0.V(),
                                                                       T, T0, strike, rate, dy);


            Console.WriteLine("Theoretical Price = " + thPrice.ToString());
            Console.WriteLine("Expected Price = " + expectedResult.ToString());

        }


    }
}
