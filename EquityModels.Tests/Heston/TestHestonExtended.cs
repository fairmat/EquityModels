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
    public class TestHestonExtended
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test, Category("BigTest")]
        public void Test()
        {
            Engine.MultiThread = true;

            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            doc.DefaultProject.NMethods.m_UseAntiteticPaths = true;

            int n_sim = 100000;
            int n_steps = 256;
            double strike = 90.0;
            double tau = 2.0;
            double rate = 0.1;
            double dy = 0.07;

            ModelParameter pStrike = new ModelParameter(strike, "strike");
            pStrike.VarName = "strike";
            rov.Symbols.Add(pStrike);

            AFunction payoff = new AFunction(rov);
            payoff.VarName = "payoff";
            payoff.m_IndependentVariables = 1;
            payoff.m_Value = (RightValue)("max(x1 - strike ; 0)");
            rov.Symbols.Add(payoff);

            AFunction zrfunc = new AFunction(rov);
            zrfunc.VarName = "zr";
            zrfunc.m_IndependentVariables = 1;
            zrfunc.m_Value = (RightValue)rate;
            rov.Symbols.Add(zrfunc);

            AFunction dyfunc = new AFunction(rov);
            dyfunc.VarName = "dy";
            dyfunc.m_IndependentVariables = 1;
            dyfunc.m_Value = (RightValue)dy;
            rov.Symbols.Add(dyfunc);

            HestonExtendedProcess process = new HestonExtendedProcess();
            process.k = (ModelParameter)2.5;
            process.theta = (ModelParameter)0.4;
            process.sigma = (ModelParameter)0.2;
            process.S0 = (ModelParameter)100.0;
            process.V0 = (ModelParameter)0.3;
            process.zrReference = (ModelParameter)"@zr";
            process.dyReference = (ModelParameter)"@dy";
            double discount = Math.Exp(-rate * tau);

            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

            // Set the discounting.
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.RiskFree;
            rfi.m_deterministicRF = 0.0;

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
            double samplePrice = discount * price.m_Value;
            double sampleDevSt = price.m_StdErr / Math.Sqrt((double)n_sim);

            // Calculates the theoretical value of the call.
            Vector param = new Vector(5);
            param[0] = process.k.V();
            param[1] = process.theta.V();
            param[2] = process.sigma.V();
            param[3] = 0.0;
            param[4] = process.V0.V();
            HestonCall hestonCall = new HestonCall();
            double theoreticalPrice = hestonCall.HestonCallPrice(param, process.S0.V(),
                                                                 tau, strike, rate, dy);
            Console.WriteLine("Theoretical Price = " + theoreticalPrice.ToString());
            Console.WriteLine("Monte Carlo Price = " + samplePrice);
            Console.WriteLine("Standard Deviation = " + sampleDevSt.ToString());
            double tol = 4.0 * sampleDevSt;
            Assert.Less(Math.Abs(theoreticalPrice - samplePrice), tol);
        }
    }
}
