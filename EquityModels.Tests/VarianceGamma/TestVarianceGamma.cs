/* Copyright (C) 2013 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Safe Khampol (safe.khampol@gmail.com)
 *            Matteo Tesser (matteo.tesser@fairmat.com)
 *            Enrico Degiuli (enrico.degiuli@fairmat.com)
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
using NUnit.Framework;

namespace VarianceGamma
{
    [TestFixture]
    public class TestVarianceGamma
    {
        [SetUp]
        public void Init()
        {
            DVPLI.PluginsManager.Init();
            Mono.Addins.AddinManager.Registry.ResetConfiguration();
            Mono.Addins.AddinManager.Registry.Update(new Mono.Addins.ConsoleProgressStatus(0));
        }

        [Test]
        public void Test()
        {
            double nu = 0.6;
            double theta = -0.2;
            double sigma = 0.2;
            double rate = 0.02;
            double dy = 0.01;
            double s0 = 1;
            double maturity = 2.0;
            double strike = 1.2;
            Vector mat = new Vector(1) + maturity;
            Vector k = new Vector(1) + strike;
            Matrix cp = new Matrix(1, 1) + 0.3;

            // Calculates the theoretical value of the call.
            double theoreticalPrice = VarianceGammaOptionsCalibration.VGCall(theta, sigma, nu,
                                                                             maturity, strike,
                                                                             dy, s0, rate);

            Engine.MultiThread = true;
            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            doc.DefaultProject.NMethods.m_UseAntiteticPaths = true;

            int n_sim = 50000;
            int n_steps = 512;

            ModelParameter paramStrike = new ModelParameter(strike, "strike");
            paramStrike.VarName = "strike";
            rov.Symbols.Add(paramStrike);

            ModelParameter paramRate = new ModelParameter(rate, "rfrate");
            paramRate.VarName = "rfrate";
            rov.Symbols.Add(paramRate);

            AFunction payoff = new AFunction(rov);
            payoff.VarName = "payoff";
            payoff.m_IndependentVariables = 1;
            payoff.m_Value = (RightValue)("max(x1 - strike ; 0)");
            rov.Symbols.Add(payoff);

            VarianceGamma process = new VarianceGamma(s0, theta, sigma, nu, rate, dy);

            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

            // Set the discounting.
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.RiskFree;
            rfi.m_deterministicRF = rate;

            OptionTree op = new OptionTree(rov);
            op.PayoffInfo.PayoffExpression = "payoff(v1)";
            op.PayoffInfo.Timing.EndingTime.m_Value = (RightValue)maturity;
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
                rov.DisplayErrors();
            }

            Assert.IsFalse(rov.HasErrors);

            ResultItem price = rov.m_ResultList[0] as ResultItem;

            double samplePrice = price.value;
            double sampleDevSt = price.stdDev / Math.Sqrt((double)n_sim);

            Console.WriteLine("Theoretical Price = " + theoreticalPrice);
            Console.WriteLine("Monte Carlo Price = " + samplePrice);
            Console.WriteLine("Standard Deviation = " + sampleDevSt.ToString());
            double tol = 4.0 * sampleDevSt;
            Assert.Less(Math.Abs(theoreticalPrice - samplePrice), tol);
        }
    }
}
