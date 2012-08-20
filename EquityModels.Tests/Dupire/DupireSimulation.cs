/* Copyright (C) 2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s):
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
using System.Collections.Generic;
using System.Text;
using DVPLDOM;
using DVPLI;
using DVPLSolver;
using DVPLUtils;
using Fairmat.Finance;
using Mono.Addins;
using NUnit.Framework;

namespace Dupire
{
    /// <summary>
    /// Tests Dupire model simulation
    /// </summary>
    [TestFixture]
    public class TestDupireSimulation
    {
        [SetUp]
        public void Init()
        {
            DVPLI.PluginsManager.Init();
            Mono.Addins.AddinManager.Registry.ResetConfiguration();
            Mono.Addins.AddinManager.Registry.Update(new Mono.Addins.ConsoleProgressStatus(6));

            //setup fake credentials for fairmat.com
            //FairmatEstimateDB.FairmatComIntegration credentials = DVPLI.UserSettings.GetSettings(typeof(FairmatEstimateDB.FairmatComIntegration)) as FairmatEstimateDB.FairmatComIntegration;
            //credentials.Username ="a@b.com";
            //credentials.Password = new DVPLI.Password("12345");
            DVPLI.Engine.Parser.NewContext();
        }

        [Test]
        public void TestSimulation()
        {
            Console.WriteLine("Inizio Test");
            // this test compares the price of a call option optained through monte carlo simulation
            // of a Dupire model with constant local volatility with the Black-Scholes price

            Engine.MultiThread = true;

            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            doc.DefaultProject.NMethods.m_UseAntiteticPaths = true;
            int n_sim = 100;
            int n_steps = 256;
            double strike = 90.0;
            double maturity = 2.0;
            double rate = 0.1;
            double dy = 0.07;
            double volatility = 0.2;
            double S0 = 100;

            ModelParameter Pstrike = new ModelParameter(strike, "","strike");
            rov.Symbols.Add(Pstrike);

            ModelParameter PS0 = new ModelParameter(S0, "","S0");
            rov.Symbols.Add(PS0);

            AFunction payoff = new AFunction(rov);
            payoff.VarName = "payoff";
            payoff.m_IndependentVariables = 1;
            payoff.m_Value =(RightValue) ("max(x1 - strike ; 0)");
            rov.Symbols.Add(payoff);

            AFunction rfunc = new AFunction(rov);
            rfunc.VarName = "r";
            rfunc.m_IndependentVariables = 1;
            rfunc.m_Value = (RightValue) rate;
            rov.Symbols.Add(rfunc);

            AFunction qfunc = new AFunction(rov);
            qfunc.VarName = "q";
            qfunc.m_IndependentVariables = 1;
            qfunc.m_Value = (RightValue) dy;
            rov.Symbols.Add(qfunc);
			
            AFunction volfunc = new AFunction(rov);
            volfunc.VarName = "localvol";
            volfunc.m_IndependentVariables = 2;
            volfunc.m_Value = (RightValue) volatility;
            rov.Symbols.Add(volfunc);

            DupireProcess process = new DupireProcess();
            process.s0 = (ModelParameter) "S0";
            process.r = (ModelParameter) "@r";
            process.q = (ModelParameter) "@q";
            process.localVol = (ModelParameter) "@localvol";

            StochasticProcessExtendible s = new StochasticProcessExtendible(rov, process);
            rov.Processes.AddProcess(s);

            //discounting
            RiskFreeInfo rfi = rov.GetDiscountingModel() as RiskFreeInfo;
            rfi.ActualizationType = EActualizationType.RiskFree;
            rfi.m_deterministicRF = rate;
            OptionTree op = new OptionTree(rov);
            op.PayoffInfo.PayoffExpression = "payoff(v1)";
            op.PayoffInfo.Timing.EndingTime.m_Value = (RightValue) maturity;
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
            double SamplePrice = price.m_Value;
            double SampleDevSt = price.m_StdErr/Math.Sqrt((double) n_sim);

            // calculation of the theoretical value of the call
            double ThPrice = BlackScholes.Call(rate, S0, strike, volatility, maturity); // todo: calcolare prezzo di BS
            Console.WriteLine("Theoretical Price  = " + ThPrice.ToString());
            Console.WriteLine("Monte Carlo Price  = " + SamplePrice);
            Console.WriteLine("Standard Deviation = " + SampleDevSt.ToString());
            double tol = 4.0 * SampleDevSt;
            Assert.LessOrEqual(Math.Abs(ThPrice - SamplePrice), tol);

        }
    }
}
