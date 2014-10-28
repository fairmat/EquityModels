/* Copyright (C) 2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Enrico Degiuli (enrico.degiuli@fairmat.com)
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
using DVPLDOM;
using DVPLI;
using DVPLSolver;
using Fairmat.Finance;
using Fairmat.MarketData;
using NUnit.Framework;

namespace Dupire
{
    /// <summary>
    /// Tests Dupire model calibration.
    /// </summary>
    [TestFixture]
    public class TestDupireCalibration
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        [Test]
        public void TestCalibrationFTSE()
        {
           CurveMarketData IData = MyObject.ReadFromXMLFile(@"c:\fermat\fairmat\prj\fairmatplugins/TestData/CallPrices/ZR20140922.xml") as CurveMarketData;
            CallPriceMarketData HData = MyObject.ReadFromXMLFile(@"c:\fermat\fairmat\prj\fairmatplugins/TestData/CallPrices/FTSE20140922.xml") as CallPriceMarketData;
           

            List<object> l = new List<object>();
            l.Add(IData);
            l.Add(HData);

            DupireEstimator DE = new DupireEstimator();
            DupireCalibrationSettings settings = new DupireCalibrationSettings();
            settings.LocalVolatilityCalculation = LocalVolatilityCalculation.Method1;
            settings.LocalVolatilityStrikes = 10;
            settings.LocalVolatilityMaturities = 10;

            EstimationResult res = DE.Estimate(l, settings);
            var localVol = res.Objects[3] as  PFunction2D.PFunction2D;
            //localVol.m
        }

        [Test]
        public void TestCalibration()
        {
            InterestRateMarketData IData = InterestRateMarketData.FromFile("../../TestData/IRMD-sample.xml");
            CallPriceMarketData HData = CallPriceMarketData.FromFile("../../TestData/CallData-sample.xml");
            //InterestRateMarketData IData = InterestRateMarketData.FromFile("../../../EquityModels.Tests/TestData/IRMD-EU-30102012-close.xml");
            //CallPriceMarketData HData = CallPriceMarketData.FromFile("../../../EquityModels.Tests/TestData/30102012-SX5E_Index-HestonData.xml");
            //CallPriceMarketData HData = ObjectSerialization.ReadFromXMLFile("../../../EquityModels.Tests/TestData/FTSE.xml") as CallPriceMarketData;


            List<object> l = new List<object>();
            l.Add(IData.DiscountingCurve);
            l.Add(HData);

            DupireEstimator DE = new DupireEstimator();
            DupireCalibrationSettings settings = new DupireCalibrationSettings();
            settings.LocalVolatilityCalculation = LocalVolatilityCalculation.Method1;
           

            //settings.LocalVolatilityCalculation = LocalVolatilityCalculation.QuantLib;
            EstimationResult res = DE.Estimate(l, settings);
            //int nmat = HData.Maturity.Length;
            //int nstrike = HData.Strike.Length;

            int i = 5; // Maturity.
            int j = 4; // Strike.

            Engine.MultiThread = true;

            Document doc = new Document();
            ProjectROV rov = new ProjectROV(doc);
            doc.Part.Add(rov);
            doc.DefaultProject.NMethods.m_UseAntiteticPaths = true;
            int n_sim = 10000;
            int n_steps = 500;
            double strike = HData.Strike[j];
            //double volatility = HData.Volatility[i, j];
            
            /*
            PFunction2D.PFunction2D impvolfunc = new PFunction2D.PFunction2D(rov);
            impvolfunc = res.Objects[3] as PFunction2D.PFunction2D;
            impvolfunc.VarName = "impvol";
            rov.Symbols.Add(impvolfunc);
            double volatility = impvolfunc.Evaluate(HData.Maturity[i], HData.Strike[j]);
             */
            double volatility = 0.2;
            double maturity = HData.Maturity[i];

            ModelParameter Pstrike = new ModelParameter(strike, string.Empty, "strike");
            rov.Symbols.Add(Pstrike);
            AFunction payoff = new AFunction(rov);
            payoff.VarName = "payoff";
            payoff.m_IndependentVariables = 1;
            payoff.m_Value = (RightValue)("max(x1 - strike ; 0)");
            rov.Symbols.Add(payoff);

            bool found;
            double S0 = PopulateHelper.GetValue("S0", res.Names, res.Values, out found);
            ModelParameter PS0 = new ModelParameter(S0, string.Empty, "S0");
            rov.Symbols.Add(PS0);
            PFunction rfunc = new PFunction(rov);
            rfunc = res.Objects[0] as PFunction;
            rfunc.VarName = "r";
            rov.Symbols.Add(rfunc);

            PFunction qfunc = new PFunction(rov);
            qfunc = res.Objects[1] as PFunction;
            qfunc.VarName = "q";
            rov.Symbols.Add(qfunc);

            PFunction2D.PFunction2D volfunc = new PFunction2D.PFunction2D(rov);
            volfunc = res.Objects[2] as PFunction2D.PFunction2D;
            volfunc.VarName = "localvol";
            rov.Symbols.Add(volfunc);
            DupireProcess process = new DupireProcess();
            process.s0 = (ModelParameter)"S0";
            process.r = (ModelParameter)"@r";
            process.q = (ModelParameter)"@q";
            process.localVol = (ModelParameter)"@localvol";
            double rate = rfunc.Evaluate(maturity);
            double dy = qfunc.Evaluate(maturity);

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

            Console.WriteLine("Surf = " + volfunc.Expr);

            // Calculation of the theoretical value of the call.
            double theoreticalPrice = BlackScholes.Call(rate, S0, strike, volatility, maturity, dy);
            Console.WriteLine("Theoretical Price  = " + theoreticalPrice.ToString());
            Console.WriteLine("Monte Carlo Price  = " + samplePrice);
            Console.WriteLine("Standard Deviation = " + sampleDevSt.ToString());
            double tol = 4.0 * sampleDevSt;
            doc.WriteToXMLFile("Dupire.fair");
            Assert.LessOrEqual(Math.Abs(theoreticalPrice - samplePrice), tol);
        }
    }
}
