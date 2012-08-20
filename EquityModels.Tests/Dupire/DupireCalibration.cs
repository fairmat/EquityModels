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
using Mono.Addins;
using NUnit.Framework;
using Fairmat.MarketData;

namespace Dupire
{
    /// <summary>
    /// Tests Dupire model calibration
    /// </summary>
    [TestFixture]
    public class TestDupireCalibration
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
        public void TestCalibration()
        {
            
            InterestRateMarketData IData = InterestRateMarketData.FromFile("../../../../../fairmatplugins/TestData/InterestRatesModels/28-06-2012-EU-close.xml");
            CallPriceMarketData HData = CallPriceMarketData.FromFile("../../../../../fairmatplugins/TestData/Heston/28062012-SX5E_Index-HestonData.xml");
            
            List<object> l = new List<object>();
            l.Add(IData);
            l.Add(HData);
   
            DupireEstimator DE = new DupireEstimator();
            Console.WriteLine("TestCalibration ok1");
            EstimationResult res = DE.Estimate(l, null);
            Console.WriteLine("TestCalibration ok2");
            Assert.AreEqual(0,0);
        }
    }
}

