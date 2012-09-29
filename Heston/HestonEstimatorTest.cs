/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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
using System.Text;
using Fairmat.Optimization;
using DVPLI;
using DVPLDOM;
using HestonExtended;
using Fairmat.Math;
using Fairmat.MarketData;

namespace HestonEstimator
{
    public class MainClass
    {
        public static void Main(string[] args)
        {
            int Caso = 1;
            if (Caso==0)
            {
                InterestRateMarketData MData = InterestRateMarketData.FromFile("../../../TestData/InterestRatesModels/05052009-EU.xml");
                CallPriceMarketData test = CallPriceMarketData.FromFile("../../../TestData/Heston/05052009-SX5E-HestonData.xml");
                EquityCalibrationData CalData = new EquityCalibrationData(test, MData);

                Matrix CallMarketPrice = (Matrix)test.CallPrice;
                Vector Maturity = (Vector)test.Maturity;
                Vector Strike = (Vector)test.Strike;
                Vector DividendYield = (Vector) test.DividendYield;
                Vector Drift = CalData.Rate - CalData.DividendYield;
                Vector Rate = CalData.Rate;

                double u, kappa, theta, sigma, rho, v0, s0, r, q, T, K, val;
                u = 1.0;
                kappa = 19.4;
                theta = 0.235;
                sigma = 0.00500999;
                rho = -0.96;
                v0 = 0.664;
                s0 = 3872.15;
                r = -0.0867303549580581;
                q = 0;
                T = 0.50;
                K = 6000;
                Vector MatBound = new Vector(2);
                Vector StrikeBound = new Vector(2);
                MatBound[0] = 0.0;
                MatBound[1] = 2.0;
                StrikeBound[0] = 0.7;
                StrikeBound[1] = 1.3;
                Matrix Volatility = new Matrix(test.CallPrice.R,test.CallPrice.C);
                HestonCallOptimizationProblem HP = new HestonCallOptimizationProblem(CallMarketPrice, Maturity, Strike, Rate, DividendYield,test.S0,MatBound,StrikeBound, Volatility);
                Complex Cval, Cu;
                Cu = u - Complex.I;
                HestonCall hc = new HestonCall(HP);

                Cval = hc.phi(u,kappa,theta,sigma,rho,s0,v0,r,T);
                Console.WriteLine("phi1 = {0}",Cval);
                Cval = hc.phi(Cu,kappa,theta,sigma,rho,s0,v0,r,T);
                Console.WriteLine("phi2 = {0}",Cval);
                val = hc.IntegrandFunc(u,kappa,theta,sigma,rho,s0,v0,r,q,T,K);
                Console.WriteLine("IntFunc = {0}",val);

                Vector x = new Vector(5);
                x[0] = kappa;
                x[1] = theta;
                x[2] = sigma;
                x[3] = rho;
                x[4] = v0;

                DateTime T1,T2;
                TimeSpan ElapsedTime;
                double Time, Time2, Time3;

                T1 = DateTime.Now;
                val = hc.HestonCallPrice(x,s0,T,K,r,q);
                T2 = DateTime.Now;
                ElapsedTime = T2-T1;
                Time = (double) ElapsedTime.Milliseconds;
                Time2 = (double) ElapsedTime.Seconds;
                Console.WriteLine("Price = {0}",val);
                Console.WriteLine("Elapsed Time = {0}",Time2+Time/1000);

                int NProve = 10;
                int NPassi = 1000;
                double val2;
                Random CasNum = new Random();
                for (int i=0; i<NProve; i++)
                {
                    for (int j=0; j<5; j++)
                    {
                        val2 = ((double) CasNum.Next(0,NPassi))/((double) NPassi);
                        x[j] = HP.Bounds.Lb[j] + (HP.Bounds.Ub[j]-HP.Bounds.Lb[j])*val2;
                    }
                    Console.Write( "Trial {0}  x = " + x.ToString(),i+1);
                    T1 = DateTime.Now;
                    val = HP.Obj(x);
                    T2 = DateTime.Now;
                    ElapsedTime = T2-T1;
                    Time = (double) ElapsedTime.Milliseconds;
                    Time2 = (double) ElapsedTime.Seconds;
                    Time3 = (double) ElapsedTime.Minutes;
                    Console.WriteLine("  Time = {0}' {1}'' Val = {2}",Time3,Time2+Time/1000,val);
                }

            }
            if (Caso==1)
            {
                TestHestonCallEstimation NewTest = new TestHestonCallEstimation();
                bool Result = NewTest.Run();
            }
        }
    }
}
