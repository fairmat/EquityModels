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
using DVPLDOM;
using DVPLI;
using Fairmat.MarketData;
using Fairmat.Math;

namespace Dupire
{
    public class DupireEstimator:IEstimator, IIntegrable
    {
        public DupireEstimator ()
        {
        }

        PFunction r;
        PFunction q;

        #region IEstimator implementation
        public Type[] GetRequirements (IEstimationSettings settings, bool multivariateRequest)
        {
            return new Type[] { typeof(InterestRateMarketData),typeof(CallPriceMarketData) };
        }

        public EstimationResult Estimate (System.Collections.Generic.List<object> MarketData, IEstimationSettings settings)
        {
            Console.WriteLine("Estimate:");
            InterestRateMarketData Mdataset = (InterestRateMarketData)MarketData[0];
            CallPriceMarketData Hdataset = (CallPriceMarketData)MarketData[1];
            EquityCalibrationData HCalData = new EquityCalibrationData(Hdataset, Mdataset);
            Console.WriteLine("Estimate ok1");
            r = new DVPLDOM.PFunction(null);
            q = new DVPLDOM.PFunction(null);
            r.Expr = (double[,]) ArrayHelper.Concat(HCalData.MaturityDY.ToArray(), HCalData.Rate.ToArray());
            q.Expr = (double[,]) ArrayHelper.Concat(HCalData.MaturityDY.ToArray(), HCalData.DividendYield.ToArray());
            PFunction2D.PFunction2D impVol = new PFunction2D.PFunction2D(Hdataset.Maturity, Hdataset.Strike, Hdataset.Volatility);
            Document doc= new Document();
            ProjectROV prj= new ProjectROV(doc);
            doc.Part.Add(prj);
            prj.Symbols.Add(impVol);
            doc.WriteToXMLFile("impVol.fair");
            
            // todo: spostare nei settings
            int nmat = 100;
            int nstrike = 100;
            Vector locVolMat = Vector.Linspace(0.0, Hdataset.Maturity[Hdataset.Maturity.Length-1], nmat);
            Vector locVolStr = Vector.Linspace(0.0, Hdataset.Strike[Hdataset.Strike.Length-1], nstrike);
            Matrix locVolMatrix = new Matrix(nmat,nstrike);
            Integrate integrate = new Integrate(this);
            double sigma, dSigmadk, num, y, den, integral;
            Vector x = new Vector(2);
            for (int i = 0; i < nmat; i++)
            {
                Console.WriteLine("Estimate ok10, i = " + i);
                integral = integrate.AdaptLobatto(0.0, locVolMat[i]);
                Console.WriteLine("Estimate ok11, integral = " + integral);
                for (int j = 0; j < nstrike; j++)
                {
                    Console.WriteLine("Estimate ok12, j = " + j);
                    x[0] = locVolMat[i];
                    Console.WriteLine("Estimate ok13");
                    x[1] = locVolStr[j];
                    Console.WriteLine("Estimate ok14");
                    Console.WriteLine("x = " + x.ToString());
                    sigma = impVol.Evaluate(x);
                    Console.WriteLine("Estimate ok15");
                    dSigmadk = impVol.Partial(x, 1);
                    Console.WriteLine("Estimate ok16");
                    num = Math.Pow(sigma, 2) + 2.0 * sigma * x[0] *
                        ( impVol.Partial(x, 0) + ( r.Evaluate(x[0]) - q.Evaluate(x[0])) * x[1] * dSigmadk );
                    Console.WriteLine("Estimate ok17");
                    y = Math.Log( locVolStr[j] / Hdataset.S0 ) + integral;
                    Console.WriteLine("Estimate ok18");
                    den = System.Math.Pow( 1.0 - x[1] * y * dSigmadk / sigma , 2) + x[1] * sigma * x[0] *
                        ( dSigmadk - 0.25 * x[1] * sigma * x[0] * dSigmadk * dSigmadk + x[1] * impVol.Partial2(x, 1) );
                    Console.WriteLine("Estimate ok19");
                    locVolMatrix[i,j] = Math.Sqrt( num / den );
                    Console.WriteLine("Estimate ok20");
                }
            }

            //create dupire outputs
            PFunction2D.PFunction2D localVol = new PFunction2D.PFunction2D();
            //localVol.


            string[] names= new string[]{"S0"};
            Vector param = new Vector(1);
            param[0] = Hdataset.S0;
            EstimationResult result = new EstimationResult(names, param);
            result.Objects = new object[3];
            result.Objects[0] = r;
            result.Objects[1] = q;
            result.Objects[2] = localVol;
            Console.WriteLine("Estimate ok30");
            return result;
        }

        public Type ProvidesTo {
            get {
                throw new System.NotImplementedException ();
            }
        }
        #endregion

        #region IIntegrable implementation
        double IIntegrable.IntegrandFunc (double x)
        {
            return (q.Evaluate(x) - r.Evaluate(x));
        }
        #endregion

    }
}

