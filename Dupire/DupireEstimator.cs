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
using DVPLDOM;
using DVPLI;
using Fairmat.MarketData;
using Fairmat.Math;
using Fairmat.Statistics;

namespace Dupire
{
    [Mono.Addins.Extension("/Fairmat/Estimator")]
    public class DupireEstimator : IEstimatorEx, IIntegrable
    {
        private PFunction r;
        private PFunction q;

        public DupireEstimator()
        {
        }

        public IEstimationSettings DefaultSettings
        {
            get
            {
                return UserSettings.GetSettings(typeof(DupireCalibrationSettings)) as DupireCalibrationSettings;
            }
        }


        #region IEstimator implementation

        /// <summary>
        /// Gets the types required by the estimator in order to work:
        /// InterestRateMarketData and CallPriceMarketData are the
        /// required types for this estimator.
        /// </summary>
        /// <param name="settings">The parameter is not used.</param>
        /// <param name="multivariateRequest">The parameter is not used.</param>
        /// <returns>
        /// An array containing the types InterestRateMarketData and CallPriceMarketData.
        /// </returns>
        public Type[] GetRequirements(IEstimationSettings settings, bool multivariateRequest)
        {
            return new Type[] { typeof(InterestRateMarketData), typeof(CallPriceMarketData) };
        }

        /// <summary>
        /// Gets the value requested by the interface ProvidesTo,
        /// returning <see cref="Dupire.DupireProcess"/> as the type.
        /// </summary>
        public Type ProvidesTo
        {
            get
            {
                return typeof(Dupire.DupireProcess);
            }
        }

        public EstimationResult Estimate(System.Collections.Generic.List<object> marketData, IEstimationSettings settings)
        {
            InterestRateMarketData Mdataset = (InterestRateMarketData)marketData[0];
            CallPriceMarketData Hdataset = (CallPriceMarketData)marketData[1];

            //gets the settings
            DupireCalibrationSettings calibrationSettings = settings as DupireCalibrationSettings;

            switch(calibrationSettings.LocalVolatilityCalculation)
            {
                case LocalVolatilityCalculation.Method1:
                    return FairmatEstimate(Mdataset, Hdataset);
                case LocalVolatilityCalculation.QuantLib:
                    return QuantLibEstimate(Mdataset, Hdataset);
                default:
                    throw new NotImplementedException("Method not implemented");
            }
        }

        #endregion

        private EstimationResult QuantLibEstimate(InterestRateMarketData Mdataset, CallPriceMarketData Hdataset)
        {
            EquityCalibrationData HCalData = new EquityCalibrationData(Hdataset, Mdataset);
            this.r = new DVPLDOM.PFunction(null);
            this.q = new DVPLDOM.PFunction(null);
            this.r.Expr = (double[,])ArrayHelper.Concat(HCalData.MaturityDY.ToArray(), HCalData.Rate.ToArray());
            this.q.Expr = (double[,])ArrayHelper.Concat(HCalData.MaturityDY.ToArray(), HCalData.DividendYield.ToArray());
            this.r.Parse(null);
            this.q.Parse(null);

            for (int i = 0; i < Hdataset.Volatility.R; i++)
            {
                for (int j = 0; j < Hdataset.Volatility.C; j++)
                {
                    Hdataset.Volatility[i,j] = Hdataset.Volatility[i,j]*Hdataset.Volatility[i,j]*Hdataset.Maturity[i];
                }
            }

            IFunction impVol = FitImplVolModel(Hdataset);

            Document doc = new Document();
            ProjectROV prj = new ProjectROV(doc);
            doc.Part.Add(prj);
            prj.Symbols.Add(impVol);
            //doc.WriteToXMLFile("impVol.fair");

            // todo: spostare nei settings
            int nmat = 100;
            int nstrike = 100;
            double lastMat = Hdataset.Maturity[Hdataset.Maturity.Length - 1];
            double lastStr = Hdataset.Strike[Hdataset.Strike.Length - 1];
            //Vector locVolMat = Vector.Linspace(lastMat / ((double) nmat), lastMat, nmat);
            //Vector locVolStr = Vector.Linspace(lastStr / ((double) nstrike), lastStr, nstrike);
            Vector locVolMat = Vector.Linspace(0.0, lastMat, nmat);
            Vector locVolStr = Vector.Linspace(0.0, lastStr, nstrike);
            Matrix locVolMatrix = new Matrix(nmat, nstrike);
            //Integrate integrate = new Integrate(this);
            //double sigma, dSigmadk, num, y, den, integral;
            double t, dt, forwardValue, y, dy, strike, strikep, strikem, w, wp, wm, dwdy;
            double d2wdy2, den1, den2, den3, strikept, strikemt, wpt, wmt, dwdt;
            Integrate integrate = new Integrate(this);

            for (int i = 0; i < nmat; i++)
            {
                t = locVolMat[i];
                forwardValue = Hdataset.S0 * Math.Exp(integrate.AdaptLobatto(0.0, t));
                for (int j = 0; j < nstrike; j++)
                {
                    strike = locVolStr[j];
                    y = Math.Log(strike / forwardValue );
                    dy = ( (Math.Abs(y) > 0.001) ? y * 0.0001 : 0.000001 );

                    // strike derivative
                    strikep = strike * Math.Exp(dy);
                    strikem = strike / Math.Exp(dy);
                    w = impVol.Evaluate(t, strike);
                    wp = impVol.Evaluate(t, strikep);
                    wm = impVol.Evaluate(t, strikem);
                    dwdy = (wp - wm) / (2.0 * dy);
                    d2wdy2 = (wp - 2.0 * w + wm) / (dy * dy);

                    // time derivative
                    if (t == 0.0)
                    {
                        dt = 0.0001;
                        strikept = strike * Math.Exp(integrate.AdaptLobatto(0.0, t + dt));
                        wpt = impVol.Evaluate(t + dt, strikept);
                        //if (wpt < w)
                        //    Console.WriteLine("Decreasing variance at strike {0} between time {1} and time {2}", strike, t, t + dt);
                        dwdt = (wpt - w) / dt;
                    }
                    else
                    {
                        dt = Math.Min(0.0001, t / 2.0);
                        strikept = strike * Math.Exp(integrate.AdaptLobatto(t, t + dt));
                        strikemt = strike * Math.Exp(-integrate.AdaptLobatto(t - dt, t));
                        wpt = impVol.Evaluate(t + dt, strikept);
                        wmt = impVol.Evaluate(t + dt, strikemt);

                        //if (wpt < w) 
                        //    Console.WriteLine("Decreasing variance at strike {0} between time {1} and time {2}", strike, t, t + dt);
                        //if (w < wmt) 
                        //    Console.WriteLine("Decreasing variance at strike {0} between time {1} and time {2}", strike, t-dt, t);
                        dwdt = (wpt - wmt) / (2.0 * dt);
                    }
                    if (dwdy == 0.0 && d2wdy2 == 0.0)
                        locVolMatrix[i, j] = Math.Sqrt(dwdt);
                    else
                    {
                        den1 = 1.0 - y / w * dwdy;
                        den2 = 0.25 * (-0.25 - 1.0 / w + y * y / w / w) * dwdy * dwdy;
                        den3 = 0.5 * d2wdy2;
                        locVolMatrix[i,j] = dwdt / (den1 + den2 + den3);
                        //if (locVolMatrix[i,j] < 0.0)
                        //    Console.WriteLine("Negative local vol^2 at strike {0} and time {1}; " +
                        //        "Black vol surface is not smooth enought.", strike, t);
                    }
                }
            }

            // Create dupire outputs.
            PFunction2D.PFunction2D localVol = new PFunction2D.PFunction2D(locVolMat, locVolStr, locVolMatrix);
            localVol.Parse(null);
            string[] names = new string[] { "S0" };
            Vector param = new Vector(1);
            param[0] = Hdataset.S0;
            EstimationResult result = new EstimationResult(names, param);
            result.Objects = new object[3];
            result.Objects[0] = this.r;
            result.Objects[1] = this.q;
            result.Objects[2] = localVol;
            Console.WriteLine("r = " + HCalData.Rate.ToString());
            Console.WriteLine("q = " + HCalData.DividendYield.ToString());
            return result;
        }


        private EstimationResult FairmatEstimate(InterestRateMarketData Mdataset, CallPriceMarketData Hdataset)
        {
            EquityCalibrationData HCalData = new EquityCalibrationData(Hdataset, Mdataset);
            this.r = new DVPLDOM.PFunction(null);
            this.q = new DVPLDOM.PFunction(null);
            this.r.Expr = (double[,])ArrayHelper.Concat(HCalData.MaturityDY.ToArray(), HCalData.Rate.ToArray());
            this.q.Expr = (double[,])ArrayHelper.Concat(HCalData.MaturityDY.ToArray(), HCalData.DividendYield.ToArray());
            this.r.Parse(null);
            this.q.Parse(null);

            IFunction impVol = FitImplVolModel(Hdataset);

            Document doc = new Document();
            ProjectROV prj = new ProjectROV(doc);
            doc.Part.Add(prj);
            prj.Symbols.Add(impVol);
            //doc.WriteToXMLFile("impVol.fair");

            // todo: spostare nei settings
            int nmat = 100;
            int nstrike = 100;
            double lastMat = Hdataset.Maturity[Hdataset.Maturity.Length - 1];
            double lastStr = Hdataset.Strike[Hdataset.Strike.Length - 1];
            Vector locVolMat = Vector.Linspace(0.0, lastMat, nmat);
            Vector locVolStr = Vector.Linspace(0.0, lastStr, nstrike);
            Matrix locVolMatrix = new Matrix(nmat, nstrike);
            Integrate integrate = new Integrate(this);
            double sigma, dSigmadk, num, y, den, integral;
            Vector x = new Vector(2);
            for (int i = 0; i < nmat; i++)
            {
                integral = integrate.AdaptLobatto(0.0, locVolMat[i]);

                int j = 0;
                x[0] = locVolMat[i];
                x[1] = locVolStr[j];
                sigma = impVol.Evaluate(x);
                dSigmadk = impVol.Partial(x, 1);
                num = Math.Pow(sigma, 2) + 2.0 * sigma * x[0] * impVol.Partial(x, 0);
                den = 1.0;
                locVolMatrix[i, j] = Math.Sqrt(num / den);

                // The rest of the cycle.
                for (j = 1; j < nstrike; j++)
                {
                    x[1] = locVolStr[j];
                    sigma = impVol.Evaluate(x);
                    dSigmadk = impVol.Partial(x, 1);
                    num = Math.Pow(sigma, 2) + 2.0 * sigma * x[0] *
                        (impVol.Partial(x, 0) + (this.r.Evaluate(x[0]) - this.q.Evaluate(x[0])) * x[1] * dSigmadk);
                    y = Math.Log(locVolStr[j] / Hdataset.S0) + integral;
                    den = System.Math.Pow(1.0 - x[1] * y * dSigmadk / sigma, 2) + x[1] * sigma * x[0] *
                        (dSigmadk - 0.25 * x[1] * sigma * x[0] * dSigmadk * dSigmadk + x[1] * impVol.Partial2(x, 1));
                    locVolMatrix[i, j] = Math.Sqrt(Math.Abs(num / den));
                    //if (double.IsNaN(locVolMatrix[i, j]))
                    //{
                    //    Console.WriteLine("num = {0}, den = {1}", num, den);
                    //    Console.WriteLine("locVolStr = {0}, locVolMat = {1}, y = {2}, dSigmadk = {3}, sigma = {4}, impvol.Partial2 = {5}", x[1], x[0], y, dSigmadk, sigma, impVol.Partial2(x,1));
                    //}
                    //Console.WriteLine("locVolMatrix[{0},{1}] = {2}, Part2 = {3}, t = {4}, S = {5}, dSigmadk = {6}, den = {7}, num = {8}, dsigma/dt = {9}", i, j, locVolMatrix[i, j], impVol.Partial2(x, 1), locVolMat[i], locVolStr[j], dSigmadk, den, num, impVol.Partial(x, 0));
                }
            }

            // Create dupire outputs.
            PFunction2D.PFunction2D localVol = new PFunction2D.PFunction2D(locVolMat, locVolStr, locVolMatrix);
            localVol.Parse(null);
            string[] names = new string[] { "S0" };
            Vector param = new Vector(1);
            param[0] = Hdataset.S0;
            EstimationResult result = new EstimationResult(names, param);
            result.Objects = new object[3];
            result.Objects[0] = this.r;
            result.Objects[1] = this.q;
            result.Objects[2] = localVol;
            Console.WriteLine("r = " + HCalData.Rate.ToString());
            Console.WriteLine("q = " + HCalData.DividendYield.ToString());
            return result;
        }

        #region IIntegrable implementation
        double IIntegrable.IntegrandFunc(double x)
        {
            return (this.q.Evaluate(x) - this.r.Evaluate(x));
        }
        #endregion

        /// <summary>
        /// This method allows to fit the implied volatility using different models.
        /// </summary>
        /// <param name="Hdataset"></param>
        /// <returns></returns>
        IFunction FitImplVolModel(CallPriceMarketData Hdataset)
        {
            int model = 0;
            switch (model)
            {
                case 0:
                    PFunction2D.PFunction2D pf = new PFunction2D.PFunction2D(Hdataset.Maturity, Hdataset.Strike, Hdataset.Volatility);
                    pf.Interpolation = DVPLUtils.EInterpolationType.LEAST_SQUARES;
                    pf.Extrapolation = DVPLUtils.ExtrapolationType.USEMODEL;
                    pf.Parse(null);
                    return pf;
                case 1:

                    //define a model for fitting the implied volatility
                    Fairmat.Statistics.LinearModel impVol = new Fairmat.Statistics.LinearModel(
                                                            new Fairmat.Statistics.Predictor[] {
                                                            delegate(Vector xx) { return 1; },
                                                            delegate(Vector xx) { return xx[0]; },
                                                            delegate(Vector xx) { return xx[1]; },
                                                            delegate(Vector xx) { return System.Math.Pow(xx[0], 2); },
                                                            delegate(Vector xx) { return System.Math.Pow(xx[1], 2); },
                                                            delegate(Vector xx) { return xx[0] * xx[1]; }, });

                    // Unroll matrix and coordinate vectors in order to make it suitable
                    // for the Quadratic model implementation.
                    int n = Hdataset.Volatility.R * Hdataset.Volatility.C;
                    Matrix xy = new Matrix(n, 2);
                    Vector z = new Vector(n);
                    int count = 0;
                    for (int x = 0; x < Hdataset.Volatility.R; x++)
                    {
                        for (int y = 0; y < Hdataset.Volatility.C; y++)
                        {
                            xy[count, Range.All] = ((Matrix)new Vector() { Hdataset.Strike[x], Hdataset.Volatility[y] }).T;
                            z[count] = Hdataset.Volatility[x, y];
                            count++;
                        }
                    }

                    impVol.Estimate(xy, z);
                    return impVol;
            }

            return null;
        }
    }
}
