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
using Fairmat.Finance;

namespace Dupire
{
    [Mono.Addins.Extension("/Fairmat/Estimator")]
    public partial class DupireEstimator : IEstimatorEx, IIntegrable
    {
        private PFunction r;
        private PFunction q;
        private DupireCalibrationSettings calibrationSettings;

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
        /// <param name="query">The parameter is not used.</param>
        /// <returns>
        /// An array containing the types InterestRateMarketData and CallPriceMarketData.
        /// </returns>
        public EstimateRequirement[] GetRequirements(IEstimationSettings settings, EstimateQuery query)
        {
            return new EstimateRequirement[] { new EstimateRequirement(typeof(DiscountingCurveMarketData),MarketRequirement.TickerMarket), 
                                               new EstimateRequirement(typeof(CallPriceMarketData)) };
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

        public EstimationResult Estimate(System.Collections.Generic.List<object> marketData, IEstimationSettings settings = null, IController controller = null, System.Collections.Generic.Dictionary<string, object> properties = null)
        {
            CurveMarketData discountingCurve = (CurveMarketData)marketData[0];
            CallPriceMarketData Hdataset = (CallPriceMarketData)marketData[1];

            //gets the settings
            calibrationSettings = settings as DupireCalibrationSettings;

            
            //return this.FairmatEstimate(discountingCurve, Hdataset);
            // Removed quantlib estimate, it does not work correctly
            switch (calibrationSettings.LocalVolatilityCalculation)
            {
                case LocalVolatilityCalculation.Method1:
                    return this.FairmatEstimate(discountingCurve, Hdataset);
                case LocalVolatilityCalculation.QuantLib:
                    return QuantLibEstimate(discountingCurve, Hdataset);
                default:
                    throw new NotImplementedException("Method not implemented");
            }
             
        }

        #endregion

        private EstimationResult FairmatEstimate(CurveMarketData discountingCurve, CallPriceMarketData Hdataset)
        {
            EquityCalibrationData HCalData = new EquityCalibrationData(Hdataset, discountingCurve);
            //HCalData.Setup(Hdataset, discountingCurve);
          
            bool hasArbitrage = HCalData.HasArbitrageOpportunity(10e-2);
            if (hasArbitrage)
                Console.WriteLine("Market data contains arbitrage opportunity");
            
            this.r = new DVPLDOM.PFunction(discountingCurve.Durations,discountingCurve.Values);
            this.q = HCalData.dyFunc as PFunction;
            
            //this.q.Expr = (double[,])ArrayHelper.Concat(HCalData.MaturityDY.ToArray(), HCalData.DividendYield.ToArray());
            this.r.Parse(null);
            this.q.Parse(null);

            Vector locVolMat, locVolStr;
            //IFunction fittedSurface = FitImplVolModel(Hdataset);
            //Matrix locVolMatrix = LocVolMatrixFromImpliedVol(Hdataset, fittedSurface, out locVolMat, out locVolStr);
            CallPriceSurface fittedSurface = CallPriceSurface.NoArbitrageSurface(HCalData);
            Matrix locVolMatrix = LocVolMatrixFromCallPrices(Hdataset, fittedSurface, out locVolMat, out locVolStr);
            Console.WriteLine(locVolMatrix);


            // Create dupire outputs.
            PFunction2D.PFunction2D localVol = new PFunction2D.PFunction2D(locVolMat, locVolStr, locVolMatrix);
            localVol.Parse(null);
            string[] names = new string[] { "S0" };
            Vector param = new Vector(1);
            param[0] = Hdataset.S0;
            EstimationResult result = new EstimationResult(names, param);
            //result.Objects = new object[3];
            result.Objects = new object[4];
            result.Objects[0] = this.r;
            result.Objects[1] = this.q;
            result.Objects[2] = localVol;
            result.Objects[3] = fittedSurface;

            //Console.WriteLine("r = " + HCalData.Rate.ToString());
            //Console.WriteLine("q = " + HCalData.DividendYield.ToString());
            return result;
        }

        private Matrix LocVolMatrixFromImpliedVol(CallPriceMarketData Hdataset, IFunction impVol, out Vector locVolMat, out Vector locVolStr)
        {

            int nmat = calibrationSettings.LocalVolatilityMaturities;
            int nstrike = calibrationSettings.LocalVolatilityStrikes;
            double lastMat = Hdataset.Maturity[Range.End];
            double lastStr = Hdataset.Strike[Range.End];
            locVolMat = Vector.Linspace(Hdataset.Maturity[0], lastMat, nmat);
            locVolStr = Vector.Linspace(Hdataset.Strike[0], lastStr, nstrike);
            Matrix locVolMatrix = new Matrix(nmat, nstrike);
           
            Integrate integrate = new Integrate(this);
            double sigma, dSigmadk, num, y, den, avgGrowthRate;
            Vector x = new Vector(2);
            for (int i = 0; i < nmat; i++)
            {
                avgGrowthRate = integrate.AdaptLobatto(0.0, locVolMat[i]);

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
                        (impVol.Partial(x, 0) + avgGrowthRate * x[1] * dSigmadk);
                    y = Math.Log(locVolStr[j] / Hdataset.S0) + avgGrowthRate;
                    den = System.Math.Pow(1.0 - x[1] * y * dSigmadk / sigma, 2) + x[1] * sigma * x[0] *
                        (dSigmadk - 0.25 * x[1] * sigma * x[0] * dSigmadk * dSigmadk + x[1] * impVol.Partial2(x, 1));
                    locVolMatrix[i, j] = Math.Sqrt(num / den);
                }
            }
            return locVolMatrix;
        }

        private Matrix LocVolMatrixFromCallPrices(CallPriceMarketData Hdataset, CallPriceSurface CallPrice, out Vector locVolMat, out Vector locVolStr)
        {
            Integrate integrate = new Integrate(this);  

            int nmat = calibrationSettings.LocalVolatilityMaturities;
            int nstrike = calibrationSettings.LocalVolatilityStrikes;

            double firstMat = CallPrice.MinMaturity;
            double lastMat = CallPrice.MaxMaturity;
           
            double firstStr =  CallPrice.MinStrike;
            double lastStr = CallPrice.MaxStrike;
            double delta = (lastStr - firstStr) / nstrike;
            locVolMat = Vector.Linspace(firstMat, lastMat, nmat);
            locVolStr = Vector.Linspace(firstStr + delta, lastStr - delta, nstrike);
            Matrix locVolMatrix = new Matrix(nmat, nstrike);

            // this next matrix is created only for debugging pourpose
            Matrix squaredLocVolMatrix = new Matrix(nmat, nstrike);
            double num, den, call, dCdt, dCdk, d2Cdk2;
            Vector x = new Vector(2);
            double h0 = 0.02*Hdataset.S0;//increment for numerical derivatives (stock)
            double ht = 0.02;// 0.25 * (lastMat - firstMat) / nmat;//increment for numerical derivatives (maturities)
            double hs = h0;
            double d2Threshold=10e-5;
            for (int i = 0; i < nmat; i++)
            { 
                double growthRate = integrate.AdaptLobatto(0.0, locVolMat[i]);
                x[0] = locVolMat[i];
                for (int j = 0; j < nstrike; j++)
                {
                    x[1] = locVolStr[j];

                    var support=CallPrice.SupportY(x[0]);
                    //if (x[1] > support[1] && j > 1)
                    if(x[1]<support[0]||x[1]>support[1])
                    {
                        //skip...fill later
                        //squaredLocVolMatrix[i, j] = squaredLocVolMatrix[i, j - 1];
                        //locVolMatrix[i, j] = locVolMatrix[i, j-1];
                    }
                    else
                    {
                        bool fail = false;
                        call = CallPrice.Evaluate(x);
                        dCdt = CallPrice.Partial(x, 0, ht);
                        //do
                       // {

                            dCdk = CallPrice.Partial(x, 1, hs);
                            d2Cdk2 = CallPrice.Partial2(x, 1, hs);
                            //if (Math.Abs(d2Cdk2) > d2Threshold || hs > 0.1 * Hdataset.S0)
                            if (Math.Abs(d2Cdk2) < d2Threshold)
                            {
                                fail = true;
                                //break;
                            }
                        //    hs *= 2;
                        //} while (true);

                        if (!fail)
                        {
                            double qq = this.q.Evaluate(x[0]);
                            num = dCdt + growthRate * x[1] * dCdk + qq * call;
                            den = x[1] * x[1] * d2Cdk2;
                            squaredLocVolMatrix[i, j] = 2.0 * num / den;
                            locVolMatrix[i, j] = Math.Sqrt(Math.Abs(2.0 * num / den));
                        }
                    }
                }
            }

            //Fill missing vol
            for (int i = 0; i < nmat; i++)
            {
                double lastVol=0;
                for (int j = 0; j < nstrike; j++)
                {
                    if (locVolMatrix[i, j] != 0)
                        lastVol = locVolMatrix[i, j];
                    else
                        locVolMatrix[i, j] = lastVol;
                }
                lastVol = 0;
                for (int j = nstrike-1; j >=0; j--)
                {
                    if (locVolMatrix[i, j] != 0)
                        lastVol = locVolMatrix[i, j];
                    else
                        locVolMatrix[i, j] = lastVol;
                }
                //

            }


            return locVolMatrix;
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
            int model = 1;
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
                            if (Hdataset.Volatility[x, y] > 0.01)
                            {
                                xy[count, Range.All] = (new Vector() { Hdataset.Maturity[x],Hdataset.Strike[y] }).T;
                                z[count] = Hdataset.Volatility[x, y];
                                count++;
                            }
                        }
                    }
                    xy.Resize(count, xy.C);
                    z.Resize(count);
                    impVol.Estimate(xy, z);
                    return impVol;
            }

            return null;
        }
    }
}
