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

/*
 * Notes
 * http://www.zeliade.com/whitepapers/zwp-0004.pdf
 */


using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using DVPLI;
using DVPLDOM;
using DVPLI.TaskScheduler;
using Fairmat.MarketData;
using Fairmat.Optimization;

namespace HestonEstimator
{
    /// <summary>
    /// Describes the time-dependent Heston optimization problem.
    /// </summary>
    public class HestonCallOptimizationProblem : IOptimizationProblem
    {
        protected CallPriceMarketData cpmd;
        /// <summary>
        /// The call market price matrix.
        /// </summary>
        internal Matrix callMarketPrice;

        /// <summary>
        /// The maturity vector relative to the callMarketPrice matrix.
        /// </summary>
        private Vector maturity;

        /// <summary>
        /// The strike vector relative to callMarketPrice matrix.
        /// </summary>
        protected Vector strike;

        /// <summary>
        /// The average rate vector up to a given maturity 
        /// </summary>
        protected Vector rate;

        /// <summary>
        /// The avereage dividend yield vector up to a given maturity
        /// </summary>
        protected Vector dividendYield;

        /// <summary>
        /// Averegate drift up to a given maturity
        /// </summary>
        //private Vector drift;

        /// <summary>
        /// Establish whether to use the boundary penalty function or not.
        /// </summary>
        public bool useBoundPenalty = false;

        /// <summary>
        /// Establish whether to use the Feller penalty function or not. 
        /// Note: it's affect a lot the multi-level single linkage algorithm performance.
        /// </summary>
        public bool useFellerPenalty = false;

        /// <summary>
        /// Weights to be used in the calibration.
        /// </summary>
        protected Matrix callWeight;
        /// <summary>
        /// Weights to be used in the calibration.
        /// </summary>
        protected Matrix putWeight;



        /// <summary>
        /// Builds objective function on relative pricing error.
        /// </summary>
        static bool optimizeRelativeError = false;
        static double pricingMin = 0.0001;
        static internal bool displayObjInfo = false;
        static internal double optionThreshold = 0.000;
        static internal int weighting = 2;  //Option weighting 0=constant w. //1=linear //2 log
        internal static bool calibrateOnCallOptions = true;
        internal static bool calibrateOnPutOptions = true;

        /// <summary>
        /// Small value used in the boundary penalty function.
        /// </summary>
        protected double smallValue = 1e-4;

        /// <summary>
        /// Value that weights the boundary penalty function.
        /// </summary>
        protected double k1 = 1e10;

        /// <summary>
        /// Value that weights the Feller inequality penalty function.
        /// </summary>
        protected double k2 = 1e2;

        /// <summary>
        /// The number of call option on which calibration is performed.
        /// </summary>
        protected internal int numCall;
        protected internal int numPut;
        /// <summary>
        /// Total volume for calls and puts.
        /// </summary>
        protected internal double totalVolume;
        /// <summary>
        /// Process starting value.
        /// </summary>
        protected double s0;


        protected Vector matBound;
        protected Vector strikeBound;

        /// <summary>
        /// Initializes a new instance of the HestonCallOptimizationProblem class using the
        /// EquityCalibrationData data structure.
        /// </summary>
        /// <param name="equityCalData">
        /// An EquityCalibrationData object containing market data for calibration.
        /// </param>
        /// <param name="matBound">
        /// A vector containing the minimum and maximum values
        /// for maturities to be used in calibration.
        /// </param>
        /// <param name="strikeBound">
        /// A vector containing the minimum and maximum values
        /// for strikes to be used in calibration.
        /// </param>
        public HestonCallOptimizationProblem(EquityCalibrationData equityCalData, Vector matBound, Vector strikeBound)
        {
            this.cpmd = equityCalData.Hdata;
            this.matBound = matBound;
            this.strikeBound = strikeBound;
            SetVariables(equityCalData.Hdata.CallPrice, equityCalData.Hdata.Maturity,
                         equityCalData.Hdata.Strike, equityCalData.CallMatrixRiskFreeRate,
                         equityCalData.CallMatrixDividendYield, equityCalData.Hdata.S0);

           

            displayObjInfo = false;
        }

        /// <summary>
        /// Initializes a new instance of the HestonCallOptimizationProblem class.
        /// </summary>
        /// <param name="callMarketPrice">A matrix containing call option market prices.</param>
        /// <param name="maturity">
        /// Vector of call option maturities relative to callMarketPrice matrix.
        /// </param>
        /// <param name="strike">
        /// Vector of call option strikes relative to callMarketPrice matrix.
        /// </param>
        /// <param name="rate">
        /// Vector of zero coupon bond rates calculated relative to maturity vector maturities.
        /// </param>
        /// <param name="dividendYield">
        /// Vector of dividend yield rates calculated relative to maturity vector maturities.
        /// </param>
        /// <param name="s0">Index/Equity value at the time of calibration.</param>
        /// <param name="matBound">
        /// A vector containing the minimum and maximum values
        /// for maturities to be used in calibration.
        /// </param>
        /// <param name="strikeBound">
        /// A vector containing the minimum and maximum values
        /// for strikes to be used in calibration.
        /// </param>
        [Obsolete]
        HestonCallOptimizationProblem(Matrix callMarketPrice, Vector maturity, Vector strike, Vector rate, Vector dividendYield, double s0, Vector matBound, Vector strikeBound)
        {
            this.matBound=matBound;
            this.strikeBound = strikeBound;
            SetVariables(callMarketPrice, maturity, strike, rate,
                         dividendYield, s0);
        }

        /// <summary>
        /// Sets several variables used to solve the optimization problem.
        /// </summary>
        /// <param name="callMarketPrice">A matrix containing call option market prices.</param>
        /// <param name="maturity">
        /// Vector of call option maturities relative to callMarketPrice matrix.
        /// </param>
        /// <param name="strike">
        /// Vector of call option strikes relative to callMarketPrice matrix.
        /// </param>
        /// <param name="rate">
        /// Vector of zero coupon bond rates calculated relative to maturity vector maturities.
        /// </param>
        /// <param name="dividendYield">
        /// Vector of dividend yield rates calculated relative to maturity vector maturities.
        /// </param>
        /// <param name="s0">Index/Equity value at the time of calibration.</param>
        /// <param name="matBound">
        /// A vector containing the minimum and maximum values
        /// for maturities to be used in calibration.
        /// </param>
        /// <param name="strikeBound">
        /// A vector containing the minimum and maximum values
        /// for strikes to be used in calibration.</param>
        private void SetVariables(Matrix callMarketPrice, Vector maturity, Vector strike, Vector rate, Vector dividendYield, double s0)
        {
            this.s0 = s0;

            //var rf = new PFunction(maturity, rate);
            var dy = new PFunction(maturity, dividendYield);
            //var rfF = new Fairmat.Math.Integrate(x => rf.Evaluate(x));
            var dyF = new Fairmat.Math.Integrate(x => dy.Evaluate(x));

            this.rate = new Vector(maturity.Length);
            this.dividendYield = new Vector(maturity.Length);
            for (int z = 0; z < maturity.Length; z++)
            {
                //this.rate[z] = rfF.AdaptLobatto(0, maturity[z]) / maturity[z];
                this.dividendYield[z] = dyF.AdaptLobatto(0, maturity[z]) / maturity[z];
            }
            this.rate = rate;

            this.maturity = maturity;
            //this.drift = this.rate - this.dividendYield;
            this.strike = strike;
            this.callMarketPrice = callMarketPrice;
            this.numCall = 0;

            callWeight = new Matrix(this.callMarketPrice.R, this.callMarketPrice.C);
            putWeight = new Matrix(this.callMarketPrice.R, this.callMarketPrice.C);

            for (int r = 0; r < this.callMarketPrice.R; r++)
            {
                if (this.maturity[r] >= matBound[0] && this.maturity[r]<= matBound[1])
                {
                    for (int c = 0; c < this.callMarketPrice.C; c++)
                    {
                        if (this.strike[c] >= s0 * strikeBound[0] && this.strike[c] <= s0 * strikeBound[1])
                        {
                            if (calibrateOnCallOptions)
                                if (this.callMarketPrice[r, c] > s0 * optionThreshold && this.cpmd.CallVolume[r, c] > 0)
                                {
                                    this.callWeight[r, c] = CalculateWeight(this.cpmd.CallVolume[r, c]);
                                    this.numCall++;
                                    this.totalVolume += CalculateWeight(this.cpmd.CallVolume[r, c]) ;
                                }
                            if (calibrateOnPutOptions)
                                if (this.cpmd.PutPrice != null && this.cpmd.PutPrice[r, c] > s0 * optionThreshold && this.cpmd.PutVolume[r, c] > 0)
                                {
                                    this.putWeight[r, c] = CalculateWeight(this.cpmd.PutVolume[r, c]);

                                    this.numPut++;
                                    this.totalVolume += CalculateWeight(this.cpmd.PutVolume[r, c]);
                                }
                            
                        }
                    }
                }
            }

            //calibrate minVolatility: actually in this.cpmd.Volatility there is sigma not sigma^2
            if (this.cpmd.Volatility != null)
            {
                //Rows maturities, columns strikes
                vLastMin = 0.5 * Math.Pow(this.cpmd.Volatility.Min().Min(), 2);

                v0Min = 0.99 * Math.Pow(this.cpmd.Volatility[0, Range.All].Min().Min(), 2);
                v0Max = 1.01 * Math.Pow(this.cpmd.Volatility[0, Range.All].Max().Max(), 2);
            }


            Console.WriteLine("Options weighting:\t" + weighting);
            Console.WriteLine("OptionThreshold:\t" + optionThreshold);
            Console.WriteLine("Strike Bounds:\t" + strikeBound);
            Console.WriteLine("Maturity Bounds:\t" + matBound);
            Console.WriteLine("Lb:\t" + Bounds.Lb);
            Console.WriteLine("Ub:\t" + Bounds.Ub);
            if(Engine.Verbose>=2)
                PutCallTest();
        }


        /// <summary>
        /// Finds the couple of integer {i,j} so that the values
        /// vector[i], vector[i+1],..., vector[j-1], vector[j]
        /// are all included in the interval (bound[0], bound[1])
        /// while vector[i-1] and vector[j+1] are not.
        /// </summary>
        /// <param name='vector'>
        /// Vector to be filtered.
        /// </param>
        /// <param name='bounds'>
        /// Bounds to be used for the filtering.
        /// </param>
        /// <returns>
        /// Integer array representing the index couple.
        /// </returns>
        private int[] FindExtremes(Vector vector, Vector bounds)
        {
            int[] result = new int[2];
            int i = 0;

            while (vector[i] < bounds[0])
                i++;
            result[0] = i;
            i = vector.Length - 1;
            while (vector[i] > bounds[1])
                i--;
            result[1] = i;
            return result;
        }

        #region IOptimizationProblem Members

        /// <summary>
        /// If volatility is observed, it may be useful to use this information
        /// </summary>
        protected double v0Min=0.001;
        protected double v0Max = 1;
        protected double vLastMin = 0.001;

        /// <summary>
        /// Gets the bounds for the optimization.
        /// </summary>
        public Bounds Bounds
        {
            get
            {
                Bounds bounds = new Bounds();

                // The order of parameters is k, theta, sigma, rho,  V0 (and optionally div.y)
                if (HestonConstantDriftEstimator.impliedDividends)
                {
                    bounds.Lb = (Vector)new double[] { 0, v0Min, 0.001, -1, vLastMin, 0 };
                    bounds.Ub = (Vector)new double[] { 15, v0Max, 2, 1, 1, .2 };
                }
                else
                {
                    bounds.Lb = (Vector)new double[] { 0, v0Min, 0.001, -1, vLastMin };
                    bounds.Ub = (Vector)new double[] { 15,v0Max, 2, 1, 1 };
                }
                return bounds;
            }
        }

        /// <summary>
        /// This method is unused but part of the interface.
        /// </summary>
        /// <param name="x">The parameter is not used.</param>
        /// <returns>Null as it's unused.</returns>
        public DVPLI.Vector G(DVPLI.Vector x)
        {
            return null;
        }

        /// <summary>
        /// This method is unused but part of the interface.
        /// </summary>
        /// <param name="x">The parameter is not used.</param>
        /// <returns>Nothing the function always throws a NotImplementedException.</returns>
        public DVPLI.Vector Grad(DVPLI.Vector x)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Gets a value indicating whether there are non linear constraints in this
        /// optimization problem. In this case there are not.
        /// </summary>
        public bool HasNonLinearConstraints
        {
            get
            {
                return false;
            }
        }

        /// <summary>
        /// Gets null as we have no linear constrains defined.
        /// </summary>
        public virtual LinearConstraints LinearIneqConstraints
        {
            get
            {
                return null;
            }
        }

        /// <summary>
        /// Calibration objective function.
        /// </summary>
        /// <param name='x'>
        /// The vector of parameters.
        /// </param>
        /// <returns>
        /// Objective function value.
        /// </returns>
        public virtual double Obj(DVPLI.Vector x)
        {
            double sum = 0;
            if(Engine.MultiThread && !displayObjInfo)
            {
                // Instantiate parallel computation if enabled.
                List<Task> tl = new List<Task>();

                // Context contains both input parameters and outputs.
                List<HestonCall> context = new List<HestonCall>();
                for (int r = 0; r < this.callMarketPrice.R; r++)
                {
                    if (this.maturity[r] >= this.matBound[0])
                    {
                        HestonCall hc = new HestonCall(this, x, this.s0);
                        context.Add(hc);

                        hc.T = this.maturity[r];
                        hc.rate = this.rate[r];
                        if (HestonConstantDriftEstimator.impliedDividends)
                            hc.dividend = x[Range.End];
                        else
                            hc.dividend = this.dividendYield[r];
                        hc.row = r;
                        tl.Add(Task.Factory.StartNew(this.CalculateSingleRow, hc));
                    }
                }

                tsScheduler.WaitTaskList(tl);
                for (int r = 0; r < tl.Count; r++)
                        sum += context[r].sum;
            }
            else
            {
                // Sequential version of the code, used when parallel computation is disabled.
                HestonCall hc = new HestonCall(this, x, this.s0);
                for (int r = 0; r < this.callMarketPrice.R; r++)
                {
                    if (this.maturity[r] >= this.matBound[0])
                    {
                        hc.T = this.maturity[r];
                        hc.rate = this.rate[r];
                        if (HestonConstantDriftEstimator.impliedDividends)
                            hc.dividend = x[Range.End];
                        else
                            hc.dividend = this.dividendYield[r];
                        hc.row = r;
                        hc.sum = 0;
                        this.CalculateSingleRow(hc);
                        sum += hc.sum;
                    }
                }

                var pricingErrors = hc.hestonCallPrice - this.callMarketPrice;
                if (displayObjInfo)
                {
                    avgPricingError = 0;
                    for (int r = 0; r < this.callMarketPrice.R; r++)
                        if (this.maturity[r] >= this.matBound[0])
                        {
                            for (int c = 0; c < this.callMarketPrice.C; c++)
                            {
                                if (this.callMarketPrice[r, c] > s0 * optionThreshold && this.cpmd.CallVolume[r, c] > 0)
                                    avgPricingError += Math.Abs(pricingErrors[r, c]);
                                //if (this.cpmd.PutPrice[r, c] > s0 * optionThreshold && this.cpmd.PutVolume[r, c] > 0)
                                //    avgPricingError += Math.Abs(pricingErrors[r, c]);
                            }
                        }
                        avgPricingError /= numCall;

                    int RR = Math.Min(12, this.callMarketPrice.R - 1);
                    int CC = Math.Min(12, this.callMarketPrice.C - 1);
                    Console.WriteLine("Mkt Price");
                    Console.WriteLine(this.callMarketPrice[Range.New(0,RR),Range.New(0,CC)]);
                    Console.WriteLine("Pricing Errors");
                    Console.WriteLine(pricingErrors[Range.New(0, RR), Range.New(0, CC)]);
                }   
                objCount++;
            }

            //Calculate average distance...
            sum = Math.Sqrt( sum /this.totalVolume);
            if (this.useBoundPenalty)
                sum += this.BoundPenalty(x);

            if (this.useFellerPenalty)
                sum +=this.FellerPenalty(x);
         
            return sum;
        }
        static int objCount = 0;
        
        //Average pricing error
        internal static double avgPricingError = 0;

        /// <summary>
        /// Calculates a single row of the objective function. Basically
        /// options with the same maturity and different strikes.
        /// </summary>
        /// <param name='context'>
        /// An object of type <see cref="HestonCall"/> containing the context.
        /// </param>
        private void CalculateSingleRow(object context)
        {
            CalculateSingleRowWithInterpolation(context);
            return;

            HestonCall hc = context as HestonCall;
          
            int r = hc.row;
            for (int c = 0; c < this.callMarketPrice.C; c++)
            {
                bool callCondition=this.callMarketPrice[r, c] > s0*optionThreshold && this.cpmd.CallVolume[r,c]>0;
                bool putCondition=this.cpmd.PutPrice[r, c] > s0*optionThreshold && this.cpmd.PutVolume[r,c]>0;
                if (callCondition)//||putCondition)
                {
                    hc.K = this.strike[c];
                    var callPut=hc.HestonCallPutPrice();
                    hc.hestonCallPrice[r, c] = callPut[0];
                    hc.hestonPutPrice[r, c] = callPut[1];

                    if (callCondition)
                    {
                        if (HestonCallOptimizationProblem.optimizeRelativeError)
                        {
                            double mkt = pricingMin + this.callMarketPrice[r, c];
                            double model = pricingMin + hc.hestonCallPrice[r, c];
                            hc.sum += callWeight[r,c] * Math.Pow((model - mkt) / mkt, 2);
                        }
                        else
                        {
                            hc.sum += callWeight[r,c] * Math.Pow(hc.hestonCallPrice[r, c] - this.callMarketPrice[r, c], 2);
                        }
                    }
                    /*
                    if (putCondition)
                    {
                        if (HestonCallOptimizationProblem.optimizeRelativeError)
                        {
                            double mkt = pricingMin + this.cpmd.PutPrice[r, c];
                            double model = pricingMin + hc.hestonPutPrice[r, c];
                            hc.sum += putWeight[r,c] * Math.Pow((model - mkt) / mkt, 2);
                        }
                        else
                        {
                            hc.sum += putWeight[r,c]* Math.Pow(hc.hestonPutPrice[r, c] - this.cpmd.PutPrice[r, c], 2);
                        }
                    }
                     */
                }

            }
            return;
        }


        /// <summary>
        /// Calculates call put prices for several strikes using controlled interpolation.
        /// </summary>
        /// <param name="context"></param>
        private void CalculateSingleRowWithInterpolation(object context)
        {
            HestonCall hc = context as HestonCall;
            int r = hc.row;
            hc.sum = 0;

            // Finds upper extreme for call and put
            int max_c =0;
            for (int c = this.callMarketPrice.C-1; c > 0; c--)
            {
                bool callCondition = this.callMarketPrice[r, c] > s0 * optionThreshold && this.cpmd.CallVolume[r, c] > 0;
                bool putCondition = this.cpmd.PutPrice!=null &&  this.cpmd.PutPrice[r, c] > s0 * optionThreshold && this.cpmd.PutVolume[r, c] > 0;
                if (callCondition || putCondition)
                {
                    max_c = c;
                    break;
                }
            }


             var strikes = new List<double>();
             var calls = new List<double>();
             var puts = new List<double>();

            //Evaluates in strategic points
            for (int c = 0; c < this.callMarketPrice.C; c++)
            {
                bool callCondition = this.callMarketPrice[r, c] > s0 * optionThreshold && this.cpmd.CallVolume[r, c] > 0;
                bool putCondition = this.cpmd.PutPrice!=null&& this.cpmd.PutPrice[r, c] > s0 * optionThreshold && this.cpmd.PutVolume[r, c] > 0;
                if (callCondition || putCondition)
                {
                    hc.K = this.strike[c];
                    var callPut = hc.HestonCallPutPrice();
                    strikes.Add(hc.K);
                    calls.Add(callPut[0]);
                    puts.Add(callPut[1]);
                    if (c == max_c)
                        break;

                    c += 1;//skip the subsequent strikes
                    
                    if (c > max_c)
                        c = max_c;
                }
            }
            
            // Builds interpolated call and put values.

            var callFun = new PFunction((Vector)strikes.ToArray(), (Vector)calls.ToArray());
            callFun.m_Function.iType = DVPLUtils.EInterpolationType.SPLINE;
           
            var putFun = new PFunction((Vector)strikes.ToArray(), (Vector)puts.ToArray());
            putFun.m_Function.iType = DVPLUtils.EInterpolationType.SPLINE;
            
            // Evaluates at the requested strikes

            for (int c = 0; c < this.callMarketPrice.C; c++)
            {
                bool callCondition = this.callMarketPrice[r, c] > s0 * optionThreshold && this.cpmd.CallVolume[r, c] > 0;
                bool putCondition = this.cpmd.PutPrice!=null && this.cpmd.PutPrice[r, c] > s0 * optionThreshold && this.cpmd.PutVolume[r, c] > 0;
                  

                if (callCondition)
                {
                    hc.hestonCallPrice[r, c] = callFun.Evaluate(this.strike[c]);
                    if (HestonCallOptimizationProblem.optimizeRelativeError)
                    {
                        double mkt = pricingMin + this.callMarketPrice[r, c];
                        double model = pricingMin + hc.hestonCallPrice[r, c];
                        hc.sum += callWeight[r,c] * Math.Pow((model - mkt) / mkt, 2);
                    }
                    else
                    {
                        hc.sum += callWeight[r, c] * Math.Pow(hc.hestonCallPrice[r, c] - this.callMarketPrice[r, c], 2);
                    }
                }
                
                if (putCondition)
                {
                    hc.hestonPutPrice[r, c] = putFun.Evaluate(this.strike[c]);
                    if (HestonCallOptimizationProblem.optimizeRelativeError)
                    {
                        double mkt = pricingMin + this.cpmd.PutPrice[r, c];
                        double model = pricingMin + hc.hestonPutPrice[r, c];
                        hc.sum +=  putWeight[r, c] * Math.Pow((model - mkt) / mkt, 2);
                    }
                    else
                    {
                        hc.sum += putWeight[r, c] * Math.Pow(hc.hestonPutPrice[r, c] - this.cpmd.PutPrice[r, c], 2);
                    }
                }

            }

            return;
        }


        #endregion

        /// <summary>
        /// Penalty function relative to bounds.
        /// </summary>
        /// <param name='x'>
        /// Vector of parameters.
        /// </param>
        /// <returns>
        /// Penalty value.
        /// </returns>
        private double BoundPenalty(Vector x)
        {
            Vector t1, t2, t3;
            t1 = Bounds.Lb + this.smallValue - x;
            t2 = -Bounds.Ub + this.smallValue + x;

            t3 = t1 * t1 * (t1 > 0) + t2 * t2 * (t2 > 0);
            return this.k1 * t3.Sum();
        }

        /// <summary>
        /// Penalty function in order to satisfy Feller condition:
        /// 2k theta >= sigma^2
        /// </summary>
        /// <param name='x'>
        /// Vector of parameters: x=[k, theta, sigma,...]
        /// </param>
        /// <returns>
        /// Penalty value.
        /// </returns>
        protected double FellerPenalty(Vector x)
        {
            double result;
            result = Math.Max(0, x[2] * x[2] - 2 * x[0] * x[1]);
            return this.s0*this.k2 * result * result;
        }

        /// <summary>
        /// Test method, displays a sensitivity on heston call and put prices.
        /// </summary>
        private void PutCallTest()
        {
            Console.WriteLine("Black-Sholes Calls Market Prices");
            Console.WriteLine(this.callMarketPrice);
            Console.WriteLine("Strikes");
            Console.WriteLine(this.strike);
            Console.WriteLine("Maturities");
            Console.WriteLine(this.maturity);

            var x = new Vector() {3.18344026504981,
                0.0427882999286046,
                0.644527074840708,
                -0.659960749691282,
                0.0150455464938991,
                0.0211747510984717};

            HestonCall hc = new HestonCall(this, x, this.s0);
            hc.T = .1;
            hc.rate = this.rate[0];
            Console.WriteLine("Strike\tCall\tPut");
            for (int z = 200; z < 6500; z += 1000)
            {
                hc.K = z;

                var call = hc.HestonCallPrice();
                var put = hc.HestonPutPrice();
                var callPut = hc.HestonCallPutPrice();
                Console.WriteLine(z + "\t" + callPut[0] + "\t" + callPut[1]);
            }
        }

        /// <summary>
        /// Allows different specification of call/put weighting
        /// </summary>
        /// <param name="volume"></param>
        /// <returns></returns>
        double CalculateWeight(double volume)
        {
            switch (weighting)
            {
                case 0://constant
                    return 1;//no-volume information
                case 1://linear
                    return volume;
                case 2://log
                    if (volume > 0)
                        return Math.Log(volume) + 1;
                    else
                        return 0;
                default:
                    throw new Exception("Not defined");

            }
            
        }

    }
}
