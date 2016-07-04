/* Copyright (C) 2009-2014 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Matteo Tesser (matteo.tesser@fairmat.com)
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
using System.Threading.Tasks;
using DVPLI;
using DVPLDOM;
using DVPLI.TaskScheduler;
using Fairmat.MarketData;
using Fairmat.Optimization;
using HestonEstimator;

namespace Heston
{
    /// <summary>
    /// This class declares an additional calibration option for heston
    /// </summary>
    [Mono.Addins.Extension("/Fairmat/Estimator")]
    public class HestonExtendesCalibrationBySimulationEstimator : CallEstimator
    {
        protected override HestonCallOptimizationProblem NewOptimizationProblem(EquityCalibrationData equityCalData, Vector matBound, Vector strikeBound)
        {
            return new HestonCallSimulationOptimizationProblem(equityCalData, matBound, strikeBound);
        }

        public override string Description
        {
            get { return "Calibrate against options (Monte Carlo simulation)"; }
        }

     

    }



    /// <summary>
    /// Calculates and minimize options prices differences by simulating Heston scenarios using Monte Carlo method.
    /// </summary>
    public class HestonCallSimulationOptimizationProblem : HestonCallOptimizationProblem
    {
        /// <summary>
        /// Number of time steps.
        /// </summary>
        int I;
        /// <summary>
        /// Number of Monte Carlo realizations
        /// </summary>
        int N =10000;
        /// <summary>
        /// Normal noise.
        /// </summary>
        Matrix epsilon;
        /// <summary>
        /// Instantaneous growth rate.
        /// </summary>
        Vector growth;
        /// <summary>
        /// Temporal discretization 
        /// </summary>
        double dt = 2*1.0 /52;

        IFunction dyFunc;
        IFunction zrFunc;

        /// <summary>
        /// Wrapper class used to pass data within task-threads.
        /// </summary>
        class SimulationContext
        {
         internal Vector[] paths;
         internal double s0;
         internal double v0;
         internal double kappa;
         internal double theta;
         internal double sigma;
         internal double rho;
         internal double horizon;
        }

        public HestonCallSimulationOptimizationProblem(EquityCalibrationData equityCalData, Vector matBound, Vector strikeBound)
        :base(equityCalData,matBound,strikeBound)
        {
            //this.dyFunc = equityCalData.dyFunc;
            this.dyFunc = CallEstimator.IstantaneousDividendYield(equityCalData);
            this.zrFunc = equityCalData.zrFunc;
            
            // Generates common random numbers
            I = (int)Math.Ceiling(cpmd.Maturity[Range.End] / dt);
            epsilon = new Matrix(I, 2*N);
            NewRandomNumbers();

            // Precalculates istantaneous growth rates
            growth = new Vector(I);
            for (int i = 0; i < I; i++)
            {
                double t = i * dt;
                double zr_t = zrFunc.Evaluate(t);
                growth[i] = zr_t + FunctionHelper.Partial(zr_t,zrFunc, new Vector() { t }, 0, dt) * t - dyFunc.Evaluate(t);
            }

            Console.WriteLine("Heston Simulation Calibration: Paths=" + N);
        }

        /// <summary>
        /// Generates (antithetic)random numbers for both components
        /// </summary>
        private void NewRandomNumbers()
        {
            for (int i = 0; i < I; i++)
            {
                int offset = 0;// component offset
                for (int c = 0; c < 2; c++)
                {
                    for (int s = 0; s < N / 2; s++)//uncomment this line for antithetic
                    //for (int s = 0; s < N; s++)//uncomment this line for non-antithetic
                    {
                        double x=Engine.Generator.Normal();
                        epsilon[i, offset + s] =x;
                        epsilon[i, offset + s + N / 2] = -x;//uncomment this line for antitethic
                    }
                    offset += N;
                }
            }
        }

        /*
        [Obsolete]
        unsafe void SimulateScenarios(double s0, double v0, double kappa, double theta, double sigma, double rho, double horizon, out Vector[] equityScenarios, out Vector[] volScenarios) 
        {
            double rdt = Math.Sqrt(dt);
            // v is the stochastic volatility
            var v = v0 * Vector.Ones(N);
            double* _v = v.Buffer;
            
            double chol22 = Math.Sqrt(1 - rho * rho);

            int I = (int)Math.Ceiling(horizon / dt);
            equityScenarios = new Vector[I + 1];
            equityScenarios[0] = s0 * Vector.Ones(N);

            volScenarios = new Vector[I + 1];
            volScenarios[0] = v0 * Vector.Ones(N);

            for (int i = 1; i <= I; i++)
            {
                double* _s = equityScenarios[i - 1].Buffer;
                equityScenarios[i] = new Vector(N, false);
                double* _sNew = equityScenarios[i].Buffer;

                double t = i * dt;
                double gr = growth[i - 1];
                double* eps = epsilon.GetRow(i - 1);
                for (int p = 0; p < N; p++)
                {
                     double sigmaV2=Math.Max(0,_v[p]);
                     double sigmaV = Math.Sqrt(sigmaV2);//stochastic volatility
                     double z1 = eps[2 * p];
                     double z2 = rho * z1 + chol22 * eps[2 * p + 1];
                  
                    _sNew[p] =_s[p]* Math.Exp((gr  - 0.5 * sigmaV2) * dt +sigmaV * z1 *rdt );
                    _v[p] += kappa * (theta - sigmaV2) * dt + sigma *sigmaV* rdt * z2;
                }
            }
            
            ///////////////////
            GC.KeepAlive(v); //Given that we only a copy of  v.Buffer, v must be kept alive
            ///////////////////
        }
        */


        /// <summary>
        /// Generate Underlying scenarios.
        /// </summary>
        /// <param name="s0"></param>
        /// <param name="v0"></param>
        /// <param name="kappa"></param>
        /// <param name="theta"></param>
        /// <param name="sigma"></param>
        /// <param name="rho"></param>
        /// <param name="horizon"></param>
        /// <returns>Scenarios: the first component is time step, the second one is the underlying path.</returns>
        unsafe Vector[] SimulateScenariosMT(double s0, double v0, double kappa, double theta, double sigma, double rho, double horizon)
        {
            int I = (int)Math.Ceiling(horizon / dt);
            var paths = new Vector[I + 1];
            for (int i = 0; i <= I;i++ )
                paths[i] = new Vector(N, false);

            List<Task> tlist = tsScheduler.AddIntervalJobs(SimulateScenariosCallback, N, new SimulationContext() {paths=paths,s0 = s0, v0 = v0, kappa=kappa,theta=theta,sigma=sigma, rho=rho,horizon=horizon });
            tsScheduler.WaitTaskList(tlist);
            return paths;
        }
        void SimulateScenariosCallback(object context)
        {
            Interval<int, object> pi = (Interval<int, object>)context;
            SimulationContext simData = pi.data as SimulationContext;
            SimulateScenariosCallback(simData.paths, simData.s0, simData.v0, simData.kappa, simData.theta, simData.sigma, simData.rho, simData.horizon, pi.sp, pi.ep);
        }

        unsafe void SimulateScenariosCallback(Vector[] paths, double s0, double v0, double kappa, double theta, double sigma, double rho, double horizon, int sp, int ep)
        {
            double rdt = Math.Sqrt(dt);
            int N=ep-sp;
           

           //v is the stochastic volatility
           var v = v0 * Vector.Ones(N);
           double* _v = v.Buffer;

            double chol22 = Math.Sqrt(1 - rho * rho);

            int I = (int)Math.Ceiling(horizon / dt);
           
            for(int p=sp;p<ep;p++)
                paths[0][p] = s0;
            
            
            for (int i = 1; i <= I; i++)
            {
                double* _s = &paths[i-1].Buffer[0];
                double* _sNew = &paths[i].Buffer[0];
            
                double t = i * dt;
                double gr = growth[i - 1];
                double* eps = epsilon.GetRow(i - 1);
                for (int p = sp; p < ep; p++)
                {
                    double sigmaV2 = Math.Max(0, _v[p-sp]);
                    double sigmaV = Math.Sqrt(sigmaV2);//stochastic volatility
                    double z1 = eps[2 * p];
                    double z2 = rho * z1 + chol22 * eps[2 * p + 1]; 

                    _sNew[p] =_s[p] * Math.Exp((gr - 0.5 * sigmaV2) * dt + sigmaV * z1 * rdt);
                    _v[p-sp] += kappa * (theta - sigmaV2) * dt + sigma * sigmaV * rdt * z2;
                }
            }

            GC.KeepAlive(v); //Given that we only a copy of  v.Buffer, v must be kept alive
            v.Dispose();
            return;
        }




        unsafe double PositivePartMean(Vector input)
        {
            int n = input.Length;
            double* buffer = input.Buffer;
            double sum = 0;
            for (int z = 0; z < n; z++)
                if (buffer[z] > 0) sum += buffer[z];
            GC.KeepAlive(input);// if input is the result of a temporay op.
            input.Dispose();
            return sum / n;
        }
        unsafe double SmoothedPositivePartMean(Vector input)
        {
            double scaling = Math.Min(1, s0 / 2000);

            int n = input.Length;
            double* buffer = input.Buffer;
            double sum = 0;
            for (int z = 0; z < n; z++)
                if (buffer[z] > scaling) sum += buffer[z];
                else
                    if (buffer[z] > -15) 
                            sum += scaling*Math.Exp(buffer[z]-scaling);

            GC.KeepAlive(input);// if input is the result of a temporay op.
            input.Dispose();
            return sum / n;
        }

        double CallPrice(Vector[] paths, double strike, double maturity)
        {
            int z = (int)(maturity / dt);
            return PositivePartMean(paths[z] - strike) * Math.Exp(-this.zrFunc.Evaluate(maturity) * maturity);
        }
        double PutPrice(Vector[] paths, double strike, double maturity)
        {
            int z = (int)(maturity / dt);
            return PositivePartMean(strike - paths[z]) * Math.Exp(-this.zrFunc.Evaluate(maturity) * maturity);
        }

       
        public override double Obj(DVPLI.Vector x)
        {
            double kappa = x[0];
            double theta = x[1];
            double sigma = x[2];
            double rho = x[3];
            double v0 = x[4];

            if (displayObjInfo)
            {
                Console.WriteLine(".");
            }

            var paths = SimulateScenariosMT(s0, v0, kappa, theta, sigma, rho,Math.Min(this.matBound[1],this.cpmd.Maturity[Range.End]));
            
                if(displayObjInfo)
                {
                    Console.WriteLine("Average Underlying Scenario");
                    int Z=paths.Length;// number of time steps
                    int J = paths[0].Length;//number of realizations.
                    //display average path
                    for (int z = 0; z < Z; z++)
                    {
                        double avg = 0;
                        for (int j = 0; j < J; j++)
                            avg += paths[z][j];
                        avg /= J;

                        Console.WriteLine((z * dt) + "\t" + avg);
                    }
                   
                }


                double sum = 0;
                int count = 0;
                

                for (int i = 0; i < this.cpmd.Maturity.Length; i++)
                {
                    double maturity = this.cpmd.Maturity[i];
                    if (maturity < matBound[0]|| maturity>matBound[1])
                        continue;

                    for (int j = 0; j < this.cpmd.Strike.Length; j++)
                    {
                        double strike = this.cpmd.Strike[j];
                        if (strike < strikeBound[0] * s0 || strike > strikeBound[1] * s0)
                            continue;

                        if(calibrateOnCallOptions)
                            if (this.cpmd.CallPrice[i, j] > s0 * optionThreshold)
                                if (cpmd.CallVolume[i, j] > 0)
                                {
                                    double call = CallPrice(paths, strike, maturity);
                                    sum += this.callWeight[i,j] * System.Math.Pow(this.cpmd.CallPrice[i, j] - call, 2.0);
                                    count++;
                                }

                        if(calibrateOnPutOptions)
                            if (cpmd.PutPrice != null)
                                if (cpmd.PutPrice[i, j] > s0 * optionThreshold)
                                    if (cpmd.PutVolume[i, j] > 0)
                                    {
                                        double put = PutPrice(paths, strike, maturity);

                                        sum += this.putWeight[i, j] * System.Math.Pow(cpmd.PutPrice[i, j] - put, 2.0);
                                        count++;
                                    }
                    }
                }
                double p = 0;

                //Frees memory
                for (int i = 0; i < paths.Length; i++)
                    paths[i].Dispose();
                

                if (this.useFellerPenalty)
                    p += this.FellerPenalty(x);

                if (double.IsNaN(sum) || double.IsInfinity(sum))
                    return p+10e5 * x.Norm();
                return p + Math.Sqrt(sum / this.totalVolume)/s0*100;
        }

 
    }
}
