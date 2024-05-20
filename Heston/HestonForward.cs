﻿using System;
using Heston;
using Fairmat.Math;
using DVPLI;

namespace HestonEstimator
{
    // this implementation is based on https://doi.org/10.1142/S0219024909005166 
    // FORWARD START OPTIONS UNDER STOCHASTIC VOLATILITY AND STOCHASTIC INTEREST RATES - REHEZ AHLIP and MAREK RUTKOWSKI
    public class HestonForwardAhlipRutkowski: HestonCall
    {

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonForward"/> class.
        /// </summary>
        public HestonForwardAhlipRutkowski() { }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonForward"/> class.
        /// </summary>
        /// <param name='problem'>
        /// HestonCallOptimizationProblem at which digital price calculations are to be linked.
        /// </param>
        public HestonForwardAhlipRutkowski(HestonCallOptimizationProblem problem) : base(problem) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonForward"/> class.
        /// </summary>
        /// <param name='process'>
        /// HestonProcess at which digital price calculations are to be linked.
        /// </param>
        /// <param name='strike'> 
        /// Strike of the digital option
        /// </param>
        /// <param name="timeToMaturity">
        /// TimeToMaturity of the digital option
        /// </param>
        public HestonForwardAhlipRutkowski(HestonProcess process, double strike, double timeToMaturity) : base(process, strike, timeToMaturity) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonForward"/> class.
        /// </summary>
        /// <param name='process'>
        /// HestonProcess at which digital price calculations are to be linked.
        /// </param>
        public HestonForwardAhlipRutkowski(HestonProcess process) : base(process) { }

        internal struct ModelParameters
        {
            public double kappa;
            public double theta;
            public double sigma;
            public double rho;
            public double v0;
            public double s0;
            public double r;
            public double q;
        }

        internal struct CallParameters
        {
            public double T;
            public double T0;
            public double K;
        }

       
        internal static Complex s1(Complex u, ModelParameters mp, CallParameters cp)
        { 
            var iu = Complex.I * u;
            var onePlusIu = Complex.One + iu;
            var result = -mp.kappa * mp.rho * onePlusIu / mp.sigma;
            result -= 0.5 * (1 - mp.rho * mp.rho) * Complex.Pow(onePlusIu, 2);
            result += onePlusIu / 2;
            return result;
        }

        internal static Complex s2(Complex u, ModelParameters mp)
        {
            var iu = Complex.I * u;
            var onePlusIu = 1.0 + iu;
            return mp.rho * onePlusIu / mp.sigma;
        }

        internal static Complex q1(Complex u, ModelParameters mp)
        {
            var iu = Complex.I * u;
            var preMult = -iu / 2;
            return preMult * (2 * mp.kappa * mp.rho / mp.sigma - iu * (1 - mp.rho * mp.rho) - 1.0);
        }

        internal static Complex q2(Complex u, ModelParameters mp)
        {
            var iu = Complex.I * u;
            return -iu * mp.rho / mp.sigma;
        }


        internal static Complex G1(double tau, Complex lambda, Complex mu, ModelParameters mp, double t)
        {
            var gamma = Complex.Sqrt(mp.kappa * mp.kappa + 2 * mp.sigma * mp.sigma * mu);
            var expGammaTau = Complex.Exp(gamma * tau);
            var numerator = lambda * ((mp.kappa + gamma) + expGammaTau * (gamma - mp.kappa));
            numerator += 2 * mu * (1 - Complex.Exp(gamma * t));

            var denominator = mp.sigma * mp.sigma * lambda * (expGammaTau - 1) + gamma - mp.kappa;
            denominator += expGammaTau * (gamma + mp.kappa);
            
            return numerator / denominator;

        }

        internal static Complex H1(double tau, Complex lambda, Complex mu, ModelParameters mp)
        {
            var gamma = Complex.Sqrt(mp.kappa * mp.kappa + 2 * mp.sigma * mp.sigma * mu);
            var sigma2 = (mp.sigma * mp.sigma);
            var preMult = -2 / sigma2;
            var argumentNumerator = 2 * gamma * Complex.Exp((gamma + mp.kappa) * tau / 2);
            var argumentDenominator = sigma2 * lambda * (Complex.Exp(gamma * tau) - 1) + gamma - mp.kappa + Complex.Exp(gamma * tau) * (gamma + mp.kappa);
            return preMult * Complex.Log(argumentNumerator / argumentDenominator);

        }

        internal static Complex H1Hat( double tau, Complex lambda, ModelParameters mp)
        {
            var sigma2 = (mp.sigma * mp.sigma);
            var preMult = -2 / sigma2;
            var argumentNumerator = 2 * mp.kappa * Math.Exp(2 * mp.kappa * tau);
            var argumentDenominator = sigma2 * lambda * (Math.Exp(mp.kappa * tau) - 1) + 2 * mp.kappa * Math.Exp(mp.kappa * tau);
            return preMult* Complex.Log(argumentNumerator / argumentDenominator);
                

        }
        
        internal static Complex G1Hat(ModelParameters mp,  double tau, Complex lambda)
        {
            var sigma2 = (mp.sigma * mp.sigma);
            var preMult = -2 / sigma2;
            var argumentNumerator = 2 * mp.kappa * Math.Exp(2 * mp.kappa * tau);
            var argumentDenominator = sigma2 * lambda * (Math.Exp(mp.kappa * tau) - 1) + 2 * mp.kappa * Math.Exp(mp.kappa * tau);
            return preMult * Complex.Log(argumentNumerator / argumentDenominator);

        }

        internal static Complex s2Hat(Complex u, double tau0, Complex s1, Complex s2, ModelParameters mp,  double t)
        {
            Complex iu = Complex.I * u;
            var g1 = G1(tau: tau0, lambda: s1, mu: s2, mp: mp, t:t);
            return g1 + (1 + iu) * (mp.rho / mp.sigma);
        }

        internal static Complex q2Hat(Complex u, double tau0, Complex q1, Complex q2, ModelParameters mp,  double t)
        {
            Complex iu = Complex.I * u;
            return -iu * mp.rho / mp.sigma * G1(tau: tau0, lambda: q1, mu: q2, mp:mp, t: t);
        }

        internal static Complex f1Hat(Complex u, ModelParameters mp, CallParameters cp)
        {

            var t = 0;
            var tau0 = cp.T - cp.T0;
            var tau = cp.T0;

            var iu = Complex.I * u;
            var _s1 = s1(u:u, mp: mp, cp:cp);
            var _s2 = s2(u:u, mp:mp);
            var s2H = s2Hat(u:u, tau0: tau0, s1: _s1, s2: _s2, t: t, mp: mp);


            var a = -(1 + iu) * mp.rho * mp.theta * tau0 / mp.sigma;
            var b = -mp.theta * H1(tau: tau0, lambda: _s1, mu: _s2, mp:mp);
            var c = -mp.v0 * G1Hat(tau: tau, lambda: s2H, mp:mp);
            var d = -mp.theta * H1Hat(tau: tau, lambda: s2H, mp:mp);

            return Complex.Exp(a + b + c + d);

        }

        internal static Complex f2Hat(Complex u, ModelParameters mp, CallParameters cp)
        {
            var t = 0;
            var tau0 = cp.T - cp.T0;
            var tau = cp.T0;


            var iu = Complex.I * u;
            var _q1 = q1(u: u, mp:mp);
            var _q2 = q2(u: u, mp:mp);
            var q2H = q2Hat(u: u, tau0: tau0,  q1: _q1, q2: _q2, mp:mp, t: t);


            var a = -iu * mp.rho * mp.theta * tau0 / mp.sigma;
            var b = -mp.theta * H1(tau: tau0, lambda: _q1, mu: _q2, mp:mp);
            var c = -mp.v0 * G1Hat(tau: tau, lambda: q2H, mp:mp);
            var d = -mp.theta * H1Hat(tau: tau, lambda: q2H, mp:mp);

            return Complex.Exp(a + b + c + d);

        }

        internal static double IntegrandP1(double u, ModelParameters mp, CallParameters cp)
        {
            var uComp = new Complex(u);
            var iu = Complex.I * u;
            var logK = Math.Log(cp.K);
            var f1 = f1Hat(u: uComp, mp:mp, cp:cp);
            return (f1 * Complex.Exp(-iu * logK) / iu).Re;


        }
        
        internal static double IntegrandP2(double u, ModelParameters mp, CallParameters cp)
        {
            var uComp = new Complex(u);

            var iu = Complex.I * u;
            var logK = Math.Log(cp.K);
            var f2 = f2Hat(u: uComp, cp:cp, mp:mp);
            return (f2 * Complex.Exp(-iu * logK) / iu).Re;

        }

        internal static double P1Hat( ModelParameters mp, CallParameters cp)
        {
            // compute the integral 
            double a = 1E-8;
            double b = 1000.0;
            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandP1(u: u, cp: cp, mp: mp);
            double part1 = PerformIntegral(a, b, functionToIntegrate);
            double integral = part1 + a * functionToIntegrate(a / 2.0);

            return 0.5 + 1.0 / Math.PI * integral;

        }

        internal static double P2Hat(ModelParameters mp, CallParameters cp)
        {
            // compute the integral 
            double a = 1E-8;
            double b = 1000.0;
            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandP2(u: u, mp: mp, cp: cp);
            double part1 = PerformIntegral(a, b, functionToIntegrate);
            double integral = part1 + a * functionToIntegrate(a / 2.0);

            return 0.5 + 1.0 / Math.PI * integral;

        }


        // based on corollary 6.2 of the paper
        public static double HestonForwardCallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            // adjust the theta value to match their theta definition
            // in ours approach the drift is kappa(theta-v)
            // in their approach the drift is (theta_start-kappa*v) 
            // theta_start-kappa*v = kappa*theta- kappa*v
            // theta_start = kappa * theta 

            double thetaAdjusted = kappa * theta;


            var cp = new CallParameters() { T = T, T0 = T0, K = K };
            var mp = new ModelParameters() { kappa = kappa, theta = thetaAdjusted, sigma = sigma, rho = rho, v0 = v0, s0 = s0, r = r, q = q };



            // use theta adjusted instead of theta 
            var p1 = P1Hat(cp:cp, mp:mp);
            var p2 = P2Hat(cp: cp, mp: mp);
            var B_T_T0 = Math.Exp(-r * (T - T0));
            return s0 * ( p1 - B_T_T0 * K * p2);
        }

        public static double HestonForwardCallPrice(Vector x, double s0, double T, double T0, double K, double r, double q)
        {
            var kappa = x[0];
            var theta = x[1];
            var sigma = x[2];
            var rho = x[3];
            var v0 = x[4];


            return HestonForwardCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                T0:T0, 
                K: K,
                r: r,
                q: q);
        }


    }

    
    // this solution follows the solution of lucic 
    // http://www.untag-smd.ac.id/files/Perpustakaan_Digital_1/FINANCE%20The%20Best%20of%20Wilmott%20Vol.1%20Incorporating%20the%20Quantitative%20Finance%20Review.pdf#page=426
    public class HestonForwardLucic
    {
        // we need to implment the formulas for V2 in the lucic paper 
        static int m = 2;

        // in lucic paper we have 
        // lambda * (vhat - v )
        // in our notation we have
        //kappa(theta-v)

        // lambda = kappa
        // vhat = theta 
        // eta = sigma 

        public static double u_j(double j)
        {
            return 1 / 2 - j;
        }


        public static double b(double lambda, double rho, double eta, double j, double T, double T0, double tau)
        {
            return lambda + rho * eta * (m - 2 + Theta_tau(tau, T0:T0, T:T) * (2 - m - j));
        }

        public static double Theta_tau(double x, double T, double T0)
        {
            var tau = T - T0;
            if (x >= 0.0 & x <= tau)
                return 1.0;
            return 0.0;
        }

        public static Complex Beta(double lambda, double rho, double eta, double j, double T, double T0, double tau, double kappa)
        {
            var theta = Theta_tau(tau, T0: T0, T: T);
            var result = -Complex.I * rho * eta * kappa * theta;
            result += b(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, tau: tau);
            return result;

        }

        public static Complex Alpha(double lambda, double rho, double eta, double j, double T, double T0, double tau, double kappa)
        {
            var i = Complex.I;
            var theta = Theta_tau(tau, T0: T0, T: T);
            return theta * (-kappa * kappa / 2 - i * kappa / 2 + i * kappa * kappa);
        }


        public static Complex d( double lambda, double rho, double eta, double j, double T, double T0, double kappa)
        {
            var eta2 = eta * eta;
            var beta = Beta(tau: 0.0, lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0,  kappa: kappa);
            var arg = beta * beta - 2 * eta2 * Alpha(tau:0.0, lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);


            return Complex.Sqrt(arg);
        }

        public static Complex rPlus(double lambda, double rho, double eta, double j, double T, double T0, double kappa)
        {
            var beta = Beta(tau: 0.0, lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var _d = d(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            return (beta + _d) / (eta * eta);
        }

        public static Complex rMinus(double lambda, double rho, double eta, double j, double T, double T0,  double kappa)
        {
            var beta = Beta(tau: 0.0, lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var _d = d(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            return (beta - _d) / (eta * eta);
        }

        public static Complex g(double lambda, double rho, double eta, double j, double T, double T0, double kappa)
        {
            var _rMinus = rMinus(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var _rPlus = rPlus(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            return _rMinus / _rPlus;
        }


    }


    public class HestonForwardApproximated
    {
        public static double HestonForwardCallPrice(Vector x, double s0, double T, double T0, double K, double r, double q)
        {
            var kappa = x[0];
            var theta = x[1];
            var sigma = x[2];
            var rho = x[3];
            var v0 = x[4];


            return HestonForwardCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: r,
                q: q);
        }
        public static double ExpectationCIRProcess(double x0, double a, double b, double t)
        {
            return x0 * Math.Exp(-a * t) + b * (1 - Math.Exp(-a * t));
        }

        public static double derivativeExpectationCIRProcess(double x0, double a, double b, double t)
        {
            return Math.Exp(-a * t);
        }

        // we use similar formulas that are used for pricing Forward call options under Black Scholes model
        public static double HestonForwardCallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{r}\t{q}");

            }

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;
            var c = HestonCall.HestonCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q);
            Engine.Verbose = verbosity;

            var price = s0 * Math.Exp(-q * T0) * c;
            
            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Call Option is {price}");
            }
            
            return price;
        }

        public static double HestonForwardPutPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{r}\t{q}");

            }

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;
            var p = HestonCall.HestonPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q);
            Engine.Verbose = verbosity;

            var price = s0 * Math.Exp(-q * T0) * p;

            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Put Option is {price}");
            }

            return price;
        }

        public static double HestonForwardDigitalPutPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {

            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS Digital Put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{r}\t{q}");

            }
            
            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;
            var dput = HestonDigital.HestonDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q);
            Engine.Verbose = verbosity;

            var digitalPrice = HestonDigital.DiscountFactor(r, T0) * dput;

            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Digital Put Option is {digitalPrice}");
            }

            return digitalPrice;
        }

        public static double HestonForwardDigitalCallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS Digital Call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{r}\t{q}");

            }

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;
            var dcall = HestonDigital.HestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q);
            Engine.Verbose = verbosity;

            var digitalPrice = HestonDigital.DiscountFactor(r, T0) * dcall;
            
            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Digital Call Option is {digitalPrice}");
            }
            
            return digitalPrice;
        }

        #region Greeks
        public static GreeksDerivatives HestonForwardCallWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {


            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price and greeks of a FS Call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{r}\t{q}");
            }

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;

            var c = HestonCall.HestonCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q);

            
            var fscallMarkToMarket =  s0 * Math.Exp(-q * T0) * c;

            var fsCallDelta = Math.Exp(-q * T0) * c;

            var fsCallGamma = 0.0;

            var fsCallTheta = s0 * Math.Exp(-q * T0) * HestonNumericalGreeks.ThetaCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            var fsCallVega = s0 * Math.Exp(-q * T0) * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            var rhoGreek = s0 * Math.Exp(-q * T0) * HestonRho.RhoCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r, q: q);

            Engine.Verbose = verbosity;


            var result = new GreeksDerivatives()
            {
                Deltas = (Vector)fsCallDelta,
                Gammas = (Vector)fsCallGamma,
                Theta = fsCallTheta,
                Vegas = (Vector)fsCallVega,
                MarkToMarket = fscallMarkToMarket, 
                Rho = rhoGreek
            };

            if (verbosity > 0)
            {
                Console.WriteLine("Price and Greeks of the Forward Starting Call Option");
                Console.WriteLine("Delta\tGamma\tTheta\tVega\tRho\tPrice");
                Console.WriteLine($"{result.Deltas[0]}\t{result.Gammas[0]}\t{result.Theta}\t{result.Vegas[0]}\t{result.Rho}\t{result.MarkToMarket}");
            }


            return result;

        }

        public static GreeksDerivatives HestonForwardPutWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price and greeks of a FS Put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{r}\t{q}");
            }

            Engine.Verbose = 0;
            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            var put = HestonCall.HestonPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q);


            var fsPutMarkToMarket = s0 * Math.Exp(-q * T0) * put;

            var fsPutDelta = Math.Exp(-q * T0) * put;

            var fsPutGamma = 0.0;

            var fsPutTheta = s0 * Math.Exp(-q * T0) * HestonNumericalGreeks.ThetaPut(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            var fsPutVega = s0 * Math.Exp(-q * T0) * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaPut(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            var rhoGreek = s0 * Math.Exp(-q * T0) * HestonRho.RhoPut(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r, q: q);

            Engine.Verbose = verbosity;

            var result = new GreeksDerivatives()
            {
                Deltas = (Vector)fsPutDelta,
                Gammas = (Vector)fsPutGamma,
                Theta = fsPutTheta,
                Vegas = (Vector)fsPutVega,
                MarkToMarket = fsPutMarkToMarket,
                Rho = rhoGreek
            };

            if (verbosity > 0)
            {
                Console.WriteLine("Price and Greeks of the Forward Starting Put Option");
                Console.WriteLine("Delta\tGamma\tTheta\tVega\tRho\tPrice");
                Console.WriteLine($"{result.Deltas[0]}\t{result.Gammas[0]}\t{result.Theta}\t{result.Vegas[0]}\t{result.Rho}\t{result.MarkToMarket}");
            }

            return result;

        }

        public static GreeksDerivatives HestonForwardDigitalPutWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price and greeks of a FS Digital Put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{r}\t{q}");
            }

            Engine.Verbose = 0;
            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            var dPut = HestonDigital.HestonDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q);


            var fsDPutMarkToMarket = HestonDigital.DiscountFactor(r, T0) * dPut;

            var fsDPutDelta = 0.0;

            var fsDPutGamma = 0.0;

            var fsDPutTheta = HestonDigital.DiscountFactor(r, T0) * HestonNumericalGreeks.ThetaDPut(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            var fsDPutVega = HestonDigital.DiscountFactor(r, T0) * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaDigitalPut(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            var fsDPutRho = HestonDigital.DiscountFactor(r, T0) * HestonRho.RhoDigitalPut(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            fsDPutRho += -T0 * fsDPutMarkToMarket;
            Engine.Verbose = verbosity;


            var result = new GreeksDerivatives()
            {
                Deltas = (Vector)fsDPutDelta,
                Gammas = (Vector)fsDPutGamma,
                Theta = fsDPutTheta,
                Vegas = (Vector)fsDPutVega,
                MarkToMarket = fsDPutMarkToMarket,
                Rho = fsDPutRho
            };

            if (verbosity > 0)
            {
                Console.WriteLine("Price and Greeks of the Forward Starting Digital Put Option");
                Console.WriteLine("Delta\tGamma\tTheta\tVega\tRho\tPrice");
                Console.WriteLine($"{result.Deltas[0]}\t{result.Gammas[0]}\t{result.Theta}\t{result.Vegas[0]}\t{result.Rho}\t{result.MarkToMarket}");
            }

            return result;

        }

        public static GreeksDerivatives HestonForwardDigitalCallWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price and greeks of a FS Digital Call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{r}\t{q}");
            }

            Engine.Verbose = 0;

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);
            
            var dCall = HestonDigital.HestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q);


            var fsDCallMarkToMarket = HestonDigital.DiscountFactor(r, T0) * dCall;

            var fsCallDelta = 0.0;

            var fsDCallGamma = 0.0;

            var fsDCallTheta = HestonDigital.DiscountFactor(r, T0) * HestonNumericalGreeks.ThetaDCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            var fsDCallVega = HestonDigital.DiscountFactor(r, T0) * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaDigitalCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            var fsDCallRho = HestonDigital.DiscountFactor(r, T0) * HestonRho.RhoDigitalCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);

            fsDCallRho += -T0 * fsDCallMarkToMarket;

            Engine.Verbose = verbosity;

            var result = new GreeksDerivatives()
            {
                Deltas = (Vector)fsCallDelta,
                Gammas = (Vector)fsDCallGamma,
                Theta = fsDCallTheta,
                Vegas = (Vector)fsDCallVega,
                MarkToMarket = fsDCallMarkToMarket,
                Rho = fsDCallRho
            };

            if (verbosity > 0)
            {
                Console.WriteLine("Price and Greeks of the Forward Starting Digital Call Option");
                Console.WriteLine("Delta\tGamma\tTheta\tVega\tRho\tPrice");
                Console.WriteLine($"{result.Deltas[0]}\t{result.Gammas[0]}\t{result.Theta}\t{result.Vegas[0]}\t{result.Rho}\t{result.MarkToMarket}");
            }

            return result;
        }

        #endregion


        }
}
