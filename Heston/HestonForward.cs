using System;
using Heston;
using Fairmat.Math;
using DVPLI;
using DVPLI.ExportableDomain;
using System.Xml.Schema;
using System.Diagnostics;
using Fairmat.Finance;

namespace HestonEstimator
{
    // this implementation is based on https://doi.org/10.1142/S0219024909005166 
    // FORWARD START OPTIONS UNDER STOCHASTIC VOLATILITY AND STOCHASTIC INTEREST RATES - REHEZ AHLIP and MAREK RUTKOWSKI
    public class HestonForwardAhlipRutkowski : HestonCall
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

        public struct ModelParameters
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

        public struct CallParameters
        {
            public double T;
            public double T0;
            public double K;
        }


        internal static Complex s1(double u, ModelParameters mp, CallParameters cp)
        {

            // s1 = a + b + c

            var iu = Complex.I * u;
            var onePlusIu = 1 + iu;

            var a = -mp.kappa * mp.rho * onePlusIu / mp.sigma;
            var b = - (1 - mp.rho * mp.rho) * (onePlusIu * onePlusIu) / 2;
            var c = onePlusIu / 2;

            return a + b + c;

        }

        internal static Complex s2(double u, ModelParameters mp)
        {
            // s2 = a / b

            var iu = Complex.I * u;

            var a = -mp.rho * (1 + iu);
            var b = mp.sigma;

            return a / b;

        }

        public static Complex q1(double u, ModelParameters mp)
        {
            // q1 = a * ( b + c + d)

            var iu = u * Complex.I;

            var a = -iu / 2;
            var b = 2 * mp.kappa * mp.rho / mp.sigma;
            var c = iu * (1 - mp.rho * mp.rho);
            var d = -1.0;

            return a * (b + c + d);

        }

        public static Complex q2(double u, ModelParameters mp)
        {
            // q2 = a / b 

            var iu = Complex.I * u;

            var a = - mp.rho * iu;
            var b = mp.sigma;

            return a / b;

        }


        internal static Complex G1(double tau, Complex lambda, Complex mu, ModelParameters mp, double t)
        {
            // g1 = ( a + b ) / ( c + d + e + f)

            var sigma2 = (mp.sigma * mp.sigma);
            var gamma = Gamma(mu: mu, kappa: mp.kappa, sigma: mp.sigma);

            var a = lambda * ( gamma + mp.kappa + Complex.Exp(gamma * tau) * (gamma - mp.kappa));
            var b = 2 * mu * (1 - Complex.Exp(t * gamma));

            var c = sigma2 * lambda * (Complex.Exp(gamma * tau) - 1);
            var d = gamma;
            var e = -mp.kappa;
            var f = Complex.Exp(gamma * tau) * (gamma + mp.kappa);

            return (a + b) / (c + d + e + f);


        }

        internal static Complex H1(double tau, Complex lambda, Complex mu, ModelParameters mp)
        {

            // H1 = a* ln( b / (c + d + e + f) )

            var sigma2 = (mp.sigma * mp.sigma);
            var gamma = Gamma(mu: mu, kappa: mp.kappa, sigma: mp.sigma);


            var a = -2 / sigma2;
            var b = 2 * gamma * Complex.Exp((gamma + mp.kappa) * tau / 2);
            var c = sigma2 * lambda * (Complex.Exp(gamma * tau) - 1);
            var d = gamma;
            var e = -mp.kappa;
            var f = Complex.Exp(gamma * tau) * (gamma + mp.kappa);

            return a * Complex.Log( b / (c + d + e + f) );

        }

        internal static Complex H1Hat(double tau, Complex lambda, ModelParameters mp)
        {
            // H1Hat = a * ln( b / c )

            var kappaHat = mp.kappa - mp.rho * mp.sigma;
            var sigma2 = (mp.sigma * mp.sigma);

            var a = -2 / sigma2;
            var b = 2 * kappaHat * Math.Exp(2 * kappaHat * tau);
            var c = sigma2 * lambda * (Math.Exp(kappaHat * tau) - 1);
            c += 2 * kappaHat * Math.Exp(kappaHat * tau);

            return a * Complex.Log(b / c);


        }


        internal static Complex Gamma(Complex mu, double kappa, double sigma)
        {
            var argument = kappa * kappa + 2 * sigma * sigma * mu;
            return Complex.Sqrt(argument);
        }

        internal static Complex Gamma(double mu, double kappa, double sigma)
        {
            var muC = new Complex(mu);
            return Gamma(muC, kappa, sigma);
        }

        internal static double GammaHat(double mu, double b, ModelParameters mp)
        {
            var argument = b * b + 2 * mp.sigma * mp.sigma * mu;
            return Math.Sqrt(argument);
        }

        internal static Complex G1Hat(ModelParameters mp, double tau, Complex lambda)
        {
            // TODO 
            // is Gamma of kappa or of kappaHat ? 
            var kappaHat = mp.kappa - mp.rho * mp.sigma;
            var sigma2 = (mp.sigma * mp.sigma);
            var gamma = Gamma(mu: 0.0, kappa: mp.kappa, sigma: mp.sigma);

            var numerator = 2 * lambda * kappaHat;
            var denominator = sigma2 * lambda * (Math.Exp(kappaHat * tau) - 1) + 2 * kappaHat * Complex.Exp(tau * gamma);

            return numerator / denominator;

        }

        internal static Complex s2Hat(double u, double tau0, Complex s1, Complex s2, ModelParameters mp, double t)
        {
            // s2hat = (1 + iu) * rho / sigma + g1(tau0, s1, s2)

            Complex iu = Complex.I * u;
            var g1 = G1(tau: tau0, lambda: s1, mu: s2, mp: mp, t: t);
            return g1 + (1 + iu) * (mp.rho / mp.sigma);
        }

        internal static Complex q2Hat(double u, double tau0, Complex q1, Complex q2, ModelParameters mp, double t)
        {
            // q2hat = - i * u * rho / sigma - G1(tau0, q1, q2) 
            // here there was a typo 
            Complex iu = Complex.I * u;
            return -iu * mp.rho / mp.sigma - G1(tau: tau0, lambda: q1, mu: q2, mp: mp, t: t);
        }

        internal static Complex f1Hat(double u, ModelParameters mp, CallParameters cp)
        {
            // f1Hat = exp (a + b + c + d)
            // a = -(1+iu) * rho * theta * tau0 / sigma
            // b = - theta H1(tau0, s1, s2)
            // c = -v0 * G1Hat(tau, s2Hat)
            // d = -theta * H1Hat(tau, s2Hat)

            var t = 0;
            var tau0 = cp.T - cp.T0;
            var tau = cp.T0;

            var iu = Complex.I * u;
            var _s1 = s1(u: u, mp: mp, cp: cp);
            var _s2 = s2(u: u, mp: mp);
            var s2H = s2Hat(u: u, tau0: tau0, s1: _s1, s2: _s2, t: t, mp: mp);


            var a = -(1 + iu) * mp.rho * mp.theta * tau0 / mp.sigma;
            var b = -mp.theta * H1(tau: tau0, lambda: _s1, mu: _s2, mp: mp);
            var c = -mp.v0 * G1Hat(tau: tau, lambda: s2H, mp: mp);
            var d = -mp.theta * H1Hat(tau: tau, lambda: s2H, mp: mp);

            return Complex.Exp(a + b + c + d);

        }

        

        internal static Complex f2Hat(double u, ModelParameters mp, CallParameters cp)
        {
            // f1Hat = exp (a + b + c + d)
            // a = - iu * rho * theta * tau0 / sigma
            // b = - theta H1(tau0, q1, q2)
            // c = -v0 * G1Hat(tau, q2Hat)
            // d = -theta * H1Hat(tau, q2Hat)

            var t = 0;
            var tau0 = cp.T - cp.T0;
            var tau = cp.T0;


            var iu = Complex.I * u;
            var _q1 = q1(u: u, mp: mp);
            var _q2 = q2(u: u, mp: mp);
            var q2H = q2Hat(u: u, tau0: tau0, q1: _q1, q2: _q2, mp: mp, t: t);


            var a = -iu * mp.rho * mp.theta * tau0 / mp.sigma;
            var b = -mp.theta * H1(tau: tau0, lambda: _q1, mu: _q2, mp: mp);
            var c = -mp.v0 * G1Hat(tau: tau, lambda: q2H, mp: mp);
            var d = -mp.theta * H1Hat(tau: tau, lambda: q2H, mp: mp);

            return Complex.Exp(a + b + c + d);

        }

        internal static double IntegrandP1(double u, ModelParameters mp, CallParameters cp)
        {
            // integrand of P1Hat = f1Hat(u) * exp(-i * u * ln(k) ) / iu 

            var iu = Complex.I * u;
            var logK = Math.Log(cp.K);
            var f1 = f1Hat(u: u, mp: mp, cp: cp);
            return (f1 * Complex.Exp(-iu * logK) / iu).Re;


        }

        internal static double IntegrandP2(double u, ModelParameters mp, CallParameters cp)
        {
            // integrand of P2Hat = f2Hat(u) * exp(-i * u * ln(k) ) / iu 

            var iu = Complex.I * u;
            var logK = Math.Log(cp.K);
            var f2 = f2Hat(u: u, cp: cp, mp: mp);
            return (f2 * Complex.Exp(-iu * logK) / iu).Re;


        }


        internal static double ComputeIntegral(TAEDelegateFunction1D f, double a = 1E-8, double b = 1000.0)
        {
            var res = PerformIntegral(a, b, f);
            res += a * f(a / 2.0);
            return res;

        }


        internal static double P1Call(ModelParameters mp, CallParameters cp)
        {
            double functionToIntegrate(double u)
            {
                // f1 = exp( a + b + c)
                var strikeMonetary = cp.K * mp.s0;
                var tau = cp.T;
                var iu = Complex.I * u;

                var _s1 = s1(u, mp, cp);
                var _s2 = s2(u, mp);

                var a = - mp.rho * (1 + iu) / mp.sigma * (mp.v0 + mp.theta * tau);
                var b = - mp.v0 * G1(tau: tau, _s1, _s2, mp: mp, t: 0.0);
                var c = - mp.theta * H1(tau: tau, _s1, _s2, mp);

                var f1 = Complex.Exp(a + b + c);

                var multiplier = (1/iu) * Complex.Exp( -iu * Math.Log(strikeMonetary));

                return (f1 * multiplier).Re;

            }

            var integral = ComputeIntegral(functionToIntegrate);
            return 0.5 + (1.0 / Math.PI) * integral;


        }

        internal static double P2Call(ModelParameters mp, CallParameters cp)
        {

            double functionToIntegrate(double u)
            {
                // f2 = exp( a + b + c)

                var strikeMonetary = cp.K * mp.s0;
                var tau = cp.T;

                var _q1 = q1(u:u, mp);
                var _q2 = q2(u:u, mp);
                var iu = Complex.I * u;


                var a = - iu * mp.rho / mp.sigma * (mp.v0 + mp.theta * tau);
                var b = -mp.v0 * G1(tau: tau, _q1, _q2, mp: mp, t: 0.0);
                var c = -mp.theta * H1(tau: tau, _q1, _q2, mp);

                var f2 = Complex.Exp(a + b + c);
                var multiplier = (1/iu) * Complex.Exp(-iu * Math.Log(strikeMonetary));

                return (f2 * multiplier).Re;
            }

            var integral = ComputeIntegral(functionToIntegrate);
            return 0.5 + (1.0 / Math.PI) * integral;

        }


        public static double CallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            // adjust the theta value to match their theta definition
            // in ours approach the drift is kappa(theta-v)
            // in their approach the drift is (theta_start-kappa*v) 
            // theta_start-kappa*v = kappa*theta- kappa*v
            // theta_start = kappa * theta 

            double thetaAdjusted = kappa * theta;
            var strikeMonetary = K * s0;


            var cp = new CallParameters() 
            { T = T, T0 = T0, K = K };
            var mp = new ModelParameters() 
            { kappa = kappa, theta = thetaAdjusted, sigma = sigma, rho = rho, v0 = v0, s0 = s0, r = r, q = q };

            
            var b_t_T = Math.Exp(-cp.T * mp.r);
            var p1 = Math.Max(P1Call(mp:mp, cp), 0);
            var p2 = Math.Max(P2Call(mp: mp, cp),0) ;

            return mp.s0 * p1 - b_t_T * strikeMonetary * p2;

        }



        internal static double P1Hat( ModelParameters mp, CallParameters cp)
        {
            // P1Hat = 1/2 + 1/pi * integral from 0 to infinity of IntegrandP1(u) du
            
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


        internal static Complex Alpha(Complex phi)
        {
            return 1 + phi * Complex.I;
        }

        internal static Complex Alpha(double phi)
        {
            Complex phiC = new Complex(phi);
            return Alpha(phiC);
        }

        internal static Complex Beta(Complex phi)
        {
            return 1 - phi * Complex.I;
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

            var stabilityCondition = (2 * mp.theta / sigma) > 1;
            if (!stabilityCondition)
            {
                throw new Exception("Stability condition is not satisfied");
            }


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

    public class HestonForwardJJ: HestonCall
    {

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonForward"/> class.
        /// </summary>
        public HestonForwardJJ() { }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonForward"/> class.
        /// </summary>
        /// <param name='problem'>
        /// HestonCallOptimizationProblem at which digital price calculations are to be linked.
        /// </param>
        public HestonForwardJJ(HestonCallOptimizationProblem problem) : base(problem) { }

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
        public HestonForwardJJ(HestonProcess process, double strike, double timeToMaturity) : base(process, strike, timeToMaturity) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonForward"/> class.
        /// </summary>
        /// <param name='process'>
        /// HestonProcess at which digital price calculations are to be linked.
        /// </param>
        public HestonForwardJJ(HestonProcess process) : base(process) { }

        public struct ModelParameters
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

        public struct CallParameters
        {
            public double T;
            public double T0;
            public double K;
        }


        public static double LittleB(ModelParameters mp, double t)
        {
        double kappa = mp.kappa;
        double xi = mp.sigma;
        return ((xi * xi) / (4.0 * kappa)) * (1.0 - Math.Exp(-kappa * t));
        }

        public static Complex D(ModelParameters mp, Complex u)
        {
            double kappa = mp.kappa;
            double xi = mp.sigma;
            double rho = mp.rho;
            return Complex.Sqrt(((kappa - rho * xi * u) * (kappa - rho * xi * u)) + u * (1.0 - u) * xi * xi);
        }

        public static Complex LittleGamma(ModelParameters mp, Complex u)
        {
            double kappa = mp.kappa;
            double xi = mp.sigma;
            double rho =mp.rho;
            Complex dd = D(mp, u);
            return (kappa - rho * xi * u - dd) / (kappa - rho * xi * u + dd);
        }

        public static Complex BigA(ModelParameters mp, Complex u, double tau)
        {
            double kappa = mp.kappa;
            double xi = mp.sigma;
            double rho = mp.rho;
            double theta = mp.theta;
            Complex ddd = D(mp, u);
            Complex g = LittleGamma(mp, u);
            Complex LL = Complex.Log((1.0 - g * Complex.Exp(-ddd * tau)) / (1.0 - g));
            return ((kappa * theta) / (xi * xi)) * ((kappa - rho * xi * u - ddd) * tau - 2.0 * LL);
        }

        public static  Complex BigB(ModelParameters mp, Complex u, double tau)
        {
            double kappa = mp.kappa;
            double xi = mp.sigma;
            double rho = mp.rho;
            Complex ddd = D(mp, u);
            Complex ee = Complex.Exp(-ddd * tau);
            return (1.0 / (xi * xi)) * (kappa - rho * xi * u - ddd) * ((1.0 - ee) / (1.0 - LittleGamma(mp,u) * ee));
        }


        public static Complex CfHestonFwd(ModelParameters mp, Complex u, double t, double tau)
        {
            double kappa = mp.kappa;
            double xi = mp.sigma;
            double v0 = mp.v0;
            Complex BBB = BigB(mp, u, tau);
            double BBb = LittleB(mp, t);
            return Complex.Exp(BigA(mp, u, tau) + (BBB / (1.0 - 2.0 * BBb * BBB)) * (v0 * (Math.Exp(-kappa * t))) - ((2.0 * kappa * mp.theta) / (xi * xi)) * Complex.Log(1.0 - 2.0 * BBb * BBB));
        }

        public static double IntPhi(ModelParameters mp, double w, double alpha, double k, double t, double tau)
        {
            return (Complex.Exp(-k * (alpha + Complex.I * w)) * CfHestonFwd(mp, alpha + Complex.I * w, t, tau) / ((alpha + Complex.I * w) * (1.0 - alpha - Complex.I * w))).Re;
        }

        public static double TrapezoidalRule( double a, double b, Func<double, double> f, int n=100000)
        {
            double h = (b - a) / n; // Step size
            double s = f(a) + f(b); // Initialize sum

            for (int i = 1; i < n; i++)
            {
                double x = a + h * i;
                s += 2 * f(x);
            }

            return (h / 2) * s;
        }


        //private double CallPrice(ModelParameters mp, double K, double t, double tau)
        public static double CallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {

            var mp = new ModelParameters{
                kappa = kappa,
                theta = theta,
                sigma= sigma,
                rho = rho,
                v0 = v0,
                s0 = s0,
                r = r,
                q = q };

            // adjust T 
            T = T - T0; 
            double alpha = 0.5;
            double pp = CfHestonFwd(mp, new Complex(1, 0), T0, T).Re;
            var integral = TrapezoidalRule(-1000.0, 1000.0, (x) => IntPhi(mp, x, alpha, Math.Log(K), T0, T));
            return pp - K / (2.0 * Math.PI) * integral;
        }

       

    }

}
