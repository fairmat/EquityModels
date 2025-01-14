using System;
using Heston;
using Fairmat.Math;
using DVPLI;



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

        protected struct ModelParameters
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

        protected struct CallParameters
        {
            public double T;
            public double T0;
            public double K;
        }


        public static Complex s1(Complex u, double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            var iu = Complex.I * u;
            var onePlusIu = Complex.One + iu;
            var result = -kappa * rho * onePlusIu / sigma;
            result -= 0.5 * (1 - rho * rho) * Complex.Pow(onePlusIu, 2);
            result += onePlusIu / 2;
            return result;
        }

        public static Complex s2(Complex u, double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            var iu = Complex.I * u;
            var onePlusIu = Complex.One + iu;
            return - rho * onePlusIu / sigma;
        }

        public static Complex q1(Complex u, double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            var iu = Complex.I * u;
            var preMult = -iu / 2;
            return preMult * (2 * kappa * rho / sigma + iu * (1 - rho * rho) - 1.0);

        }

        public static Complex q2(Complex u, double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            var iu = Complex.I * u;
            return -iu * rho / sigma;
        }


        public static Complex G1(double tau, Complex lambda, Complex mu, double kappa, double sigma, double t)
        {
            var gamma = Complex.Sqrt(kappa * kappa + 2 * sigma * sigma * mu);
            var expGammaTau = Complex.Exp(gamma * tau);
            var numerator = lambda * ((kappa + gamma) + expGammaTau * (gamma - kappa));
            numerator += 2 * mu * (1 - Complex.Exp(gamma * t));

            var denominator = sigma * sigma * lambda * (expGammaTau - 1) + gamma - kappa;
            denominator += expGammaTau * (gamma + kappa);

            return numerator / denominator;

        }

        public static Complex H1(double tau, Complex lambda, Complex mu, double kappa, double sigma)
        {
            var gamma = Complex.Sqrt(kappa * kappa + 2 * sigma * sigma * mu);
            var sigma2 = (sigma * sigma);
            var preMult = -2 / sigma2;
            var argumentNumerator = 2 * gamma * Complex.Exp((gamma + kappa) * tau / 2);
            var argumentDenominator = sigma2 * lambda * (Complex.Exp(gamma * tau) - 1) + gamma - kappa + Complex.Exp(gamma * tau) * (gamma + kappa);
            return preMult * Complex.Log(argumentNumerator / argumentDenominator);

        }

        //public static Complex H1Hat(double tau, Complex lambda, double sigma, double kappa)
        public static Complex H1Hat(double tau, Complex lambda, double sigma, double kappa, double rho)
        {
            var kappaHat = kappa - sigma * rho;
            var sigma2 = (sigma * sigma);
            var preMult = -2 / sigma2;
            var argumentNumerator = 2 * kappaHat * Math.Exp(kappaHat * tau);    // the multiplication by two is intentionally omitted by the exponent for calculations that differ from the paper.

            var argumentDenominator = sigma2 * lambda * (Math.Exp(kappaHat * tau) - 1) + 2 * kappaHat * Math.Exp(kappaHat * tau);
            return preMult * Complex.Log(argumentNumerator / argumentDenominator);

        }

        public static Complex G1Hat(double sigma, double kappa, double tau, Complex lambda, double rho)
        {
            var kappaHat = kappa - sigma * rho;
            var sigma2 = (sigma * sigma);
            var argumentNumerator = 2 * kappaHat * lambda;
            var argumentDenominator = sigma2 * lambda * (Math.Exp(kappaHat * tau) - 1) + 2 * kappaHat * Math.Exp(kappaHat * tau);
            return argumentNumerator / argumentDenominator;
        }

        public static Complex s2Hat(Complex u, double tau0, double rho, double sigma, Complex s1, Complex s2, double kappa, double t)
        {
            Complex iu = Complex.I * u;
            var g1 = G1(tau: tau0, lambda: s1, mu: s2, kappa: kappa, sigma: sigma, t: t);
            return g1 + (1 + iu) * (rho / sigma);
        }

        public static Complex q2Hat(Complex u, double tau0, double rho, double sigma, Complex q1, Complex q2, double kappa, double t)
        {
            Complex iu = Complex.I * u;
            return -iu * rho / sigma - G1(tau: tau0, lambda: q1, mu: q2, kappa: kappa, sigma: sigma, t: t);
        }

        public static Complex f1Hat(Complex u, double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            var tau0 = T - T0;
            var tau = T0;
            var t = 0;

            var iu = Complex.I * u;
            var _s1 = s1(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0);
            var _s2 = s2(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0);
            var s2H = s2Hat(u: u, tau0: tau0, rho: rho, sigma: sigma, s1: _s1, s2: _s2, kappa: kappa, t: t);


            var a = -(1 + iu) * rho * theta * tau0 / sigma;
            var b = -theta * H1(tau: tau0, lambda: _s1, mu: _s2, kappa: kappa, sigma: sigma);
            var c = -v0 * G1Hat(tau: tau, lambda: s2H, sigma: sigma, kappa: kappa, rho: rho);
            var d = -theta * H1Hat(tau: tau, lambda: s2H, sigma: sigma, kappa: kappa, rho: rho);
            var info = Complex.Exp(a + b + c + d);

            return Complex.Exp(a + b + c + d);

        }

        public static Complex f2Hat(Complex u, double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            var tau0 = T - T0;
            var tau = T0;
            var t = 0;

            var iu = Complex.I * u;
            var _q1 = q1(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0);
            var _q2 = q2(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0);
            var q2H = q2Hat(u: u, tau0: tau0, rho: rho, sigma: sigma, q1: _q1, q2: _q2, kappa: kappa, t: t);


            var a = -iu * rho * theta * tau0 / sigma;
            var b = -theta * H1(tau: tau0, lambda: _q1, mu: _q2, kappa: kappa, sigma: sigma);
            var c = -v0 * G1Hat(tau: tau, lambda: q2H, sigma: sigma, kappa: kappa, rho: rho);
            var d = -theta * H1Hat(tau: tau, lambda: q2H, sigma: sigma, kappa: kappa, rho: rho);

            return Complex.Exp(a + b + c + d);

        }

        public static double IntegrandP1(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            var uComp = new Complex(u);
            var iu = Complex.I * u;
            var logK = Math.Log(K);
            var f1 = f1Hat(u: uComp, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0);
            var integrand_f1 = (f1 * Complex.Exp(-iu * logK) / iu).Re;
            return (f1 * Complex.Exp(-iu * logK) / iu).Re;


        }

        public static double IntegrandP2(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            var uComp = new Complex(u);

            var iu = Complex.I * u;
            var logK = Math.Log(K);
            var f2 = f2Hat(u: uComp, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0);
            return (f2 * Complex.Exp(-iu * logK) / iu).Re;

        }

        public static double P1Hat(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {  
            // compute the integral 
            double a = 1E-8;
            double b = 1000.0;
            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandP1(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, T0: T0, K: K);
            double part1 = PerformIntegral(a, b, functionToIntegrate);
            double integral = part1 + a * functionToIntegrate(a / 2.0);
            return 0.5 + 1.0 / Math.PI * integral;

        }

        public static double P2Hat(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            // compute the integral 
            double a = 1E-8;
            double b = 10000.0;
            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandP2(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, T0: T0, K: K);
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

            // use theta adjusted instead of theta 
            var p1 = P1Hat(kappa: kappa, theta: thetaAdjusted, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0);
            var p2 = P2Hat(kappa: kappa, theta: thetaAdjusted, rho: rho, v0: v0, sigma: sigma, s0: s0, K: K, r: r, q: q, T: T, T0: T0);
            var B_T_T0 = Math.Exp(-r * (T - T0));
            return s0 * (p1 - B_T_T0 * K * p2);
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
                T0: T0,
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
            return lambda + rho * eta * (m - 2 + Theta_tau(tau, T0: T0, T: T) * (2 - m - j));
        }

        public static double Theta_tau(double tau, double T, double T0)
        {
            var x = T - T0;
            if (tau >= 0.0 & tau <= x)
                return 1.0;
            return 0.0;
        }

        public static Complex Beta(double lambda, double rho, double eta, double j, double T, double T0, double tau, double kappa)
        {
            var theta = Theta_tau(tau: tau, T0: T0, T: T);
            var result = -Complex.I * rho * eta * kappa * theta;
            result += b(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, tau: tau);
            return result;

        }

        public static Complex Alpha(double lambda, double rho, double eta, double j, double T, double T0, double tau, double kappa)
        {
            var i = Complex.I;
            var theta = Theta_tau(tau: tau, T0: T0, T: T);
            return theta * (-kappa * kappa / 2 - i * kappa / 2 + i * j * kappa);
        }


        public static Complex d(double lambda, double rho, double eta, double j, double T, double T0, double kappa)
        {
            var eta2 = eta * eta;
            var beta0 = Beta(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, tau: 0.0, kappa: kappa);
            var arg = beta0 * beta0 - 2 * eta2 * Alpha(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, tau: 0.0, kappa: kappa);


            return Complex.Sqrt(arg);
        }

        public static Complex rPlus(double lambda, double rho, double eta, double j, double T, double T0, double kappa)
        {
            var beta0 = Beta(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, tau: 0.0, kappa: kappa);
            var _d = d(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            return (beta0 + _d) / (eta * eta);
        }

        public static Complex rMinus(double lambda, double rho, double eta, double j, double T, double T0, double kappa)
        {
            var beta = Beta(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, tau: 0.0, kappa: kappa);
            var _d = d(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            return (beta - _d) / (eta * eta);
        }

        public static Complex g(double lambda, double rho, double eta, double j, double T, double T0, double kappa)
        {
            var _rMinus = rMinus(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var _rPlus = rPlus(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            return _rMinus / _rPlus;
        }


        public static Complex D_small_tau(double lambda, double rho, double eta, double j, double T, double T0, double kappa, double tau)
        {
            var _d = d(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var _rMinus = rMinus(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var _g = g(lambda, rho, eta, j, T, T0, kappa);
            var numerator = 1 - Complex.Exp(-_d * tau);
            var denominator = 1 - _g * Complex.Exp(-_d * tau);

            return _rMinus * numerator / denominator;

        }

        

        public static Complex C_small_tau(double lambda, double rho, double eta, double j, double T, double T0, double kappa, double tau)
        {
            var _d = d(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var _rMinus = rMinus(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var _g = g(lambda, rho, eta, j, T, T0, kappa);
            var log_numerator = 1 - _g * Complex.Exp(-_d * tau);
            var log_denominator = 1 - _g;
            var eta2 = eta * eta;

            return lambda * (_rMinus * tau - 2 * Complex.Log(log_numerator / log_denominator) / (eta2));

        }


        public static Complex c(double lambda, double rho, double eta, double j, double T, double T0, double kappa)
        {
            var betaT = Beta(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, tau: T, kappa: kappa);
            var tau = (T - T0);
            var _D = D_small_tau(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa, tau: tau);
            var eta2 = eta * eta;
            return 2 * betaT / (eta2 * _D) - 1;

        }


        public static Complex D_big_tau(double lambda, double rho, double eta, double j, double T, double T0, double kappa, double tau)
        {
            var betaT = Beta(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, tau: T, kappa: kappa);
            var _c = c(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var numerator = 2 * betaT;
            var denominator = eta * eta * (1 + _c * Complex.Exp(betaT * (tau - T + T0)));

            return    numerator / denominator;

        }
        


        public static Complex C_big_tau(double lambda, double rho, double eta, double j, double T, double T0, double kappa, double tau)
        {
            var betaT = Beta(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, tau: T, kappa: kappa);
            var _c = c(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa);
            var tauC = (T - T0);
            var _C = C_small_tau(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa, tau: tauC);
            var part2 = (2 * betaT * lambda * (tau - T + T0)) / (eta * eta);
            var log_numerator = 1 + _c * Complex.Exp(betaT * (tau - T + T0));
            var log_denominator = 1 + _c;
            var eta2 = eta * eta;
            var premult = -2 * lambda / eta2;
            var part3 = premult * Complex.Log(log_numerator / log_denominator);
            return _C + part2 + part3;

        }


        public static Complex P_tilda(double lambda, double rho, double eta, double j, double T, double T0, double kappa, double tau, double vhat, double v)
        {
            var i = Complex.I;
            var C = C_big_tau(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa, tau: tau);
            var D = D_big_tau(lambda: lambda, rho: rho, eta: eta, j: j, T: T, T0: T0, kappa: kappa, tau:tau);
            var arg_exp = C * vhat + D * v;
            var pre_mult = 1 / (i * kappa);  
            var result = pre_mult * Complex.Exp(arg_exp);
            
            return result;

        }

        public static double x(double T, double T0, double rate, double K)
        {
            var P_0_T0 = Math.Exp(-rate * T0);
            var P_0_T = Math.Exp(-rate * T);

            return Math.Log(P_0_T0 / P_0_T / K);
        }


        public static Complex PerformIntegral(double kMin, double kMax, Func<double, Complex> F)
        {
             
            int steps = 11;
            double deltak = (kMax - kMin) / steps;
            Complex sum = Complex.Zero;

            for (int i = 0; i <= steps; i++)
            {
                double kappa = kMin + i * deltak;
                Complex integrand = F (kappa) ; 
                sum += integrand;
            }

            // Approssimazione dell'integrale usando la regola dei trapezi
            Complex result = deltak * sum ;
            return result;
        }



        public static Complex P(double lambda, double rho, double eta, double j, double T, double T0, double tau, double vhat, double v, double rate, double K)
        {
            var i = Complex.I;
            var _x = x(T: T, T0: T0, rate: rate, K: K);
            Func<double, Complex> F = (kappa) => (P_tilda(lambda: lambda,rho: rho,eta: eta,j: j,T: T,T0: T0, kappa: kappa, tau: tau, vhat: vhat, v: v) * Complex.Exp(Complex.I * kappa * _x));
            Complex result = PerformIntegral(-1, 1, F) / (2.0 * Math.PI);


            return result;

        }


        

        public static Complex HestonForwardLucicCallPrice(double lambda, double rho, double eta, double T, double T0, double tau, double vhat, double v, double rate, double K)
        {
            var _x = x(T: T, T0: T0, rate: rate, K: K);
            var P1 = P(lambda: lambda, rho: rho, eta: eta, j: 1, T: T, T0: T0, tau: tau, vhat: vhat, v: v, rate: rate, K: K);
            var P0 = P(lambda: lambda, rho: rho, eta: eta, j: 0, T: T, T0: T0, tau: tau, vhat: vhat, v: v, rate: rate, K: K);
            //P0.Re = Math.Max(P0.Re, 0); 
            var exp_x = Math.Exp(_x);
            return K * (exp_x * P1 - P0) * Math.Exp(-rate * T);

        }




    }


    public class HestonForwardApproximated
    {
        public static double HestonForwardCallPrice(Vector x, double s0, double T, double T0, double K, double r, double q, double? Tp = null)
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
                q: q,
                timeToPaymentDate: Tp?? T
                );
        }
        public static double HestonForwardPutPrice(Vector x, double s0, double T, double T0, double K, double r, double q, double? Tp = null)
        {
            var kappa = x[0];
            var theta = x[1];
            var sigma = x[2];
            var rho = x[3];
            var v0 = x[4];


            return HestonForwardPutPrice(
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
                q: q,
                paymentDate: Tp ?? T
                );
        }

        public static double HestonForwardPercentageCallPrice(Vector x, double s0, double T, double T0, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            var kappa = x[0];
            var theta = x[1];
            var sigma = x[2];
            var rho = x[3];
            var v0 = x[4];


            return HestonForwardPercentageCallPrice(
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
                q: q,
                discountingFactorFunction: discountingFactorFunction,
                paymentTime: paymentTime
                );
        }
        public static double HestonForwardPercentagePutPrice(Vector x, double s0, double T, double T0, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            var kappa = x[0];
            var theta = x[1];
            var sigma = x[2];
            var rho = x[3];
            var v0 = x[4];


            return HestonForwardPercentagePutPrice(
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
                q: q,
                discountingFactorFunction: discountingFactorFunction,
                paymentTime: paymentTime
                );
        }

        public static double ExpectationCIRProcess(double x0, double a, double b, double t)
        {
            return x0 * Math.Exp(-a * t) + b * (1 - Math.Exp(-a * t));
        }

        public static double derivativeExpectationCIRProcess(double x0, double a, double b, double t)
        {
            return Math.Exp(-a * t);
        }

        public static double tDerivativeExpectationCIRProcess(double x0, double a, double b, double t)
        {
            return -a * x0 * Math.Exp(-a * t) + a * b * Math.Exp(-a * t);
        }

        // we use similar formulas that are used for pricing Forward call options under Black Scholes model
        public static double HestonForwardCallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction=null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tTimeToPaymentDate\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{timeToPaymentDate}\t{r}\t{q}");

            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            var Tp = timeToPaymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);


            Engine.Verbose = 0;

            var c = HestonCall.HestonUndiscountedCallPrice(
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

            var forwardPrice = s0 * Math.Exp((r-q) * T0);
            var price = df * forwardPrice * c;

            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Call Option is {price}");
            }

            return price;
        }

        // we use similar formulas that are used for pricing Forward call options under Black Scholes model
        public static double HestonForwardPercentageCallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            int verbosity = Engine.Verbose;
            double Tp = paymentTime ?? T;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS percentage call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tTimeToPaymentDate\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");

            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;
            var c = HestonCall.HestonUndiscountedCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );
            Engine.Verbose = verbosity;

            // call(K, T-S) --> Df(0,T-S) * undiscountedCall

            // this is DF(0,Tp) * undiscountedCall
            var discountingFactor = discountingFactorFunction(0, Tp);
            var forwardStartingCall = discountingFactor * c;


            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Percentage Call Option is {forwardStartingCall}");
            }

            return forwardStartingCall;
        }

        public static double HestonForwardPutPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? paymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            int verbosity = Engine.Verbose;

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentDate ?? T;


            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\ttimeToPaymentDate\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");

            }

            var df = discountingFactorFunction(0, Tp);

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;
            var p = HestonCall.HestonUndiscountedPutPrice(
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

            var forwardPrice = s0 * Math.Exp((r-q) * T0);
            var price = df * forwardPrice * p;

            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Put Option is {price}");
            }

            return price;
        }

        public static double HestonForwardPercentagePutPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            int verbosity = Engine.Verbose;
            double Tp = paymentTime ?? T;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS percentage put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Put Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tTimeToPaymentDate\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");

            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;
            var c = HestonCall.HestonUndiscountedPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );
            Engine.Verbose = verbosity;

            // call(K, T-S) --> Df(0,T-S) * undiscountedCall

            // this is DF(0,Tp) * undiscountedCall
            var discountingFactor = discountingFactorFunction(0, Tp);
            var forwardStartingPut = discountingFactor * c;


            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Percentage Put Option is {forwardStartingPut}");
            }

            return forwardStartingPut;
        }

        public static double HestonForwardDigitalPutPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? Tp = null, Func<double, double, double> discountingFactorFunction = null)
        {

            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS Digital Put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\ttimeToPaymentDate\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");

            }

            if (discountingFactorFunction is null)
            {
                discountingFactorFunction = (t, T) => Math.Exp(-r * (T - t));
            }

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;

            // undiscounted digital put 
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
                q: q,
                discountingFactorFunction: (_, _) => 1.0);

            Engine.Verbose = verbosity;

            var tp = Tp ?? T;
            var digitalPrice = discountingFactorFunction(0, tp) * dput;

            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Digital Put Option is {digitalPrice}");
            }

            return digitalPrice;
        }

        public static double HestonForwardDigitalCallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? Tp = null, Func<double, double, double> discountingFactorFunction = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price of a FS Digital Call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\ttimeToPaymentDate\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");

            }

            if (discountingFactorFunction is null)
            {
                discountingFactorFunction = (t, T) => Math.Exp(-r * (T - t));
            }

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;

            // undiscounted digital call
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
                q: q,
                discountingFactorFunction: (_, _) => 1.0
                );
            Engine.Verbose = verbosity;

            var digitalPrice = discountingFactorFunction(0, Tp ?? T) * dcall;

            if (verbosity > 0)
            {
                Console.WriteLine($"Price of the Forward Starting Digital Call Option is {digitalPrice}");
            }

            return digitalPrice;
        }

        #region Greeks

        public static double FSCallCalculateDelta(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {

            var fsCallPrice = HestonForwardCallPrice(kappa: kappa, K: K, theta: theta,
                T: T,
                T0: T0,
                q: q,
                rho: rho,
                sigma: sigma,
                v0: v0,
                s0: s0,
                r: r,
                timeToPaymentDate: timeToPaymentDate,
                discountingFactorFunction: discountingFactorFunction);

            return fsCallPrice / s0;
        }

        public static double FSCallCalculateGamma(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0)
        {
            return 0.0;
        }

        public static double FSCallCalculateTheta(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            return HestonNumericalGreeks.ThetaFSCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);
        }

        public static double FSCallCalculateVega(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            if (discountingFactorFunction is null)
            {
                discountingFactorFunction = (t, T) => Math.Exp(-r * (T - t));
            }

            var Tp = timeToPaymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);
            var forwardPrice = s0 * Math.Exp((r - q) * T0);
            
            return df * forwardPrice * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaUndiscountedCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: ExpectationCIRProcess(v0, kappa, theta, T0),
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);
        }

        public static double FSCallCalculateRho(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            if (discountingFactorFunction is null)
            {
                discountingFactorFunction = (t, T) => Math.Exp(-r * (T - t));
            }

            var Tp = timeToPaymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);
            var forwardPrice = s0 * Math.Exp((r - q) * T0);
            var derivativeDf = -Tp * df;
            var derivativeForwardPrice = T0 * forwardPrice;
            var rhoUndiscounted = HestonRho.RhoUndiscountedCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: ExpectationCIRProcess(v0, kappa, theta, T0),
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r, q: q);

            var undiscountedCallPrice = HestonCall.HestonUndiscountedCallPrice(kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: ExpectationCIRProcess(v0, kappa, theta, T0),
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r, q: q);

            var rhoGreek = derivativeDf * forwardPrice * undiscountedCallPrice;
            rhoGreek += df * derivativeForwardPrice * undiscountedCallPrice;
            rhoGreek += df * forwardPrice * rhoUndiscounted;

            return rhoGreek;
        }

        public static GreeksDerivatives HestonForwardCallWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
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

            if (discountingFactorFunction is null)
            {
                discountingFactorFunction = (t, T) => Math.Exp(-r * (T - t));
            }

            var Tp = timeToPaymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);
            var forwardPrice = s0 * Math.Exp((r - q) * T0);
            

            var undiscountedCallPrice = HestonCall.HestonUndiscountedCallPrice(kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: ExpectationCIRProcess(v0, kappa, theta, T0),
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r, q: q);

            Engine.Verbose = 0;



            var fscallMarkToMarket = df * forwardPrice * undiscountedCallPrice;

            /*
            var fsCallDelta = Math.Exp(-q * T0) * c;

            var fsCallGamma = 0.0;

            var fsCallTheta = HestonNumericalGreeks.ThetaFSCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q);


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
            */
            var fsCallDelta = HestonForwardApproximated.FSCallCalculateDelta(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);

            var fsCallGamma = HestonForwardApproximated.FSCallCalculateGamma(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q);

            var fsCallTheta = HestonForwardApproximated.FSCallCalculateTheta(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);


            var fsCallVega = HestonForwardApproximated.FSCallCalculateVega(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);

            var rhoGreek = HestonForwardApproximated.FSCallCalculateRho(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);


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


        public static double FSPCallCalculateDelta(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null)
        {
            return 0.0;
        }

        public static double FSPCallCalculateGamma(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null)
        {
            return 0.0;
        }

        public static double FSPCallCalculateTheta(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null,double? paymentTime = null)
        {
            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);
            var discountingRate = -System.Math.Log(discountingFactor) / Tp;
            var relativeBump = 0.01/100;
            var derivativeDiscountFactor = (discountingFactorFunction(0, Tp + Tp * relativeBump) - discountingFactor) / (Tp * relativeBump);
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);
            var undiscountedCall = HestonCall.HestonUndiscountedCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );
            var hestonVega = HestonVega.VegaUndiscountedCall(
                        kappa: kappa,
                        theta: theta,
                        rho: rho,
                        v0: v_T0,
                        sigma: sigma,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q
                        );

            return - derivativeDiscountFactor * undiscountedCall - tDerivativeExpectationCIRProcess(v0, kappa, theta, T0) * hestonVega * discountingFactor;
        }

        public static double FSPCallCalculateVega(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);

            return discountingFactor * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaUndiscountedCall(
                        kappa: kappa,
                        theta: theta,                        
                        rho: rho,
                        v0: ExpectationCIRProcess(v0, kappa, theta, T0),
                        sigma: sigma,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q
                        );
        }

        public static double FSPCallCalculateRho(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega
            
            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);
            var derivativeDiscountFactor = - Tp * discountingFactor;
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);
            var undiscountedCall = HestonCall.HestonUndiscountedCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );
            var hestonRho = HestonRho.RhoUndiscountedCall(
                        kappa: kappa,
                        theta: theta,
                        rho: rho,
                        v0: v_T0,
                        sigma: sigma,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r, q: q
                        );

            return derivativeDiscountFactor * undiscountedCall + discountingFactor * hestonRho;
        }

        public static GreeksDerivatives HestonForwardPercentageCallWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
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

            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;

            var c = HestonCall.HestonUndiscountedCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );

            var fscallMarkToMarket = s0 * Math.Exp(-q * T0) * c;

            var fsCallDelta = 0.0;

            var fsCallGamma = 0.0;

            var fsCallTheta = HestonForwardApproximated.FSPCallCalculateTheta(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, discountingFactorFunction: discountingFactorFunction, paymentTime: paymentTime);


            var fsCallVega = HestonForwardApproximated.FSPCallCalculateVega(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, discountingFactorFunction: discountingFactorFunction, paymentTime: paymentTime);

            var rhoGreek = HestonForwardApproximated.FSPCallCalculateRho(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, discountingFactorFunction: discountingFactorFunction, paymentTime: paymentTime);
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


        public static double FSPutCalculateDelta(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = timeToPaymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);

            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            var put = HestonCall.HestonUndiscountedPutPrice(
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

            return Math.Exp((r-q) * T0) * put * df;
        }

        public static double FSPutCalculateGamma(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0)
        {
            return 0.0;
        }

        public static double FSPutCalculateTheta(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            return HestonNumericalGreeks.ThetaFSPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);
        }

        public static double FSPutCalculateVega(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = timeToPaymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);
            var forwardPrice = s0 * Math.Exp((r - q) * T0);

            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);
            

            return df * forwardPrice * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaUndiscountedPut(
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
        }

        public static double FSPutCalculateRho(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = timeToPaymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);
            var forwardPrice = s0 * Math.Exp((r - q) * T0);
            var derivativeFwdPrice = forwardPrice * T0;
            var derivativeDf = -Tp * df;
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            var undiscountedRho = HestonRho.RhoUndiscountedPut(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: v_T0,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r, q: q);


            var undiscountedPut = HestonCall.HestonUndiscountedPutPrice(kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: ExpectationCIRProcess(v0, kappa, theta, T0),
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r, q: q);


            var rhoGreek = derivativeDf * forwardPrice * undiscountedPut;
            rhoGreek += df * derivativeFwdPrice * undiscountedPut;
            rhoGreek += df * forwardPrice * undiscountedRho;

            return rhoGreek;
        }

        public static GreeksDerivatives HestonForwardPutWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
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

            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = timeToPaymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);

            var put = HestonCall.HestonUndiscountedPutPrice(
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

            var fwdPrice = s0 * Math.Exp((r - q) * T0);
            var fsPutMarkToMarket = df * fwdPrice * put;

            var fsPutDelta = HestonForwardApproximated.FSPutCalculateDelta(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);

            var fsPutGamma = HestonForwardApproximated.FSPutCalculateGamma(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q);

            var fsPutTheta = HestonForwardApproximated.FSPutCalculateTheta(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);


            var fsPutVega = HestonForwardApproximated.FSPutCalculateVega(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);

            var rhoGreek = HestonForwardApproximated.FSPutCalculateRho(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, timeToPaymentDate: timeToPaymentDate, discountingFactorFunction: discountingFactorFunction);


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

        public static double FSPPutCalculateDelta(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null)
        {
            return 0.0;
        }

        public static double FSPPutCalculateGamma(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null)
        {
            return 0.0;
        }

        public static double FSPPutCalculateTheta(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);
            var discountingRate = -System.Math.Log(discountingFactor) / Tp;
            var relativeBump = 0.01 / 100;
            var derivativeDiscountFactor = (discountingFactorFunction(0, Tp + Tp * relativeBump) - discountingFactor) / (Tp * relativeBump);
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);
            var undiscountedPut = HestonCall.HestonUndiscountedPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );
            var hestonVega = HestonVega.VegaUndiscountedPut(
                        kappa: kappa,
                        theta: theta,
                        rho: rho,
                        v0: v_T0,
                        sigma: sigma,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q
                        );

            return -derivativeDiscountFactor * undiscountedPut - tDerivativeExpectationCIRProcess(v0, kappa, theta, T0) * hestonVega * discountingFactor;
        }

        public static double FSPPutCalculateVega(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);

            return discountingFactor * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaUndiscountedPut(
                        kappa: kappa,
                        theta: theta,
                        rho: rho,
                        v0: ExpectationCIRProcess(v0, kappa, theta, T0),
                        sigma: sigma,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q
                        );
        }

        public static double FSPPutCalculateRho(double s0, double K, double T, double T0, double r, double q, double kappa, double theta, double sigma, double rho, double v0, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);
            var derivativeDiscountFactor = -Tp * discountingFactor;
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);
            var undiscountedCall = HestonCall.HestonUndiscountedPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );
            var hestonRho = HestonRho.RhoUndiscountedPut(
                        kappa: kappa,
                        theta: theta,
                        rho: rho,
                        v0: v_T0,
                        sigma: sigma,
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r, q: q
                        );

            return derivativeDiscountFactor * undiscountedCall + discountingFactor * hestonRho;
        }


        public static double FSDPutCalculateDelta(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            return 0.0;
        }

        public static double FSDPutCalculateGamma(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            return 0.0;
        }

        public static double FSDPutCalculateTheta(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? paymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            return HestonNumericalGreeks.ThetaFSDPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, paymentDate:paymentDate, discountingFactorFunction: discountingFactorFunction);
        }

        public static double FSDPutCalculateVega(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? paymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            if (discountingFactorFunction==null)
            {
                 discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);

            return df * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * -1.0 * HestonVega.VegaUndiscountedDigitalCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: ExpectationCIRProcess(v0, kappa, theta, T0),
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);
        }

        public static double FSDPutCalculateRho(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? paymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);
            var derivativeDiscountingFactor = -Tp * df;

            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            var res = df * HestonRho.RhoUndiscountedPut(
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

            var undiscountedPutPrice = 1.0 - HestonDigital.UndiscountedHestonDigitalCallPrice(
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

            res += derivativeDiscountingFactor * undiscountedPutPrice;

            return res;
        }

        public static GreeksDerivatives HestonForwardDigitalPutWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? paymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price and greeks of a FS Digital Put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\ttimeToPaymentDate\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{paymentDate}\t{r}\t{q}");
            }

            Engine.Verbose = 0;

            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentDate ?? T;

            var df = discountingFactorFunction(0, Tp); 

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            var undiscountedPutOption = 1.0 - HestonDigital.UndiscountedHestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );


            var fsDPutMarkToMarket = df * undiscountedPutOption;

            var fsDPutDelta = 0.0;

            var fsDPutGamma = 0.0;

            // theta 
            var fsDPutTheta = HestonNumericalGreeks.ThetaFSDPut(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, discountingFactorFunction:discountingFactorFunction, paymentDate:paymentDate);

            // vega 
            var fsDPutVega = - df * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaUndiscountedDigitalCall(
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

            // rho 
            var fsDPutRho = df * HestonRho.RhoUndiscountedDigitalPut(
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

            fsDPutRho += -Tp * df * undiscountedPutOption;



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

        public static GreeksDerivatives HestonForwardPercentagePutWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
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

            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            Engine.Verbose = 0;

            var c = HestonCall.HestonUndiscountedPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );

            var fscallMarkToMarket = s0 * Math.Exp(-q * T0) * c;

            var fsCallDelta = 0.0;

            var fsCallGamma = 0.0;

            var fsCallTheta = HestonForwardApproximated.FSPPutCalculateTheta(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, discountingFactorFunction: discountingFactorFunction, paymentTime:Tp);


            var fsCallVega = HestonForwardApproximated.FSPPutCalculateVega(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, discountingFactorFunction: discountingFactorFunction, paymentTime: Tp);

            var rhoGreek = HestonForwardApproximated.FSPPutCalculateRho(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, discountingFactorFunction: discountingFactorFunction, paymentTime: Tp);
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
        public static double FSDCallCalculateDelta(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            return 0.0;
        }

        public static double FSDCallCalculateGamma(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0)
        {
            return 0.0;
        }

        public static double FSDCallCalculateTheta(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? paymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            return HestonNumericalGreeks.ThetaFSDCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, paymentDate:paymentDate, discountingFactorFunction:discountingFactorFunction);

        }

        public static double FSDCallCalculateVega(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? paymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            var Tp = paymentDate ?? T;

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            var df = discountingFactorFunction(0, Tp);

            return df * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaUndiscountedDigitalCall(
                        kappa: kappa,
                        theta: theta,
                        sigma: sigma,
                        rho: rho,
                        v0: ExpectationCIRProcess(v0, kappa, theta, T0),
                        s0: 1.0,
                        T: T - T0,
                        K: K,
                        r: r,
                        q: q);
        }

        public static double FSDCallCalculateRho(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? paymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            if (discountingFactorFunction == null)
            {
               // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);
            var derivativeDiscountingFactor = -Tp * df;

            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            var res = df * HestonRho.RhoUndiscountedDigitalCall(
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

            res += derivativeDiscountingFactor * HestonDigital.UndiscountedHestonDigitalCallPrice(
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

            return res; 

        }

        public static GreeksDerivatives HestonForwardDigitalCallWithGreeks(double kappa, double theta, double rho, double v0, double sigma, double s0, double K, double r, double q, double T, double T0, double? paymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating price and greeks of a FS Digital Call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\ttimeToPaymentDate\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{paymentDate}\t{r}\t{q}");
            }

            Engine.Verbose = 0;

            if (discountingFactorFunction == null)
            { 
               // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentDate ?? T;
            var df = discountingFactorFunction(0, Tp);

            // v follows a CIR process so we take its expectations 
            var v_T0 = ExpectationCIRProcess(v0, kappa, theta, T0);

            var undiscountedDCall = HestonDigital.UndiscountedHestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v_T0,
                s0: 1.0,
                T: T - T0,
                K: K,
                r: r,
                q: q
                );


            var fsDCallMarkToMarket = df * undiscountedDCall;

            var fsCallDelta = 0.0;

            var fsDCallGamma = 0.0;

            var fsDCallTheta = HestonNumericalGreeks.ThetaFSDCall(kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, T: T, T0: T0, K: K, r: r, q: q, discountingFactorFunction:discountingFactorFunction, paymentDate:paymentDate);

            var fsDCallVega = df * derivativeExpectationCIRProcess(v0, kappa, theta, T0) * HestonVega.VegaUndiscountedDigitalCall(
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

            var fsDCallRho = df * HestonRho.RhoUndiscountedDigitalCall(
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

            fsDCallRho += -Tp * df * undiscountedDCall;

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
