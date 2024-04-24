using System;
using Heston;
using Fairmat.Math;
using DVPLI;

namespace HestonEstimator
{
    public class HestonDigital: HestonCall
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonDigital"/> class.
        /// </summary>
        public HestonDigital() { }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonDigital"/> class.
        /// </summary>
        /// <param name='problem'>
        /// HestonCallOptimizationProblem at which digital price calculations are to be linked.
        /// </param>
        public HestonDigital(HestonCallOptimizationProblem problem) : base(problem) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonDigital"/> class.
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
        public HestonDigital(HestonProcess process, double strike, double timeToMaturity) : base(process, strike, timeToMaturity) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="HestonEstimator.HestonDigital"/> class.
        /// </summary>
        /// <param name='process'>
        /// HestonProcess at which digital price calculations are to be linked.
        /// </param>
        public HestonDigital(HestonProcess process) : base(process) { }


        /// <summary>
        /// Calculates the price of a digital call option using the Heston model.
        /// </summary>
        /// <param name="strike">The strike price of the option.</param>
        /// <param name="timeToMaturity">The time to maturity of the option.</param>
        /// <returns>The price of the digital call option.</returns>
        public double HestonDigitalCallPrice(double strike, double timeToMaturity)
        {
            this.T = timeToMaturity;
            this.K = strike;

            return  HestonDigitalCallPrice(
                kappa: this.kappa,
                theta: this.theta,
                sigma: this.sigma,
                rho: this.rho,
                v0: this.v0,
                s0: this.s0,
                T: this.T,
                K: this.K,
                r: this.rate,
                q: this.dividend);
        }

        /// <summary>
        /// Calculates the price of a digital call option using the Heston model.
        /// </summary>
        /// <param name="x">A vector containing the model parameters [kappa, theta, sigma, rho, v0].</param>
        /// <param name="s0">The initial stock price.</param>
        /// <param name="T">Time to maturity of the option.</param>
        /// <param name="K">The strike price of the option.</param>
        /// <param name="r">The risk-free interest rate.</param>
        /// <param name="q">The dividend yield.</param>
        /// <returns>The price of the digital call option.</returns>
        public static double HestonDigitalCallPrice(Vector x, double s0, double T, double K, double r, double q)
        {
            return HestonDigitalCallPrice(
                kappa: x[0],
                theta: x[1],
                rho: x[3],
                v0: x[4],
                sigma: x[2],
                s0: s0,
                T: T,
                K: K,
                r: r,
                q: q
                );
        }

        /// <summary>
        /// Calculates the price of a digital put option using the Heston model.
        /// </summary>
        /// <param name="x">A vector containing the model parameters [kappa, theta, sigma, rho, v0].</param>
        /// <param name="s0">The initial stock price.</param>
        /// <param name="T">Time to maturity of the option.</param>
        /// <param name="K">The strike price of the option.</param>
        /// <param name="r">The risk-free interest rate.</param>
        /// <param name="q">The dividend yield.</param>
        /// <returns>The price of the digital put option.</returns>
        public static double HestonDigitalPutPrice(Vector x, double s0, double T, double K, double r, double q)
        {
            return HestonDigitalPutPrice(
                kappa: x[0],
                theta: x[1],
                rho: x[3],
                v0: x[4],
                sigma: x[2],
                s0: s0,
                T: T,
                K: K,
                r: r,
                q: q
                );
        }

        /// <summary>
        /// Calculates the price of a digital call option using the Heston model.
        /// </summary>
        /// <returns>The price of the digital call option.</returns>
        public static double HestonDigitalCallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            var undiscountedCallPrice = UndiscountedHestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                rho: rho,
                v0: v0,
                sigma: sigma,
                s0: s0,
                T: T,
                K: K,
                r: r,
                q: q
                );

            return DiscountFactor(rate: r, T: T) * undiscountedCallPrice;
        }

        /// <summary>
        /// Calculates the price of a digital put option using the Heston model.
        /// </summary>
        /// <param name="strike">The strike price of the option.</param>
        /// <param name="timeToMaturity">The time to maturity of the option.</param>
        /// <returns>The price of the digital put option.</returns>
        public double HestonDigitalPutPrice(double strike, double timeToMaturity)
        {
            this.T = timeToMaturity;
            this.K = strike;

            return HestonDigitalPutPrice(
                kappa: this.kappa,
                theta: this.theta,
                sigma: this.sigma,
                rho: this.rho,
                v0: this.v0,
                s0: this.s0,
                T: this.T,
                K: this.K,
                r: this.rate,
                q: this.dividend);
        }

        /// <summary>
        /// Calculates the price of a digital put option using the Heston model.
        /// </summary>
        /// <returns>The price of the digital put option.</returns>
        public static double HestonDigitalPutPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            // Price of the digital put is discountFactor * P(St<K) = discountFactor * (1 - P(St>K))
            var undiscountedCallPrice = UndiscountedHestonDigitalCallPrice(
                kappa:kappa, 
                theta:theta,
                rho:rho,
                v0:v0,
                sigma:sigma,
                s0:s0,
                T:T,
                K:K,
                r:r,
                q:q
                );

            var undiscountedPutPrice = 1.0 - undiscountedCallPrice;
            return DiscountFactor(rate: r, T:T) * undiscountedPutPrice;
        }

        /// <summary>
        /// Calculates the undiscounted price of a digital call option using the Heston model.
        /// </summary>
        /// <returns>The undiscounted price of the digital call option.</returns>
        private static double UndiscountedHestonDigitalCallPrice(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            // the formula to price a digital call is given by 
            //  exp(-r*T) (0.5 +  1/pi * integral from 0 to infinity of Re((exp(-i*ln(K)*u)*phi(u)/(i*u))du )

            // we can reuse a some functions define in the Heston Call class 
            double a = 1E-8;
            double b = 1000.0;

            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandFunc(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K);

            double part1 = PerformIntegral(a, b, functionToIntegrate);
            double integral = part1 + a * functionToIntegrate(a / 2.0);

            return (0.5 + 1 / Math.PI * integral);

        }

        /// <summary>
        /// Calculates the discount factor.
        /// </summary>
        /// <returns>The discount factor.</returns>
        private static double DiscountFactor(double rate, double T ) { return Math.Exp(-rate * T); }

        /// <summary>
        /// Calculates the integrand function for the digital option.
        /// </summary>
        /// <param name="u">The integration variable.</param>
        /// <returns>The value of the integrand function at the given u.</returns>
        private double IntegrandFunc(double u)
        {
            return IntegrandFunc(
               u: u,
               kappa: this.kappa,
               theta: this.theta,
               sigma: this.sigma,
               rho: this.rho,
               s0: this.s0,
               v0: this.v0,
               r: this.rate,
               q: this.dividend,
               T: this.T,
               K: this.K);
        }

        public static double IntegrandFunc(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            Complex Iu = Complex.I * u;
            Complex A = Complex.Exp(-Iu * Math.Log(K));
            Complex Cu = new Complex(u);
 
            var f2 = Phi(
               u: Cu,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r - q,
               T: T
               );


            Complex complexVal = A * f2 / Iu;
            return complexVal.Re;

        }





    }
}
