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
            return  HestonDigitalCallPrice();
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
        public double HestonDigitalCallPrice(Vector x, double s0, double T, double K, double r, double q)
        {
            this.kappa = x[0];
            this.theta = x[1];
            this.sigma = x[2];
            this.rho = x[3];
            this.v0 = x[4];
            this.s0 = s0;
            this.T = T;
            this.K = K;
            this.rate = r;
            this.dividend = q;

            return HestonDigitalCallPrice();
        }

        /// <summary>
        /// Calculates the price of a digital call option using the Heston model.
        /// </summary>
        /// <returns>The price of the digital call option.</returns>
        public double HestonDigitalCallPrice()
        {
            return DiscountFactor() * UndiscountedHestonDigitalCallPrice();
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
            return HestonDigitalPutPrice();
        }

        /// <summary>
        /// Calculates the price of a digital put option using the Heston model.
        /// </summary>
        /// <returns>The price of the digital put option.</returns>
        public double HestonDigitalPutPrice()
        {
            // Price of the digital put is discountFactor * P(St<K) = discountFactor * (1 - P(St>K))
            var undiscountedPutPrice = 1.0 - UndiscountedHestonDigitalCallPrice();
            return DiscountFactor() * undiscountedPutPrice;
        }

        /// <summary>
        /// Calculates the undiscounted price of a digital call option using the Heston model.
        /// </summary>
        /// <returns>The undiscounted price of the digital call option.</returns>
        private double UndiscountedHestonDigitalCallPrice()
        {
            // the formula to price a digital call is given by 
            //  exp(-r*T) (0.5 +  1/pi * integral from 0 to infinity of Re((exp(-i*ln(K)*u)*phi(u)/(i*u))du )

            // we can reuse a some functions define in the Heston Call class 
            double a = 1E-8;
            double b = 1000.0;

            double part1 = HestonCall.PerformIntegral(a, b, DigitalIntegrandFunc);
            double integral = part1 + a * DigitalIntegrandFunc(a / 2.0);

            return (0.5 + 1 / Math.PI * integral);

        }

        /// <summary>
        /// Calculates the discount factor.
        /// </summary>
        /// <returns>The discount factor.</returns>
        private double DiscountFactor() { return Math.Exp(-this.rate * this.T); }

        /// <summary>
        /// Calculates the integrand function for the digital option.
        /// </summary>
        /// <param name="u">The integration variable.</param>
        /// <returns>The value of the integrand function at the given u.</returns>
        private double DigitalIntegrandFunc(double u)
        {
            Complex Iu = Complex.I * u;
            Complex A = Complex.Exp(-Iu * Math.Log(this.K));

            // f2 comes from the call framework 
            var complexValue = A * f2(u) / Iu;
            // return the real part of the integrand 
            return complexValue.Re;

        }





    }
}
