using System;
using HestonEstimator;
using Fairmat.Math;


namespace Heston
{
    internal class HestonDigital: HestonCall
    {
        // define the constructors from the base class 
        public HestonDigital() { }

        public HestonDigital(HestonCallOptimizationProblem problem) : base(problem) { }

        public HestonDigital(HestonProcess process, double strike, double timeToMaturity) : base(process, strike, timeToMaturity) { } 

        public double HestonDigitalCallPrice(double strike, double timeToMaturity)
        {
            this.T = timeToMaturity;
            this.K = strike;
            return  HestonDigitalCallPrice();
        }
        
        public double HestonDigitalCallPrice()
        {
            return DiscountFactor() * UndiscountedHestonDigitalCallPrice();
        }

        public double HestonDigitalPutPrice()
        {
            // Price of the digital put is discountFactor * P(St<K) = discountFactor * (1 - P(St>K))
            var undiscountedPutPrice = 1.0 - UndiscountedHestonDigitalCallPrice();
            return DiscountFactor() * undiscountedPutPrice;
        }


        private double UndiscountedHestonDigitalCallPrice()
        {
            // the formula to price a digital call is given by 
            //  exp(-r*T) (0.5 +  1/pi * integral from 0 to infinity of Re((exp(-i*ln(K)*u)*phi(u)/(i*u))du )

            // we can reuse a some functions define in the Heston Call class 
            double a = 1E-8;
            double b = 1000.0;

            double part1 = HestonCall.PerformIntegral(a, b, DigitalIntegrandFunc);
            double integral = part1 + a * DigitalIntegrandFunc(a / 2.0);

            return (0.5 *  + 1 / Math.PI * integral);

        }
        private double DiscountFactor() { return Math.Exp(-this.rate * this.T); }

        private double DigitalIntegrandFunc(double u)
        {
            Complex Iu = Complex.I * u;
            Complex A = Complex.Exp(-Iu * Math.Log(this.K));

            // f2 comes for the call framework 
            var complexValue = A * f2(u) / Iu;
            // return the real part of the integrand 
            return complexValue.Re;

        }





    }
}
