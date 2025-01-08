using System;
using System.Runtime.ExceptionServices;
using DVPLDOM;
using DVPLI;
using Fairmat.Math;
using Heston;
using TAEvaluator;

namespace HestonEstimator
{
    // class that calculates the Heston greeks numerically
    public class HestonNumericalGreeks
    {
        /// <summary>
        /// Calculates the numerical approximation of the derivative of a given function at a specified point using the finite difference method, considering a percentage bump around the point.
        /// </summary>
        /// <param name="bumpPercentage">
        /// The percentage by which to bump the initial value for finite differencing.
        /// </param>
        /// <param name="init">
        /// The initial value at which to compute the derivative.
        /// </param>
        /// <param name="fun">
        /// The function for which to compute the derivative.
        /// </param>
        /// <returns>
        /// The numerical approximation of the derivative.
        /// </returns>
        public static double GreeksBumper(double bumpPercentage, double init, Func<double, double> fun, double? unbumpedPrice = null)
        {
            var deltaX = bumpPercentage * init;

            if (unbumpedPrice.HasValue)
            {
                var bumpedPrice = fun(init + deltaX);
                return (bumpedPrice - unbumpedPrice.Value) / (deltaX);
            }
            else
                return (fun(init + deltaX) - fun(init - deltaX)) / (2 * deltaX);
        }

        #region call options numerical greeks

        public static (double, double) DeltaGammaCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, double? Tp = null, Func<double, double, double> discountingFactorFunction = null)
        {
            if (discountingFactorFunction == null)
            {
               discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonCall.HestonCallPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    q: q,
                    timeToPaymentDate: Tp, 
                    discountingFactorFunction: discountingFactorFunction
                    );
            }

            double deltaS = bumpPercentage * s0;

            double callPriceBumpUp = HestonCall.HestonCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                r: r,
                q: q, 
                timeToPaymentDate: Tp, 
                discountingFactorFunction: discountingFactorFunction
                );

            double callPriceBumpDown = HestonCall.HestonCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: Tp, 
                discountingFactorFunction:discountingFactorFunction
                );

            var delta = (callPriceBumpUp - callPriceBumpDown) / (2 * deltaS);
            var gamma = (callPriceBumpUp - 2 * unBumpedPrice.Value + callPriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }

        public static double VegaCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            int verbosity = Engine.Verbose;

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T, timeToPaymentDate, r, q);
            }

            
            Func<double, double> callPrice = (double initialVariance) => HestonCall.HestonCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: timeToPaymentDate, 
                discountingFactorFunction: discountingFactorFunction
                );

            Engine.Verbose = 0;
            var vega = GreeksBumper(bumpPercentage, v0, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vega);
            }
            return vega;
        }

        public static double RhoCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? Tp = null, Func<double, double, double> discountingFactorFunction = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, Tp, r, q);
            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            double DiscountingFunctionAsFunctionOfIR(double IR, double t, double T)
            {
                Func<double, double, double> discountingRate = (t, T) => -System.Math.Log(discountingFactorFunction(t, T)) / (T - t);
                var spread = discountingRate(t, T) - r;
                // r = discountingRate(t,T) - spread
                // discountingRate(t,T) = r + spread
                return System.Math.Exp(-(IR + spread) * (T - t));
            }

            Func<double, double> callPrice = (double interestRate) => HestonCall.HestonCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                K: K,
                r: interestRate,
                q: q,
                timeToPaymentDate:Tp, 
                discountingFactorFunction: (t,T) => DiscountingFunctionAsFunctionOfIR(interestRate, t, T)
                );


            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;
        }

        public static double ThetaCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T, paymentTime, r, q);
            }

            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);
            var s = Tp - T;

            Func<double, double> callPrice = (double timeToMat) => HestonCall.HestonCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: timeToMat,
                K: K,
                r: r,
                q: q,
                discountingFactorFunction: discountingFactorFunction,
                timeToPaymentDate: s + timeToMat

                );

            Engine.Verbose = 0;
            var thetaGreek = (-1.0)*GreeksBumper(bumpPercentage, T, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Theta");
                Console.WriteLine(thetaGreek);
            }

            return thetaGreek;
        }

        #endregion 

        #region put options numerical greeks

        public static (double, double) DeltaGammaPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {

            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonCall.HestonPutPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    q: q,
                    timeToPaymentDate: paymentTime,
                    discountingFactorFunction: discountingFactorFunction
                    );
            }

            double deltaS = bumpPercentage * s0;

            double putPriceBumpUp = HestonCall.HestonPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: paymentTime,
                discountingFactorFunction: discountingFactorFunction);

            double putPriceBumpDown = HestonCall.HestonPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: paymentTime,
                discountingFactorFunction: discountingFactorFunction);

            var delta = (putPriceBumpUp - putPriceBumpDown) / (2 * deltaS);
            var gamma = (putPriceBumpUp - 2 * unBumpedPrice.Value + putPriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }

        public static double VegaPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);
            }

            Func<double, double> price = (double initialVariance) => HestonCall.HestonPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                K: K,
                r: r,
                q: q);


            Engine.Verbose = 0;
            var vega = GreeksBumper(bumpPercentage, v0, price);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vega);
            }
            
            return vega;
        }

        public static double RhoPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);
            }

            double DiscountingFunctionAsFunctionOfIR(double IR, double t, double T)
            {
                Func<double, double, double> discountingRate = (t, T) => -System.Math.Log(discountingFactorFunction(t, T)) / (T - t);
                var spread = discountingRate(t, T) - r;
                // r = discountingRate(t,T) - spread
                // discountingRate(t,T) = r + spread
                return System.Math.Exp(-(IR + spread) * (T - t));
            }

            Func<double, double> price = (double interestRate) => HestonCall.HestonPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                K: K,
                r: interestRate,
                q: q,
                timeToPaymentDate: paymentTime,
                discountingFactorFunction: (t, T) => DiscountingFunctionAsFunctionOfIR(interestRate, t, T));

            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, price);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;
        }

        public static double ThetaPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? timeToPaymentDate = null, Func<double, double, double> discountingFactorFunction = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T,timeToPaymentDate, r, q);

            }
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }
            // derive the formula for the vega

            var Tp = timeToPaymentDate ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);
            var s = Tp - T;

            Func<double, double> price = (double timeToMat) => HestonCall.HestonPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: timeToMat,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: s + timeToMat,
                discountingFactorFunction: discountingFactorFunction
                );


            Engine.Verbose = 0;
            var thetaGreek = (-1.0)*GreeksBumper(bumpPercentage, T, price);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
                Console.WriteLine("Theta");
                Console.WriteLine(thetaGreek);

            return thetaGreek;
        }

        #endregion

        #region digital options numerical greeks

        public static (double, double) DeltaGammaDCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, double? timeToPaymentDate = null)
        {

            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonDigital.HestonDigitalCallPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    q: q,
                    timeToPaymentDate: timeToPaymentDate
                    );
            }

            double deltaS = bumpPercentage * s0;

            double PriceBumpUp = HestonDigital.HestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );

            double PriceBumpDown = HestonDigital.HestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );

            var delta = (PriceBumpUp - PriceBumpDown) / (2 * deltaS);
            var gamma = (PriceBumpUp - 2 * unBumpedPrice.Value + PriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }

        public static (double, double) DeltaGammaDPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, double? timeToPaymentDate = null)
        {
            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonDigital.HestonDigitalPutPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    q: q,
                    timeToPaymentDate:timeToPaymentDate);
            }

            double deltaS = bumpPercentage * s0;

            double PriceBumpUp = HestonDigital.HestonDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );

            double PriceBumpDown = HestonDigital.HestonDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );

            var delta = (PriceBumpUp - PriceBumpDown) / (2 * deltaS);
            var gamma = (PriceBumpUp - 2 * unBumpedPrice.Value + PriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }

        public static double VegaDCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? timeToPaymentDate = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);
            }

            Func<double, double> price = (double initialVariance) => HestonDigital.HestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );

            Engine.Verbose = 0;
            var vega = GreeksBumper(bumpPercentage, v0, price);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vega);
            }

            return vega;
        }

        public static double VegaDPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? timeToPaymentDate = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T,timeToPaymentDate, r, q);
            }

            Func<double, double> price = (double initialVariance) => HestonDigital.HestonDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );

            Engine.Verbose = 0;
            var vega = GreeksBumper(bumpPercentage, v0, price);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vega);
            }

            return vega;
        }

        public static double RhoDCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? timeToPaymentDate = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T,timeToPaymentDate, r, q);
            }

            Func<double, double> price = (double interestRate) => HestonDigital.HestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                K: K,
                r: interestRate,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );

            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, price);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;
        }

        public static double RhoDPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? timeToPaymentDate = null)
        {
            int verbosity = Engine.Verbose;

            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T,timeToPaymentDate, r, q);
            }

            Func<double, double> price = (double interestRate) => HestonDigital.HestonDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                K: K,
                r: interestRate,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );

            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, price);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;
        }

        public static double ThetaDCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? timeToPaymentDate = null)
        {

            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T,timeToPaymentDate, r, q);
            }

            Func<double, double> price = (double timeToMat) => HestonDigital.HestonDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: timeToMat,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );


            Engine.Verbose = 0;
            var thetaGreek = (- 1.0) * GreeksBumper(bumpPercentage, T, price);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Theta");
                Console.WriteLine(thetaGreek);
            }

            return thetaGreek;
        }

        public static double ThetaDPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double bumpPercentage = 0.01, double? timeToPaymentDate = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T,timeToPaymentDate, r, q);
            }

            Func<double, double> price = (double timeToMat) => HestonDigital.HestonDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: timeToMat,
                K: K,
                r: r,
                q: q,
                timeToPaymentDate: timeToPaymentDate
                );

            Engine.Verbose = 0;
            var thetaGreek = (-1.0) * GreeksBumper(bumpPercentage, T, price);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Theta");
                Console.WriteLine(thetaGreek);
            }

            return thetaGreek;
        }

        #endregion

        #region forward starting options numerical greeks

        public static (double, double) DeltaGammaFSPCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double T0, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, Func<double, double, double> discountingFactorFunction = null, double? timeToPaymentDate = null)
        {

            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonForwardApproximated.HestonForwardPercentageCallPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    T0: T0,
                    q: q,
                    paymentTime: timeToPaymentDate,
                    discountingFactorFunction: discountingFactorFunction
                    );
            }

            double deltaS = bumpPercentage * s0;

            double callPriceBumpUp = HestonForwardApproximated.HestonForwardPercentageCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                T0: T0,
                r: r,
                q: q,
                paymentTime: timeToPaymentDate,
                discountingFactorFunction: discountingFactorFunction
                );

            double callPriceBumpDown = HestonForwardApproximated.HestonForwardPercentageCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                T0: T0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                paymentTime: timeToPaymentDate,
                discountingFactorFunction: discountingFactorFunction
                );

            var delta = (callPriceBumpUp - callPriceBumpDown) / (2 * deltaS);
            var gamma = (callPriceBumpUp - 2 * unBumpedPrice.Value + callPriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }


        public static double VegaFSPCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, Func<double, double, double> discountingFactorFunction = null, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, Tp, r, q);
            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            Func<double, double> callPrice = (double initialVariance) => HestonForwardApproximated.HestonForwardPercentageCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: r,
                q: q,
                paymentTime: Tp,
                discountingFactorFunction: discountingFactorFunction
                );

            Engine.Verbose = 0;
            var vegaGreek = GreeksBumper(bumpPercentage, v0, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vegaGreek);
            }

            return vegaGreek;
        }

        public static double ThetaFSPCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, Func<double, double, double> discountingFactorFunction = null, double? Tp = null)
        {

            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T, Tp, r, q);
            }


            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            if (Tp is null)
            {

                var M = T - T0;
                Func<double, double> callPrice = (double strikeDate) => HestonForwardApproximated.HestonForwardPercentageCallPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: strikeDate + M,
                    T0: strikeDate,
                    K: K,
                    r: r,
                    q: q,
                    discountingFactorFunction: discountingFactorFunction
                    );

                Engine.Verbose = 0;
                var thetaGreek = -GreeksBumper(bumpPercentage, T0, callPrice);
                Engine.Verbose = verbosity;

                if (verbosity > 0)
                {
                    Console.WriteLine("Theta");
                    Console.WriteLine(thetaGreek);
                }

                return thetaGreek;

            }
            else
            {
                var TminusT0 = T - T0;
                var TpminusT0 = Tp.Value - T0;

                Func<double, double> callPrice = (double strikeDate) => HestonForwardApproximated.HestonForwardPercentageCallPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: strikeDate + TminusT0,
                    T0: strikeDate,
                    K: K,
                    r: r,
                    q: q,
                    discountingFactorFunction: discountingFactorFunction,
                    paymentTime: strikeDate + TpminusT0
                    );

                Engine.Verbose = 0;
                var thetaGreek = -GreeksBumper(bumpPercentage, T0, callPrice);
                Engine.Verbose = verbosity;

                if (verbosity > 0)
                {
                    Console.WriteLine("Theta");
                    Console.WriteLine(thetaGreek);
                }

                return thetaGreek;


            }


        }

        public static double RhoFSPCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.0001, Func<double, double, double> discountingFactorFunction = null, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T, Tp, r, q);
            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }


            double DiscountingFunctionAsFunctionOfIR(double IR, double t, double T)
            {
                Func<double, double, double> discountingRate = (t, T) => -System.Math.Log(discountingFactorFunction(t, T)) / (T - t);
                var spread = discountingRate(t, T) - r;
                // r = discountingRate(t,T) - spread
                // discountingRate(t,T) = r + spread
                return System.Math.Exp(-(IR + spread) * (T - t));
            }


            Func<double, double> callPrice = (double interestRate) => HestonForwardApproximated.HestonForwardPercentageCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: interestRate,
                q: q,
                paymentTime: Tp,
                discountingFactorFunction: (x, y) => DiscountingFunctionAsFunctionOfIR(interestRate, x, y)
                );

            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;
        }

        public static (double, double) DeltaGammaFSCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double T0, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, double? timeToPaymentDate = null)
        {

            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonForwardApproximated.HestonForwardCallPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    T0:T0,
                    q: q,
                    Tp: timeToPaymentDate
                    );
            }

            double deltaS = bumpPercentage * s0;

            double callPriceBumpUp = HestonForwardApproximated.HestonForwardCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                T0:T0,
                r: r,
                q: q,
                Tp: timeToPaymentDate
                );

            double callPriceBumpDown = HestonForwardApproximated.HestonForwardCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                T0:T0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                Tp: timeToPaymentDate
                );

            var delta = (callPriceBumpUp - callPriceBumpDown) / (2 * deltaS);
            var gamma = (callPriceBumpUp - 2 * unBumpedPrice.Value + callPriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }

        public static double VegaFSCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, Tp, r, q);
            }

            Func<double, double> callPrice = (double initialVariance) => HestonForwardApproximated.HestonForwardCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                T0:T0,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            Engine.Verbose = 0;
            var vegaGreek = GreeksBumper(bumpPercentage, v0, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vegaGreek);
            }

            return vegaGreek;
        }

        public static double ThetaFSCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T,Tp, r, q);
            }

            var M = T - T0;
            Func<double, double> callPrice = (double strikeDate) => HestonForwardApproximated.HestonForwardCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: strikeDate + M,
                T0: strikeDate,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            Engine.Verbose = 0;
            var thetaGreek = - GreeksBumper(bumpPercentage, T0, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Theta");
                Console.WriteLine(thetaGreek);
            }

            return thetaGreek;
        }

        public static double RhoFSCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T,Tp, r, q);
            }

            Func<double, double> callPrice = (double interestRate) => HestonForwardApproximated.HestonForwardCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: interestRate,
                q: q,
                Tp: Tp
                );

            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;
        }

        public static (double, double) DeltaGammaFSPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double T0, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, double? Tp = null)
        {

            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonForwardApproximated.HestonForwardPutPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    T0: T0,
                    q: q,
                    Tp: Tp
                    );
            }

            double deltaS = bumpPercentage * s0;

            double putPriceBumpUp = HestonForwardApproximated.HestonForwardPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                T0: T0,
                r: r,
                q: q,
                Tp: Tp
                );

            double putPriceBumpDown = HestonForwardApproximated.HestonForwardPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                T0: T0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            var delta = (putPriceBumpUp - putPriceBumpDown) / (2 * deltaS);
            var gamma = (putPriceBumpUp - 2 * unBumpedPrice.Value + putPriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }

        public static double VegaFSPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T,Tp, r, q);
            }

            Func<double, double> putPrice = (double initialVariance) => HestonForwardApproximated.HestonForwardPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            Engine.Verbose = 0;
            var vegaGreek = GreeksBumper(bumpPercentage, v0, putPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vegaGreek);
            }

            return vegaGreek;

        }

        public static double ThetaFSPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T,Tp, r, q);
            }

            var M = T - T0;
            Func<double, double> putPrice = (double strikeDate) => HestonForwardApproximated.HestonForwardPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: strikeDate + M,
                T0: strikeDate,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            Engine.Verbose = 0;
            var thetaGreek = - GreeksBumper(bumpPercentage, T0, putPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Theta");
                Console.WriteLine(thetaGreek);
            }

            return thetaGreek;
        }

        public static double RhoFSPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T,Tp, r, q);
            }

            Func<double, double> callPrice = (double interestRate) => HestonForwardApproximated.HestonForwardPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: interestRate,
                q: q,
                Tp: Tp
                );

            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;
        }

        public static (double, double) DeltaGammaFSPPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double T0, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, Func<double, double, double> discountingFactorFunction = null, double? timeToPaymentDate = null)
        {

            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonForwardApproximated.HestonForwardPercentagePutPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    T0: T0,
                    q: q,
                    paymentTime: timeToPaymentDate,
                    discountingFactorFunction: discountingFactorFunction
                    );
            }

            double deltaS = bumpPercentage * s0;

            double callPriceBumpUp = HestonForwardApproximated.HestonForwardPercentagePutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                T0: T0,
                r: r,
                q: q,
                paymentTime: timeToPaymentDate,
                discountingFactorFunction: discountingFactorFunction
                );

            double callPriceBumpDown = HestonForwardApproximated.HestonForwardPercentagePutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                T0: T0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                paymentTime: timeToPaymentDate,
                discountingFactorFunction: discountingFactorFunction
                );

            var delta = (callPriceBumpUp - callPriceBumpDown) / (2 * deltaS);
            var gamma = (callPriceBumpUp - 2 * unBumpedPrice.Value + callPriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }


        public static double VegaFSPPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, Func<double, double, double> discountingFactorFunction = null, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, Tp, r, q);
            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            Func<double, double> callPrice = (double initialVariance) => HestonForwardApproximated.HestonForwardPercentagePutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: r,
                q: q,
                paymentTime: Tp,
                discountingFactorFunction: discountingFactorFunction
                );

            Engine.Verbose = 0;
            var vegaGreek = GreeksBumper(bumpPercentage, v0, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vegaGreek);
            }

            return vegaGreek;
        }

        public static double ThetaFSPPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, Func<double, double, double> discountingFactorFunction = null, double? Tp = null)
        {

            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T, Tp, r, q);
            }


            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            if (Tp is null)
            {

                var M = T - T0;
                Func<double, double> callPrice = (double strikeDate) => HestonForwardApproximated.HestonForwardPercentagePutPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: strikeDate + M,
                    T0: strikeDate,
                    K: K,
                    r: r,
                    q: q,
                    discountingFactorFunction: discountingFactorFunction
                    );

                Engine.Verbose = 0;
                var thetaGreek = -GreeksBumper(bumpPercentage, T0, callPrice);
                Engine.Verbose = verbosity;

                if (verbosity > 0)
                {
                    Console.WriteLine("Theta");
                    Console.WriteLine(thetaGreek);
                }

                return thetaGreek;

            }
            else
            {
                var TminusT0 = T - T0;
                var TpminusT0 = Tp.Value - T0;

                Func<double, double> callPrice = (double strikeDate) => HestonForwardApproximated.HestonForwardPercentagePutPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: strikeDate + TminusT0,
                    T0: strikeDate,
                    K: K,
                    r: r,
                    q: q,
                    discountingFactorFunction: discountingFactorFunction,
                    paymentTime: strikeDate + TpminusT0
                    );

                Engine.Verbose = 0;
                var thetaGreek = -GreeksBumper(bumpPercentage, T0, callPrice);
                Engine.Verbose = verbosity;

                if (verbosity > 0)
                {
                    Console.WriteLine("Theta");
                    Console.WriteLine(thetaGreek);
                }

                return thetaGreek;


            }


        }

        public static double RhoFSPPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, Func<double, double, double> discountingFactorFunction = null, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tTp\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", s0, K, T, Tp, r, q);
            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }


            double DiscountingFunctionAsFunctionOfIR(double IR, double t, double T)
            {
                Func<double, double, double> discountingRate = (t, T) => -System.Math.Log(discountingFactorFunction(t, T)) / (T - t);
                var spread = discountingRate(t, T) - r;
                // r = discountingRate(t,T) - spread
                // discountingRate(t,T) = r + spread
                return System.Math.Exp(-(IR + spread) * (T - t));
            }


            Func<double, double> callPrice = (double interestRate) => HestonForwardApproximated.HestonForwardPercentagePutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: interestRate,
                q: q,
                paymentTime: Tp,
                discountingFactorFunction: (x, y) => DiscountingFunctionAsFunctionOfIR(interestRate, x, y)
                );

            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;
        }

        public static (double, double) DeltaGammaFSDCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double T0, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, double? Tp = null)
        {

            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonForwardApproximated.HestonForwardDigitalCallPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    T0: T0,
                    q: q,
                    Tp: Tp
                    );
            }

            

            double deltaS = bumpPercentage * s0;

            double callPriceBumpUp = HestonForwardApproximated.HestonForwardDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                T0: T0,
                r: r,
                q: q,
                Tp: Tp
                );

            double callPriceBumpDown = HestonForwardApproximated.HestonForwardDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                T0: T0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            var delta = (callPriceBumpUp - callPriceBumpDown) / (2 * deltaS);
            var gamma = (callPriceBumpUp - 2 * unBumpedPrice.Value + callPriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }

        public static double VegaFSDCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tTp\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");
            }

            Func<double, double> callPrice = (double initialVariance) => HestonForwardApproximated.HestonForwardDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            Engine.Verbose = 0;
            var vega = GreeksBumper(bumpPercentage, v0, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vega);
            }

            return vega;
        }

        public static double ThetaFSDCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\t{Tp}\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");
            }

            var M = T - T0;
            Func<double, double> callPrice = (double strikeDate) => HestonForwardApproximated.HestonForwardDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: strikeDate+M,
                T0: strikeDate,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            Engine.Verbose = 0;
            var thetaGreek = - GreeksBumper(bumpPercentage, T0, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Theta");
                Console.WriteLine(thetaGreek);
            }

            return thetaGreek;
        }

        public static double RhoFSDCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.001, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tTp\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");
            }

            Func<double, double> callPrice = (double interestRate) => HestonForwardApproximated.HestonForwardDigitalCallPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: interestRate,
                q: q,
                Tp:Tp
                );

            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, callPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }
            
            return rhoGreek;
        }

        public static (double, double) DeltaGammaFSDPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double T0, double K, double r, double q, double bumpPercentage = 0.01, double? unBumpedPrice = null, double? Tp = null)
        {

            if (!unBumpedPrice.HasValue)
            {
                unBumpedPrice = HestonForwardApproximated.HestonForwardDigitalPutPrice(
                    kappa: kappa,
                    theta: theta,
                    sigma: sigma,
                    rho: rho,
                    v0: v0,
                    s0: s0,
                    T: T,
                    K: K,
                    r: r,
                    T0: T0,
                    q: q,
                    Tp: Tp);
            }

            double deltaS = bumpPercentage * s0;

            double putPriceBumpUp = HestonForwardApproximated.HestonForwardDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0 + deltaS,
                T: T,
                K: K,
                T0: T0,
                r: r,
                q: q,
                Tp: Tp
                );

            double putPriceBumpDown = HestonForwardApproximated.HestonForwardDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                T0: T0,
                s0: s0 - deltaS,
                T: T,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            var delta = (putPriceBumpUp - putPriceBumpDown) / (2 * deltaS);
            var gamma = (putPriceBumpUp - 2 * unBumpedPrice.Value + putPriceBumpDown) / (deltaS * deltaS);
            return (delta, gamma);
        }

        public static double VegaFSDPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tTp\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");
            }

            Func<double, double> putPrice = (double initialVariance) => HestonForwardApproximated.HestonForwardDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: initialVariance,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: r,
                q: q,
                Tp:Tp
                );

            Engine.Verbose = 0;
            var vegaGreek = GreeksBumper(bumpPercentage, v0, putPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vegaGreek);
            }

            return vegaGreek;
        }

        public static double ThetaFSDPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.01, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating theta with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tTp\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");
            }

            var M = T - T0;
            Func<double, double> putPrice = (double strikeDate) => HestonForwardApproximated.HestonForwardDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: strikeDate+M,
                T0: strikeDate,
                K: K,
                r: r,
                q: q,
                Tp: Tp
                );

            Engine.Verbose = 0;
            var thetaGreek = - GreeksBumper(bumpPercentage, T0, putPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Theta");
                Console.WriteLine(thetaGreek);
            }

            return thetaGreek;
        }

        public static double RhoFSDPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, double T0, double bumpPercentage = 0.001, double? Tp = null)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tT0\tTp\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{T0}\t{Tp}\t{r}\t{q}");
            }

            Func<double, double> putPrice = (double interestRate) => HestonForwardApproximated.HestonForwardDigitalPutPrice(
                kappa: kappa,
                theta: theta,
                sigma: sigma,
                rho: rho,
                v0: v0,
                s0: s0,
                T: T,
                T0: T0,
                K: K,
                r: interestRate,
                q: q,
                Tp: Tp);

            Engine.Verbose = 0;
            var rhoGreek = GreeksBumper(bumpPercentage, r, putPrice);
            Engine.Verbose = verbosity;

            if (verbosity > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;
        }

        #endregion

    }


    // class that calculates the analytical delta of options with the Heston model
    public class HestonDelta : HestonCall
    {
        public static Complex Phi(Complex u, double kappa, double theta, double sigma, double rho, double v0, double s0, double r, double T)
        {
            Complex I = Complex.I;
            var addittionalTem = I * u / s0;
            var val = HestonCall.Phi(
               u: u,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r,
               T: T
               ) * addittionalTem;
            return val;
        }

        public static double IntegrandFunc(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            var int1 = IntegrandFunc1(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T: T, K: K, r: r, q: q);
            var int2 = IntegrandFunc2(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T: T, K: K, r: r, q: q);
            return int1 - K * int2;
        }

        public static double IntegrandFunc1(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            Complex Iu = Complex.I * u;
            Complex I = Complex.I;
            Complex A = Complex.Exp(-Iu * Math.Log(K));

            var f1 = Phi(
               u: u - I,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r - q,
               T: T
               );


            Complex complexVal1 = A * f1 / Iu;
            return complexVal1.Re;

        }

        internal static Complex Phi(double u, double kappa, double theta, double sigma, double rho, double v0, double s0, double r, double T)
        {
            Complex Cu = new Complex(u);
            return Phi(Cu, kappa, theta, sigma, rho, v0, s0, r, T);
        }

        public static double IntegrandFunc2(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            Complex Iu = Complex.I * u;
            Complex A = Complex.Exp(-Iu * Math.Log(K));


            var f2 = Phi(
               u: u,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r - q,
               T: T
               );


            Complex complexVal2 = A * f2 / Iu;
            return complexVal2.Re;

        }

        public static double DeltaCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {

            if (Engine.Verbose > 0)
            {
                Console.WriteLine("Calculating delta of a call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentTime ?? T;
            double F = 1.0 * Math.Exp((r - q) * T);
            double firstTerm = 0.5 * F;
            double a = 1E-8;
            double b = 1000.0;


            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandFunc(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K);
            double part1 = PerformIntegral(a, b, functionToIntegrate);

            double integral = part1 + a * functionToIntegrate(a / 2.0);

            double delta = (firstTerm + integral / Math.PI);

            double adjustedDelta = discountingFactorFunction(0,Tp) * delta;

            if (Engine.Verbose > 0)
                Console.WriteLine("Delta: {0}", adjustedDelta);


            return adjustedDelta;
        }

        public static double DeltaPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {

            var verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating delta of a put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }
            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentTime ?? T;
            Engine.Verbose = 0;
            var deltaCall = DeltaCall(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q, paymentTime: paymentTime, discountingFactorFunction: discountingFactorFunction);
            Engine.Verbose = verbosity;


            return deltaCall  - discountingFactorFunction(0, Tp) * Math.Exp((r-q) * T);
        }

        private static double DeltaDigital(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            // the formula to price a digital call is given by 
            //  exp(-r*T) (0.5 +  1/pi * integral from 0 to infinity of Re((exp(-i*ln(K)*u)*phi(u)/(i*u))du )

            // we can reuse a some functions define in the Heston Call class 
            double a = 1E-8;
            double b = 1000.0;

            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandFunc2(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K);

            double part1 = PerformIntegral(a, b, functionToIntegrate);
            double integral = part1 + a * functionToIntegrate(a / 2.0);

            return (1 / Math.PI * integral);
        }

        public static double DeltaDigitalCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            var delta = HestonDigital.DiscountFactor(rate: r, T: T) * DeltaDigital(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
            return delta;
        }

        public static double DeltaDigitalPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            return -DeltaDigitalCall(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
        }

    }


    // class that calculates the analytical gamma of options with the Heston model
    public class HestonGamma : HestonCall
    {

        public static Complex Phi(Complex u, double kappa, double theta, double sigma, double rho, double v0, double s0, double r, double T)
        {
            Complex I = Complex.I;
            var addittionalTerm = I * u / s0;
            var additionalTerm2 = -I * u / (s0 * s0);


            var phi = HestonCall.Phi(
               u: u,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r,
               T: T
               );

            var val1 = phi * additionalTerm2;
            var val2 = addittionalTerm * addittionalTerm * phi;
            return val1 + val2;
        }

        public static double IntegrandFunc(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            var int1 = IntegrandFunc1(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T: T, K: K, r: r, q: q);
            var int2 = IntegrandFunc2(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T: T, K: K, r: r, q: q);
            return int1 - K * int2;
        }

        public static double IntegrandFunc1(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            Complex Iu = Complex.I * u;
            Complex I = Complex.I;
            Complex A = Complex.Exp(-Iu * Math.Log(K));

            var f1 = Phi(
               u: u - I,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r - q,
               T: T
               );


            Complex complexVal1 = A * f1 / Iu;
            return complexVal1.Re;

        }

        internal static Complex Phi(double u, double kappa, double theta, double sigma, double rho, double v0, double s0, double r, double T)
        {
            Complex Cu = new Complex(u);
            return Phi(Cu, kappa, theta, sigma, rho, v0, s0, r, T);
        }

        public static double IntegrandFunc2(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            Complex Iu = Complex.I * u;
            Complex A = Complex.Exp(-Iu * Math.Log(K));


            var f2 = Phi(
               u: u,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r - q,
               T: T
               );


            Complex complexVal2 = A * f2 / Iu;
            return complexVal2.Re;

        }

        public static double GammaCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {

            if (Engine.Verbose > 0)
            {
                Console.WriteLine("Calculating gamma with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }

            if (discountingFactorFunction == null)
            {
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentTime ?? T;
            double firstTerm = 0.0;
            double a = 1E-8;
            double b = 1000.0;

           
            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandFunc(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K);
            double part1 = PerformIntegral(a, b, functionToIntegrate);

            double integral = part1 + a * functionToIntegrate(a / 2.0);

            double gamma = (firstTerm + integral / Math.PI);

            double adjustedGamma = discountingFactorFunction(0,Tp) * gamma;

            if (Engine.Verbose > 0)
                Console.WriteLine("Gamma: {0}", adjustedGamma);


            return adjustedGamma;
        }

        public static double GammaPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            return GammaCall(
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
        }

        private static double GammaDigital(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            // the formula to price a digital call is given by 
            //  exp(-r*T) (0.5 +  1/pi * integral from 0 to infinity of Re((exp(-i*ln(K)*u)*phi(u)/(i*u))du )

            // we can reuse a some functions define in the Heston Call class 
            double a = 1E-8;
            double b = 1000.0;

            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandFunc2(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K);

            double part1 = PerformIntegral(a, b, functionToIntegrate);
            double integral = part1 + a * functionToIntegrate(a / 2.0);

            return (1 / Math.PI * integral);
        }

        public static double GammaDigitalCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            int verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating gamma with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine($"{kappa}\t{theta}\t{sigma}\t{rho}\t{v0}");
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine($"{s0}\t{K}\t{T}\t{r}\t{q}");
            }

            var gamma = HestonDigital.DiscountFactor(rate: r, T: T) * GammaDigital(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
            
            if (verbosity > 0)
            {
                Console.WriteLine("Gamma");
                Console.WriteLine(gamma);
            }
            return gamma;
        }

        public static double GammaDigitalPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            return - GammaDigitalCall(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
        }

    }


    // class that calculates the analytical vega of options with the Heston model
    public class HestonVega : HestonCall
    {
        internal static Complex Phi(double u, double kappa, double theta, double sigma, double rho, double v0, double s0, double r, double T)
        {
            Complex Cu = new Complex(u);
            return Phi(Cu, kappa, theta, sigma, rho, v0, s0, r, T);
        }

        /// <summary>
        /// Calculates the Heston characteristic function.
        /// </summary>
        /// <remarks>
        /// The calculation does not follow the original paper by Heston
        /// but uses the form found in the paper
        /// <br />
        /// Schoutens, W., Simons E. and Tistaert, J. (2004)
        /// A perfect calibration! Now what?, Wilmott Magazine, March 2005, 66–78.
        /// <br />
        /// As stated in:
        /// Albrecher, H., Mayer, Ph., Schoutens, W and Tistaert, J. (2007)
        /// The Little Heston Trap. Wilmott Magazine, January Issue, 83-92.
        /// <br />
        /// This form does not have the discontinuity issue that affects the form
        /// found in the original Heston paper.
        /// </remarks>
        /// <param name="u">
        /// Complex value at which the characteristic function is to be calculated.
        /// </param>
        /// <param name="kappa">Heston volatility mean reversion speed parameter.</param>
        /// <param name="theta">Heston volatility mean reversion level parameter.</param>
        /// <param name="sigma">Heston volatility of volatility parameter.</param>
        /// <param name="rho">
        /// Correlation between the two Wiener processes in the Heston dynamics.
        /// </param>
        /// <param name="s0">Starting value for the stock process.</param>
        /// <param name="v0">Starting value for the volatility process.</param>
        /// <param name="r">Risk free rate.</param>
        /// <param name="T">Time at which the characteristic function is to be calculated.</param>
        /// <returns>The value of the characteristic function.</returns>
        public static Complex Phi(Complex u, double kappa, double theta, double sigma, double rho, double v0, double s0, double r, double T)
        {
            Complex d, g, A, B, val;
            Complex I = Complex.I;
            double ss = sigma * sigma;
            Complex tmp1 = I * rho * sigma * u;
            d = Complex.Sqrt(Complex.Pow(tmp1 - kappa, 2.0) + ss * (I * u + u * u));
            Complex tmp2 = kappa - tmp1;
            Complex par = tmp2 - d;
            g = par / (tmp2 + d);

            Complex edT = Complex.Exp(-d * T);
            Complex numArg = 1.0 - g * edT;
            A = (theta * kappa) * (par * T - 2.0 * Complex.Log(numArg / (1.0 - g))) / (ss);
            B = (par * (1.0 - edT) / numArg) / ss;

            val = Complex.Exp(I * u * (Math.Log(s0) + r * T) + A + B * v0) * B;
            
            return val;
        }

        public static double IntegrandFunc(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            var int1 = IntegrandFunc1(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T: T, K: K, r: r, q: q);
            var int2 = IntegrandFunc2(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T: T, K: K, r: r, q: q);
            return int1 - int2;
        }

        public static double IntegrandFunc1(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            Complex Iu = Complex.I * u;
            Complex I = Complex.I;
            Complex A = Complex.Exp(-Iu * Math.Log(K));

            var f1 = Phi(
               u: u - I,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r - q,
               T: T
               );


            Complex complexVal1 = A * f1 / Iu;
            return complexVal1.Re;

        }

        public static double IntegrandFunc2(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            Complex Iu = Complex.I * u;
            Complex A = Complex.Exp(-Iu * Math.Log(K));


            var f2 = Phi(
               u: u,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r - q,
               T: T
               );


            Complex complexVal2 = A * f2 / Iu;
            return K * complexVal2.Re;

        }

        public static double VegaCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentTime ?? T;
            var v = VegaUndiscountedCall(kappa, theta, rho, v0, sigma, s0, T, K, r, q);

            double adjustedVega = discountingFactorFunction(0, Tp) * v;

            if (Engine.Verbose > 0)
                Console.WriteLine("Vega: {0}", adjustedVega);


            return adjustedVega;
        }

        public static double VegaUndiscountedCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            if (Engine.Verbose > 0)
            {
                Console.WriteLine("Calculating vega of a call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }

            double a = 1E-8;
            double b = 1000.0;

            TAEDelegateFunction1D functionToIntegrate = (double u) => IntegrandFunc(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K);
            double part1 = PerformIntegral(a, b, functionToIntegrate);

            double integral = part1 + a * functionToIntegrate(a / 2.0);

            double adjustedVega = integral / Math.PI;

            if (Engine.Verbose > 0)
                Console.WriteLine("Vega: {0}", adjustedVega);


            return adjustedVega;
        }

        public static double VegaPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            var verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega of a put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }

            Engine.Verbose = 0;
            var vegaCall = VegaCall(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
            Engine.Verbose = verbosity;

            return vegaCall;
            
        }

        public static double VegaUndiscountedPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            var verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating vega of a put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }

            Engine.Verbose = 0;
            var vegaPut = VegaUndiscountedCall(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
            Engine.Verbose = verbosity;

            return vegaPut;

        }

        private static double VegaDigital(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            // the formula to price a digital call is given by 
            //  exp(-r*T) (0.5 +  1/pi * integral from 0 to infinity of Re((exp(-i*ln(K)*u)*phi(u)/(i*u))du )

            // we can reuse a some functions define in the Heston Call class 
            double a = 1E-8;
            double b = 1000.0;

            // here integrand func 2 must be divided by K because it is in the formulation of the call option  
            TAEDelegateFunction1D functionToIntegrate = 
                (double u) => IntegrandFunc2(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K) / K;

            double part1 = PerformIntegral(a, b, functionToIntegrate);
            double integral = part1 + a * functionToIntegrate(a / 2.0);

            return (1 / Math.PI * integral);
        }

        public static double VegaDigitalCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            if (Engine.Verbose > 0)
            {
                Console.WriteLine("Calculating vega with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }

            var vega = HestonDigital.DiscountFactor(rate: r, T: T) * VegaDigital(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
            
            if (Engine.Verbose > 0)
            {
                Console.WriteLine("Vega");
                Console.WriteLine(vega);
            }
            
            return vega;
        }

        public static double VegaDigitalPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            return -VegaDigitalCall(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
        }


    }


    // class that calculates the analytical rho of options with the Heston model
    public class HestonRho : HestonCall
    {
        public static Complex Phi(Complex u, double kappa, double theta, double sigma, double rho, double v0, double s0, double r, double T)
        {
            Complex I = Complex.I;
            var additionalTerm = I * u * T;
            var val = HestonCall.Phi(
               u: u,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r,
               T: T
               ) * additionalTerm;
            return val;
        }

        public static double IntegrandFunc(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            var int1 = IntegrandFunc1(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T: T, K: K, r: r, q: q);
            var int2 = IntegrandFunc2(u: u, kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T: T, K: K, r: r, q: q);
            return int1 - K * int2;
        }

        public static double IntegrandFunc1(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            Complex Iu = Complex.I * u;
            Complex I = Complex.I;
            Complex A = Complex.Exp(-Iu * Math.Log(K));

            var f1 = Phi(
               u: u - I,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r - q,
               T: T
               );


            Complex complexVal1 = A * f1 / Iu;
            return complexVal1.Re;

        }

        internal static Complex Phi(double u, double kappa, double theta, double sigma, double rho, double v0, double s0, double r, double T)
        {
            Complex Cu = new Complex(u);
            return Phi(Cu, kappa, theta, sigma, rho, v0, s0, r, T);
        }

        public static double IntegrandFunc2(double u, double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {

            Complex Iu = Complex.I * u;
            Complex A = Complex.Exp(-Iu * Math.Log(K));


            var f2 = Phi(
               u: u,
               kappa: kappa,
               theta: theta,
               sigma: sigma,
               rho: rho,
               v0: v0,
               s0: s0,
               r: r - q,
               T: T
               );


            Complex complexVal2 = A * f2 / Iu;
            return complexVal2.Re;

        }

        public static double RhoCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentTime ?? T;
            var discountingFactor = discountingFactorFunction(0, Tp);
            var derivativeDiscontingFactor = - Tp * discountingFactor;

            var rhoUndiscountedCall = RhoUndiscountedCall(kappa, theta, rho, v0, sigma, s0, T, K, r, q);
            var undiscountedCall = HestonCall.HestonUndiscountedCallPrice(kappa, theta, rho, v0, sigma, s0, T, K, r, q);
            return discountingFactor * rhoUndiscountedCall + derivativeDiscontingFactor * undiscountedCall;
            
        }

        public static double RhoUndiscountedCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            
            if (Engine.Verbose > 0)
            {
                Console.WriteLine("Calculating rho of a call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }
            double a = 1E-8;
            double b = 1000.0;
            var firstTerm = (0.5 * (s0 * Math.Exp((r - q) * T) - K));
            var derivativeFirstTerm = (0.5 * s0 * T * Math.Exp((r - q) * T));

            TAEDelegateFunction1D functionToIntegrate = (double u) =>
                 IntegrandFunc(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K);

            double part1 = PerformIntegral(a, b, functionToIntegrate);

            double integral = part1 + a * functionToIntegrate(a / 2.0);

            double rhoGreek = (integral / Math.PI + derivativeFirstTerm);


            if (Engine.Verbose > 0)
                Console.WriteLine("Rho Undiscounted Call: {0}", rhoGreek);


            return rhoGreek;
        }

        public static double RhoPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q, Func<double, double, double> discountingFactorFunction = null, double? paymentTime = null)
        {
            var verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho of a put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }

            if (discountingFactorFunction == null)
            {
                // setting discounting Rate equal to the risk free rate
                discountingFactorFunction = (t, T) => System.Math.Exp(-r * (T - t));
            }

            var Tp = paymentTime ?? T;
            Engine.Verbose = 0;
            var rhoCall = RhoCall(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q, discountingFactorFunction: discountingFactorFunction, paymentTime: paymentTime);
            Engine.Verbose = verbosity;

            return rhoCall + Tp * discountingFactorFunction(0, Tp) * (s0 * Math.Exp((r - q) * T) - K) - discountingFactorFunction(0, Tp) * T * s0 * Math.Exp((r - q) * T);
         
        }
        public static double RhoUndiscountedPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            var verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho of a put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }


            Engine.Verbose = 0;
            var rhoPut = RhoUndiscountedCall(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
            Engine.Verbose = verbosity;

            return rhoPut - s0 * T * Math.Exp((r - q) * T);
            //return rhoPut;
        }

        private static double RhoDigital(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            // the formula to price a digital call is given by 
            //  exp(-r*T) (0.5 +  1/pi * integral from 0 to infinity of Re((exp(-i*ln(K)*u)*phi(u)/(i*u))du )

            // we can reuse a some functions define in the Heston Call class 
            double a = 1E-8;
            double b = 1000.0;

            TAEDelegateFunction1D functionToIntegrate = (double u) =>
                IntegrandFunc2(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K)
                -
                T * HestonDigital.IntegrandFunc2(u: u, kappa: kappa, theta: theta, sigma: sigma, rho: rho, v0: v0, s0: s0, r: r, q: q, T: T, K: K);

            double part1 = PerformIntegral(a, b, functionToIntegrate);
            double integral = part1 + a * functionToIntegrate(a / 2.0);

            return (1 / Math.PI * integral);
        }

        public static double RhoDigitalCall(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            if (Engine.Verbose > 0)
            {
                Console.WriteLine("Calculating rho of a digital call with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }

            var ds = HestonDigital.DiscountFactor(rate: r, T: T);
            var rho1 = RhoDigital(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
            var rhoGreek = ds * (rho1 - T / 2);

            if (Engine.Verbose > 0)
            {
                Console.WriteLine("Rho");
                Console.WriteLine(rhoGreek);
            }

            return rhoGreek;

        }
        
        public static double RhoDigitalPut(double kappa, double theta, double rho, double v0, double sigma, double s0, double T, double K, double r, double q)
        {
            var verbosity = Engine.Verbose;
            if (verbosity > 0)
            {
                Console.WriteLine("Calculating rho of a digital put with Heston model");
                Console.WriteLine("Heston Parameters");
                Console.WriteLine("kappa\ttheta\tsigma\trho\tv0");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", kappa, theta, sigma, rho, v0);
                Console.WriteLine("Call Option Information");
                Console.WriteLine("s0\tK\tT\tr\tq");
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", s0, K, T, r, q);

            }

            Engine.Verbose = 0;
            var rhoDigitalCall = RhoDigitalCall(kappa: kappa, theta: theta, rho: rho, v0: v0, sigma: sigma, s0: s0, T, K, r, q);
            Engine.Verbose = verbosity;

            return - rhoDigitalCall - T * Math.Exp(-r * T);

        }

    }
}
