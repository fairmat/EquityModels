/* Copyright (C) 2013 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Safe Khampol (safe.khampol@gmail.com)
 *            Matteo Tesser (matteo.tesser@fairmat.com)
 *            Enrico Degiuli (enrico.degiuli@fairmat.com)
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
using System.Linq;
using DVPLI;
using Fairmat.MarketData;
using NUnit.Framework;

namespace VarianceGamma
{
    [TestFixture]
    public class TestVarianceGammaCalibration
    {
        [Test]
        public void Test()
        {
            // Yield (dividend).
            double q = 0;

            // Exercise price of the option (Strike price).
            double k = 100;
            Random rand = new Random();

            // Number of observed options samples.
            int nm = 9;

            // Number of observed options samples.
            int nk = 11;

            // Drift theta of VG model.
            double theta = -0.02 + rand.NextDouble() * 0.04;

            // Volatility for VG model.
            double sigma = 0.01 + rand.NextDouble() * 0.49;

            // Nu for VG model.
            double nu = 0.01 + rand.NextDouble() * 1.99;

            // Risk free rate of return.
            double r = rand.NextDouble() * 0.4;

            // Current price of the underlying asset.
            double s0 = 60 + rand.NextDouble() * 90;

            // Time to expiration of the option.
            double t = 0.05 + rand.NextDouble() * 0.95;
            Vector x = Test(nm, nk, q, s0, r, t, theta, sigma, nu);
            Vector benchmark = (Vector)new double[] { theta, sigma, nu };
            double distance = Math.Sqrt((x - benchmark).Scalar(x - benchmark));
            Console.WriteLine(distance);
            double thetaa = x[0];
            double sigmaa = x[1];
            double nuu = x[2];
            Graph(thetaa, sigmaa, nuu, t, k, q, r);
            VGSimulation(thetaa, sigmaa, nuu, t, k, q, s0, r);
        }

        /// <summary>
        /// Computation of Variance Gamma and Black-Scholes European Option Pricing .
        /// </summary>
        /// <param name="theta"></param>
        /// <param name="sigma"></param>
        /// <param name="nu"></param>
        /// <param name="t"></param>
        /// <param name="k"></param>
        /// <param name="q"></param>
        /// <param name="r"></param>
        private static void Graph(double theta, double sigma, double nu, double t, double k, double q, double r)
        {
            int num = 200;

            double[] cbls = new double[num];
            double[] pbls = new double[num];
            double[] callvgp = new double[num];
            double[] putvgp = new double[num];

            // The spot price goes from 0 to 200.
            for (int s = 0; s < num; s++)
            {
                Vector bs_var = BlackSholes(sigma, t, k, q, s, r);
                Vector vgopt_var = VGOpt(theta, sigma, nu, t, k, q, s, r);
                cbls[s] = bs_var[0];
                pbls[s] = bs_var[1];
                callvgp[s] = vgopt_var[0];
                putvgp[s] = vgopt_var[1];
            }

            // Print results.
            Console.WriteLine(" Call Option Prices in the European Black-Scholes");
            for (int s = 0; s < num; s++)
            {
                Console.Write(cbls[s] + " ");
            }

            Console.WriteLine("\n Call Option Prices in the Variance Gamma Model");
            for (int s = 0; s < num; s++)
            {
                Console.Write(callvgp[s] + " ");
            }

            Console.WriteLine("\n Put Option Prices in the European Black-Scholes");
            for (int s = 0; s < num; s++)
            {
                Console.Write(pbls[s] + " ");
            }

            Console.WriteLine("\n Put Option Prices in the Variance Gamma Model");
            for (int s = 0; s < num; s++)
            {
                Console.Write(putvgp[s] + " ");
            }
        }

        /// <summary>
        /// Variance Gamma Simulation of Price Path.
        /// </summary>
        /// <param name="theta"></param>
        /// <param name="sigma"></param>
        /// <param name="nu"></param>
        /// <param name="t"></param>
        /// <param name="k"></param>
        /// <param name="q"></param>
        /// <param name="s0"></param>
        /// <param name="r"></param>
        private static void VGSimulation(double theta, double sigma, double nu, double t, double k, double q, double s0, double r)
        {
            int n = 1000;
            double omega = Math.Log(1 - sigma * sigma * nu * 0.5 - theta * nu) / nu;
            double dt = t / n;

            // Simulating VG as Gamma time-changed Brownian Motion.
            double[] s = new double[n];
            double[] vg = new double[n];
            double[] tt = new double[n];
            double[] dg = new double[n];
            double[] z = new double[n];

            vg[0] = 0;
            s[0] = 0;
            tt[0] = 0;

            for (int i = 1; i < n; i++)
            {
                tt[i] = tt[i - 1] + dt;
                dg[i] = GammaRnd(dt / nu, nu);

                z[i] = Engine.Generator.Normal();
                vg[i] = vg[i - 1] + theta * dg[i] + sigma * Math.Sqrt(dg[i]) * z[i];
                s[i] = s0 * Math.Exp((r - q + omega) * tt[i] + vg[i]);
            }

            Console.WriteLine("Simulating VG as Gamma time-changed Brownian Motion");
            for (int i = 0; i < n; i++)
            {
                Console.Write(tt[i] + "   ");
                Console.WriteLine(s[i]);
            }

            // Simulating VG as difference of Gamma.
            double[] ttt = new double[n];
            s = new double[n];
            vg = new double[n];

            vg[0] = 0;
            s[0] = s0;
            ttt[0] = 0;
            double uq = (Math.Sqrt(theta * theta + (2 * sigma * sigma) / nu) - theta / 2) / 2;
            double up = (Math.Sqrt(theta * theta + (2 * sigma * sigma) / nu) + theta / 2) / 2;
            for (int i = 1; i < n; i++)
            {
                ttt[i] = ttt[i - 1] + dt;
                double g1 = GammaRnd(dt / nu, up * nu);
                double g2 = GammaRnd(dt / nu, uq * nu);
                vg[i] = vg[i - 1] + g1 - g2;
                s[i] = s0 * Math.Exp((r - q + omega) * ttt[i] + vg[i]);
            }

            Console.WriteLine("Simulating VG as difference of Gamma");
            for (int i = 0; i < n; i++)
            {
                Console.Write(ttt[i] + "   ");
                Console.WriteLine(s[i]);
            }
        }

        private static Vector Test(int nm, int nk, double q, double s0, double r, double t, double theta, double sigma, double nu)
        {
            // Simulate synthetic data.
            Vector m = new Vector(nm);
            Vector k = new Vector(nk);
            Matrix cp = new Matrix(nm, nk);
            Random rand = new Random();
            for (int i = 0; i < nm; i++)
            {
                m[i] = 0.01 + rand.NextDouble() * 0.99;
            }

            for (int i = 0; i < nk; i++)
            {
                k[i] = 60 + rand.NextDouble() * 90;
            }

            for (int i = 0; i < nm; i++)
            {
                for (int j = 0; j < nk; j++)
                {
                    cp[i, j] = VarianceGammaOptionsCalibration.VGCall(theta, sigma, nu, m[i], k[j], q, s0, r);
                }
            }

            Console.WriteLine("Benchmark value");
            Console.WriteLine(new Vector() { theta, sigma, nu });

            Console.WriteLine("Call prices");
            Console.WriteLine((Vector)cp);

            // VGDiff at optimum.
            double fopt = VarianceGammaOptimizationProblem.VGDiff(new Vector() { theta, sigma, nu }, q, s0, k, r, cp, m);
            Console.WriteLine("fopt");
            Console.WriteLine(fopt);

            VarianceGammaOptionsCalibration c = new VarianceGammaOptionsCalibration();
            List<object> marketData = new List<object>();

            EquitySpotMarketData espmd = new EquitySpotMarketData();
            CallPriceMarketData cpmd = new CallPriceMarketData();
            espmd.Price = s0;
            espmd.RiskFreeRate = r;
            espmd.DividendYield = q;
            cpmd.Strike = k;
            cpmd.Maturity = m;
            cpmd.CallPrice = cp;

            marketData.Add(espmd);
            marketData.Add(cpmd);

            EstimationResult res = c.Estimate(marketData, null);
            return (Vector)res.Values;
        }

        /// <summary>
        /// Black Scholes European Call and Put Option.
        /// </summary>
        /// <param name="sigma"></param>
        /// <param name="t"></param>
        /// <param name="k"></param>
        /// <param name="q"></param>
        /// <param name="s0"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        private static Vector BlackSholes(double sigma, double t, double k, double q, double s0, double r)
        {
            // Compute d1 in BS model.
            double d1 = (Math.Log(s0 / k) + (r + 0.5 * sigma * sigma) * t) / (sigma * Math.Sqrt(t));

            // Compute d2 in BS model.
            double d2 = d1 - (sigma * Math.Sqrt(t));
            Fairmat.Statistics.Normal nom = new Fairmat.Statistics.Normal(0, 1);

            // The cumulative normal distribution functions for calls
            double ncd1 = nom.Cdf(d1);
            double ncd2 = nom.Cdf(d2);

            // The cumulative normal distribution functions for puts
            double npd1 = nom.Cdf(-d1);
            double npd2 = nom.Cdf(-d2);
            Vector bs_var = (Vector)new double[] { 0, 0 };

            // Calculate the cbls.
            bs_var[0] = s0 * ncd1 - k * Math.Exp(-r * t) * ncd2;

            // Calculate the pbls.
            bs_var[1] = k * Math.Exp(-r * t) * npd2 - s0 * npd1;
            return bs_var;
        }

        /// <summary>
        /// Variance Gamma Monte Carlo European Call and Put Option.
        /// </summary>
        /// <param name="theta"></param>
        /// <param name="sigma"></param>
        /// <param name="nu"></param>
        /// <param name="t"></param>
        /// <param name="k"></param>
        /// <param name="q"></param>
        /// <param name="s0"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        private static Vector VGOpt(double theta, double sigma, double nu, double t, double k, double q, double s0, double r)
        {
            int nn = 100000;
            double omega = Math.Log(1 - (sigma * sigma) * nu * 0.5 - theta * nu) / nu;
            double c = t / nu;
            double[] uni, eps, g, x, st, payoff;
            uni = new double[nn];
            eps = new double[nn];
            g = new double[nn];
            x = new double[nn];
            st = new double[nn];
            payoff = new double[nn];
            double a = (r - q) * t;
            double b = omega * t;
            for (int n = 0; n < nn; n++)
            {
                Random rand = new Random();
                uni[n] = rand.NextDouble();
                eps[n] = Engine.Generator.Normal();
                g[n] = GammaInv(uni[n], c, nu);
                x[n] = g[n] * theta + Math.Sqrt(g[n]) * sigma * eps[n];
                st[n] = s0 * Math.Exp(a + x[n] + b);
                payoff[n] = Math.Max(st[n] - k, 0);
            }

            double callvgp = Math.Exp(-r * t) * payoff.Sum() / nn;
            double putvgp = callvgp - s0 * Math.Exp(-q * t) + k * Math.Exp(-r * t);
            Vector vgopt_var = (Vector)new double[] { 0, 0 };
            vgopt_var[0] = callvgp;
            vgopt_var[1] = putvgp;
            return vgopt_var;
        }

        /// <summary>
        /// Gamma Random Number.
        /// </summary>
        /// <param name="shape"></param>
        /// <param name="scale"></param>
        /// <returns></returns>
        private static double GammaRnd(double shape, double scale)
        {
            double d, c, x, xsquared, v, u;

            if (shape >= 1.0)
            {
                d = shape - 1.0 / 3.0;
                c = 1.0 / Math.Sqrt(9.0 * d);

                while (true)
                {
                    do
                    {
                        x = Engine.Generator.Normal();
                        v = 1.0 + c * x;
                    }
                    while (v <= 0.0);
                    v = v * v * v;
                    Random rand = new Random();
                    u = rand.NextDouble();
                    xsquared = x * x;
                    if (u < 1.0 - .0331 * xsquared * xsquared || Math.Log(u) < 0.5 * xsquared + d * (1.0 - v + Math.Log(v)))
                        return scale * d * v;
                }
            }
            else if (shape <= 0.0)
            {
                string msg = string.Format("Shape must be positive. Received {0}.", shape);
                throw new ArgumentOutOfRangeException(msg);
            }
            else
            {
                double g = GammaRnd(shape + 1.0, 1.0);
                Random rand = new Random();
                double w = rand.NextDouble();
                return scale * g * Math.Pow(w, 1.0 / shape);
            }
        }

        private static double GammaInv(double p, double a, double b)
        {
            return (new Fairmat.Statistics.Gamma(a, b)).InvCdf(p);
        }
    }
}
