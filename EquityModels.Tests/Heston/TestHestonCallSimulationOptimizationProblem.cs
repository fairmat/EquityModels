/* Copyright (C) 2009-2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
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
using DVPLI;
using Fairmat.MarketData;
using HestonEstimator;
using NUnit.Framework;

namespace Heston
{
    /// <summary>
    /// Covers <see cref="HestonCallSimulationOptimizationProblem"/> and its factory wrapper
    /// <see cref="HestonExtendesCalibrationBySimulationEstimator"/>, which previously had zero
    /// test coverage: the Monte-Carlo objective's happy path, its NaN/Infinity guard, and the
    /// Feller-penalty toggle.
    /// </summary>
    [TestFixture]
    public class TestHestonCallSimulationOptimizationProblem
    {
        [SetUp]
        public void Init()
        {
            TestCommon.TestInitialization.CommonInitialization();
        }

        private static EquityCalibrationData LoadEquityCalibrationData()
        {
            var (iData, hData) = TestCommon.TestInitialization.LoadHestonSampleMarketData();
            return new EquityCalibrationData(hData, iData.DiscountingCurve);
        }

        private static void DefaultBounds(out Vector matBound, out Vector strikeBound)
        {
            matBound = new Vector(2);
            matBound[0] = 1.0 / 12;
            matBound[1] = 6;
            strikeBound = new Vector(2);
            strikeBound[0] = 0.4;
            strikeBound[1] = 1.6;
        }

        [Test, Category("BigTest")]
        public void TestObjReturnsFiniteValue()
        {
            Engine.MultiThread = true;
            var equityCalData = LoadEquityCalibrationData();
            DefaultBounds(out Vector matBound, out Vector strikeBound);

            var problem = new HestonCallSimulationOptimizationProblem(equityCalData, matBound, strikeBound);

            Vector x = new Vector(5);
            x[0] = 1.5;   // kappa
            x[1] = 0.4;   // theta
            x[2] = 0.3;   // sigma
            x[3] = -0.5;  // rho
            x[4] = 0.3;   // v0

            double obj = problem.Obj(x);

            Console.WriteLine("Obj = " + obj);

            Assert.IsFalse(double.IsNaN(obj));
            Assert.IsFalse(double.IsInfinity(obj));
            Assert.GreaterOrEqual(obj, 0.0);
        }

        [Test, Category("BigTest")]
        public void TestObjHandlesExtremeParametersWithoutPropagatingNaN()
        {
            Engine.MultiThread = true;
            var equityCalData = LoadEquityCalibrationData();
            DefaultBounds(out Vector matBound, out Vector strikeBound);

            var problem = new HestonCallSimulationOptimizationProblem(equityCalData, matBound, strikeBound);

            // An unstable (anti-mean-reverting, since kappa is negative) parameter set with a
            // large vol-of-vol: this drives the simulated variance/price paths towards
            // overflow, exercising the same code path as the NaN/Infinity guard in Obj
            // (`if (double.IsNaN(sum) || double.IsInfinity(sum)) return p + 10e5 * x.Norm();`)
            // regardless of which side of the check it lands on.
            Vector x = new Vector(5);
            x[0] = -2000; // kappa (unstable: pushes variance away from theta once it drifts above it)
            x[1] = 1;     // theta
            x[2] = 50;    // sigma
            x[3] = -0.9;  // rho
            x[4] = 5;     // v0 (starts above theta, so the instability kicks in immediately)

            double obj = problem.Obj(x);

            Console.WriteLine("Obj (extreme parameters) = " + obj);

            // Whether or not the guard's true-branch fires, Obj must never itself return NaN.
            Assert.IsFalse(double.IsNaN(obj));
        }

        [Test, Category("BigTest")]
        public void TestFellerPenaltyIncreasesObjective()
        {
            Engine.MultiThread = true;
            var equityCalData = LoadEquityCalibrationData();
            DefaultBounds(out Vector matBound, out Vector strikeBound);

            // A parameter set that clearly violates the Feller condition sigma^2 > 2*kappa*theta.
            Vector x = new Vector(5);
            x[0] = 0.5;   // kappa
            x[1] = 0.1;   // theta
            x[2] = 1.5;   // sigma (sigma^2 = 2.25 >> 2*kappa*theta = 0.1)
            x[3] = -0.5;  // rho
            x[4] = 0.2;   // v0

            var problemWithoutPenalty = new HestonCallSimulationOptimizationProblem(equityCalData, matBound, strikeBound);
            problemWithoutPenalty.useFellerPenalty = false;
            double objWithoutPenalty = problemWithoutPenalty.Obj(x);

            var problemWithPenalty = new HestonCallSimulationOptimizationProblem(equityCalData, matBound, strikeBound);
            problemWithPenalty.useFellerPenalty = true;
            double objWithPenalty = problemWithPenalty.Obj(x);

            Console.WriteLine("Obj without Feller penalty = " + objWithoutPenalty);
            Console.WriteLine("Obj with Feller penalty    = " + objWithPenalty);

            Assert.Greater(objWithPenalty, objWithoutPenalty);
        }

        [Test]
        public void TestSimulationEstimatorFactoryOverride()
        {
            var equityCalData = LoadEquityCalibrationData();
            DefaultBounds(out Vector matBound, out Vector strikeBound);

            var estimator = new HestonExtendesCalibrationBySimulationEstimator();

            Assert.AreEqual("Calibrate against options (Monte Carlo simulation)", estimator.Description);

            var problem = TestNewOptimizationProblemViaReflection(estimator, equityCalData, matBound, strikeBound);
            Assert.IsInstanceOf<HestonCallSimulationOptimizationProblem>(problem);
        }

        private static object TestNewOptimizationProblemViaReflection(CallEstimator estimator, EquityCalibrationData equityCalData, Vector matBound, Vector strikeBound)
        {
            var method = typeof(CallEstimator).GetMethod("NewOptimizationProblem", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            return method.Invoke(estimator, new object[] { equityCalData, matBound, strikeBound });
        }
    }
}
