using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using HestonEstimator;

namespace Heston
{
    internal class HestonDigital: HestonCall
    {
        // define the constructors from the base class 
        public HestonDigital() { }

        public HestonDigital(HestonCallOptimizationProblem problem) : base(problem) { }

        public HestonDigital(HestonProcess process, double strike, double timeToMaturity) : base (process, strike, timeToMaturity) { } 

    }
}
