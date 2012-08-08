using System;
using DVPLI;
using DVPLDOM;

namespace Dupire
{
	public class DupireEstimator:IEstimator
	{
		public DupireEstimator ()
		{
		}	

		#region IEstimator implementation
		public Type[] GetRequirements (IEstimationSettings settings, bool multivariateRequest)
		{
			throw new System.NotImplementedException ();
		}

		public EstimationResult Estimate (System.Collections.Generic.List<object> data, IEstimationSettings settings)
		{
			IFunction impliedVol=  new PFunction2D.PFunction2D(null);

			Vector x= new Vector(){1,2};
			impliedVol.Evaluate(x);
			impliedVol.Partial(x,0);


			IFunction r = new DVPLDOM.PFunction(null);
			IFunction q = new DVPLDOM.PFunction(null);

			//create dupire outputs
			PFunction2D.PFunction2D localVol= new PFunction2D.PFunction2D(null);
			//localVol.
			return null;
		}

		public Type ProvidesTo {
			get {
				throw new System.NotImplementedException ();
			}
		}
		#endregion




	}
}

