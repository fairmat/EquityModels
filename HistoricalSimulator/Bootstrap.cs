using System;
using DVPLI;

namespace HistoricalSimulator
{
    public class Bootstrap
    {
        //this can be done in setup
       Vector[] CreateReturnSamples(Tuple<DateTime,Vector>[] Samples)
       {
            Vector[] returns = new Vector[Samples.Length-1];
            for(int z=1;z<Samples.Length;z++)
                returns[z]=  Vector.Log( Samples[z].Item2/Samples[z-1].Item2);
            return returns;
       }

       void Simulate(double[] dates,IMatrixSlice outDynamic)
       {
            Vector[] returns= null;//todo: get returns

            int C= returns[0].Length;//number of components
            for(int c=0;c<C;c++)
                outDynamic[0,c]=1;

            //assumes equispaced dates, otherwise we can normalize returns
            for(int i=1;i<dates.Length;i++)
            {
                //select a return using bootstrap
                double u=Engine.Generator.Uniform();
                int z =  (int)(u*returns.Length);
                for(int c=0;c<C;c++)
                    outDynamic[i,c]=outDynamic[i-1,c]*Math.Exp (returns[z][c]);
            }
       }

    }

}

