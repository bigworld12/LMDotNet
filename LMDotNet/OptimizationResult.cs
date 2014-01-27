using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace LMDotNet
{
    public class OptimizationResult
    {
        /// <summary>
        /// determined optimia
        /// </summary>
        public double[] optimizedParameters;

        /// <summary>
        /// norm of the residue vector fvec
        /// </summary> 
        public double fnorm;

        /// <summary>
        /// actual number of iterations
        /// </summary>
        public int nfev;

        /// <summary>
        /// Status indicator (converged, failed, ...)
        /// </summary>
        public SolverStatus outcome;

        /// <summary>
        /// Status message
        /// </summary>
        public string outcomeMessage;

        /// <summary>
        /// Set when function evaluation requests termination
        /// </summary>
        public int userbreak;
    }    
}
