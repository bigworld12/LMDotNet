namespace LMDotNet
{
    /// <summary>
    /// Information about the optimization outcome; if the
    /// optimization procedure converged, optimizedParameters
    /// contains the parameters, which minimize the user
    /// supplied function.
    /// </summary>
    public class OptimizationResult
    {
        /// <summary>
        /// Determined optimial parameters
        /// </summary>
        public double[] optimizedParameters;

        /// <summary>
        /// Norm of the residue vector fvec
        /// </summary> 
        public double errorNorm;

        /// <summary>
        /// actual number of iterations
        /// </summary>
        public int iterations;

        /// <summary>
        /// Status indicator (converged, failed, ...)
        /// </summary>
        public SolverStatus outcome;

        /// <summary>
        /// Status message
        /// </summary>
        public string message;

        /// <summary>
        /// Set when function evaluation requests termination
        /// </summary>
        public bool terminatedByUserRequest;
    }    
}
