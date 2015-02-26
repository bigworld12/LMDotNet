namespace LMDotNet
{
    /// <summary>
    /// Information about the optimization outcome; if the
    /// optimization procedure converged, OptimizedParameters
    /// contains the parameters, which minimize the user
    /// supplied function.
    /// </summary>
    public sealed class OptimizationResult
    {
        /// <summary>
        /// Determined optimial parameters
        /// </summary>
        public double[] OptimizedParameters { get; private set; }

        /// <summary>
        /// Norm of the residue vector fvec
        /// </summary> 
        public double ErrorNorm { get; private set; }

        /// <summary>
        /// actual number of iterations
        /// </summary>
        public int Iterations { get; private set; }

        /// <summary>
        /// Status indicator (converged, failed, ...)
        /// </summary>
        public SolverStatus Outcome { get; private set; }

        /// <summary>
        /// Status message
        /// </summary>
        public string Message { get; private set; }

        /// <summary>
        /// Set when function evaluation requests termination
        /// </summary>
        public bool TerminatedByUserRequest { get; private set; }

        public OptimizationResult(
            double[] optimizedParameters, 
            double errorNorm, 
            int iterations, 
            SolverStatus outcome, 
            string message,
            bool terminatedByUser) 
        {
            this.OptimizedParameters = optimizedParameters;
            this.ErrorNorm = errorNorm;
            this.Iterations = iterations;
            this.Outcome = outcome;
            this.Message = message;
            this.TerminatedByUserRequest = terminatedByUser;
        }
    }    
}
