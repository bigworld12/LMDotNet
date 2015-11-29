namespace LMDotNet
{
    /// <summary>
    /// Encodes possible result conditions for the LM solver
    /// </summary>
    public enum SolverStatus
    {
        /// <summary>
        /// sum of squares below underflow limit
        /// </summary>
        Underflow = 0,

        /// <summary>
        /// the relative error in the sum of squares is at most tol
        /// </summary>
        ConvergedSumSq = 1,

        /// <summary>
        /// the relative error of the parameter vector is at most tol
        /// </summary>
        ConvergedParam = 2,

        /// <summary>
        /// both errors are at most tol
        /// </summary>
        ConvergedBoth = 3,

        /// <summary>
        /// trapped by degeneracy; increasing epsilon might help
        /// </summary>
        Trapped = 4,

        /// <summary>
        /// number of function calls exceeding preset patience/maxIterations
        /// </summary>
        Exhausted = 5,

        /// <summary>
        /// ftol &lt; tol: cannot reduce sum of squares any further
        /// </summary>
        FailedFTOL = 6,

        /// <summary>
        /// xtol &lt; tol: cannot improve approximate solution any further
        /// </summary>
        FailedXTOL = 7,

        /// <summary>
        /// gtol &lt; tol: cannot improve approximate solution any further
        /// </summary>
        FailedGTOL = 8,

        /// <summary>
        /// not enough memory
        /// </summary>
        Crashed = 9,

        /// <summary>
        /// fatal coding error: improper input parameters
        /// </summary>
        Exploded = 10,

        /// <summary>
        /// break requested within function evaluation
        /// </summary>
        Stopped = 11
    }    
}
