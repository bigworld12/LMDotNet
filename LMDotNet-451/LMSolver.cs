using LMDotNet.Native;
using System;
using System.Diagnostics;

namespace LMDotNet
{
    /// <summary>
    /// Levenberg-Marquardt non-linear least squares solver
    /// based on lmfit
    /// </summary>
    public sealed class LMSolver
    {
        /// <summary>Relative error desired in the sum of squares.
        /// Termination occurs when both the actual and
        /// predicted relative reductions in the sum of squares
        /// are at most ftol.</summary>
        public double Ftol { get; set; }

        /// <summary> Relative error between last two approximations.
        /// Termination occurs when the relative error between
        /// two consecutive iterates is at most xtol.</summary>
        public double Xtol { get; set; }

        /// <summary> Orthogonality desired between fvec and its derivs.
        /// Termination occurs when the cosine of the angle
        /// between fvec and any column of the Jacobian is at
        /// most gtol in absolute value. (measure of degeneracy) </summary>
        public double Gtol { get; set; }

        /// <summary> Step used to calculate the Jacobian, should be
        /// slightly larger than the relative error in the
        /// user - supplied functions.</summary>
        public double Epsilon { get; set; }

        /// <summary> Used in determining the initial step bound. This
        /// bound is set to the product of stepbound and the
        /// Euclidean norm of diag*x if nonzero, or else to
        /// stepbound itself. In most cases stepbound should lie
        /// in the interval (0.1,100.0). Generally, the value
        /// 100.0 is recommended.</summary>         
        public double InitialStepbound { get; set; }

        /// <summary>
        /// Used to set the maximum number of function evaluations
        /// to patience*(number_of_parameters+1)
        /// </summary>
        public int MaxIterations { get; set; }

        /// <summary>
        /// If true, the variables will be rescaled internally. Recommended value is 1.
        /// </summary>
        public bool ScaleDiagonal { get; set; }

        /// <summary>
        /// true: print status messages to stdout
        /// </summary>
        public bool VerboseOutput { get; set; }

        // from lmmin.c
        private const double defaultTolerance = 1.0e-14;

        private static readonly string[] outcomeMessages =
          { "found zero (sum of squares below underflow limit)",
            "converged  (the relative error in the sum of squares is at most tol)",
            "converged  (the relative error of the parameter vector is at most tol)",
            "converged  (both errors are at most tol)",
            "trapped    (by degeneracy; increasing epsilon might help)",
            "exhausted  (number of function calls exceeding preset patience)",
            "failed     (ftol<tol: cannot reduce sum of squares any further)",
            "failed     (xtol<tol: cannot improve approximate solution any further)",
            "failed     (gtol<tol: cannot improve approximate solution any further)",
            "crashed    (not enough memory)",
            "exploded   (fatal coding error: improper input parameters)",
            "stopped    (break requested within function evaluation)" };

        /// <param name="ftol">Relative error desired in the sum of squares.
        /// Termination occurs when both the actual and
        /// predicted relative reductions in the sum of squares
        /// are at most ftol.</param>
        /// <param name="xtol">Relative error between last two approximations.
        /// Termination occurs when the relative error between
        /// two consecutive iterates is at most xtol.</param>
        /// <param name="gtol">Orthogonality desired between fvec and its derivs.
        /// Termination occurs when the cosine of the angle
        /// between fvec and any column of the Jacobian is at
        /// most gtol in absolute value.</param>
        /// <param name="epsilon">Step used to calculate the Jacobian, should be
        /// slightly larger than the relative error in the
        /// user - supplied functions.</param>
        /// <param name="stepbound">Used in determining the initial step bound. This
        /// bound is set to the product of stepbound and the
        /// Euclidean norm of diag*x if nonzero, or else to
        /// stepbound itself. In most cases stepbound should lie
        /// in the interval (0.1,100.0). Generally, the value
        /// 100.0 is recommended.</param>
        /// <param name="patience">Used to set the maximum number of function evaluations
        /// to patience*(number_of_parameters+1)</param>
        /// <param name="scaleDiagonal">If 1, the variables will be rescaled internally. Recommended value is 1.</param>
        /// <param name="verbose">true: print status messages to stdout</param>
        public LMSolver(double ftol = defaultTolerance,
                        double xtol = defaultTolerance,
                        double gtol = defaultTolerance,
                        double epsilon = defaultTolerance,
                        double stepbound = 100.0,
                        int patience = 100,
                        bool scaleDiagonal = true,
                        bool verbose = false) {
            this.Ftol = ftol;
            this.Xtol = xtol;
            this.Gtol = gtol;
            this.Epsilon = epsilon;
            this.InitialStepbound = stepbound;
            this.MaxIterations = patience;
            this.ScaleDiagonal = scaleDiagonal;
            this.VerboseOutput = verbose;
        }               

        /// <summary>
        /// Calls the native lmmin API from the lmfit package
        /// </summary>
        /// <param name="fun">The user supplied function to update the residue vector</param>
        /// <param name="parameters">initial guess for the parameters</param>
        /// <param name="allocate">Memory allocator</param>
        /// <param name="deallocate">Memory deallocator</param>
        /// <param name="mData">Number of data points == number of equations == length of the residue vector</param>
        /// <returns>Optimization outcome and optimal paramters, if successful</returns>
        private OptimizationResult CallNativeSolver(
            LMDelegate fun, double[] parameters, 
            AllocatorDelegate allocate, DeallocatorDelegate deallocate,
            int mData)
        {
            // build control structure understood by lmmin
            LMControlStruct ctrl = new LMControlStruct {
                ftol = this.Ftol,
                gtol = this.Gtol,
                xtol = this.Xtol,
                patience = this.MaxIterations,
                epsilon = this.Epsilon,
                msgfile = IntPtr.Zero,
                m_maxpri = -1,
                n_maxpri = -1,
                scale_diag = this.ScaleDiagonal ? 1 : 0,
                stepbound = this.InitialStepbound,
                verbosity = this.VerboseOutput ? 3 : 0
            };

            LMStatusStruct stat = new LMStatusStruct();
            // call native lmmin from lmfit package
            LMFit.lmmin(parameters.Length, parameters, mData, IntPtr.Zero, fun, ref ctrl, ref stat, allocate, deallocate);

            // extract results from lmmin's result data struct
            OptimizationResult result = new OptimizationResult(
                parameters, 
                stat.fnorm, 
                stat.nfev, 
                (SolverStatus)stat.outcome, 
                LMSolver.outcomeMessages[stat.outcome], 
                stat.userbreak > 0 ? true : false
            );

            return result;
        }
        
       /// <summary>
       /// Determines the vector x that minimizes the squared L2-norm of a user-supplied
       /// function f, i.e. it determines x_opt = argmin_x ||f(x)||²
       /// </summary>
       /// <param name="f">Evaluates the system at the point x;<br/>
       /// first parameter:  x    (IN;  length = length(initialGuess));
       /// second parameter: f(x) (OUT; length = nDataPoints)</param>
       /// <param name="x0">Initial guess for x_opt; length determines the length of x</param>
       /// <param name="nDataPoints">Length of f(x) 
       /// == number of datapoints (for regression)        
       /// == number of equations (for solving NLS)
       /// == length of the residue vector; invariant: nDataPoints &gt;= length(x0)</param>
       /// <returns>Optimum x_opt (if successful) and solution status</returns>
       public OptimizationResult Minimize(Action<double[], double[]> f, double[] x0, int nDataPoints) {
            var pool = new PinnedArrayPool<double>();

            // optimizedPars must be allocated via allocator, because
            // the first callback-call (== call to nativeFun) pases a 
            // pointer to optimizedPars in the "par" parameter
            var pOptimizedPars = pool.AllocatePinnedArray(x0.Length);
            double[] optimizedPars = pool[pOptimizedPars];
            x0.CopyTo(optimizedPars, 0);

            var result = CallNativeSolver(
                // translate Action<double[], double[]> to LMDelegate
                (par, m_dat, dataPtr, fvec, userbreak) => f(pool[par], pool[fvec]),
                optimizedPars,
                pool.AllocatePinnedArray,
                pool.UnpinArray, 
                nDataPoints);

            // pinned managed arrays allocated by lmmin may be garbage collected 
            // starting from here (if not referenced anymore)
            pool.UnpinArray(pOptimizedPars);
            GC.KeepAlive(pool);

            return result;
        }

        /// <summary>
        /// Solve a system of non-linear equations (in a least-squares sense,
        /// i.e. by optimizing parameters to minimize a residue vector)
        /// </summary>
        /// <param name="f">Computes the residuals based on the current parameters;
        /// first parameter: current parameter vector (IN; length = length(initialGuess));
        /// second parameter: residuals (OUT; length = length(initialGuess))</param>
        /// <param name="x0">Initial guess for the free variables; length determines
        /// the number of free variables and the number of equations, and thus, residuals</param>
        /// <returns>Optimized parameters and status</returns>
        public OptimizationResult Solve(Action<double[], double[]> f, double[] x0) {
            return Minimize(f, x0, x0.Length);
        }

        /// <summary>
        /// 1D-curve fitting (non-linear regression): optimize a parameter vector beta
        /// for a model equation to minimize the sum of squared residuals, i.e.:
        /// x: sample point, y: measured data, y': predicted value by the model;
        /// Residual for datapoint i: eps_i = y_i - model(x_i, beta);
        /// find beta_opt = argmin_beta ||eps_i||²
        /// </summary>
        /// <param name="model">Regression model to fit to the data points; 
        /// first parameter: sample location x_i;
        /// second parameter: parameter vector beta;
        /// result: model prediction y' = model(x_i, beta)
        /// </param>
        /// <param name="beta0">Initial guess for the model parameter vector</param>
        /// <param name="xs">sampling locations</param>
        /// <param name="ys">samples (data points)</param>
        /// <returns>Optimized model parameters and status</returns>
        public OptimizationResult FitCurve(Func<double, double[], double> model, double[] beta0, double[] xs, double[] ys) {
            Debug.Assert(xs.Length == ys.Length);

            Action<double[], double[]> fun = (parameters, residuals) => {
                for (int i = 0; i < xs.Length; ++i) {
                    residuals[i] = ys[i] - model(xs[i], parameters);
                }
            };
            return Minimize(fun, beta0, xs.Length);
        }

        /// <summary>
        /// 2D-curve fitting (non-linear regression): optimize a parameter vector beta
        /// for a model equation to minimize the sum of squared residuals, i.e.:
        /// (x, y): sample point, z: measured data, z': predicted value by the model;
        /// Residual for datapoint i: eps_i = z_i - model(x_i, y_i, beta);
        /// find beta_opt = argmin_beta ||eps_i||²
        /// </summary>
        /// <param name="model">Regression model to fit to the data points; 
        /// first parameter: sample location x_i;
        /// second parameter: sample location y_i;
        /// third parameter: parameter vector beta;
        /// result: model prediction z' = model(x_i, y_i, beta)
        /// </param>
        /// <param name="beta0">Initial guess for the model parameter vector</param>
        /// <param name="xs">First coordinate of the sampling locations</param>
        /// <param name="ys">Second coordinate of the sampling locations</param>
        /// <param name="zs">samples (data points) at (x, y) locations</param>
        /// <returns>Optimized model parameters and status</returns>
        public OptimizationResult FitSurface(Func<double, double, double[], double> model, double[] beta0, double[] xs, double[] ys, double[] zs) {
            Debug.Assert(xs.Length == ys.Length && ys.Length == zs.Length);
            
            Action<double[], double[]> fun = (parameters, residuals) => {
                for (int i = 0; i < xs.Length; ++i) {
                    residuals[i] = zs[i] - model(xs[i], ys[i], parameters);
                }
            };
            return Minimize(fun, beta0, xs.Length);
        }

        /// <summary>
        /// nD-curve fitting (non-linear regression): optimize a parameter vector beta
        /// for a model equation to minimize the sum of squared residuals, i.e.:
        /// x: n-dim sample point, z: measured data, z': predicted value by the model;
        /// Residual for datapoint i: eps_i = z_i - model(x_i, beta);
        /// find beta_opt = argmin_beta ||eps_i||²
        /// </summary>
        /// <param name="model">Regression model to fit to the data points; 
        /// first parameter: sample location x_i;
        /// second parameter: parameter vector beta;
        /// result: model prediction z' = model(x_i, beta)
        /// </param>
        /// <param name="beta0">Initial guess for the model parameter vector</param>
        /// <param name="samplePoints">Sample locations; first index: coordinate; second index: sample number</param>
        /// <param name="samples">Samples (data points) for each sample locations</param>
        /// <returns>Optimized model parameters and status</returns>
        public OptimizationResult Fit(Func<double[], double[], double> model, double[] beta0, double[][] samplePoints, double[] samples) {
            Action<double[], double[]> fun = (parameters, residuals) => {
                for (int i = 0; i < samples.Length; ++i) {
                    residuals[i] = samples[i] - model(samplePoints[i], parameters);
                }
            };
            return Minimize(fun, beta0, samples.Length);
        }
    }
}
