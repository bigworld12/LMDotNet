using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LMDotNet.Native;
using System.Runtime.InteropServices;

namespace LMDotNet
{
    public class LMSolver
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

        private OptimizationResult SolveNative(
            LMDelegate fun, double[] parameters, 
            AllocaterDelegate allocate, DeallocatorDelegate deallocate,
            int mData)
        {
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

            OptimizationResult result = new OptimizationResult {
                errorNorm = stat.fnorm,
                iterations = stat.nfev,
                optimizedParameters = parameters,
                message = LMSolver.outcomeMessages[stat.outcome],
                terminatedByUserRequest = stat.userbreak > 0 ? true : false,
                outcome = (SolverStatus)stat.outcome
            };

            return result;
        }
        
       /// <summary>
       /// Determines a parameter vector that minimizes the sum of squared elements 
       /// of a residue vector that is computed by a user-supplied function fun. 
       /// </summary>
       /// <param name="fun">Updates the residuals based on the current parameters;
       /// first parameter: current parameter vector (IN; length = length(initialGuess));
       /// second parameter: residuals (OUT; length = nDataPoints)</param>
       /// <param name="initialGuess">Initial guess for the free variables; length determines
       /// the number of free variables</param>
       /// <param name="nDataPoints">Length of the residuals vector == number of datapoints == number of equations</param>
       /// <returns>Optimized parameters and status</returns>
       public OptimizationResult Minimize(Action<double[], double[]> fun, double[] initialGuess, int nDataPoints) {
            var allocator = new PinnedArrayAllocator<double>();

            // optimizedPars must be allocated via allocator, because
            // the first callback-call (== call to nativeFun) pases a 
            // pointer to optimizedPars in the "par" parameter
            var pOptimizedPars = allocator.AllocatePinnedArray(initialGuess.Length);
            double[] optimizedPars = allocator[pOptimizedPars];
            initialGuess.CopyTo(optimizedPars, 0);

            var result = SolveNative(
                // translate Action<double[], double[]> to LMDelegate
                (par, m_dat, dataPtr, fvec, userbreak) => fun(allocator[par], allocator[fvec]),
                optimizedPars,
                allocator.AllocatePinnedArray,
                allocator.UnpinArray, 
                nDataPoints);

            // pinned managed arrays allocated by lmmin may be garbage collected 
            // starting from here (if not referenced anymore)
            allocator.UnpinArray(pOptimizedPars);
            GC.KeepAlive(allocator);

            return result;
        }

        /// <summary>
        /// Solve a system of nonlinear equations (in a least-squares sense,
        /// i.e. by fitting parameters to minimize a residue vector)
        /// </summary>
        /// <param name="fun">Updates the residuals based on the current parameters;
        /// first parameter: current parameter vector (IN; length = length(initialGuess));
        /// second parameter: residuals (OUT; length = length(initialGuess))</param>
        /// <param name="initialGuess">Initial guess for the free variables; length determines
        /// the number of free variables and the number of equations, and thus, residues</param>
        /// <returns>Optimized parameters and status</returns>
        public OptimizationResult Solve(Action<double[], double[]> fun, double[] initialGuess) {
            return Minimize(fun, initialGuess, initialGuess.Length);
        }

        // convenient 1D-curve fitting à la lmcurve; user only needs to
        // supply samples and model function;
        // don't wrap lmcurve, but call SolveNative/MinimizeSumOfSquares/lmmin
        public OptimizationResult FitCurve(Func<double, double[], double> model, double[] initialGuess, Tuple<double, double>[] samples) {            
            throw new NotImplementedException();
        }
        /*typedef struct {
            const double *t;
            const double *y;
            double (*f) (double t, const double *par);
        } lmcurve_data_struct;

        void lmcurve_evaluate( const double *par, int m_dat, const void *data,
                               double *fvec, int *info )
        {
            int i;
            for ( i = 0; i < m_dat; i++ )
                fvec[i] =
                    ((lmcurve_data_struct*)data)->y[i] -
                    ((lmcurve_data_struct*)data)->f(
                        ((lmcurve_data_struct*)data)->t[i], par );
        }

        void lmcurve( int n_par, double *par, int m_dat, 
                      const double *t, const double *y,
                      double (*f)( double t, const double *par ),
                      const lm_control_struct *control,
                      lm_status_struct *status )
        {
            lmcurve_data_struct data;
            data.t = t;
            data.y = y;
            data.f = f;

            lmmin( n_par, par, m_dat, (const void*) &data,
                   lmcurve_evaluate, control, status, 
                   &malloc_array_allocator, &free_deallocator );
        }*/
    }
}
