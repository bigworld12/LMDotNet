using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LMDotNet.Native;

namespace LMDotNet
{
    public class LMSolver
    {
        /// <summary>Relative error desired in the sum of squares.
        /// Termination occurs when both the actual and
        /// predicted relative reductions in the sum of squares
        /// are at most ftol.</summary>
        double ftol;

        /// <summary> Relative error between last two approximations.
        /// Termination occurs when the relative error between
        /// two consecutive iterates is at most xtol.</summary>
        double xtol;

        /// <summary> Orthogonality desired between fvec and its derivs.
        /// Termination occurs when the cosine of the angle
        /// between fvec and any column of the Jacobian is at
        /// most gtol in absolute value. </summary>
        double gtol;

        /// <summary> Step used to calculate the Jacobian, should be
        /// slightly larger than the relative error in the
        /// user - supplied functions.</summary>
        double epsilon;

        /// <summary> Used in determining the initial step bound. This
        /// bound is set to the product of stepbound and the
        /// Euclidean norm of diag*x if nonzero, or else to
        /// stepbound itself. In most cases stepbound should lie
        /// in the interval (0.1,100.0). Generally, the value
        /// 100.0 is recommended.</summary>         
        double stepbound;

        /// <summary>
        /// Used to set the maximum number of function evaluations
        /// to patience*(number_of_parameters+1)
        /// </summary>
        int patience;

        /// <summary>
        /// If 1, the variables will be rescaled internally. Recommended value is 1.
        /// </summary>
        int scale_diag;

        /// <summary>
        /// true: print status messages to stdout
        /// </summary>
        bool verbose;

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

        /// <summary>
        /// 
        /// </summary>
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
        /// <param name="scale_diag">If 1, the variables will be rescaled internally. Recommended value is 1.</param>
        /// <param name="verbose">true: print status messages to stdout</param>
        public LMSolver(double ftol = defaultTolerance,
                        double xtol = defaultTolerance,
                        double gtol = defaultTolerance,
                        double epsilon = defaultTolerance,
                        double stepbound = 100.0,
                        int patience = 100,
                        int scale_diag = 1,
                        bool verbose = false) {
            this.ftol = ftol;
            this.xtol = xtol;
            this.gtol = gtol;
            this.epsilon = epsilon;
            this.stepbound = stepbound;
            this.patience = patience;
            this.scale_diag = scale_diag;
            this.verbose = verbose;
        }

        public OptimizationResult Solve(LMDelegate fun, double[] initialGuess) {
            LMStatusStruct stat = new LMStatusStruct();

            LMControlStruct ctrl = new LMControlStruct {
                ftol = this.ftol,
                gtol = this.gtol,
                xtol = this.xtol,
                patience = this.patience,
                epsilon = this.epsilon,
                msgfile = IntPtr.Zero,
                m_maxpri = -1,
                n_maxpri = -1,
                scale_diag = this.scale_diag,
                stepbound = this.stepbound,
                verbosity = this.verbose ? 31 : 0
            };
            
            double[] optimizedPars = new double[initialGuess.Length];
            initialGuess.CopyTo(optimizedPars, 0);
            
            LMFit.lmmin(initialGuess.Length, optimizedPars, initialGuess.Length, IntPtr.Zero, fun, ref ctrl, ref stat);

            OptimizationResult result = new OptimizationResult {
                fnorm = stat.fnorm,
                nfev = stat.nfev,
                optimizedParameters = optimizedPars,
                outcomeMessage = LMSolver.outcomeMessages[stat.outcome],
                userbreak = stat.userbreak,
                outcome = (SolverStatus)stat.outcome
            };

            return result;
        }
    }
}
