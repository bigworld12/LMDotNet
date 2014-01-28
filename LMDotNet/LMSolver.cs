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

        public static OptimizationResult Solve(LMDelegate fun, double[] initialGuess, SolverSettings configuration) {
            LMStatusStruct stat = new LMStatusStruct();

            LMControlStruct ctrl = new LMControlStruct {
                ftol = configuration.ftol,
                gtol = configuration.gtol,
                xtol = configuration.xtol,
                patience = configuration.patience,
                epsilon = configuration.epsilon,
                msgfile = IntPtr.Zero,
                m_maxpri = -1,
                n_maxpri = -1,
                scale_diag = configuration.scale_diag,
                stepbound = configuration.stepbound,
                verbosity = configuration.verbose ? 31 : 0
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
