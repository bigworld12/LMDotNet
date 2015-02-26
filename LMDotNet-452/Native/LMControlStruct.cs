using System;
using System.Runtime.InteropServices;

namespace LMDotNet.Native
{
    /// <summary>
    /// Controls the LMA solver options of LMFit.lmmin
    /// (corresponds to lm_control_struct from lmfit)
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    struct LMControlStruct
    {
        /// <summary>
        /// Relative error desired in the sum of squares.
        /// Termination occurs when both the actual and
        /// predicted relative reductions in the sum of squares
        /// are at most ftol.
        /// </summary>
        public double ftol;

        /// <summary>
        /// Relative error between last two approximations.
        /// Termination occurs when the relative error between
        /// two consecutive iterates is at most xtol.
        /// </summary>
        public double xtol;

        /// <summary>
        /// Orthogonality desired between fvec and its derivs.
        /// Termination occurs when the cosine of the angle
        /// between fvec and any column of the Jacobian is at
        /// most gtol in absolute value.
        /// </summary>
        public double gtol;

        /// <summary>
        /// Step used to calculate the Jacobian, should be
        /// slightly larger than the relative error in the
        /// user-supplied functions.
        /// </summary>
        public double epsilon;

        /// <summary>
        /// Used in determining the initial step bound. This
        /// bound is set to the product of stepbound and the
        /// Euclidean norm of diag*x if nonzero, or else to
        /// stepbound itself. In most cases stepbound should lie
        /// in the interval (0.1,100.0). Generally, the value
        /// 100.0 is recommended.
        /// </summary>
        public double stepbound;

        /// <summary>
        /// Used to set the maximum number of function evaluations
        /// to patience*(number_of_parameters+1).
        /// </summary>
        public int patience;

        /// <summary>
        /// If 1, the variables will be rescaled internally.
        /// Recommended value is 1.
        /// </summary>
        public int scale_diag;

        /// <summary>
        /// Progress messages will be written to this file.
        /// </summary>
        public IntPtr msgfile;
        
        /// <summary>
        /// OR'ed: 1: print some messages; 2: print Jacobian.
        /// </summary>
        public int verbosity;
        
        /// <summary>
        /// -1, or max number of parameters to print.
        /// </summary>
        public int n_maxpri;
        
        /// <summary>
        /// -1, or max number of residuals to print.
        /// </summary>
        public int m_maxpri;
    }
}
