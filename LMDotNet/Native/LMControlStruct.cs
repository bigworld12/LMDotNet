using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace LMDotNet.Native
{
    /// <summary>
    /// Controls the LMA solver options of LMFit.lmmin
    /// (corresponds lm_control_struct from lmfit)
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    struct LMControlStruct
    {
        public double ftol;      /* Relative error desired in the sum of squares.
                         Termination occurs when both the actual and
                         predicted relative reductions in the sum of squares
                         are at most ftol. */
        public double xtol;      /* Relative error between last two approximations.
                         Termination occurs when the relative error between
                         two consecutive iterates is at most xtol. */
        public double gtol;      /* Orthogonality desired between fvec and its derivs.
                         Termination occurs when the cosine of the angle
                         between fvec and any column of the Jacobian is at
                         most gtol in absolute value. */
        public double epsilon;   /* Step used to calculate the Jacobian, should be
                         slightly larger than the relative error in the
                         user-supplied functions. */
        public double stepbound; /* Used in determining the initial step bound. This
                         bound is set to the product of stepbound and the
                         Euclidean norm of diag*x if nonzero, or else to
                         stepbound itself. In most cases stepbound should lie
                         in the interval (0.1,100.0). Generally, the value
                         100.0 is recommended. */
        public int patience;     /* Used to set the maximum number of function evaluations
                         to patience*(number_of_parameters+1). */
        public int scale_diag;   /* If 1, the variables will be rescaled internally.
                         Recommended value is 1. */
        public IntPtr msgfile;    /* Progress messages will be written to this file. */
        public int verbosity;    /* OR'ed: 1: print some messages; 2: print Jacobian. */
        public int n_maxpri;     /* -1, or max number of parameters to print. */
        public int m_maxpri;     /* -1, or max number of residuals to print. */
    }
}
