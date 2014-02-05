using System.Runtime.InteropServices;

namespace LMDotNet.Native
{
    /// <summary>
    /// Contains the status of the optimization process
    /// (corresponds to lm_status_struct)
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    struct LMStatusStruct
    {
        public double fnorm;     /* norm of the residue vector fvec. */
        public int nfev;         /* actual number of iterations. */
        public int outcome;      /* Status indicator. Nonnegative values are used as index
                                    for the message text lm_infmsg, set in lmmin.c. */
        public int userbreak;    /* Set when function evaluation requests termination. */
    }
}
