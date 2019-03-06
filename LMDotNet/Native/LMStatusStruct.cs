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
        /// <summary>
        /// norm of the residue vector fvec.
        /// </summary>
        public double fnorm;

        /// <summary>
        /// actual number of iterations.
        /// </summary>
        public int nfev;

        /// <summary>
        /// Status indicator. Nonnegative values are used as index
        /// for the message text lm_infmsg, set in lmmin.c.
        /// </summary>
        public int outcome;

        /// <summary>
        /// Set when function evaluation requests termination.
        /// </summary>
        public int userbreak;
    }
}
