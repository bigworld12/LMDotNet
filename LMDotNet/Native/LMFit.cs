using System;
using System.Runtime.InteropServices;

namespace LMDotNet.Native
{
    /// <summary>
    /// APIs of lmfit
    /// </summary>
    static class LMFit
    {
        /// <summary>
        /// Signature of the lmmin function, the core API of lmfit: performs
        /// generic non-linear least-squares minimization using the 
        /// Levenberg-Marquardt algorithm
        /// </summary>
        /// <param name="n_par">Number of free variables/parameters</param>
        /// <param name="par">Initial guess for the parameters; contains the optimized parameters afterwards</param>
        /// <param name="m_dat">Number of data points/equations</param>
        /// <param name="data">Additional data passed to the evaluate callback (void*)</param>
        /// <param name="evaluate">User-supplied callback, which evaluates the system with the current parameters and updates the residue vector</param>
        /// <param name="control">Settings for the solver</param>
        /// <param name="status">Result/status of the optimzation process</param>
        /// <param name="arrayAllocator">Allocator to use for allocating arrays</param>
        /// <param name="arrayDeallocator">Deallocator for freeing memory allocated using arrayAllocator</param>
        [DllImport("lmfit.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern void lmmin(
            int n_par, double[] par, int m_dat, IntPtr data, LMDelegate evaluate, 
            ref LMControlStruct control, ref LMStatusStruct status, 
            AllocaterDelegate arrayAllocator, DeallocatorDelegate arrayDeallocator);        
    }
}
