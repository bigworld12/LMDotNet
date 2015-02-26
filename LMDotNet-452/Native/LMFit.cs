using System;
using System.IO;
using System.Runtime.InteropServices;

namespace LMDotNet.Native
{
    /// <summary>
    /// APIs of lmfit
    /// </summary>
    static class LMFit
    {
        /// <summary>
        /// File name of the native lmfit (shared) library
        /// </summary>
        const string dllName = "lmfit.dll";

        /// <summary>
        /// Pointer to the dynamically (pre-)loaded lmfit library
        /// </summary>
        static IntPtr lmfitDllHandle;

        [DllImport("kernel32", SetLastError = true, CharSet = CharSet.Unicode)]
        private static extern IntPtr LoadLibrary(string lpFileName);
                
        // Preload the native lmfit library with the correct bitness;
        // => calls to lmmin are automatically resolved to this already loaded library
        static LMFit() {
            string dllDir;
            if (IntPtr.Size == 8) {
                dllDir = "lmfit64";
            }
            else {
                dllDir = "lmfit32";
            }

            lmfitDllHandle = LoadLibrary(Path.Combine(Path.Combine(Path.GetFullPath("."), dllDir), dllName));
        }

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
        [DllImport(dllName, CallingConvention = CallingConvention.Cdecl)]
        internal static extern void lmmin(
            int n_par, double[] par, int m_dat, IntPtr data, LMDelegate evaluate, 
            ref LMControlStruct control, ref LMStatusStruct status, 
            AllocatorDelegate arrayAllocator, DeallocatorDelegate arrayDeallocator);        
    }
}
