using System;
using System.Runtime.InteropServices;

namespace LMDotNet.Native
{
    /// <summary>
    /// Signature of the allocator passed into LMFit.lmmin
    /// </summary>
    /// <param name="count">Number of doubles to allocate</param>
    /// <returns>Pointer to the base address</returns>
    [UnmanagedFunctionPointer(CallingConvention.StdCall)]
    delegate IntPtr AllocaterDelegate(int count);

    /// <summary>
    /// Signature of the deallocator ("free") passed into LMFit.lmmin
    /// </summary>
    /// <param name="ptr">Managed array to "free" (unpin)</param>
    [UnmanagedFunctionPointer(CallingConvention.StdCall)]
    delegate void DeallocatorDelegate(IntPtr ptr);

    /// <summary>
    /// Signature of the user-defined equation (system) that 
    /// is to be evaluated
    /// </summary>
    /// <param name="par">[In] Current values of the free variables/parameters (double*)</param>
    /// <param name="m_dat">Number of equations/data points</param>
    /// <param name="data">[In] Auxilliary data (void*) (for curve fitting; otherwise usually 0/null)</param>
    /// <param name="fvec">[Out] Residue vector resulting from evaluating the system using the parameters in par (double*)</param>
    /// <param name="userbreak">[Out] Request termination if *userbreak == 1 (int*) (usually 0)</param>
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    delegate void LMDelegate(IntPtr par, int m_dat, IntPtr data, IntPtr fvec, IntPtr userbreak);
}
