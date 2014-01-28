using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace LMDotNet
{
    /// <summary>
    /// Signature of the user-defined equation (system) that 
    /// is to be evaluated
    /// </summary>
    /// <param name="par">[In] Current values of the free variables/parameters (double*)</param>
    /// <param name="m_dat">Number of equations/data points</param>
    /// <param name="data">[In] Auxilliary data (void*) (usually 0/null)</param>
    /// <param name="fvec">[Out] Residue vector resulting from evaluating the system using the parameters in par (double*)</param>
    /// <param name="userbreak">[Out] Request termination if *userbreak == 1 (int*) (usually 0)</param>
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    public delegate void LMDelegate([In] IntPtr par, int m_dat, [In] IntPtr data, [Out] IntPtr fvec, [Out] IntPtr userbreak);
}
