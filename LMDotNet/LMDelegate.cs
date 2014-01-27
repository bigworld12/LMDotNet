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
    /// <param name="m_dat">Number of equations/data points (double*)</param>
    /// <param name="data">[In] Auxilliary data (void*) (usually 0/null)</param>
    /// <param name="fvec">[Out] Residue vector resulting from evaluating the system using the parameters in par</param>
    /// <param name="userbreak">[Out] Request termination (usually 0)</param>
    public delegate void LMDelegate([In] IntPtr par, int m_dat, [In] IntPtr data, [Out] IntPtr fvec, [Out] IntPtr userbreak);
}
