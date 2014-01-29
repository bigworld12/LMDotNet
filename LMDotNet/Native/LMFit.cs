using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace LMDotNet.Native
{
    public static class LMFit
    {
        [DllImport("lmfit.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void lmmin(int n_par, double[] par, int m_dat, IntPtr data, LMDelegate evaluate, ref LMControlStruct control, ref LMStatusStruct status, DoubleArrayAllocatorDelegate arrayAllocator);
        /*void lmmin( int n_par, double *par, int m_dat, const void *data, void (*evaluate) (const double *par, int m_dat, const void *data, double *fvec, int *userbreak), const lm_control_struct *control, lm_status_struct *status );*/
    }
}
