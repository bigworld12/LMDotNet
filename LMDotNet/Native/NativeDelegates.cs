using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace LMDotNet.Native
{
    [UnmanagedFunctionPointer(CallingConvention.StdCall)]
    public delegate IntPtr DoubleArrayAllocatorDelegate(int count); 
}
