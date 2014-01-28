﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace LMDotNet.Native
{
    [StructLayout(LayoutKind.Sequential)]
    public struct LMStatusStruct
    {
        public double fnorm;     /* norm of the residue vector fvec. */
        public int nfev;         /* actual number of iterations. */
        public int outcome;      /* Status indicator. Nonnegative values are used as index
                         for the message text lm_infmsg, set in lmmin.c. */
        public int userbreak;    /* Set when function evaluation requests termination. */
    }
}