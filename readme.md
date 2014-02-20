# LM.NET — Levenberg-Marquardt algorithm (LMA) for .NET

LM.NET provides an implementation of the Levenberg-Marquardt algorithm (LMA) for solving non-linear least-squares problems, e.g. for curve fitting (regression) or solving non-linear systems (in a least-squares sense).

LM.NET is a light-weight wrapper around the existing LMA package [lmfit](http://apps.jcns.fz-juelich.de/doku/sc/lmfit) written in pure C, itself derived from the original Fortran implementation [lmdif.f](http://www.netlib.org/minpack/lmdif.f) found in [MINPACK-1](http://www.netlib.org/minpack/). In the repository, the subdirectory **lmfit** contains a Visual C++ project that compiles lmfit to a DLL, which exports (a customized version of) lmfit's core API [lmmin](http://apps.jcns.fz-juelich.de/man/lmmin.html). The managed class library `LMDotNet.dll` selects the correct native DLL to load depending on the bitness of the .NET process at runtime (cf. `LMDotNet-451\Native\LMFit.cs`).

## Installation
There are two ways to get LM.NET:

### Via NuGet
This is most likely the easiest way to set up LM.NET. Just type

    PM> Install-Package LMDotNet 
    
in the NuGet Package Manager Console (or use the built-in NuGet Package Manager GUI in Visual Studio). The package automatically chooses the right version of the class library depending on the target framework version of the current project (v3.5 up to v4.5.1). The packge also automatically installs the native libraries and corresponding debug symbols into the `lmfit32` (x86) and `lmfit64` (x64) project subdirectories. The correct version is then chosen automatically at runtime by LM.NET.

*Important note:* The `lmfit` binaries are built using Visual Studio 2013; you thus need to make sure that you have installed the [Visual C++ Redistributable Packages for Visual Studio 2013](http://www.microsoft.com/en-us/download/details.aspx?id=40784).

### Building from Source
Download the current version of the LM.NET source using Mercurial

    hg clone https://bitbucket.org/frank_niemeyer/lmdotnet

... and build the Visual Studio solution. The folders `LMDotNet-[version]\bin\[Release|Debug]` will contain the class library for the corresponding .NET framework version.

## Redistribution
If you want to redistribute software that uses LM.NET, you need to include the managed class library `LMDotNet.dll` as well as the native libraries it depends on in your application package. Assuming that your main binary is `.\myapp.exe`, assure the following layout:

 - `.\LMDotNet.dll` - the managed class library
 - `.\lmfit32\lmfit.dll` - x86 version of lmfit
 - `.\lmfit64\lmfit.dll` - x64 version of lmfit

If you know that your .NET application will always run in either x86 or x64 mode (i.e. the platform target is **not** "AnyCPU") you can optionally omit the unneeded native library version.

## Quickstart: Solving a system of non-linear equations

The following example shows how to implement the [nls solving demo problem](http://apps.jcns.fz-juelich.de/doku/sc/lmfit) using the LM.NET library:

```
#!c#
    var lmaSolver = new LMSolver(verbose: true);
    
    // find pairs (x, y) for which x² + y² = 1 and y = x²
    var res1 = lmaSolver.Solve((p, r) => {
        r[0] = p[0] * p[0] + p[1] * p[1] - 1.0; // unit circle
        r[1] = p[1] - p[0] * p[0]; },           // parabola
        new[] { 1.0, 1.0 });                    // initial guess for (x, y)
    
    var res2 = lmaSolver.Solve((p, r) => {
        r[0] = p[0] * p[0] + p[1] * p[1] - 1.0; // unit circle
        r[1] = p[1] - p[0] * p[0]; },           // parabola
        new[] { -1.0, 1.0 });                   // initial guess for (x, y)
    
    // (x, y) pairs for which the squared L2 norm of F -> min
    Console.WriteLine("1st solution: x = {0}, y = {1}", 
        res1.optimizedParameters[0], res1.optimizedParameters[1]);            
    Console.WriteLine("2nd solution: x = {0}, y = {1}", 
        res2.optimizedParameters[0], res2.optimizedParameters[1]);            
```

## License
LM.NET is distributed under the terms of the FreeBSD license (BSD 2-clause license):

Copyright (c) 2013, 2014, Frank Niemeyer
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
