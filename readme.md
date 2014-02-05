# LM.NET — Levenberg-Marquardt algorithm (LMA) for .NET

LM.NET provides an implementation of the Levenberg-Marquardt algorithm (LMA) for solving non-linear least-squares problems, e.g. for curve fitting (regression) or solving non-linear systems (in a least-squares sense).

LM.NET is a light-weight wrapper around the existing LMA package [lmfit](http://apps.jcns.fz-juelich.de/doku/sc/lmfit) written in pure C, itself derived from the original Fortran implementation [lmdif.f](http://www.netlib.org/minpack/lmdif.f) found in [MINPACK-1](http://www.netlib.org/minpack/). In the repository, the subdirectory lmfit contains a Visual C++ project that compiles lmfit as a DLL, which currently only exports the [lmmin API](http://apps.jcns.fz-juelich.de/man/lmmin.html). Note that due to its dependency on the native implementation, the resulting managed class library LMDotNet.dll is not platform agnostic anymore (e.g. x86 vs. x64).

## Quickstart: Solving a system of non-linear equations

The following example shows how to implement the [nls solving demo problem](http://apps.jcns.fz-juelich.de/doku/sc/lmfit:nonlinear-equations-example) using the LM.NET library:

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

