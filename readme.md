# LM.NET — Levenberg-Marquardt algorithm (LMA) for .NET

LMDotNet provides an implementation of the Levenberg-Marquardt algorithm (LMA) for solving non-linear least-squares problems, e.g. for curve fitting (regression) or solving non-linear systems (in a least-squares sense).

## Quickstart

The following example shows how to implement a [nls solving demo problem](http://apps.jcns.fz-juelich.de/doku/sc/lmfit:nonlinear-equations-example) using the LMDotNet library (this project is part of the repository: LMDotNet.CSharpTest)

```
#!c#
using System;
using LMDotNet;

class LMDemo
{
    // demo: intersect parabola with unit circle
    static void IntersectUnitCircleParabola(double[] parameters, double[] residuals) {
        // evaluate non-linear system/residuals based on
        // parameter vector p
        // unit circle        x^2 + y^2 = 1
        residuals[0] = parameters[0] * parameters[0] + parameters[1] * parameters[1] - 1.0;
        // standard parabola  y = x^2
        residuals[1] = parameters[1] - parameters[0] * parameters[0];                  
    }

    static void Main(string[] args) {
        var lmaSolver = new LMSolver(verbose: true);

        // First solution: start at (1.0, 1.0)
        var res1 = lmaSolver.Solve(LMDemo.IntersectUnitCircleParabola, new[] { 1.0, 1.0 });
        // second solution: start at (-1.0, 1.0)
        var res2 = lmaSolver.Solve(LMDemo.IntersectUnitCircleParabola, new[] { -1.0, 1.0 });

        Console.WriteLine();            
        Console.WriteLine("=============================================================");
        if (res1.outcome == SolverStatus.ConvergedBoth ||
            res1.outcome == SolverStatus.ConvergedParam ||
            res1.outcome == SolverStatus.ConvergedSumSq) {
            Console.WriteLine("1st solution: x = {0}, y = {1}", res1.optimizedParameters[0], res1.optimizedParameters[1]);
        }
        if (res2.outcome == SolverStatus.ConvergedBoth ||
            res2.outcome == SolverStatus.ConvergedParam ||
            res2.outcome == SolverStatus.ConvergedSumSq) {
            Console.WriteLine("2nd solution: x = {0}, y = {1}", res2.optimizedParameters[0], res2.optimizedParameters[1]);
        }

        Console.WriteLine();
        Console.WriteLine();

        // Alternative: use a lambda expression:
        var res3 = lmaSolver.Solve((p, r) => {
            r[0] = p[0] * p[0] + p[1] * p[1] - 1.0;            
            r[1] = p[1] - p[0] * p[0]; },
            new[] { 1.0, 1.0 });
                    
        Console.WriteLine("=============================================================");
        if (res3.outcome == SolverStatus.ConvergedBoth ||
            res3.outcome == SolverStatus.ConvergedParam ||
            res3.outcome == SolverStatus.ConvergedSumSq) {
            Console.WriteLine("1st solution: x = {0}, y = {1}", res3.optimizedParameters[0], res3.optimizedParameters[1]);
        }
    }
}
```

## Implementation details

LMDotNet is a light-weight wrapper around the existing LMA package [lmfit](http://apps.jcns.fz-juelich.de/doku/sc/lmfit) written in pure C, itself derived from the original Fortran implementation found in [MINPACK-1](http://www.netlib.org/minpack/). In the repository, the subdirectory lmfit contains a Visual C++ project that compiles lmfit as a DLL, which currently only exports the [lmmin API](http://apps.jcns.fz-juelich.de/man/lmmin.html). Note that due to its dependency on the native implementation, the resulting managed class library LMDotNet.dll is not platform agnostic anymore (e.g. x86 vs. x64).