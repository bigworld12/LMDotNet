// LMDotNet.h

#pragma once

using namespace System;
using namespace System::Runtime::InteropServices;

namespace LMDotNet {
    public ref struct LMSolver
    {
        static OptimizationResult^ Solve(LMDelegate^ fun, array<double>^ initialGuess, SolverSettings^ configuration);

    private:
        static array<String^>^ outcomeMessages =
          { "found zero (sum of squares below underflow limit)",
            "converged  (the relative error in the sum of squares is at most tol)",
            "converged  (the relative error of the parameter vector is at most tol)",
            "converged  (both errors are at most tol)",
            "trapped    (by degeneracy; increasing epsilon might help)",
            "exhausted  (number of function calls exceeding preset patience)",
            "failed     (ftol<tol: cannot reduce sum of squares any further)",
            "failed     (xtol<tol: cannot improve approximate solution any further)",
            "failed     (gtol<tol: cannot improve approximate solution any further)",
            "crashed    (not enough memory)",
            "exploded   (fatal coding error: improper input parameters)",
            "stopped    (break requested within function evaluation)" };
    };
}
