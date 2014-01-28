using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LMDotNet;

namespace LMDotNet.CSharpTest
{

    class Program
    {
        // Function to be evaluated must have a signature according to
        // LMDotNet.LMDelegate (native callback, called from lmmin C code), e.g.:
        //      demo: intersect parabola with unit circle
        static unsafe void IntersectUnitCircleParabola(IntPtr par, int m_dat, IntPtr data, IntPtr fvec, IntPtr userbreak) {
            // par and fvec are native double* (double arrays), passed in
            // from the native C implementation of lmmin
            var p = (double*)par.ToPointer();  // input: the current parameter vector (here: x and y)
            var f = (double*)fvec.ToPointer(); // output: the new residuals (to be computed)
            
            // evaluate non-linear system/residuals based on
            // parameter vector p
            f[0] = p[0] * p[0] + p[1] * p[1] - 1.0; // unit circle        x^2 + y^2 = 1
            f[1] = p[1] - p[0] * p[0];            // standard parabola  y = x^2
        }

        static void Main(string[] args) {
            var lmaSolver = new LMSolver(verbose: true);
            
            // First solution: start at (1.0, 1.0)
            var res1 = lmaSolver.Solve(Program.IntersectUnitCircleParabola, new[] { 1.0, 1.0 });            
            // second solution: start at (-1.0, 1.0)
            var res2 = lmaSolver.Solve(Program.IntersectUnitCircleParabola, new[] { -1.0, 1.0 });
            
            Console.WriteLine();
            Console.WriteLine("=============================================================");
            if (res1.outcome == SolverStatus.ConvergedBoth ||
                res1.outcome == SolverStatus.ConvergedParam ||
                res1.outcome == SolverStatus.ConvergedSumSq) {
                Console.WriteLine("1st solution: x = {0}, y = {1}", res1.optimizedParameters[0], res1.optimizedParameters[1]);
            }
            if (res2.outcome == SolverStatus.ConvergedBoth ||
                res2.outcome == SolverStatus.ConvergedParam ||
                res2.outcome == SolverStatus.ConvergedSumSq) 
            {
                Console.WriteLine("2nd solution: x = {0}, y = {1}", res2.optimizedParameters[0], res2.optimizedParameters[1]);
            }


            /////////// New in v1.2.0: Using the "safe" API ///////////////////
            ///// at each iteration: //////////////////////////////////////////
            /////// copies parameter vector from unmanaged to manged heap /////
            /////// copies residual vector from managed to unmanaged heap /////

            lmaSolver.VerboseOutput = false; // switch off solver msg

            // Func-based: requires a double[] allocation (= new residuals) for each iteration            
            var res3 = lmaSolver.Solve(p => new[] { p[0] * p[0] + p[1] * p[1] - 1.0, // unit circle 0 = x^2 + y^2 - 1
                                                    p[1] - p[0] * p[0] },            // parabola    0 = y - x^2
                                       initialGuess: new[] { 1.0, 1.0 });
            if (res3.outcome == SolverStatus.ConvergedBoth ||
                res3.outcome == SolverStatus.ConvergedParam ||
                res3.outcome == SolverStatus.ConvergedSumSq) 
            {
                Console.WriteLine("Func-based -- 1st solution: x = {0}, y = {1}", res3.optimizedParameters[0], res3.optimizedParameters[1]);
            }

            // Action-based: requires no additional allocations, updates old residual array r
            var res4 = lmaSolver.Solve((p, r) => { r[0] = p[0] * p[0] + p[1] * p[1] - 1.0;  // unit circle 0 = x^2 + y^2 - 1
                                                   r[1] = p[1] - p[0] * p[0]; },            // parabola    0 = y - x^2
                                       initialGuess: new[] { 1.0, 1.0 });
            if (res4.outcome == SolverStatus.ConvergedBoth ||
                res4.outcome == SolverStatus.ConvergedParam ||
                res4.outcome == SolverStatus.ConvergedSumSq) {
                Console.WriteLine("Action-based -- 1st solution: x = {0}, y = {1}", res4.optimizedParameters[0], res4.optimizedParameters[1]);
            }
        }
    }
}
