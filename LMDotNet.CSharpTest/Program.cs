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
        // LMDotNet.LMDelegate, e.g.:
        //      demo: intersect parabola with unit circle
        static unsafe void IntersectUnitCircleParabola(IntPtr par, int m_dat, IntPtr data, IntPtr fvec, IntPtr userbreak) {
            // par and fvec are native double* (double arrays), passed in
            // from the native C implementation of lmmin
            var p = (double*)par.ToPointer();  // input: the current parameter vector (here: x and y)
            var f = (double*)fvec.ToPointer(); // output: the new residuals (to be computed)
            
            // evaluate non-linear system/residuals based on
            // parameter vector p
            f[0] = p[0] * p[0] + p[1] * p[1] - 1; // unit circle        x^2 + y^2 = 1
            f[1] = p[1] - p[0] * p[0];            // standard parabola  y = x^2
        }

        static void Main(string[] args) {
            // solver configuration
            var ctrl = SolverSettings.defaultSettings;
            ctrl.verbose = true; // print status messages
            
            // First solution: start at (1.0, 1.0)
            var res1 = LMSolver.Solve(Program.IntersectUnitCircleParabola, new[] { 1.0, 1.0 }, ctrl);            
            // second solution: start at (-1.0, 1.0)
            var res2 = LMSolver.Solve(Program.IntersectUnitCircleParabola, new[] { -1.0, 1.0 }, ctrl);

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
        }
    }
}
