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
        private const bool native = true;

        /// <summary>
        /// Demonstrates how to solve systems of non-linear
        /// equations using the LMA
        /// </summary>        
        static void SolveExamples() {
            var solver = new LMSolver(useNativeBackend: native);

            // 1st example:
            // Aim: Find (x, y) for which
            //      y = -0.5x² + 3.0 &&
            //      y = exp(-x)
            // is true.
            // Transform to a minimization problem:
            //      r[0] = 0 = y + 0.5x² - 3.0
            //      r[1] = 0 = y - exp(-x)            
            Action<double[], double[]> f = (p, r) => {
                // Here: p[0] = x; p[1] = y
                var x = p[0];
                var y = p[1];
                r[0] = y + 0.5 * x * x - 3.0;
                r[1] = y - Math.Exp(-x);
            };

            // -> Find (x, y) that minimizes sum(r²) using LMA
            // first solution
            var result1 = solver.Solve(f, new[] { -2.0, 0.0 });
            // second solution
            var result2 = solver.Solve(f, new[] {  2.0, 0.0 });

            Console.WriteLine("Find (x, y) for which");
            Console.WriteLine("      y = -0.5x² + 3.0 &&");
            Console.WriteLine("      y = exp(-x)");
            Console.WriteLine("is true.");
            Console.WriteLine("1st solution: x = {0}, y = {1}", result1.OptimizedParameters[0], result1.OptimizedParameters[1]);
            Console.WriteLine("2nd solution: x = {0}, y = {1}", result2.OptimizedParameters[0], result2.OptimizedParameters[1]);
            Console.WriteLine();

            // 2nd Example:
            // Solve
            //      y = -x² + 6 &&
            //      y = -2x - 2            
            result1 = solver.Solve((p, r) => {
                r[0] = p[1] + p[0] * p[0] - 6.0;   // 0 = y + x² - 6
                r[1] = p[1] + 2.0 * p[0] + 2.0;    // 0 = y + 2x + 2
            }, new[] { 0.0, 0.0 });

            result2 = solver.Solve((p, r) => {
                r[0] = p[1] + p[0] * p[0] - 6.0;
                r[1] = p[1] + 2.0 * p[0] + 2.0;
            }, new[] { 10.0, 0.0 });
                        
            Console.WriteLine("Solve");
            Console.WriteLine("      y = -x² + 6 &&");
            Console.WriteLine("      y = -2x - 2");
            Console.WriteLine("1st solution: x = {0}, y = {1}", result1.OptimizedParameters[0], result1.OptimizedParameters[1]);
            Console.WriteLine("2nd solution: x = {0}, y = {1}", result2.OptimizedParameters[0], result2.OptimizedParameters[1]);
        }

        /// <summary>
        /// Demonstrantes the usage of the generic minimization API
        /// LMSolver.Minimize
        /// </summary>
        static void GenericMinimizationExamples() {
            var solver = new LMSolver(useNativeBackend: native);

            // grid points
            var xs = new[] { -1.0, -1.0,  1.0, 1.0 };
            var zs = new[] { -1.0,  1.0, -1.0, 1.0 };
            // data points/samples
            var ys = new[] { 0.0, 1.0, 1.0, 2.0 };

            // Aim: fit a regression model to the samples
            //    Model: y' = β0 + β1 * x + β2 * z
            //    Regressors: x and z (provided as grid points)
            //    Unknown paramters: β0, β1, β2
            //    Design matrix X = [1 xs zs] (4x3 matrix)
            //    Parameter vector ε = [β0 β1 β2]^T
            //    Residue: ε = y - y' = y - X β
            //    Find β' = argmin_β sum(ε²)
            var fit = solver.Minimize((p, r) => {
                // compute the residue vector ε = y - y' based on the
                // current parameters p (== β0, β1, β2)
                var β0 = p[0];
                var β1 = p[1];
                var β2 = p[2];                
                for (int i = 0; i < ys.Length; ++i)
                    r[i] = ys[i] - (β0 + β1 * xs[i] + β2 * zs[i]);
            }, new[] { 0.0, 0.0, 0.0 }, // initial guess for β'
            ys.Length);                 // number of samples

            Console.WriteLine();
            Console.WriteLine("Aim: fit a linear 2D regression model to four samples");
            Console.WriteLine("   Model: y' = β0 + β1 * x + β2 * z");
            Console.WriteLine("   Regressors: x and z (grid points)");
            Console.WriteLine("   Unknown paramters: β0, β1, β2");
            Console.WriteLine("   Design matrix X = [1 xs zs] (4x3 matrix)");
            Console.WriteLine("   Parameter vector ε = [β0 β1 β2]^T");
            Console.WriteLine("   Residue: ε = y - y' = y - X β");
            Console.WriteLine("   Find β' = argmin_β sum(ε²)");
            Console.WriteLine("Fit: y' = {0} + {1} x + {2} z", fit.OptimizedParameters[0], fit.OptimizedParameters[1], fit.OptimizedParameters[2]);
        }

        /// <summary>
        /// Curve fitting using the Fit... convenience functions
        /// </summary>
        static void CurveFittingExamples() {
            var solver = new LMSolver(useNativeBackend: native);

            // 1st example: fit a model that is non-linear in its parameters
            //      y' = p1 * sin(p2 * t + p3) + p4
            // to some noisy data (t, y)
            // => minimize sum((y - y')²)
            //
            // for generating some noise
            var rand = new Random(12345);
            // sample points (here: equidistant)
            var ts   = Enumerable.Range(0, 100).Select(i => i * 0.01).ToArray();
            // signal: "measured" values at sample points
            var sins = Enumerable.Range(0, 100)
                .Select(i => 2.0 * Math.Sin(3.0 * ts[i] + 0.25) + rand.NextDouble() * 0.02)
                .ToArray();
            // FitCurve evaluates the model p[0] * Math.Sin(p[1] * t + p[2]) + p[3]
            // for each t and computes the residual based on the corresponding y value
            var fit1D = solver.FitCurve((t, p) => p[0] * Math.Sin(p[1] * t + p[2]) + p[3], new[] { 1.0, 1.0, 1.0, 1.0 }, ts, sins);
            Console.WriteLine();
            Console.WriteLine("Aim: fit a sine wave to noisy data");
            Console.WriteLine("Fit: y' = {0} sin({1} x + {2}) + {3}", fit1D.OptimizedParameters[0], fit1D.OptimizedParameters[1], fit1D.OptimizedParameters[2], fit1D.OptimizedParameters[3]);

            // 2nd example: fit a linear plane to some sample points 
            // (cf. GenericMinimizationExamples for doing the same with
            // the generic Minimization API)
            // grid points
            var xs = new[] { -1.0, -1.0, 1.0, 1.0 };
            var zs = new[] { -1.0, 1.0, -1.0, 1.0 };
            // data points/samples
            var ys = new[] { 0.0, 1.0, 1.0, 2.0 };            
            var fitPlane = solver.FitSurface((x, z, p) => p[0] + p[1] * x + p[2] * z, new[] { 0.0, 0.0, 0.0 }, xs, zs, ys);
            Console.WriteLine();
            Console.WriteLine("Aim: fit a plane to four data points");
            Console.WriteLine("Fit: y = {0} + {1} x + {2} z", fitPlane.OptimizedParameters[0], fitPlane.OptimizedParameters[1], fitPlane.OptimizedParameters[2]);
        }
        
        static void Main(string[] args) {
            Program.SolveExamples();            
            Program.CurveFittingExamples();
            Program.GenericMinimizationExamples();
        }
    }
}
