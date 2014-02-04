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
        /// <summary>
        /// Demonstrates how to solve systems of non-linear
        /// equations using the LMA
        /// </summary>        
        static void SolveExamples() {
            var solver = new LMSolver();

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
            Console.WriteLine("1st solution: x = {0}, y = {1}", result1.optimizedParameters[0], result1.optimizedParameters[1]);
            Console.WriteLine("2nd solution: x = {0}, y = {1}", result2.optimizedParameters[0], result2.optimizedParameters[1]);
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
            Console.WriteLine("is true.");
            Console.WriteLine("1st solution: x = {0}, y = {1}", result1.optimizedParameters[0], result1.optimizedParameters[1]);
            Console.WriteLine("2nd solution: x = {0}, y = {1}", result2.optimizedParameters[0], result2.optimizedParameters[1]);
            Console.WriteLine();
        }

        /// <summary>
        /// Demonstrantes the usage of the generic minimization API
        /// LMSolver.Minimize
        /// </summary>
        static void GenericMinimizationExamples() {
            var solver = new LMSolver();

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

            Console.WriteLine("Aim: fit a 2D regression model to the 4 samples");
            Console.WriteLine("   Model: y' = β0 + β1 * x + β2 * z");
            Console.WriteLine("   Regressors: x and z (grid points)");
            Console.WriteLine("   Unknown paramters: β0, β1, β2");
            Console.WriteLine("   Design matrix X = [1 xs zs] (4x3 matrix)");
            Console.WriteLine("   Parameter vector ε = [β0 β1 β2]^T");
            Console.WriteLine("   Residue: ε = y - y' = y - X β");
            Console.WriteLine("   Find β' = argmin_β sum(ε²)");
            Console.WriteLine("Fit: y = {0} + {1} x + {2} z", fit.optimizedParameters[0], fit.optimizedParameters[1], fit.optimizedParameters[2]);
        }

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
            Program.SolveExamples();
            Program.GenericMinimizationExamples();

            var lmaSolver = new LMSolver(verbose: true);

            // First solution: start at (1.0, 1.0)
            var res1 = lmaSolver.Solve(Program.IntersectUnitCircleParabola, new[] { 1.0, 1.0 });
            // second solution: start at (-1.0, 1.0)
            var res2 = lmaSolver.Solve(Program.IntersectUnitCircleParabola, new[] { -1.0, 1.0 });
            
            Console.WriteLine();            
            Console.WriteLine("=== Static method as callback ===============================");
            Console.WriteLine("1st solution: {0}", res1.message);
            Console.WriteLine("1st solution: x = {0}, y = {1}", res1.optimizedParameters[0], res1.optimizedParameters[1]);
            Console.WriteLine("2nd solution: {0}", res2.message);
            Console.WriteLine("2nd solution: x = {0}, y = {1}", res2.optimizedParameters[0], res2.optimizedParameters[1]);
        
            Console.WriteLine();

            ///////////// Alternative: use a lambda expression: ////////////
            lmaSolver.VerboseOutput = false;

            var res3 = lmaSolver.Solve((p, r) => {
                r[0] = p[0] * p[0] + p[1] * p[1] - 1.0;            
                r[1] = p[1] - p[0] * p[0]; },
                new[] { 1.0, 1.0 });                        
            Console.WriteLine("=== Lambda as callback ======================================");
            Console.WriteLine("1st solution: x = {0}, y = {1}", res3.optimizedParameters[0], res3.optimizedParameters[1]);

            ///////////// 2nd example ////////////
            Console.WriteLine();
            Console.WriteLine("=== 2nd Example ======================================");
            
            var res4 = lmaSolver.Solve((p, r) => 
            { 
                r[0] = p[1] + p[0] * p[0] - 6.0;   // y = -x² + 6
                r[1] = p[1] +  2.0 * p[0] + 2.0;   // y = -2x - 2
            }, new[] { 0.0, 0.0 });
            
            var res5 = lmaSolver.Solve((p, r) =>
            { 
                r[0] = p[1] + p[0] * p[0] - 6.0;
                r[1] = p[1] +  2.0 * p[0] + 2.0; 
            },  new[] { 10.0, 0.0 });            

            Console.WriteLine("Solution 1: x = {0}, y = {1}", res4.optimizedParameters[0], res4.optimizedParameters[1]);
            Console.WriteLine("Solution 2: x = {0}, y = {1}", res5.optimizedParameters[0], res5.optimizedParameters[1]);

            ///////////// Generic least-squares minimization ////////////
            // Example: surface fitting

            // grid points
            var xs = new[] { -1.0, -1.0,  1.0,  1.0 };
            var zs = new[] { -1.0,  1.0, -1.0,  1.0 };
            // data points
            var ys = new[] {  0.0,  1.0,  1.0,  2.0 };

            var fit = lmaSolver.Minimize((p, r) => {
                for (int i = 0; i < ys.Length; ++i)
                    // residual for data point i, 
                    // assuming a model y = p0 + p1 * x + p2 * z
                    r[i] = ys[i] - (p[0] + p[1] * xs[i] + p[2] * zs[i]); }, 
                new[] { 0.0, 0.0, 0.0 }, 
                ys.Length);

            Console.WriteLine();
            Console.WriteLine("=== Fit surface via lmmin ============================");
            Console.WriteLine("Fit: y = {0} + {1} x + {2} z", fit.optimizedParameters[0], fit.optimizedParameters[1], fit.optimizedParameters[2]);

            ///////////// convenient 1D curve fitting ////////////
            var rand = new Random();
            var ts   = Enumerable.Range(0, 100).Select(i => i * 0.01).ToArray();
            var sins = Enumerable.Range(0, 100).Select(i => 2.0 * Math.Sin(3.0 * ts[i] + 0.1) /*+ rand.NextDouble() * 0.02*/).ToArray();
            
            var sinFit = lmaSolver.FitCurve((x, p) => p[0] * Math.Sin(p[1] * x + p[2]), new[] { 1.0, 1.0, 1.0 }, ts, sins);
            Console.WriteLine();
            Console.WriteLine("=== Fit sine to noisy data ===========================");
            Console.WriteLine("Fit: y = {0} sin({1} x + {2})", sinFit.optimizedParameters[0], sinFit.optimizedParameters[1], sinFit.optimizedParameters[2]);

            ///////////// convenient 2D surface fitting ////////////
            var surfsamples = Enumerable
                .Range(0, ys.Length)
                .Select(i => Tuple.Create(xs[i], zs[i], ys[i]))
                .ToArray();
            
            var surfFit = lmaSolver.FitSurface((x, y, p) => p[0] + p[1] * x + p[2] * y, new[] { 0.0, 0.0, 0.0 }, xs, zs, ys);
            Console.WriteLine();
            Console.WriteLine("=== Fit surface ======================================");
            Console.WriteLine("Fit: y = {0} + {1} x + {2} z", surfFit.optimizedParameters[0], surfFit.optimizedParameters[1], surfFit.optimizedParameters[2]);
        }
    }
}
