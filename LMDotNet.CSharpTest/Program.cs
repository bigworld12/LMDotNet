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
            // Example 1: surface fitting

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
            Console.WriteLine("=== Fit surface ======================================");
            Console.WriteLine("Fit: y = {0} + {1} x + {2} z", fit.optimizedParameters[0], fit.optimizedParameters[1], fit.optimizedParameters[2]);
        }
    }
}
