namespace LMDotNet

type SolverSettings = 
    { /// <summary>Relative error desired in the sum of squares.
      /// Termination occurs when both the actual and
      /// predicted relative reductions in the sum of squares
      /// are at most ftol.</summary>
      ftol: float

      /// <summary> Relative error between last two approximations.
      /// Termination occurs when the relative error between
      /// two consecutive iterates is at most xtol.</summary>
      xtol: float

      /// <summary> Orthogonality desired between fvec and its derivs.
      /// Termination occurs when the cosine of the angle
      /// between fvec and any column of the Jacobian is at
      /// most gtol in absolute value. (measure of degeneracy) </summary>
      gtol: float

      /// <summary> Step used to calculate the Jacobian, should be
      /// slightly larger than the relative error in the
      /// user - supplied functions.</summary>
      epsilon: float

      /// <summary> Used in determining the initial step bound. This
      /// bound is set to the product of stepbound and the
      /// Euclidean norm of diag*x if nonzero, or else to
      /// stepbound itself. In most cases stepbound should lie
      /// in the interval (0.1,100.0). Generally, the value
      /// 100.0 is recommended.</summary>         
      initialStepbound: float

      /// <summary>
      /// Used to set the maximum number of function evaluations
      /// to patience * (number_of_parameters + 1)
      /// </summary>
      patience: int

      /// <summary>
      /// If true, the variables will be rescaled internally. Recommended value is 1.
      /// </summary>
      scaleDiagonal: bool

      /// <summary>
      /// true: print status messages to stdout
      /// </summary>
      verboseOutput: bool
      
      /// <summary>
      /// choose solver implementation (native or managed)
      /// </summary>
      optimizerBackend: LMBackend }

module LMA =
    open System
    open LMDotNet
    
    let defaultSettings =
        let machineEpsilon = 2.2204460492503131e-16; // smallest normalized dp value
        let defaultTolerance = 30.0 * machineEpsilon;
        { ftol = defaultTolerance
          xtol = defaultTolerance
          gtol = defaultTolerance
          epsilon = defaultTolerance
          initialStepbound = 100.0
          patience = 100
          scaleDiagonal = true
          verboseOutput = false
          optimizerBackend = LMBackend.NativeLmmin }

    let init settings =
        let solver = LMSolver()
        solver.Ftol <- settings.ftol
        solver.Xtol <- settings.xtol
        solver.Gtol <- settings.gtol
        solver.Epsilon <- settings.epsilon
        solver.InitialStepbound <- settings.initialStepbound
        solver.ScaleDiagonal <- settings.scaleDiagonal
        solver.VerboseOutput <- settings.verboseOutput
        solver.Patience <- settings.patience
        solver.OptimizerBackend <- settings.optimizerBackend
        solver
    
    /// <summary>
    /// Determines the vector x that minimizes the squared L2-norm of a user-supplied
    /// function f, i.e. it determines x_opt = argmin_x ||f(x)||²
    /// </summary>
    /// <param name="solver">Configured solver<br/>
    /// <param name="f">Evaluates the system at the point x;<br/>
    /// first parameter:  x    (IN;  length = length(initialGuess));
    /// second parameter: f(x) (OUT; length = nDataPoints)</param>
    /// <param name="x0">Initial guess for x_opt; length determines the length of x</param>
    /// <param name="n">Length of f(x) 
    /// == number of datapoints (for regression)        
    /// == number of equations (for solving NLS)
    /// == length of the residue vector; invariant: n &gt;= length(x0)</param>
    /// <returns>Optimum x_opt (if successful) and solution status</returns>
    let inline minimize (solver: LMSolver) x0 n f = solver.Minimize(Action<_,_>(f), x0, n)

    let inline solve (solver: LMSolver) x0 f = solver.Solve(Action<_,_>(f), x0)

    let inline fitCurve (solver: LMSolver) beta0 xs ys f = solver.FitCurve(Func<_,_,_>(f), beta0, xs, ys)

    let inline fitSurface (solver: LMSolver) beta0 xs ys zs f = solver.FitSurface(Func<_,_,_,_>(f), beta0, xs, ys, zs)

    let inline fit (solver: LMSolver) beta0 grid data f = solver.Fit(Func<_,_,_>(f), beta0, grid, data)
