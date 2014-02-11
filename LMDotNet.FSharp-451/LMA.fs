namespace LMDotNet.FSharp

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
      /// to patience*(number_of_parameters+1)
      /// </summary>
      maxIterations: int

      /// <summary>
      /// If true, the variables will be rescaled internally. Recommended value is 1.
      /// </summary>
      scaleDiagonal: bool

      /// <summary>
      /// true: print status messages to stdout
      /// </summary>
      verboseOutput: bool }

module LMA =
    open System
    open LMDotNet
    
    let defaultSettings =
        let defaultTolerance = 1.0e-14;
        { ftol = defaultTolerance
          xtol = defaultTolerance
          gtol = defaultTolerance
          epsilon = defaultTolerance
          initialStepbound = 100.0
          maxIterations = 100
          scaleDiagonal = true
          verboseOutput = false }

    let init settings =
        let solver = LMSolver()
        solver.Ftol <- settings.ftol
        solver.Xtol <- settings.xtol
        solver.Gtol <- settings.gtol
        solver.Epsilon <- settings.epsilon
        solver.InitialStepbound <- settings.initialStepbound
        solver.ScaleDiagonal <- settings.scaleDiagonal
        solver.VerboseOutput <- settings.verboseOutput
        solver

    let minimize settings x0 n f =
        let solver = init settings        
        solver.Minimize(Action<_,_>(f), x0, n)
