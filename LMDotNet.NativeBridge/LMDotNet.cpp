// This is the main DLL file.

#include "stdafx.h"

#include "LMDotNet.h"

namespace LMDotNet {
    typedef void (*lm_callback)(const double*, int, const void*, double*, int*);
    
    OptimizationResult^ LMSolver::Solve(LMDelegate^ fun, array<double>^ initialGuess, SolverSettings^ configuration)
    {
        lm_status_struct stat;

        lm_control_struct ctrl;
        ctrl.ftol = configuration->ftol;
        ctrl.gtol = configuration->gtol;
        ctrl.xtol = configuration->xtol;
        ctrl.patience = configuration->patience;
        ctrl.epsilon = configuration->epsilon;
        ctrl.msgfile = nullptr;
        ctrl.m_maxpri = -1;
        ctrl.n_maxpri = -1;
        ctrl.scale_diag = configuration->scale_diag;
        ctrl.stepbound = configuration->stepbound;
        ctrl.verbosity = configuration->verbose ? 31 : 0; // 31: print jacobian as well
        
        // copy initial guess to output array (will be updated inplace)
        auto optimizedParameters = gcnew array<double>(initialGuess->Length);
        Array::Copy(initialGuess, optimizedParameters, initialGuess->Length);
        pin_ptr<double> p_optimizedParameters = &optimizedParameters[0];

        // make sure delegate doesn't get collected during the optimization
        GCHandle fh = GCHandle::Alloc(fun);
        auto p_fun = static_cast<lm_callback>(Marshal::GetFunctionPointerForDelegate(fun).ToPointer());
        
        lmmin(initialGuess->Length, p_optimizedParameters, initialGuess->Length, nullptr, p_fun, &ctrl, &stat);

        // copy/translate results to an OptimizationResult
        OptimizationResult^ result = gcnew OptimizationResult();
        result->fnorm = stat.fnorm;
        result->nfev = stat.nfev;
        result->optimizedParameters = optimizedParameters;
        result->outcomeMessage = LMSolver::outcomeMessages[stat.outcome];
        result->userbreak = stat.userbreak;
        result->outcome = (SolverStatus)stat.outcome;

        return result;
    }
}