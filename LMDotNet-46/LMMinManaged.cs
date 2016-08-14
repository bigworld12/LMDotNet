using System;

namespace LMDotNet
{
    /// <summary>
    /// Direct C# port (copy & paste, as far as possible) of lmmin.c
    /// </summary>
    /// <remarks>
    /// /*
    /// * Library:   lmfit(Levenberg-Marquardt least squares fitting)
    /// *
    /// * File:      lmmin.c
    /// *
    /// * Contents:  Levenberg-Marquardt minimization.
    /// *
    /// * Copyright: MINPACK authors, The University of Chicago (1980-1999)
    /// * Joachim Wuttke, Forschungszentrum Juelich GmbH(2004-2013)
    /// *
    /// * License:   see../COPYING(FreeBSD)
    /// *
    /// * Homepage:  apps.jcns.fz-juelich.de/lmfit
    /// */
    /// </remarks>
    static class LMMinManaged
    {
        private const double LM_DWARF = double.MinValue;
        private static readonly double LM_SQRT_DWARF = Math.Sqrt(double.MinValue);
        private static readonly double LM_SQRT_GIANT = Math.Sqrt(double.MaxValue);
        private const double LM_MACHEP = 2.2204460492503131e-16;
        private const double LM_USERTOL = 30.0 * LM_MACHEP;
        private const double p1 = 0.1;
        private const double p0001 = 1.0e-4;

        private static void PrintPars(int nout, double[] par, double fnorm) {
            int i;
            for (i = 0; i < nout; ++i)
                Console.Write($"{par[i]:G16} ");
            Console.WriteLine($" {fnorm:G16}");
        }

        private static double Square(double x) {
            return x * x;
        }

    /*****************************************************************************/
    /*  lmmin (main minimization routine)                                        */
    /*****************************************************************************/
    /*
     *   This routine contains the core algorithm of our library.
     *
     *   It minimizes the sum of the squares of m nonlinear functions
     *   in n variables by a modified Levenberg-Marquardt algorithm.
     *   The function evaluation is done by the user-provided routine 'evaluate'.
     *   The Jacobian is then calculated by a forward-difference approximation.
     *
     *   Parameters:
     *
     *      n is the number of variables (INPUT, positive integer).
     *
     *      x is the solution vector (INPUT/OUTPUT, array of length n).
     *        On input it must be set to an estimated solution.
     *        On output it yields the final estimate of the solution.
     *
     *      m is the number of functions to be minimized (INPUT, positive integer).
     *        It must fulfill m>=n.
     *
     *      data is a pointer that is ignored by lmmin; it is however forwarded
     *        to the user-supplied functions evaluate and printout.
     *        In a typical application, it contains experimental data to be fitted.
     *
     *      evaluate is a user-supplied function that calculates the m functions.
     *        Parameters:
     *          n, x, m, data as above.
     *          fvec is an array of length m; on OUTPUT, it must contain the
     *            m function values for the parameter vector x.
     *          userbreak is an integer pointer. When *userbreak is set to a 
     *            nonzero value, lmmin will terminate.
     *
     *      control contains INPUT variables that control the fit algorithm,
     *        as declared and explained in lmstruct.h
     *
     *      status contains OUTPUT variables that inform about the fit result,
     *        as declared and explained in lmstruct.h
     */
    public static void LMMin(int n, double[] x, int m, /*const void* data,*/
            Action<double[], double[]> evaluate,
            ref Native.LMControlStruct C,
            ref Native.LMStatusStruct S) 
        {
            double[] fvec, diag, fjac, qtf, wa1, wa2, wa3, wf;
            int[] ipvt;
            int j, i;
            double actred, dirder, fnorm, fnorm1, gnorm, pnorm, prered, ratio, step, sum, temp, temp1, temp2, temp3;
            int maxfev = C.patience * (n + 1);
            int outer, inner;  /* loop counters, for monitoring */
            bool inner_success; /* flag for loop control */
            double lmpar = 0;  /* Levenberg-Marquardt parameter */
            double delta = 0;
            double xnorm = 0;
            double eps = Math.Sqrt(Math.Max(C.epsilon, LM_MACHEP)); /* for forward differences */
            int nout = C.n_maxpri == -1 ? n : Math.Min(C.n_maxpri, n);
            
            /* Default status info; must be set ahead of first return statements */
            S.outcome = 0;    /* status code */
            S.userbreak = 0;
            S.nfev = 0;      /* function evaluation counter */

            /***  Check input parameters for errors.  ***/
            if (n <= 0) {
                Console.Error.WriteLine($"lmmin: invalid number of parameters {n}");
                S.outcome = 10; /* invalid parameter */
                return;
            }
            if (m < n) {
                Console.Error.WriteLine($"lmmin: number of data points {m} smaller than number of parameters {n}");
                S.outcome = 10;
                return;
            }
            if (C.ftol < 0.0 || C.xtol < 0.0 || C.gtol < 0.0) {
                Console.Error.WriteLine($"lmmin: negative tolerance (at least one of {C.ftol} {C.xtol} {C.gtol}");
                S.outcome = 10;
                return;
            }
            if (maxfev <= 0) {
                Console.Error.WriteLine($"lmmin: nonpositive function evaluations limit {maxfev}");
                S.outcome = 10;
                return;
            }
            if (C.stepbound <= 0.0) {
                Console.Error.WriteLine($"lmmin: nonpositive stepbound {C.stepbound}");
                S.outcome = 10;
                return;
            }
            if (C.scale_diag != 0 && C.scale_diag != 1) {
                Console.Error.WriteLine($"lmmin: logical variable scale_diag={C.scale_diag} should be 0 or 1");
                S.outcome = 10;
                return;
            }

            /***  Allocate work space.  ***/
            fvec = new double[m];
            diag = new double[n];
            qtf = new double[n];
            fjac = new double[n * m];
            wa1 = new double[n];
            wa2 = new double[n];
            wa3 = new double[n];
            wf = new double[m];
            ipvt = new int[n];
            
            if (C.scale_diag != 1) {
                for (j = 0; j < n; j++)
                    diag[j] = 1.0;
            }

            /***  Evaluate function at starting point and calculate norm.  ***/
            evaluate(x, fvec);
            S.nfev = 1;
            if (S.userbreak != 0)
                goto terminate;
            fnorm = EuclideanNorm(m, 0, fvec);
            if (C.verbosity > 0) {
                Console.WriteLine("lmmin start ");
                PrintPars(nout, x, fnorm);
            }
            if (fnorm <= LM_DWARF) {
                S.outcome = 0; /* sum of squares almost zero, nothing to do */
                goto terminate;
            }

            /***  The outer loop: compute gradient, then descend.  ***/
            for (outer = 0; ; ++outer) {
                
                /***  [outer]  Calculate the Jacobian.  ***/
                for (j = 0; j < n; j++) {
                    temp = x[j];
                    step = Math.Max(eps* eps, eps * Math.Abs(temp));
                    x[j] += step; /* replace temporarily */
                    evaluate(x, wf);
                    ++S.nfev;
                    if (S.userbreak != 0)
                        goto terminate;
                    for (i = 0; i<m; i++)
                        fjac[j * m + i] = (wf[i] - fvec[i]) / step;
                    x[j] = temp; /* restore */
                }

                /***  [outer]  Compute the QR factorization of the Jacobian.  ***/
                /*      fjac is an m by n array. The upper n by n submatrix of fjac 
                 *        is made to contain an upper triangular matrix r with diagonal
                 *        elements of nonincreasing magnitude such that
                 *
                 *              p^T*(jac^T*jac)*p = r^T*r
                 *
                 *              (NOTE: ^T stands for matrix transposition),
                 *
                 *        where p is a permutation matrix and jac is the final calculated
                 *        Jacobian. Column j of p is column ipvt(j) of the identity matrix.
                 *        The lower trapezoidal part of fjac contains information generated
                 *        during the computation of r.
                 *
                 *      ipvt is an integer array of length n. It defines a permutation
                 *        matrix p such that jac*p = q*r, where jac is the final calculated
                 *        Jacobian, q is orthogonal (not stored), and r is upper triangular
                 *        with diagonal elements of nonincreasing magnitude. Column j of p
                 *        is column ipvt(j) of the identity matrix.
                 */
                QRFactorization(m, n, fjac, ipvt, wa1, wa2, wa3);                       
                /* return values are ipvt, wa1=rdiag, wa2=acnorm */

                /***  [outer]  Form q^T * fvec and store first n components in qtf.  ***/
                for (i = 0; i < m; i++)
                    wf[i] = fvec[i];

                for (j = 0; j<n; j++) {
                    temp3 = fjac[j * m + j];
                    if (temp3 != 0.0) {
                        sum = 0;
                        for (i = j; i<m; i++)
                            sum += fjac[j * m + i] * wf[i];
                        temp = -sum / temp3;
                        for (i = j; i<m; i++)
                            wf[i] += fjac[j * m + i] * temp;
                    }
                    fjac[j * m + j] = wa1[j];
                    qtf[j] = wf[j];
                }

                /***  [outer]  Compute norm of scaled gradient and detect degeneracy.  ***/
                gnorm = 0;
                for (j = 0; j<n; j++) {
                    if (wa2[ipvt[j]] == 0)
                        continue;
                    sum = 0.0;
                    for (i = 0; i <= j; i++)
                        sum += fjac[j * m + i] * qtf[i];
                    gnorm = Math.Max(gnorm, Math.Abs(sum / wa2[ipvt[j]] / fnorm));
                }

                if (gnorm <= C.gtol) {
                    S.outcome = 4;
                    goto terminate;
                }

                /***  [outer]  Initialize / update diag and delta. ***/
                if (outer == 0) { 
                    /* first iteration only */
                    if (C.scale_diag != 0) {
                        /* diag := norms of the columns of the initial Jacobian */
                        for (j = 0; j<n; j++)
                            diag[j] = wa2[j] != 0.0 ? wa2[j] : 1;
                        /* xnorm := || D x || */
                        for (j = 0; j<n; j++)
                            wa3[j] = diag[j] * x[j];
                        xnorm = EuclideanNorm(n, 0, wa3);
                        if (C.verbosity >= 2) {
                            Console.Write("lmmin diag  ");
                            PrintPars(nout, x, xnorm);
                        }
                        /* only now print the header for the loop table */
                        if (C.verbosity >= 3) {
                            Console.Write("  o  i     lmpar    prered          ratio    dirder      delta      pnorm                 fnorm");
                            for (i = 0; i < nout; ++i)
                                Console.Write($"               p{i}");
                            Console.WriteLine();
                        }
                    }
                    else {
                        xnorm = EuclideanNorm(n, 0, x);
                    }
                    /* initialize the step bound delta. */
                    if (xnorm != 0.0)
                        delta = C.stepbound * xnorm;
                    else
                        delta = C.stepbound;
                }
                else {
                    if (C.scale_diag != 0) {
                        for (j = 0; j < n; j++)
                            diag[j] = Math.Max(diag[j], wa2[j] );
                    }
                }

                /***  The inner loop. ***/
                inner = 0;
                do {
                    /***  [inner]  Determine the Levenberg-Marquardt parameter.  ***/
                    LMParameter(n, fjac, m, ipvt, diag, qtf, delta, ref lmpar, wa1, wa2, wf, wa3);                    
                    /* used return values are fjac (partly), lmpar, wa1=x, wa3=diag*x */

                    /* predict scaled reduction */
                    pnorm = EuclideanNorm(n, 0, wa3);
                    var pdf = pnorm / fnorm;
                    temp2 = lmpar * pdf * pdf;
                    for (j = 0; j<n; j++) {
                        wa3[j] = 0;
                        for (i = 0; i <= j; i++)
                            wa3[i] -= fjac[j * m + i] * wa1[ipvt[j]];
                    }
                    temp1 = Square(EuclideanNorm(n, 0, wa3) / fnorm);
                    prered = temp1 + 2 * temp2;
                    dirder = -temp1 + temp2; /* scaled directional derivative */

                    /* at first call, adjust the initial step bound. */
                    if (outer == 0 && pnorm < delta)
                        delta = pnorm;

                    /***  [inner]  Evaluate the function at x + p.  ***/
                    for (j = 0; j < n; j++)
                        wa2[j] = x[j] - wa1[j];

                    evaluate(wa2, wf);
                    ++S.nfev;
                    if (S.userbreak != 0)
                        goto terminate;
                    fnorm1 = EuclideanNorm(m, 0, wf);

                    /***  [inner]  Evaluate the scaled reduction.  ***/

                    /* actual scaled reduction */
                    actred = 1 - Square(fnorm1 / fnorm);

                    /* ratio of actual to predicted reduction */
                    ratio = prered != 0.0? actred/prered : 0;

                    if (C.verbosity == 2 ) {
                        Console.Write($"lmmin ({outer}:{inner}) ");
                        PrintPars(nout, wa2, fnorm1);
                    }
                    else if (C.verbosity >= 3) {
                        Console.Write($"{outer} {inner} {lmpar} {prered} {ratio} {dirder} {delta} {pnorm} {fnorm1}");
                        for (i = 0; i < nout; ++i)
                            Console.Write($" {wa2[i]}");
                        Console.WriteLine();
                    }

                    /* update the step bound */
                    if (ratio <= 0.25) {
                        if (actred >= 0)
                            temp = 0.5;
                        else if (actred > -99) /* -99 = 1-1/0.1^2 */
                            temp = Math.Max(dirder / (2 * dirder + actred), 0.1);
                        else
                            temp = 0.1;
                        delta = temp * Math.Min(delta, pnorm / 0.1);
                        lmpar /= temp;
                    }
                    else if (ratio >= 0.75) {
                        delta = 2 * pnorm;
                        lmpar *= 0.5;
                    }
                    else if (lmpar == 0.0) {
                        delta = 2 * pnorm;
                    }

                    /***  [inner]  On success, update solution, and test for convergence.  ***/
                    inner_success = ratio >= p0001;
                    if (inner_success) {
                        /* update x, fvec, and their norms */
                        if (C.scale_diag != 0) {
                            for (j = 0; j<n; j++) {
                                x[j] = wa2[j];
                                wa2[j] = diag[j] * x[j];
                            }
                        }
                        else {
                            for (j = 0; j<n; j++)
                                x[j] = wa2[j];
                        }
                        for (i = 0; i<m; i++)
                            fvec[i] = wf[i];
                        xnorm = EuclideanNorm(n, 0, wa2);
                        fnorm = fnorm1;
                    }

                    /* convergence tests */ 
                    S.outcome = 0;
                    if (fnorm <= LM_DWARF)
                        goto terminate;  /* success: sum of squares almost zero */
                    /* test two criteria (both may be fulfilled) */
                    if (Math.Abs(actred) <= C.ftol && prered <= C.ftol && ratio <= 2)
                        S.outcome = 1;  /* success: x almost stable */
                    if (delta <= C.xtol * xnorm)
                        S.outcome += 2; /* success: sum of squares almost stable */
                    if (S.outcome != 0) {
                        goto terminate;
                    }

                    /***  [inner]  Tests for termination and stringent tolerances.  ***/
                    if (S.nfev >= maxfev) {
                        S.outcome = 5;
                        goto terminate;
                    }
                    if (Math.Abs(actred) <= LM_MACHEP &&
                        prered <= LM_MACHEP && ratio <= 2) {
                        S.outcome = 6;
                        goto terminate;
                    }
                    if (delta <= LM_MACHEP * xnorm) {
                        S.outcome = 7;
                        goto terminate;
                    }
                    if (gnorm <= LM_MACHEP) {
                        S.outcome = 8;
                        goto terminate;
                    }

                    /***  [inner]  End of the loop. Repeat if iteration unsuccessful.  ***/
                    ++inner;
                } while (!inner_success);
            /***  [outer]  End of the loop. ***/
            };

        terminate:
            S.fnorm = EuclideanNorm(m, 0, fvec);
            if (C.verbosity >= 2)
                Console.WriteLine($"lmmin outcome ({S.outcome}) xnorm {xnorm} ftol {C.ftol} xtol {C.xtol}");                
            if (C.verbosity % 2 != 0) {
                Console.Write("lmmin final ");
                PrintPars(nout, x, S.fnorm);
            }
            if (S.userbreak == 1) /* user-requested break */
                S.outcome = 11;
        } /*** lmmin. ***/

        /*****************************************************************************/
        /*  lm_lmpar (determine Levenberg-Marquardt parameter)                       */
        /*****************************************************************************/
        /*     Given an m by n matrix a, an n by n nonsingular diagonal
         *     matrix d, an m-vector b, and a positive number delta,
         *     the problem is to determine a value for the parameter
         *     par such that if x solves the system
         *
         *          a*x = b  and  sqrt(par)*d*x = 0
         *
         *     in the least squares sense, and dxnorm is the euclidean
         *     norm of d*x, then either par=0 and (dxnorm-delta) < 0.1*delta,
         *     or par>0 and abs(dxnorm-delta) < 0.1*delta.
         *
         *     Using lm_qrsolv, this subroutine completes the solution of the problem
         *     if it is provided with the necessary information from the
         *     qr factorization, with column pivoting, of a. That is, if
         *     a*p = q*r, where p is a permutation matrix, q has orthogonal
         *     columns, and r is an upper triangular matrix with diagonal
         *     elements of nonincreasing magnitude, then lmpar expects
         *     the full upper triangle of r, the permutation matrix p,
         *     and the first n components of qT*b. On output
         *     lmpar also provides an upper triangular matrix s such that
         *
         *          p^T*(a^T*a + par*d*d)*p = s^T*s.
         *
         *     s is employed within lmpar and may be of separate interest.
         *
         *     Only a few iterations are generally needed for convergence
         *     of the algorithm. If, however, the limit of 10 iterations
         *     is reached, then the output par will contain the best
         *     value obtained so far.
         *
         *     parameters:
         *
         *      n is a positive integer input variable set to the order of r.
         *
         *      r is an n by n array. on input the full upper triangle
         *        must contain the full upper triangle of the matrix r.
         *        on OUTPUT the full upper triangle is unaltered, and the
         *        strict lower triangle contains the strict upper triangle
         *        (transposed) of the upper triangular matrix s.
         *
         *      ldr is a positive integer input variable not less than n
         *        which specifies the leading dimension of the array r.
         *
         *      ipvt is an integer input array of length n which defines the
         *        permutation matrix p such that a*p = q*r. column j of p
         *        is column ipvt(j) of the identity matrix.
         *
         *      diag is an input array of length n which must contain the
         *        diagonal elements of the matrix d.
         *
         *      qtb is an input array of length n which must contain the first
         *        n elements of the vector (q transpose)*b.
         *
         *      delta is a positive input variable which specifies an upper
         *        bound on the euclidean norm of d*x.
         *
         *      par is a nonnegative variable. on input par contains an
         *        initial estimate of the levenberg-marquardt parameter.
         *        on OUTPUT par contains the final estimate.
         *
         *      x is an OUTPUT array of length n which contains the least
         *        squares solution of the system a*x = b, sqrt(par)*d*x = 0,
         *        for the output par.
         *
         *      sdiag is an array of length n needed as workspace; on OUTPUT
         *        it contains the diagonal elements of the upper triangular matrix s.
         *
         *      aux is a multi-purpose work array of length n.
         *
         *      xdi is a work array of length n. On OUTPUT: diag[j] * x[j].
         *
         */
        private static void LMParameter(int n, double[] r, int ldr, int[] ipvt, double[] diag,
                      double[] qtb, double delta, ref double par, double[] x,
                      double[] sdiag, double[] aux, double[] xdi) 
        {
            int i, iter, j, nsing;
            double dxnorm, fp, fp_old, gnorm, parc, parl, paru;
            double sum, temp;
            
            /*** lmpar: compute and store in x the gauss-newton direction. if the
                 jacobian is rank-deficient, obtain a least squares solution. ***/
            nsing = n;
            for (j = 0; j < n; j++) {
                aux[j] = qtb[j];
                if (r[j * ldr + j] == 0 && nsing == n)
                    nsing = j;
                if (nsing < n)
                    aux[j] = 0;
            }
            for (j = nsing - 1; j >= 0; j--) {
                aux[j] = aux[j] / r[j + ldr * j];
                temp = aux[j];
                for (i = 0; i < j; i++)
                    aux[i] -= r[j * ldr + i] * temp;
            }

            for (j = 0; j < n; j++)
                x[ipvt[j]] = aux[j];

            /*** lmpar: initialize the iteration counter, evaluate the function at the
                 origin, and test for acceptance of the gauss-newton direction. ***/
            for (j = 0; j < n; j++)
                xdi[j] = diag[j] * x[j];
            dxnorm = EuclideanNorm(n, 0, xdi);
            fp = dxnorm - delta;
            if (fp <= p1 * delta) {
                par = 0;
                return;
            }

            /*** lmpar: if the jacobian is not rank deficient, the newton
                 step provides a lower bound, parl, for the 0. of
                 the function. otherwise set this bound to 0.. ***/
            parl = 0;
            if (nsing >= n) {
                for (j = 0; j < n; j++)
                    aux[j] = diag[ipvt[j]] * xdi[ipvt[j]] / dxnorm;

                for (j = 0; j < n; j++) {
                    sum = 0.0;
                    for (i = 0; i < j; i++)
                        sum += r[j * ldr + i] * aux[i];
                    aux[j] = (aux[j] - sum) / r[j + ldr * j];
                }
                temp = EuclideanNorm(n, 0, aux);
                parl = fp / delta / temp / temp;
            }

            /*** lmpar: calculate an upper bound, paru, for the 0. of the function. ***/
            for (j = 0; j < n; j++) {
                sum = 0;
                for (i = 0; i <= j; i++)
                    sum += r[j * ldr + i] * qtb[i];
                aux[j] = sum / diag[ipvt[j]];
            }
            gnorm = EuclideanNorm(n, 0, aux);
            paru = gnorm / delta;
            if (paru == 0.0)
                paru = LM_DWARF / Math.Min(delta, p1);

            /*** lmpar: if the input par lies outside of the interval (parl,paru),
                 set par to the closer endpoint. ***/
            par = Math.Max(par, parl);
            par = Math.Min(par, paru);
            if (par == 0.0)
                par = gnorm / dxnorm;

            /*** lmpar: iterate. ***/
            for (iter = 0; ; iter++) {
                /** evaluate the function at the current value of par. **/
                if (par == 0.0)
                    par = Math.Max(LM_DWARF, 0.001 * paru);
                temp = Math.Sqrt(par);
                for (j = 0; j < n; j++)
                    aux[j] = temp * diag[j];

                QRSolve(r, ldr, ipvt, aux, qtb, x, sdiag, xdi);
                /* return values are r, x, sdiag */

                for (j = 0; j < n; j++)
                    xdi[j] = diag[j] * x[j]; /* used as output */
                dxnorm = EuclideanNorm(n, 0, xdi);
                fp_old = fp;
                fp = dxnorm - delta;

                /** if the function is small enough, accept the current value
                    of par. Also test for the exceptional cases where parl
                    is zero or the number of iterations has reached 10. **/
                if (Math.Abs(fp) <= p1 * delta
                    || (parl == 0.0 && fp <= fp_old && fp_old < 0.0)
                    || iter == 10) {
                    break; /* the only exit from the iteration. */
                }

                /** compute the Newton correction. **/
                for (j = 0; j < n; j++)
                    aux[j] = diag[ipvt[j]] * xdi[ipvt[j]] / dxnorm;

                for (j = 0; j < n; j++) {
                    aux[j] = aux[j] / sdiag[j];
                    for (i = j + 1; i < n; i++)
                        aux[i] -= r[j * ldr + i] * aux[j];
                }
                temp = EuclideanNorm(n, 0, aux);
                parc = fp / delta / temp / temp;

                /** depending on the sign of the function, update parl or paru. **/
                if (fp > 0)
                    parl = Math.Max(parl, par);
                else if (fp < 0)
                    paru = Math.Min(paru, par);
                /* the case fp==0 is precluded by the break condition  */

                /** compute an improved estimate for par. **/
                par = Math.Max(parl, par + parc);
            }
        }

        /*****************************************************************************/
        /*  lm_qrfac (QR factorization, from lapack)                                 */
        /*****************************************************************************/
        /*
         *     This subroutine uses Householder transformations with column
         *     pivoting (optional) to compute a qr factorization of the
         *     m by n matrix a. That is, qrfac determines an orthogonal
         *     matrix q, a permutation matrix p, and an upper trapezoidal
         *     matrix r with diagonal elements of nonincreasing magnitude,
         *     such that a*p = q*r. The Householder transformation for
         *     column k, k = 1,2,...,min(m,n), is of the form
         *
         *          i - (1/u(k))*u*uT
         *
         *     where u has zeroes in the first k-1 positions. The form of
         *     this transformation and the method of pivoting first
         *     appeared in the corresponding linpack subroutine.
         *
         *     Parameters:
         *
         *      m is a positive integer input variable set to the number
         *        of rows of a.
         *
         *      n is a positive integer input variable set to the number
         *        of columns of a.
         *
         *      a is an m by n array. On input a contains the matrix for
         *        which the qr factorization is to be computed. On OUTPUT
         *        the strict upper trapezoidal part of a contains the strict
         *        upper trapezoidal part of r, and the lower trapezoidal
         *        part of a contains a factored form of q (the non-trivial
         *        elements of the u vectors described above).
         *
         *      ipvt is an integer OUTPUT array of length lipvt. This array
         *        defines the permutation matrix p such that a*p = q*r.
         *        Column j of p is column ipvt(j) of the identity matrix.
         *
         *      rdiag is an OUTPUT array of length n which contains the
         *        diagonal elements of r.
         *
         *      acnorm is an OUTPUT array of length n which contains the
         *        norms of the corresponding columns of the input matrix a.
         *        If this information is not needed, then acnorm can coincide
         *        with rdiag.
         *
         *      wa is a work array of length n.
         *
         */
        private static void QRFactorization(int m, int n, double[] a, int[] ipvt, double[] rdiag, double[] acnorm, double[] wa) {
            int k, kmax, minmn;
            double ajnorm, sum, temp;

            /*** qrfac: compute initial column norms and initialize several arrays. ***/
            for (var j = 0; j < n; j++) {
                acnorm[j] = EuclideanNorm(m, j * m, a);
                rdiag[j] = acnorm[j];
                wa[j] = rdiag[j];
                ipvt[j] = j;
            }

            /*** qrfac: reduce a to r with Householder transformations. ***/
            minmn = Math.Min(m, n);
            for (var j = 0; j < minmn; j++) {
                /** bring the column of largest norm into the pivot position. **/
                kmax = j;
                for (k = j + 1; k < n; k++)
                    if (rdiag[k] > rdiag[kmax])
                        kmax = k;
                if (kmax == j)
                    goto pivot_ok;

                for (var i = 0; i < m; i++) {
                    temp = a[j * m + i];
                    a[j * m + i] = a[kmax * m + i];
                    a[kmax * m + i] = temp;
                }
                rdiag[kmax] = rdiag[j];
                wa[kmax] = wa[j];
                k = ipvt[j];
                ipvt[j] = ipvt[kmax];
                ipvt[kmax] = k;

            pivot_ok:
                /** compute the Householder transformation to reduce the
                    j-th column of a to a multiple of the j-th unit vector. **/
                ajnorm = EuclideanNorm(m - j, j * m + j, a);
                if (ajnorm == 0.0) {
                    rdiag[j] = 0;
                    continue;
                }

                if (a[j * m + j] < 0.0)
                    ajnorm = -ajnorm;
                for (var i = j; i < m; i++)
                    a[j * m + i] /= ajnorm;
                a[j * m + j] += 1;

                /** apply the transformation to the remaining columns
                    and update the norms. **/
                for (k = j + 1; k < n; k++) {
                    sum = 0;

                    for (var i = j; i < m; i++)
                        sum += a[j * m + i] * a[k * m + i];

                    temp = sum / a[j + m * j];

                    for (var i = j; i < m; i++)
                        a[k * m + i] -= temp * a[j * m + i];

                    if (rdiag[k] != 0.0) {
                        temp = a[m * k + j] / rdiag[k];
                        temp = Math.Max(0.0, 1 - temp * temp);
                        rdiag[k] *= Math.Sqrt(temp);
                        temp = rdiag[k] / wa[k];
                        if (0.05 * temp * temp <= LM_MACHEP) {
                            rdiag[k] = EuclideanNorm(m - j - 1, m * k + j + 1, a);
                            wa[k] = rdiag[k];
                        }
                    }
                }

                rdiag[j] = -ajnorm;
            }
        }

        /*
         *     Given an m by n matrix a, an n by n diagonal matrix d,
         *     and an m-vector b, the problem is to determine an x which
         *     solves the system
         *
         *          a*x = b  and  d*x = 0
         *
         *     in the least squares sense.
         *
         *     This subroutine completes the solution of the problem
         *     if it is provided with the necessary information from the
         *     qr factorization, with column pivoting, of a. That is, if
         *     a*p = q*r, where p is a permutation matrix, q has orthogonal
         *     columns, and r is an upper triangular matrix with diagonal
         *     elements of nonincreasing magnitude, then qrsolv expects
         *     the full upper triangle of r, the permutation matrix p,
         *     and the first n components of (q transpose)*b. The system
         *     a*x = b, d*x = 0, is then equivalent to
         *
         *          r*z = q^T*b,  p^T*d*p*z = 0,
         *
         *     where x = p*z. If this system does not have full rank,
         *     then a least squares solution is obtained. On output qrsolv
         *     also provides an upper triangular matrix s such that
         *
         *          p^T *(a^T *a + d*d)*p = s^T *s.
         *
         *     s is computed within qrsolv and may be of separate interest.
         *
         *     Parameters
         *
         *      n is a positive integer input variable set to the order of r.
         *
         *      r is an n by n array. On input the full upper triangle
         *        must contain the full upper triangle of the matrix r.
         *        On OUTPUT the full upper triangle is unaltered, and the
         *        strict lower triangle contains the strict upper triangle
         *        (transposed) of the upper triangular matrix s.
         *
         *      ldr is a positive integer input variable not less than n
         *        which specifies the leading dimension of the array r.
         *
         *      ipvt is an integer input array of length n which defines the
         *        permutation matrix p such that a*p = q*r. Column j of p
         *        is column ipvt(j) of the identity matrix.
         *
         *      diag is an input array of length n which must contain the
         *        diagonal elements of the matrix d.
         *
         *      qtb is an input array of length n which must contain the first
         *        n elements of the vector (q transpose)*b.
         *
         *      x is an OUTPUT array of length n which contains the least
         *        squares solution of the system a*x = b, d*x = 0.
         *
         *      sdiag is an OUTPUT array of length n which contains the
         *        diagonal elements of the upper triangular matrix s.
         *
         *      wa is a work array of length n.
         *
         */
        private static void QRSolve(double[] r, int ldr, int[] ipvt, double[] diag, double[] qtb, double[] x, double[] sdiag, double[] wa) {
            var n = x.Length;

            int kk;
            double qtbpj, sum, temp;
            double _sin, _cos, _tan, _cot; /* local variables, not functions */

            /*** qrsolv: copy r and q^T*b to preserve input and initialize s.
                 in particular, save the diagonal elements of r in x. ***/
            for (var j = 0; j < n; j++) {
                for (var i = j; i < n; i++)
                    r[j * ldr + i] = r[i * ldr + j];
                x[j] = r[j * ldr + j];
                wa[j] = qtb[j];
            }

            /*** qrsolv: eliminate the diagonal matrix d using a Givens rotation. ***/
            for (var j = 0; j < n; j++) {
                /*** qrsolv: prepare the row of d to be eliminated, locating the
                     diagonal element using p from the qr factorization. ***/
                if (diag[ipvt[j]] == 0.0)
                    goto L90;
                for (var k = j; k < n; k++)
                    sdiag[k] = 0.0;
                sdiag[j] = diag[ipvt[j]];
                /*** qrsolv: the transformations to eliminate the row of d modify only 
                     a single element of qT*b beyond the first n, which is initially 0. ***/
                qtbpj = 0.0;
                for (var k = j; k < n; k++) {
                    /** determine a Givens rotation which eliminates the
                        appropriate element in the current row of d. **/
                    if (sdiag[k] == 0.0)
                        continue;
                    kk = k + ldr * k;
                    if (Math.Abs(r[kk]) < Math.Abs(sdiag[k])) {
                        _cot = r[kk] / sdiag[k];
                        _sin = 1 / Math.Sqrt(1 + _cot * _cot);
                        _cos = _sin * _cot;
                    }
                    else {
                        _tan = sdiag[k] / r[kk];
                        _cos = 1 / Math.Sqrt(1 + _tan * _tan);
                        _sin = _cos * _tan;
                    }
                    /** compute the modified diagonal element of r and
                        the modified element of ((q^T)*b,0). **/
                    r[kk] = _cos * r[kk] + _sin * sdiag[k];
                    temp = _cos * wa[k] + _sin * qtbpj;
                    qtbpj = -_sin * wa[k] + _cos * qtbpj;
                    wa[k] = temp;
                    /** accumulate the tranformation in the row of s. **/
                    for (var i = k + 1; i < n; i++) {
                        temp = _cos * r[k * ldr + i] + _sin * sdiag[i];
                        sdiag[i] = -_sin * r[k * ldr + i] + _cos * sdiag[i];
                        r[k * ldr + i] = temp;
                    }
                }

            L90:
                /** store the diagonal element of s and restore
                    the corresponding diagonal element of r. **/
                sdiag[j] = r[j * ldr + j];
                r[j * ldr + j] = x[j];
            }

            /*** qrsolv: solve the triangular system for z. if the system is
                singular, then obtain a least squares solution. ***/
            var nsing = n;
            for (var j = 0; j < n; j++) {
                if (sdiag[j] == 0.0 && nsing == n)
                    nsing = j;
                if (nsing < n)
                    wa[j] = 0;
            }

            for (var j = nsing - 1; j >= 0; j--) {
                sum = 0;
                for (var i = j + 1; i < nsing; i++)
                    sum += r[j * ldr + i] * wa[i];
                wa[j] = (wa[j] - sum) / sdiag[j];
            }

            /*** qrsolv: permute the components of z back to components of x. ***/
            for (var j = 0; j < n; j++)
                x[ipvt[j]] = wa[j];
        }

        /*     Given an n-vector x, this function calculates the
         *     euclidean norm of x.
         *
         *     The euclidean norm is computed by accumulating the sum of
         *     squares in three different sums. The sums of squares for the
         *     small and large components are scaled so that no overflows
         *     occur. Non-destructive underflows are permitted. Underflows
         *     and overflows do not occur in the computation of the unscaled
         *     sum of squares for the intermediate components.
         *     The definitions of small, intermediate and large components
         *     depend on two constants, LM_SQRT_DWARF and LM_SQRT_GIANT. The main
         *     restrictions on these constants are that LM_SQRT_DWARF**2 not
         *     underflow and LM_SQRT_GIANT**2 not overflow.
         */
        private static double EuclideanNorm(int n, int offset, double[] x) {
            double agiant = LM_SQRT_GIANT / n;
            double s1 = 0.0;
            double s2 = 0.0;
            double s3 = 0.0;
            double x1max = 0.0;
            double x3max = 0.0;

            double xabs;
            double temp;
            /** sum squares. **/
            for (var i = offset; i < offset + n; ++i) {
                xabs = Math.Abs(x[i]);
                if (xabs > LM_SQRT_DWARF) {
                    if (xabs < agiant) {
                        s2 += xabs * xabs;
                    }
                    else if (xabs > x1max) {
                        temp = x1max / xabs;
                        s1 = 1 + s1 * temp * temp;
                        x1max = xabs;
                    }
                    else {
                        temp = xabs / x1max;
                        s1 += temp * temp;
                    }
                }
                else if (xabs > x3max) {
                    temp = x3max / xabs;
                    s3 = 1 + s3 * temp * temp;
                    x3max = xabs;
                }
                else if (xabs != 0.0) {
                    temp = xabs / x3max;
                    s3 += temp * temp;
                }
            }

            /** calculation of norm. **/
            if (s1 != 0)
                return x1max * Math.Sqrt(s1 + (s2 / x1max) / x1max);
            else if (s2 != 0)
                if (s2 >= x3max)
                    return Math.Sqrt(s2 * (1 + (x3max / s2) * (x3max * s3)));
                else
                    return Math.Sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
            else
                return x3max * Math.Sqrt(s3);
        }
    }
}
