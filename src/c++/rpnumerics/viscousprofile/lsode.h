#ifndef _LSODE_
#define _LSODE_

// Return codes. Taken from lsode.f:
//
//  on output, istate has the following values and meanings.
//   1  means nothing was done, as tout was equal to t with
//  istate = 1 on input.  (however, an internal counter was
//  set to detect and prevent repeated calls of this type.)
//   2  means the integration was performed successfully.
//  -1  means an excessive amount of work (more than mxstep
//  steps) was done on this call, before completing the
//  requested task, but the integration was otherwise
//  successful as far as t.  (mxstep is an optional input
//  and is normally 500.)  to continue, the user may
//  simply reset istate to a value .gt. 1 and call again
//  (the excess work step counter will be reset to 0).
//  in addition, the user may increase mxstep to avoid
//  this error return (see below on optional inputs).
//  -2  means too much accuracy was requested for the precision
//  of the machine being used.  this was detected before
//  completing the requested task, but the integration
//  was successful as far as t.  to continue, the tolerance
//  parameters must be reset, and istate must be set
//  to 3.  the optional output tolsf may be used for this
//  purpose.  (note.. if this condition is detected before
//  taking any steps, then an illegal input return
//  (istate = -3) occurs instead.)
//  -3  means illegal input was detected, before taking any
//  integration steps.  see written message for details.
//  note..  if the solver detects an infinite loop of calls
//  to the solver with illegal input, it will cause
//  the run to stop.
//  -4  means there were repeated error test failures on
//  one attempted step, before completing the requested
//  task, but the integration was successful as far as t.
//  the problem may have a singularity, or the input
//  may be inappropriate.
//  -5  means there were repeated convergence test failures on
//  one attempted step, before completing the requested
//  task, but the integration was successful as far as t.
//  this may be caused by an inaccurate jacobian matrix,
//  if one is being used.
//  -6  means ewt(i) became zero for some i during the
//  integration.  pure relative error control (atol(i)=0.0)
//  was requested on a variable which has now vanished.
//  the integration was successful as far as t.
//

#ifndef LSODE_NOTHING_DONE
#define LSODE_NOTHING_DONE 1
#endif

#ifndef LSODE_OK
#define LSODE_OK 2
#endif

#ifndef LSODE_EXCESSIVE_WORK
#define LSODE_EXCESSIVE_WORK (-1)
#endif

#ifndef LSODE_TOO_MUCH_ACCURACY_REQUESTED
#define LSODE_TOO_MUCH_ACCURACY_REQUESTED (-2)
#endif

#ifndef LSODE_ILLEGAL_INPUT
#define LSODE_ILLEGAL_INPUT (-3)
#endif

#ifndef LSODE_SINGULARITY
#define LSODE_SINGULARITY (-4)
#endif

#ifndef LSODE_CONVERGENCE_FAILURE
#define LSODE_CONVERGENCE_FAILURE (-5)
#endif

#ifndef LSODE_VARIABLE_VANISHED
#define LSODE_VARIABLE_VANISHED (-6)
#endif

#ifdef __cplusplus
extern "C"{
#endif
    int lsode_(int (*)(int *, double *, double *, double *, int *, double *), int *, double *, double *, double *,
            int *, double *, double *, int *, int *, int *, double *, int *,
            int *, int *, int(*)(int *, double *, double *, int *, int *, double *, int *), int *, int*, double*);
#ifdef __cplusplus
}
#endif

#endif // _LSODE_

