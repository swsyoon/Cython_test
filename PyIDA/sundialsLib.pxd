# -*- coding: utf-8 -*-
cdef extern from "sundials/sundials_types.h":
    ctypedef double realtype
    ctypedef bint booleantype

cdef extern from "sundials/sundials_nvector.h":
    cdef struct _generic_N_Vector:
        void *content
        
    ctypedef _generic_N_Vector *N_Vector
    N_Vector N_VNew_Serial(long int vec_length)
    void N_VDestroy_Serial(N_Vector v)
    void N_VPrint_Serial(N_Vector v)

cdef extern from "nvector/nvector_serial.h":
    cdef N_Vector N_VMake_Serial(long int vec_length, realtype *v_data)
    
    cdef struct _N_VectorContent_Serial:
        long int length
        realtype *data
        
    ctypedef _N_VectorContent_Serial *N_VectorContent_Serial
    

cdef extern from "sundials/sundials_direct.h":
    cdef struct _DlsMat:
        pass
        
    ctypedef _DlsMat *DlsMat
    cdef realtype* DENSE_COL(DlsMat A, int j)

cdef extern from "cvode/cvode_dense.h":
    int CVDense(void *cvode_mem, int N)
    ctypedef int (*CVDlsDenseJacFn)(int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn djac)
    
cdef extern from "ida/ida.h":
    int  IDA_NORMAL
    int  IDA_ONE_STEP
    int  IDA_YA_YDP_INIT
    int  IDA_Y_INIT
    int  IDA_SUCCESS
    int  IDA_TSTOP_RETURN
    int  IDA_ROOT_RETURN
    int  IDA_WARNING
    int  IDA_MEM_NULL
    int  IDA_ILL_INPUT
    int  IDA_NO_MALLOC
    int  IDA_TOO_MUCH_WORK
    int  IDA_TOO_MUCH_ACC
    int  IDA_ERR_FAIL
    int  IDA_CONV_FAIL
    int  IDA_LINIT_FAIL
    int  IDA_LSETUP_FAIL
    int  IDA_LSOLVE_FAIL
    int  IDA_RES_FAIL
    int  IDA_CONSTR_FAIL
    int  IDA_REP_RES_ERR
    int  IDA_MEM_FAIL
    int  IDA_BAD_T
    int  IDA_BAD_EWT
    int  IDA_FIRST_RES_FAIL
    int  IDA_LINESEARCH_FAIL
    int  IDA_NO_RECOVERY
    int  IDA_RTFUNC_FAIL
    
    ctypedef int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
    ctypedef int (*IDARootFn)(realtype t, N_Vector y, N_Vector yp, realtype *gout, void *user_data)
    
    void *IDACreate()
    int IDASetUserData(void *ida_mem, void *user_data)
    int IDASetMaxOrd(void *ida_mem, int maxord)
    int IDASetMaxNumSteps(void *ida_mem, long int mxsteps)
    int IDASetInitStep(void *ida_mem, realtype hin)
    int IDASetMaxStep(void *ida_mem, realtype hmax)
    int IDASetStopTime(void *ida_mem, realtype tstop)
    int IDASetNonlinConvCoef(void *ida_mem, realtype epcon)
    int IDASetMaxErrTestFails(void *ida_mem, int maxnef)
    int IDASetMaxNonlinIters(void *ida_mem, int maxcor)
    int IDASetMaxConvFails(void *ida_mem, int maxncf)
    int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg)
    int IDASetId(void *ida_mem, N_Vector id)
    int IDASetConstraints(void *ida_mem, N_Vector constraints)
    int IDASetRootDirection(void *ida_mem, int *rootdir)
    int IDASetNoInactiveRootWarn(void *ida_mem)
    int IDAInit(void *ida_mem, IDAResFn res, realtype t0, N_Vector yy0, N_Vector yp0)
    int IDAReInit(void *ida_mem, realtype t0, N_Vector yy0, N_Vector yp0)
    int IDASStolerances(void *ida_mem, realtype reltol, realtype abstol)
    int IDASVtolerances(void *ida_mem, realtype reltol, N_Vector abstol)
    int IDASetNonlinConvCoefIC(void *ida_mem, realtype epiccon)
    int IDASetMaxNumStepsIC(void *ida_mem, int maxnh)
    int IDASetMaxNumJacsIC(void *ida_mem, int maxnj)
    int IDASetMaxNumItersIC(void *ida_mem, int maxnit)
    int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff)
    int IDASetStepToleranceIC(void *ida_mem, realtype steptol)
    int IDARootInit(void *ida_mem, int nrtfn, IDARootFn g)
    int IDACalcIC(void *ida_mem, int icopt, realtype tout1)
    int IDASolve(void *ida_mem, realtype tout, realtype *tret, N_Vector yret, N_Vector ypret, int itask)
    int IDAGetSolution(void *ida_mem, realtype t,  N_Vector yret, N_Vector ypret)
    int IDAGetWorkSpace(void *ida_mem, long int *lenrw, long int *leniw)
    int IDAGetNumSteps(void *ida_mem, long int *nsteps)
    int IDAGetNumResEvals(void *ida_mem, long int *nrevals)
    int IDAGetNumLinSolvSetups(void *ida_mem, long int *nlinsetups)
    int IDAGetNumErrTestFails(void *ida_mem, long int *netfails)
    int IDAGetNumBacktrackOps(void *ida_mem, long int *nbacktr)
    int IDAGetConsistentIC(void *ida_mem, N_Vector yy0_mod, N_Vector yp0_mod)
    int IDAGetLastOrder(void *ida_mem, int *klast)
    int IDAGetCurrentOrder(void *ida_mem, int *kcur)
    int IDAGetActualInitStep(void *ida_mem, realtype *hinused)
    int IDAGetLastStep(void *ida_mem, realtype *hlast)
    int IDAGetCurrentStep(void *ida_mem, realtype *hcur)
    int IDAGetCurrentTime(void *ida_mem, realtype *tcur)
    int IDAGetTolScaleFactor(void *ida_mem, realtype *tolsfact)
    int IDAGetErrWeights(void *ida_mem, N_Vector eweight)
    int IDAGetEstLocalErrors(void *ida_mem, N_Vector ele)
    int IDAGetNumGEvals(void *ida_mem, long int *ngevals)
    int IDADlsGetNumJacEvals(void *ida_mem, long int *njevals)
    int IDADlsGetNumResEvals(void *ida_mem, long int *nrevalsLS)
    int IDAGetRootInfo(void *ida_mem, int *rootsfound)
    int IDAGetIntegratorStats(void *ida_mem, long int *nsteps, 
                                          long int *nrevals, long int *nlinsetups, 
                                          long int *netfails, int *qlast, int *qcur, 
                                          realtype *hinused, realtype *hlast, realtype *hcur, 
                                          realtype *tcur)
    int IDAGetNumNonlinSolvIters(void *ida_mem, long int *nniters)
    int IDAGetNumNonlinSolvConvFails(void *ida_mem, long int *nncfails)
    int IDAGetNonlinSolvStats(void *ida_mem, long int *nniters,  long int *nncfails)
    char *IDAGetReturnFlagName(int flag)
    void IDAFree(void **ida_mem)

cdef extern from "ida/ida_dense.h":
    int IDADense(void *ida_mem, int Neq)
    
    ctypedef int (*IDADlsDenseJacFn)(int Neq, realtype tt, realtype cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int IDADlsSetDenseJacFn(void *ida_mem, IDADlsDenseJacFn djac)
