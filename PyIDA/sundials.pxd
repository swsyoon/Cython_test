# -*- coding: utf-8 -*-
ctypedef double realtype

cdef class N_Vector_Serial:
    cdef void *thisptr
    cdef Py_ssize_t __shape[1]
    
    cdef int getsize(N_Vector_Serial self)
    cdef double *getdata(N_Vector_Serial self)
    cdef double getitem(N_Vector_Serial self, int index)
    cdef setitem(N_Vector_Serial self, int index, double value)
    cpdef copy(N_Vector_Serial self)

cdef double *N_Vector_data(void *vector)
cdef int N_Vector_length(void *vector)

cdef struct cv_userdata:
    void *RHSF
    void *ROOTF
    void *JACF
    void *SW
    void *Y
    
cdef class SolverStatistics:
    cdef readonly long int nsteps
    cdef readonly long int nfevals
    cdef readonly long int njevals
    cdef readonly long int ngevals
    cdef readonly long int netfails
    cdef readonly long int nniters
    cdef readonly long int nncfails


cdef struct ida_userdata:
    void *RESF
    void *ROOTF
    void *JACF
    void *SW
    void *Y
    void *YDOT

cdef class IDASettings:
    cdef public int maxord
    cdef public long int mxsteps
    cdef public realtype tstop
    cdef public realtype hmax
    cdef public bint suppressalg
    cdef public bint lsoff
    cdef public icopt
    cdef public object idv
    cdef public realtype reltol
    cdef public realtype abstol
    cdef public object abstolv
    
    cdef public object RES
    cdef public object ROOT
    cdef public object SW
    cdef public object JAC

cdef class IDASolver

cdef class IDAIterator:
    cdef IDASolver solver
    
    cdef public realtype dt
    cdef readonly realtype t
    
    cdef realtype tret
    cdef realtype troot
    cdef N_Vector_Serial yroot
    cdef N_Vector_Serial ydotroot
    cdef realtype hlast
    
    cdef bint root
    cdef bint stop
    cdef bint last
    
    cdef Next(self)
    
cdef class IDASolver:
    cdef void *thisptr
    cdef void *ResFn
    cdef void *RootFn
    cdef void *JacFn
    
    cdef realtype t0
    
    cdef N_Vector_Serial y
    cdef N_Vector_Serial yv
    cdef N_Vector_Serial ydot
    cdef N_Vector_Serial ydotv
    
    cdef object RESF
    cdef object ROOTF
    cdef object JACF
    cdef public object SW
    
    cdef ida_userdata userdata
    cdef public IDASettings settings
    
    cpdef GetConsistentIC(self, realtype direction)
    cpdef __handleRoot(self, realtype event_time, N_Vector_Serial yroot, N_Vector_Serial ydotroot)
    cpdef IDAIterator iter(self, realtype t0, realtype dt)
    cpdef step(self, realtype tf)
