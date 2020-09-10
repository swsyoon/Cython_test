# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 19:05:36 2020

@author: swsyo
"""
cdef extern: 
    int IDACalcIC(void *ida_mem, int icopt, real tout1, real epicfac, 
               int maxnh, int maxnj, int maxnit, int lsoff, real steptol)

    
cpdef CyIDACalcIC(void *ida_mem, icopt, tout1, epicfac, maxnh, maxnj, maxnit, lsoff, steptol):
    
    cdef int icopt
    cdef int maxnh
    cdef int maxnj
    cdef int maxnit
    cdef int lsoff
    cdef real tout1
    cdef real epicfac
    cdef real steptol
    
    IDACalcIC(void *ida_mem, icopt,  tout1,  epicfac, maxnh,  maxnj,  maxnit, lsoff, steptol)