# Include the ctridiag() function as declared in ctridiag.h
# We use a cdef since this function will only be called from Cython
cdef extern from "cssor_src.h":
    void cssor(double* U, int m, int n, double omega, double tol, int maxiters, int* info);  # This is C function calling
# Define a Cython wrapper for cssor().
# Accept four NumPy arrays of doubles
# This will not check for out of bounds memory accesses.
cpdef cycssor(double[:,::1] U,  double omega): # [:,:] 2D matrix, This is Cython function
    cdef int n = U.shape[1]  # # of columns
    cdef int m = U.shape[0]# # of roclws
    cdef double tol = 10e-8
    cdef int maxiters = 10000
    cdef int *info = NULL
    cssor(&U[0,0], m, n, omega, tol, maxiters, info) # calling C function inside python function
    