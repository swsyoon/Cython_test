# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 19:05:36 2020

@author: swsyo
"""
#https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html?highlight=pointer#pass-data-from-a-c-function-via-pointer

from libc.stdlib cimport malloc, free
from sundialsLib cimport *
from sundials cimport *

class SundialsError(Exception):
    pass

cdef class SolverStatistics:
    """Class to encapsulate statistics for solver"""
    def __repr__(self):
        s = []
        s.append("nsteps = %d"      % self.nsteps)
        s.append("nfevals = %d"     % self.nfevals)
        s.append("njevals = %d"     % self.njevals)
        s.append("ngevals = %d"     % self.ngevals)
        s.append("netfails = %d"    % self.netfails)
        s.append("nniters = %d"     % self.nniters)
        s.append("nncfails = %d"    % self.nncfails)
    
        return "(%s)" % (", ".join(s))
    
    def __str__(self):
        return '%s%s' % (self.__class__.__name__, self.__repr__())


cdef int ida_res(realtype t, N_Vector yv, N_Vector yvdot, N_Vector residual, void* user_data) with gil:
    """
    Wraps  Python res-callback function to obtain IDA required interface
    see also ctypedef statement above
    """
    cdef object[double, ndim=1, mode="c"] buf
    cdef ida_userdata *ud = <ida_userdata *>user_data
    
    cdef double *data = (<N_VectorContent_Serial>residual.content).data
    cdef int i, size = (<N_VectorContent_Serial>residual.content).length
    
    # copy data
    cdef N_Vector_Serial Y = <N_Vector_Serial>ud.Y
    for i in range(Y.getsize()):
        Y.setitem(i, (<N_VectorContent_Serial>yv.content).data[i])
    
    cdef N_Vector_Serial YDOT = <N_Vector_Serial>ud.YDOT
    for i in range(YDOT.getsize()):
        YDOT.setitem(i, (<N_VectorContent_Serial>yvdot.content).data[i])
        
    if not ud.SW == NULL:
        res = (<object>ud.RESF)(t, Y, YDOT, <object>ud.SW)
    else:
        res = (<object>ud.RESF)(t, Y, YDOT)
    
    try:
        buf = res
    except TypeError:
        for i in range(size):
            data[i] = res[i]
    else:
        for i in range(size):
            data[i] = buf[i]
    
    return 0

cdef int ida_root(realtype t, N_Vector yv, N_Vector yvdot, realtype *gout,
                  void* user_data) with gil:
    """
    Wraps  Python root-callback function to obtain IDA required interface
    see also ctypedef statement above
    """
    cdef object[double, ndim=1, mode="c"] buf
    cdef ida_userdata *ud = <ida_userdata *>user_data
    cdef int i
    
    # copy data
    cdef N_Vector_Serial Y = <N_Vector_Serial>ud.Y
    for i in range(Y.getsize()):
        Y.setitem(i, (<N_VectorContent_Serial>yv.content).data[i])
    
    cdef N_Vector_Serial YDOT = <N_Vector_Serial>ud.YDOT
    for i in range(YDOT.getsize()):
        YDOT.setitem(i, (<N_VectorContent_Serial>yvdot.content).data[i])
        
    out = (<object>ud.ROOTF)(t, Y, YDOT, <object>ud.SW)
    
    try:
        buf = out
    except TypeError:
        for i in range(len(out)):
            gout[i] = out[i]
    else:
        for i in range(len(out)):
            gout[i] = buf[i]
            
    return 0

cdef int ida_jac(int Neq, realtype t, realtype c, N_Vector yv, N_Vector yvdot,
                 N_Vector residual, DlsMat Jac, void* user_data, N_Vector tmp1,
                 N_Vector tmp2, N_Vector tmp3) with gil:
    """
    Wraps Python jacobian-callback function to obtain IDA required interface.
    """
    cdef object[double, ndim=2, mode="c"] buf
    cdef ida_userdata *ud = <ida_userdata *>user_data
    cdef realtype* col_i
    cdef int i, j
    
    # copy data
    cdef N_Vector_Serial Y = <N_Vector_Serial>ud.Y
    for i in range(Y.getsize()):
        Y.setitem(i, (<N_VectorContent_Serial>yv.content).data[i])
    
    cdef N_Vector_Serial YDOT = <N_Vector_Serial>ud.YDOT
    for i in range(YDOT.getsize()):
        YDOT.setitem(i, (<N_VectorContent_Serial>yvdot.content).data[i])
        
    if not ud.SW == NULL:
        jac = (<object>ud.JACF)(c, t, Y, YDOT, <object>ud.SW)
    else:
        jac = (<object>ud.JACF)(c, t, Y, YDOT)
    
    try:
        buf = jac
    except TypeError:
        # list/tuple etc.
        for i in range(Neq):
            col_i = DENSE_COL(Jac, i)
            for j in range(Neq):
                col_i[j] = jac[j][i]
    else:
        # fast c array
        for i in range(Neq):
            col_i = DENSE_COL(Jac, i)
            for j in range(Neq):
                col_i[j] = buf[j,i]
                
    return 0
    
class IDAError(SundialsError):
    def __init__(self, value = 0, msg = None):
        self.value = value
        self.msg = msg
        
    def __str__(self):
        if self.value >= 0:
            return repr(self.msg)
        else:
            return IDAGetReturnFlagName(self.value)

class IDARootException(IDAError):
    def __init__(self, realtype t, y, ydot, SW):
        self.t = t
        self.y = y
        self.ydot = ydot
        self.SW = SW
        
    def __str__(self):
        args = self.t, repr(self.SW)
        return "CVodeRoot(t = %f, SW = %s)" % args
        
cdef class IDASettings:
    """Class to encapsulate settings for IDA"""
    def __init__(self, **kwargs):
        # set default
        self.maxord = -1
        self.mxsteps = 500
        self.tstop = float('+inf')
        self.hmax = float('+inf')
        
        self.suppressalg = False
        
        self.lsoff = False
        self.icopt = 'IDA_Y_INIT'
        
        self.reltol = -1.0
        self.abstol = -1.0
        self.abstolv = None
        
        self.RES = None
        self.ROOT = None
        self.SW = None
        self.JAC = None
        
        # process arguments
        self.update(**kwargs)
            
        if self.maxord == -1:
            self.maxord = 5
                            
        if self.reltol == -1.0:
            raise IDAError(0, "IDASettings: 'reltol' parameter must be given")
        
        if self.abstol == -1.0 and self.abstolv is None:
            raise IDAError(0, "IDASettings: 'abstol' parameter must be given")
                          
        if not self.SW is None:
            self.SW = [bool(value) for value in self.SW]
    
    def __repr__(self):
        s = []
        s.append("maxord = %d" % self.maxord)
        s.append("mxsteps = %d" % self.mxsteps)
        s.append("tstop = %g" % self.tstop)
        s.append("hmax = %g" % self.hmax)
        s.append("suppressalg = %s" % self.suppressalg)
        s.append("lsoff = %s" % self.lsoff)
        s.append("icopt = %s" % self.icopt)
        s.append("reltol = %g" % self.reltol)
        s.append("abstol = %s" % (self.abstolv or self.abstol))
        
        return "(%s)" % (", ".join(s))
    
    def __str__(self):
        return '%s%s' % (self.__class__.__name__, self.__repr__())
        
    def update(self, **kwargs):
        cdef int size, i
        
        for name, value in kwargs.iteritems():
            if name == 'abstol':
                try:
                    size = len(value)
                except TypeError:
                    self.abstol = value
                else:
                    self.abstolv = N_Vector_Serial(size)
                    for i in range(size):
                        self.abstolv[i] = value[i]
                        
            elif name == "id":
                size = len(value)
                self.idv = N_Vector_Serial(size)
                for i in range(size):
                    self.idv[i] = value[i]
            else:
                try:
                    setattr(self, name, value)
                except AttributeError:
                    raise IDAError(0, "CVodeSettings: Unknown parameter: '%s'" % name)

cdef class IDAIterator:
    """
    Single step iterator
    """
    def __init__(self, solver, realtype t0, realtype dt):
        self.solver = solver
        self.dt = dt
        self.t = t0 + dt
        
        # ensure first step is run
        self.tret = t0
        self.hlast = 0.
        
        self.root = False
        self.stop = False
        self.last = False
        
    def __iter__(self):
        return self
    
    cdef Next(self):
        """
        Interpolate result at time t and return
        solution vector. The time must not be
        past the current internal time of the solver.
        """
        cdef N_Vector_Serial y, ydot
        cdef realtype tcur
        cdef int flag
        
        y = N_Vector_Serial(self.solver.y.getsize())
        ydot = N_Vector_Serial(self.solver.ydot.getsize())
        
        while True:
            if self.root and self.t >= self.troot:
                self.root = False
                self.solver.__handleRoot(self.troot, self.yroot, self.ydotroot)
                
            elif (self.tret - self.hlast) <= self.t <= self.tret:
                tcur = self.t
                self.t += self.dt
                
                if self.stop:
                    # last step
                    self.last = True
                
                flag = IDAGetSolution(self.solver.thisptr, tcur, <N_Vector>y.thisptr, <N_Vector>ydot.thisptr)
                if flag != IDA_SUCCESS:
                    raise IDAError(flag)
                
                break
                
            else:
                y = self.solver.y
                ydot = self.solver.ydot
                flag = IDASolve(self.solver.thisptr, self.t, &self.tret,
                                <N_Vector>y.thisptr,
                                <N_Vector>ydot.thisptr, IDA_ONE_STEP)
            
                if flag < 0:
                    raise IDAError(flag)
                
                IDAGetLastStep(self.solver.thisptr, &self.hlast)
                
                if flag == IDA_ROOT_RETURN:
                    self.root = True
                    self.troot = self.tret
                    self.yroot = y.copy()
                    self.ydotroot = ydot.copy()
    
                elif flag == IDA_TSTOP_RETURN:
                    self.stop = True
                    
                continue
        
        return y, ydot
        
    def __next__(self):
        cdef N_Vector_Serial y, ydot
        cdef realtype tcur
        
        if self.last:
            raise StopIteration()
        
        # save current step
        tcur = self.t
        
        y, ydot = self.Next()
        
        return tcur, y, ydot
        
cdef class IDASolver:
    """Class to wrap IDA"""
    def __init__(self, settings = None, **kwargs):
        if settings is None:
            self.settings = IDASettings(**kwargs)
        else:
            self.settings = settings
            
            if kwargs:
                settings.update(**kwargs)
                
        self.ResFn = <void *>ida_res
        self.RootFn = <void *>ida_root
        self.JacFn = <void *>ida_jac

    def __dealloc__(self):
        if self.thisptr != NULL:
            IDAFree(&self.thisptr)
    
    def stat(self):
        """
        Return solver statistics
        """
        cdef SolverStatistics stat = SolverStatistics()
    
        IDAGetNumSteps(self.thisptr, &stat.nsteps)
        IDAGetNumResEvals(self.thisptr, &stat.nfevals)
        IDADlsGetNumJacEvals(self.thisptr, &stat.njevals)
        IDAGetNumGEvals(self.thisptr, &stat.ngevals)
        IDAGetNumErrTestFails(self.thisptr, &stat.netfails)
        IDAGetNumNonlinSolvIters(self.thisptr, &stat.nniters)
        IDAGetNumNonlinSolvConvFails(self.thisptr, &stat.nncfails)
        
        return stat
    
    def init(self, realtype t0, y0, yd0, **kwargs):
        """
        Initialize the solver or restart the solver
        """
        cdef IDASettings settings
        cdef N_Vector_Serial abstolv, idv
        cdef int i, flag
        
        # update settings
        settings = self.settings
        if kwargs:
            settings.update(**kwargs)
        
        # store initial values
        self.t0 = t0
        
        try:
            # create from buffer type: numpy/array etc.
            self.y = N_Vector_Serial_Copy_Array(y0)
        except TypeError:
            # create from list/tuple etc.
            y0 = [float(value) for value in y0]
            
            self.y = N_Vector_Serial(len(y0))
            for i in range(self.y.size):
                self.y[i] = y0[i]
                
        try:
            # create from buffer type: numpy/array etc.
            self.ydot = N_Vector_Serial_Copy_Array(yd0)
        except TypeError:
            # create from list/tuple etc.
            yd0 = [float(value) for value in yd0]
            
            self.ydot = N_Vector_Serial(len(yd0))
            for i in range(self.ydot.size):
                self.ydot[i] = yd0[i]
        
        # Initialization
        if self.thisptr == NULL:
            self.thisptr = IDACreate()
            if self.thisptr == NULL:
                raise IDAError(IDA_MEM_NULL)
            
            flag = IDAInit(self.thisptr, <IDAResFn>self.ResFn, self.t0, <N_Vector>self.y.thisptr,
                           <N_Vector>self.ydot.thisptr)
            if flag != IDA_SUCCESS:
                raise IDAError(flag)
            
            flag = IDASetMaxOrd(self.thisptr, settings.maxord)
            if flag != IDA_SUCCESS:
                raise IDAError(flag)
                
            flag = IDADense(self.thisptr, self.y.size)
            if flag != IDA_SUCCESS:
                raise IDAError(flag)
                
            self.RESF = settings.RES
            self.userdata = ida_userdata(<void*>settings.RES, <void*>NULL,
                                         <void*>NULL, <void*>NULL, <void*>NULL,
                                         <void*>NULL)
                                        
            if not settings.ROOT is None:
                self.ROOTF = settings.ROOT
                self.userdata.ROOTF = <void*>settings.ROOT
            
                if not settings.SW is None:
                    self.SW = settings.SW
                    self.userdata.SW = <void*>settings.SW
                else:
                    raise IDAError(0, "Expected switch vector")
                
                flag = IDARootInit(self.thisptr, len(settings.SW), <IDARootFn>self.RootFn)
                if flag != IDA_SUCCESS:
                    raise IDAError(flag)
            
            if not settings.JAC is None:
                self.JACF = settings.JAC
                self.userdata.JACF = <void*>settings.JAC
                
                flag = IDADlsSetDenseJacFn(self.thisptr, <IDADlsDenseJacFn>self.JacFn)
                if flag != IDA_SUCCESS:
                    raise IDAError(flag)
                    
            self.yv = self.y.copy()
            self.userdata.Y = <void*>self.yv
            
            self.ydotv = self.ydot.copy()
            self.userdata.YDOT = <void*>self.ydotv
        
        # Reinitialization
        else:
            flag = IDAReInit(self.thisptr, self.t0, <N_Vector>self.y.thisptr,
                             <N_Vector>self.ydot.thisptr)
            if flag != IDA_SUCCESS:
                raise IDAError(flag)
            
            if self.RESF != settings.RES:
                self.RESF = settings.RES
                self.userdata.RESF = <void*>settings.RES
            
            if self.ROOTF != settings.ROOT:
                self.ROOTF = settings.ROOT
                self.userdata.ROOTF = <void*>settings.ROOT
            
            if self.SW != settings.SW:
                self.SW = settings.SW
                self.userdata.SW = <void*>settings.SW
            
            if self.JACF != settings.JAC:
                self.JACF = settings.JAC
                self.userdata.JACF = <void*>settings.JAC
        
        if settings.abstolv is None:
            flag = IDASStolerances(self.thisptr, settings.reltol, settings.abstol)
        else:
            abstolv = settings.abstolv
            flag = IDASVtolerances(self.thisptr, settings.reltol, <N_Vector>abstolv.thisptr)
            
        if flag != IDA_SUCCESS:
            raise IDAError(flag)
        
        flag = IDASetMaxNumSteps(self.thisptr, settings.mxsteps)
        if flag != IDA_SUCCESS:
            raise IDAError(flag)
                
        flag = IDASetMaxStep(self.thisptr, settings.hmax)
        if flag != IDA_SUCCESS:
            raise IDAError(flag)
        
        flag = IDASetSuppressAlg(self.thisptr, settings.suppressalg)
        if flag != IDA_SUCCESS:
            raise IDAError(flag)
        
        if not settings.idv is None:
            idv = settings.idv
            flag = IDASetId(self.thisptr, <N_Vector>idv.thisptr)
            if flag != IDA_SUCCESS:
                raise IDAError(flag)
                
        if settings.tstop != float('+inf'):
            flag = IDASetStopTime(self.thisptr, settings.tstop)
            if flag != IDA_SUCCESS:
                raise IDAError(flag)
        
        flag = IDASetUserData(self.thisptr, <void*>&self.userdata)
        if flag != IDA_SUCCESS:
            raise IDAError(flag)
            
    cpdef GetConsistentIC(self, realtype direction):
        """
        Ensure the initial conditions are consistent.
        """
        cdef realtype tout = self.t0 + direction
        cdef int icopt, flag
        
        if self.settings.icopt == 'IDA_Y_INIT':
            icopt = IDA_Y_INIT
        elif self.settings.icopt == 'IDA_YA_YDP_INIT':
            icopt = IDA_YA_YDP_INIT
        else:
            msg = "Expected either 'IDA_Y_INIT' or 'IDA_YA_YDP_INIT' for icopt"
            raise IDAError(0, msg)
        
        flag = IDASetLineSearchOffIC(self.thisptr, self.settings.lsoff)
        if flag != IDA_SUCCESS:
            raise IDAError(flag)
            
        flag = IDACalcIC(self.thisptr, icopt, tout)
        if flag != IDA_SUCCESS:
            raise IDAError(flag)
        
        flag = IDAGetConsistentIC(self.thisptr, <N_Vector>self.y.thisptr, <N_Vector>self.ydot.thisptr)
        if flag != IDA_SUCCESS:
            raise IDAError(flag)
            
    cpdef __handleRoot(self, realtype event_time, N_Vector_Serial yroot, N_Vector_Serial ydotroot):
        """
        Create root expection with info on time
        and solution vector.
        """
        cdef int *event_info
        cdef int i, size = len(self.settings.SW)
        
        # Allocate memory for the event_info
        event_info = <int *>malloc(size * sizeof(int))
        IDAGetRootInfo(self.thisptr, event_info)
        
        SW = [bool(event_info[i]) for i in range(size)]
        
        raise IDARootException(event_time, yroot, ydotroot, SW)
    
    cpdef IDAIterator iter(self, realtype t0, realtype dt):
        """
        Run solver to next root or to 
        tstop if set
        """
        return IDAIterator(self, t0, dt)
        
    cpdef step(self, realtype tf):
        """
        Run solver to time tf, next root change
        or to tstop if set
        """
        cdef int flag
        cdef realtype tret
        
        flag = IDASolve(self.thisptr, tf, &tret, <N_Vector>self.y.thisptr, <N_Vector>self.ydot.thisptr, IDA_NORMAL)
        
        if flag < 0:
            raise IDAError(flag)
        
        if flag == IDA_ROOT_RETURN:
            self.__handleRoot(tret, self.y.copy(), self.ydot.copy())
        
        return self.y.copy(), self.ydot.copy()

################################### N_Vector_Serial##########################

cdef class N_Vector_Serial:
    """
    Wrapper arround sundials N_Vector_Serial
    """
    # buffer interface
    cdef __cythonbufferdefaults__ = {"ndim": 1, "mode": "c"}
    
    def __init__(N_Vector_Serial self, int size):
        assert size > 0
        self.thisptr = N_VNew_Serial(size)
        
    def __dealloc__(N_Vector_Serial self):
        if not self.thisptr == NULL:
            N_VDestroy_Serial(<N_Vector>self.thisptr)
    
    cdef int getsize(N_Vector_Serial self):
        cdef N_Vector v = <N_Vector>self.thisptr
        return (<N_VectorContent_Serial>v.content).length
    
    cdef double *getdata(N_Vector_Serial self):
        cdef N_Vector v = <N_Vector>self.thisptr
        return (<N_VectorContent_Serial>v.content).data
        
    cdef double getitem(N_Vector_Serial self, int index):
        cdef N_Vector v = <N_Vector>self.thisptr
        return (<N_VectorContent_Serial>v.content).data[index]
    
    cdef setitem(N_Vector_Serial self, int index, double value):
        cdef N_Vector v = <N_Vector>self.thisptr
        (<N_VectorContent_Serial>v.content).data[index] = value
        
    def __getbuffer__(N_Vector_Serial self, Py_buffer* buffer, int flags):
        self.__shape[0] = self.getsize()
        
        buffer.buf = <void *>self.getdata()
        buffer.obj = self
        buffer.len = self.getsize() * sizeof(double)
        buffer.readonly = 0
        buffer.format = <char*>"d"
        buffer.ndim = 1
        buffer.shape = <Py_ssize_t *>&self.__shape
        buffer.strides = NULL
        buffer.suboffsets = NULL
        buffer.itemsize = sizeof(double)
        buffer.internal = NULL
        
    def __releasebuffer__(N_Vector_Serial self, Py_buffer* buffer):
        pass
    
    property dtype:
        def __get__(N_Vector_Serial self):
            return "d"
    
    property size:
        def __get__(N_Vector_Serial self):
            return self.getsize()
    
    property itemsize:
        def __get__(N_Vector_Serial self):
            return sizeof(double)
    
    property nbytes:
        def __get__(N_Vector_Serial self):
            return sizeof(double) * self.getsize()
    
    def __repr__(N_Vector_Serial self):
        cdef double *data = self.getdata()
        cdef int size = self.getsize()
        cdef int index
        cdef char *typestr = "%g"
            
        ret = [None] * size
        for index in range(size):
            ret[index] = typestr % data[index]
         
        return "({})".format( (b",".join(ret)).decode('utf-8') )
    
    def __str__(N_Vector_Serial self):
        return '%s%s' % (self.__class__.__name__, self.__repr__())
        
    def __len__(N_Vector_Serial self):
        return self.getsize()
    
    def __getitem__(N_Vector_Serial self, int index):
        assert index < self.getsize()
        return self.getitem(index)
    
    def __setitem__(N_Vector_Serial self, int index, double value):
        assert index < self.getsize()
        self.setitem(index, value)
    
    cpdef copy(N_Vector_Serial self):
        cdef N_Vector_Serial ret = N_Vector_Serial.__new__(N_Vector_Serial, 0)
        cdef double *data = self.getdata()
        cdef int i, size = self.getsize()
        
        ret.thisptr = N_VNew_Serial(size)
        
        for i in range(size):
            ret.setitem(i, data[i])
        
        return ret
    
cdef double *N_Vector_data(void *vector):
    cdef N_Vector v = <N_Vector>vector
    return (<N_VectorContent_Serial>v.content).data

cdef int N_Vector_length(void *vector):
    cdef N_Vector v = <N_Vector>vector
    return (<N_VectorContent_Serial>v.content).length
    
def N_Vector_Serial_Copy_Array(object[double, ndim=1] arg):
    # fast copy construct bypassing __init__
    cdef N_Vector_Serial ret = N_Vector_Serial.__new__(N_Vector_Serial, 0)
    cdef int i, size
    
    size = len(arg)
    ret.thisptr = N_VNew_Serial(size)
    
    for i in range(size):
        ret.setitem(i, arg[i])
    
    return ret     