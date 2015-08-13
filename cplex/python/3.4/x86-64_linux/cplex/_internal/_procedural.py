# --------------------------------------------------------------------------
# File: _procedural.py 
# ---------------------------------------------------------------------------
# Licensed Materials - Property of IBM
# 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
# Copyright IBM Corporation 2008, 2014. All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with
# IBM Corp.
# ------------------------------------------------------------------------

from __future__ import print_function

from . import _pycplex as CR

from ._constants import *

from . import _list_array_utils as LAU
from ..exceptions import CplexSolverError, CplexError

from random import randrange
from .. import six
from ..six.moves import map
from ..six.moves import zip

default_encoding     = "UTF-8"
cpx_default_encoding = "UTF-8"

def cpx_decode(my_str, enc):
    if isinstance(my_str, six.binary_type):
        if enc == cpx_default_encoding:
            return my_str
        else:
            return six.text_type(my_str, enc).encode(cpx_default_encoding)
    else:
        assert isinstance(my_str, six.text_type)
        return my_str.encode(cpx_default_encoding)

def cpx_decode_noop3(my_str, enc):
    if six.PY2:
        return cpx_decode(my_str, enc)
    else:
        return my_str

def _cpx_encode_py2(my_str, enc):
    assert six.PY2
    if enc == cpx_default_encoding:
        return my_str
    else:
        return six.text_type(my_str, cpx_default_encoding).encode(enc)

def _cpx_encode_py3(my_str, enc):
    assert six.PY3
    if isinstance(my_str, six.binary_type):
        return my_str.decode(enc)
    else:
        assert isinstance(my_str, six.text_type)
        if enc == cpx_default_encoding:
            return my_str
        else:
            return my_str.encode(cpx_default_encoding).decode(enc)

def cpx_encode_noop3(my_str, enc):
    if six.PY2:
        return cpx_encode(my_str, enc)
    else:
        return my_str

def cpx_encode(my_str, enc):
    if six.PY2:
        return _cpx_encode_py2(my_str, enc)
    else:
        return _cpx_encode_py3(my_str, enc)

def cpx_transcode(name, enc):
    if isinstance(name, six.text_type):
        return name.encode("utf-8")
    else:
        return six.text_type(name, enc).encode("utf-8")

def _getstrbfr():
    return "?" + str(randrange(0, 1000000000, 1))

def _safeDoubleArray(arraylen):
    # Make sure that we never request a zero-length array.  This results in
    # a malloc(0) call in the SWIG layer.  On AIX this returns NULL which
    # causes problems.  By ensuring that the array is at least size 1, we
    # avoid these problems and the overhead should be negligable.
    if arraylen <= 0:
        arraylen = 1
    return CR.doubleArray(arraylen)

def _safeIntArray(arraylen):
    # See comment for _safeDoubleArray above.
    if arraylen <= 0:
        arraylen = 1
    return CR.intArray(arraylen)

def getstatstring(env, statind, enc=default_encoding):
    output = []
    CR.CPXXgetstatstring(env, statind, output)
    return cpx_encode(output[0], enc)

def geterrorstring(env, errcode):
    output = []
    CR.CPXXgeterrorstring(env, errcode, output)
    return output[0]

def cb_geterrorstring(env, status):
    return CR.cb_geterrorstring(env, status)

def unset_py_terminator():
    CR.unset_py_terminator()

def set_py_terminator():
    CR.set_py_terminator()

def sigint_swap():
    CR.sigint_swap()

def pack_env_lp_ptr(env, lp):
    return CR.pack_env_lp_ptr(env, lp)

def Pylolmat_to_CHBmat(lolmat, get_indices, r_c, enc):
    return CR.Pylolmat_to_CHBmat(lolmat, get_indices, r_c, cpx_transcode, enc)

def free_CHBmat(lolmat):
    CR.free_CHBmat(lolmat)

class StatusChecker:
    def __init__(self):
        class NoOp:
            pass
        self._pyenv = NoOp()
        self._pyenv._callback_exception = None
    def __call__(self, env, status, from_cb=0):
        error_string = None
        try:
            if self._pyenv._callback_exception is not None:
                callback_exception = self._pyenv._callback_exception
                self._pyenv._callback_exception = None
                unset_py_terminator()
                if isinstance(callback_exception, Exception) and \
                       len(callback_exception.args) == 2 and \
                       callback_exception.args[0] == "ERROR":
                    error_string = callback_exception.args[1]
                else:
                    try:
                        if isinstance(callback_exception, Exception):
                            raise callback_exception
                        if hasattr(callback_exception[1], "args"):
                            callback_exception = callback_exception[0](
                                *callback_exception[1].args)
                        else:
                            callback_exception = callback_exception[0](
                                callback_exception[1])
                        raise callback_exception
                    except KeyboardInterrupt:
                        print("KeyboardInterrupt")
                        return
        except ReferenceError:
            pass
        if status != 0:
            if error_string is None:
                if from_cb == 1:
                    error_string = cb_geterrorstring(env, status)
                else:
                    error_string = geterrorstring(env, status)
            raise CplexSolverError(error_string, env, status)

check_status = StatusChecker()

def set_status_checker():
    CR.set_status_checker(check_status)

# Environment

def version(env):
    return CR.CPXXversion(env)

def versionnumber(env):
    ver = CR.intPtr()
    status = CR.CPXXversionnumber(env, ver)
    return ver.value()

def openCPLEX():
    status = CR.intPtr()
    env = CR.CPXXopenCPLEX(status)
    check_status(env, status.value())
    return env

def closeCPLEX(env):
    envp = CR.CPXENVptrPtr()
    envp.assign(env)
    status = CR.CPXXcloseCPLEX(envp)
    check_status(env, status)
    return

def getchannels(env):
    results = CR.CPXCHANNELptrPtr()
    warning = CR.CPXCHANNELptrPtr()
    error   = CR.CPXCHANNELptrPtr()
    log     = CR.CPXCHANNELptrPtr()
    status = CR.CPXXgetchannels(env, results, warning, error, log)
    check_status(env, status)
    return (results.value(), warning.value(), error.value(), log.value())

def addfuncdest(env, channel, fileobj):
    status = CR.CPXXaddfuncdest(env, channel, fileobj)
    check_status(env, status)
    return

def delfuncdest(env, channel, fileobj):
    status = CR.CPXXdelfuncdest(env, channel, fileobj)
    check_status(env, status)
    return

def setlpcallbackfunc(env, cbhandle):
    status = CR.CPXXsetlpcallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setnetcallbackfunc(env, cbhandle):
    status = CR.CPXXsetnetcallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def settuningcallbackfunc(env, cbhandle):
    status = CR.CPXXsettuningcallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setheuristiccallbackfunc(env, cbhandle):
    status = CR.CPXXsetheuristiccallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setlazyconstraintcallbackfunc(env, cbhandle):
    status = CR.CPXXsetlazyconstraintcallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setusercutcallbackfunc(env, cbhandle):
    status = CR.CPXXsetusercutcallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setincumbentcallbackfunc(env, cbhandle):
    status = CR.CPXXsetincumbentcallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setnodecallbackfunc(env, cbhandle):
    status = CR.CPXXsetnodecallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setbranchcallbackfunc(env, cbhandle):
    status = CR.CPXXsetbranchcallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setbranchnosolncallbackfunc(env, cbhandle):
    status = CR.CPXXsetbranchnosolncallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setsolvecallbackfunc(env, cbhandle):
    status = CR.CPXXsetsolvecallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setinfocallbackfunc(env, cbhandle):
    status = CR.CPXXsetinfocallbackfunc(env, cbhandle)
    check_status(env, status)
    return

def setmipcallbackfunc(env, cbhandle):
    status = CR.CPXXsetmipcallbackfunc(env, cbhandle)
    check_status(env, status)
    return

# Parameters

def setintparam(env, whichparam, newvalue):
    status = CR.CPXXsetintparam(env, whichparam, newvalue)
    check_status(env, status)
    return

def setlongparam(env, whichparam, newvalue):
    status = CR.CPXXsetlongparam(env, whichparam, newvalue)
    check_status(env, status)
    return

def setdblparam(env, whichparam, newvalue):
    status = CR.CPXXsetdblparam(env, whichparam, newvalue)
    check_status(env, status)
    return

def setstrparam(env, whichparam, newvalue):
    status = CR.CPXXsetstrparam(env, whichparam, newvalue)
    check_status(env, status)
    return

def getintparam(env, whichparam):
    curval = CR.intPtr()
    status = CR.CPXXgetintparam(env, whichparam, curval)
    check_status(env, status)
    return curval.value()

def getlongparam(env, whichparam):
    curval = CR.cpxlongPtr()
    status = CR.CPXXgetlongparam(env, whichparam, curval)
    check_status(env, status)
    return curval.value()

def getdblparam(env, whichparam):
    curval = CR.doublePtr()
    status = CR.CPXXgetdblparam(env, whichparam, curval)
    check_status(env, status)
    return curval.value()

def getstrparam(env, whichparam):
    output = []
    status = CR.CPXXgetstrparam(env, whichparam, output)
    check_status(env, status)
    return output[0]

def infointparam(env, whichparam):
    default = CR.intPtr()
    minimum = CR.intPtr()
    maximum = CR.intPtr()
    status  = CR.CPXXinfointparam(env, whichparam, default, minimum, maximum)
    check_status(env, status)
    return (default.value(), minimum.value(), maximum.value())

def infolongparam(env, whichparam):
    default = CR.cpxlongPtr()
    minimum = CR.cpxlongPtr()
    maximum = CR.cpxlongPtr()
    status  = CR.CPXXinfolongparam(env, whichparam, default, minimum, maximum)
    check_status(env, status)
    return (default.value(), minimum.value(), maximum.value())

def infodblparam(env, whichparam):
    default = CR.doublePtr()
    minimum = CR.doublePtr()
    maximum = CR.doublePtr()
    status  = CR.CPXXinfodblparam(env, whichparam, default, minimum, maximum)
    check_status(env, status)
    return (default.value(), minimum.value(), maximum.value())

def infostrparam(env, whichparam):
    output = []
    status = CR.CPXXinfostrparam(env, whichparam, output)
    check_status(env, status)
    return output[0]

def getparamname(env, whichparam, enc=default_encoding):
    output = []
    status = CR.CPXXgetparamname(env, whichparam, output)
    check_status(env, status)
    return cpx_encode(output[0]. enc)
    
def getparamnum(env, param_name):
    output = CR.intPtr()
    status = CR.CPXXgetparamnum(env, param_name, output)
    check_status(env, status)
    return output.value()
    
def getparamtype(env, param_name):
    output = CR.intPtr()
    status = CR.CPXXgetparamtype(env, param_name, output)
    check_status(env, status)
    return output.value()

def readcopyparam(env, filename):
    status = CR.CPXXreadcopyparam(env, filename)
    check_status(env, status)
    return 

def writeparam(env, filename):
    status = CR.CPXXwriteparam(env, filename)
    check_status(env, status)
    return 
    
def getchgparam(env):
    count    = CR.intPtr()
    surplus  = CR.intPtr()
    paramnum = LAU.int_list_to_array([])
    pspace   = 0
    status   = CR.CPXXgetchgparam(env, count, paramnum, pspace, surplus)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if surplus.value() == 0:
        return []
    pspace   = -surplus.value()
    paramnum = _safeIntArray(pspace)
    status   = CR.CPXXgetchgparam(env, count, paramnum, pspace, surplus)
    check_status(env, status)
    return LAU.int_array_to_list(paramnum, pspace)
    
def tuneparam(env, lp, int_param_values, dbl_param_values, str_param_values):
    tuning_status = CR.intPtr()
    intcnt = len(int_param_values)
    dblcnt = len(dbl_param_values)
    strcnt = len(str_param_values)
    intnum = [int_param_values[i][0] for i in range(intcnt)]
    intval = [int_param_values[i][1] for i in range(intcnt)]
    dblnum = [dbl_param_values[i][0] for i in range(dblcnt)]
    dblval = [dbl_param_values[i][1] for i in range(dblcnt)]
    strnum = [str_param_values[i][0] for i in range(strcnt)]
    strval = [str_param_values[i][1] for i in range(strcnt)]
    sigint_swap()
    status = CR.CPXXtuneparam(env, lp,
                              intcnt,
                              LAU.int_list_to_array(intnum),
                              LAU.int_list_to_array_trunc_int32(intval),
                              dblcnt,
                              LAU.int_list_to_array(dblnum),
                              LAU.double_list_to_array(dblval),
                              strcnt,
                              LAU.int_list_to_array(strnum),
                              [cpx_decode(x, default_encoding) for x in strval],
                              tuning_status)
    sigint_swap()
    check_status(env, status)
    return tuning_status.value()

def tuneparamprobset(env, filenames, filetypes, int_param_values,
                     dbl_param_values, str_param_values):
    tuning_status = CR.intPtr()
    intcnt = len(int_param_values)
    dblcnt = len(dbl_param_values)
    strcnt = len(str_param_values)
    intnum = [int_param_values[i][0] for i in range(intcnt)]
    intval = [int_param_values[i][1] for i in range(intcnt)]
    dblnum = [dbl_param_values[i][0] for i in range(dblcnt)]
    dblval = [dbl_param_values[i][1] for i in range(dblcnt)]
    strnum = [str_param_values[i][0] for i in range(strcnt)]
    strval = [str_param_values[i][1] for i in range(strcnt)]
    sigint_swap()
    status = CR.CPXXtuneparamprobset(
        env, len(filenames),
        [cpx_decode(x, default_encoding) for x in filenames],
        [cpx_decode(x, default_encoding) for x in filetypes],
        intcnt, LAU.int_list_to_array(intnum),
        LAU.int_list_to_array_trunc_int32(intval),
        dblcnt, LAU.int_list_to_array(dblnum),
        LAU.double_list_to_array(dblval),
        strcnt, LAU.int_list_to_array(strnum),
        [cpx_decode(x, default_encoding) for x in strval],
        tuning_status)
    sigint_swap()
    check_status(env, status)
    return tuning_status.value()



# Cplex

def createprob(env, probname, enc=default_encoding):
    status = CR.intPtr()
    lp = CR.CPXXcreateprob(env, status, cpx_decode_noop3(probname, enc))
    check_status(env, status.value())
    return lp

def readcopyprob(env, lp, filename, filetype=""):
    if filetype == "":
        status = CR.CPXXreadcopyprob(env, lp, filename)
    else:
        status = CR.CPXXreadcopyprob(env, lp, filename, filetype)
    check_status(env, status)
    return

def cloneprob(env, lp):
    status = CR.intPtr()
    lp     = CR.CPXXcloneprob(env, lp, status)
    check_status(env, status.value())
    return lp

def freeprob(env, lp):
    lpp = CR.CPXLPptrPtr()
    lpp.assign(lp)
    status = CR.CPXXfreeprob(env, lpp)
    check_status(env, status)
    return
    
def mipopt(env, lp):
    sigint_swap()
    status = CR.CPXXmipopt(env, lp)
    sigint_swap()
    check_status(env, status)
    return

def distmipopt(env, lp):
    sigint_swap()
    status = CR.CPXXdistmipopt(env, lp)
    sigint_swap()
    check_status(env, status)
    return

def copyvmconfig(env, xmlstring):
    status = CR.CPXXcopyvmconfig(env, xmlstring)
    check_status(env, status)

def readcopyvmconfig(env, file):
    status = CR.CPXXreadcopyvmconfig(env, file)
    check_status(env, status)

def delvmconfig(env):
    status = CR.CPXXdelvmconfig(env)
    check_status(env, status)

def hasvmconfig(env):
    hasvmconfig_p = CR.intPtr()
    status = CR.CPXEhasvmconfig(env, hasvmconfig_p)
    check_status(env, status)
    return hasvmconfig_p.value() != 0

def qpopt(env, lp):
    sigint_swap()
    status = CR.CPXXqpopt(env, lp)
    sigint_swap()
    check_status(env, status)
    return

def baropt(env, lp):
    sigint_swap()
    status = CR.CPXXbaropt(env, lp)
    sigint_swap()
    check_status(env, status)
    return

def hybbaropt(env, lp, method):
    sigint_swap()
    status = CR.CPXXhybbaropt(env, lp, method)
    sigint_swap()
    check_status(env, status)
    return

def hybnetopt(env, lp, method):
    sigint_swap()
    status = CR.CPXXhybnetopt(env, lp, method)
    sigint_swap()
    check_status(env, status)
    return
    
def lpopt(env, lp):
    sigint_swap()
    status = CR.CPXXlpopt(env, lp)
    sigint_swap()
    check_status(env, status)
    return

def primopt(env, lp):
    status = CR.CPXXprimopt(env, lp)
    check_status(env, status)
    return

def dualopt(env, lp):
    status = CR.CPXXdualopt(env, lp)
    check_status(env, status)
    return

def siftopt(env, lp):
    status = CR.CPXXsiftopt(env, lp)
    check_status(env, status)
    return

def feasopt(env, lp, rhs, rng, lb, ub):
    sigint_swap()
    status = CR.CPXXfeasopt(env, lp, LAU.double_list_to_array(rhs),
                           LAU.double_list_to_array(rng),
                           LAU.double_list_to_array(lb),
                           LAU.double_list_to_array(ub))
    sigint_swap()
    check_status(env, status)
    return

def feasoptext(env, lp, grppref, grpbeg, grpind, grptype):
    grpcnt = len(grppref)
    concnt = len(grpind)
    sigint_swap()
    status = CR.CPXXfeasoptext(env, lp, grpcnt, concnt,
                              LAU.double_list_to_array(grppref), LAU.int_list_to_array(grpbeg),
                              LAU.int_list_to_array(grpind), grptype)
    sigint_swap()
    check_status(env, status)
    return

def delnames(env, lp):
    status = CR.CPXXdelnames(env, lp)
    check_status(env, status)
    return

def writeprob(env, lp, filename, filetype = ""):
    if filetype == "":
        status = CR.CPXXwriteprob(env, lp, filename)
    else:
        status = CR.CPXXwriteprob(env, lp, filename, filetype)
    check_status(env, status)
    return

def embwrite(env, lp, filename):
    status = CR.CPXXembwrite(env, lp, filename)
    check_status(env, status)
    return

def dperwrite(env, lp, filename, epsilon):
    status = CR.CPXXdperwrite(env, lp, filename, epsilon)
    check_status(env, status)
    return
    
def pperwrite(env, lp, filename, epsilon):
    status = CR.CPXXpperwrite(env, lp, filename, epsilon)
    check_status(env, status)
    return
    
def preslvwrite(env, lp, filename):
    objoff = CR.doublePtr()
    status = CR.CPXXpreslvwrite(env, lp, filename, objoff)
    check_status(env, status)
    return objoff.value()
    
def dualwrite(env, lp, filename):
    objshift = CR.doublePtr()
    status   = CR.CPXXdualwrite(env, lp, filename, objshift)
    check_status(env, status)
    return objshift.value()

def chgprobtype(env, lp, type):
    status = CR.CPXXchgprobtype(env, lp, type)
    check_status(env, status)
    return

def chgprobtypesolnpool(env, lp, type, soln):
    status = CR.CPXXchgprobtypesolnpool(env, lp, type, soln)
    check_status(env, status)
    return
    
def getprobtype(env, lp):
    return CR.CPXXgetprobtype(env, lp)

def chgprobname(env, lp, probname, enc=default_encoding):
    status = CR.CPXXchgprobname(env, lp, cpx_decode_noop3(probname, enc))
    check_status(env, status)

def getprobname(env, lp, enc=default_encoding):
    inoutlist = [0]
    status = CR.CPXXgetprobname(env, lp, inoutlist)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inoutlist == [0]:
        return cpx_encode_noop3("", enc)
    status = CR.CPXXgetprobname(env, lp, inoutlist)
    check_status(env, status)
    return cpx_encode_noop3(inoutlist[0], enc)

def getnumcols(env, lp):
    return CR.CPXXgetnumcols(env, lp)

def getnumint(env, lp):
    return CR.CPXXgetnumint(env, lp)

def getnumbin(env, lp):
    return CR.CPXXgetnumbin(env, lp)

def getnumsemicont(env, lp):
    return CR.CPXXgetnumsemicont(env, lp)

def getnumsemiint(env, lp):
    return CR.CPXXgetnumsemiint(env, lp)

def getnumrows(env, lp):
    return CR.CPXXgetnumrows(env, lp)

def populate(env, lp):
    sigint_swap()
    status = CR.CPXXpopulate(env, lp)
    sigint_swap()
    check_status(env, status)
    return

def _getnumusercuts(env, lp):
    return CR.CPXEgetnumusercuts(env, lp)

def _getnumlazyconstraints(env, lp):
    return CR.CPXEgetnumlazyconstraints(env, lp)

def _hasgeneralconstraints(env, lp):
    for which in range(CPX_CON_SOS + 1, CPX_CON_LAST_CONTYPE):
        if CR.CPXEgetnumgconstrs(env, lp, which) > 0:
            return True
    return False

def getnummipstarts(env, lp):
    return CR.CPXXgetnummipstarts(env, lp)

def getnumqconstrs(env, lp):
    return CR.CPXXgetnumqconstrs(env, lp)

def getnumindconstrs(env, lp):
    return CR.CPXXgetnumindconstrs(env, lp)

def getnumsos(env, lp):
    return CR.CPXXgetnumsos(env, lp)

def cleanup(env, lp, eps):
    status = CR.CPXXcleanup(env, lp, eps)
    check_status(env, status)
    return

def basicpresolve(env, lp):
    numcols = CR.CPXXgetnumcols(env, lp)
    numrows = CR.CPXXgetnumcols(env, lp)
    redlb   = _safeDoubleArray(numcols)
    redub   = _safeDoubleArray(numcols)
    rstat   = _safeIntArray(numrows)
    status  = CR.CPXXbasicpresolve(env, lp, redlb, redub, rstat)
    check_status(env, status)
    return (LAU.double_array_to_list(redlb, numcols),
            LAU.double_array_to_list(redub, numcols),
            LAU.int_array_to_list(rstat, numrows))

def pivotin(env, lp, rlist):
    status = CR.CPXXpivotin(env, lp, LAU.int_list_to_array(rlist), len(rlist))
    check_status(env, status)
    return

def pivotout(env, lp, clist):
    status = CR.CPXXpivotout(env, lp, LAU.int_list_to_array(clist), len(clist))
    check_status(env, status)
    return
    
def pivot(env, lp, jenter, jleave, leavestat):
    status = CR.CPXXpivot(env, lp, jenter, jleave, leavestat)
    check_status(env, status)
    return

def strongbranch(env, lp, goodlist, itlim):
    goodlen = len(goodlist)
    downpen = _safeDoubleArray(goodlen)
    uppen   = _safeDoubleArray(goodlen)
    sigint_swap()
    status = CR.CPXXstrongbranch(env, lp, goodlen, goodlist, downpen, uppen, itlim)
    sigint_swap()
    check_status(env, status)
    return (LAU.double_array_to_list(downpen, goodlen),
            LAU.double_array_to_list(uppen, goodlen))

def completelp(env, lp):
    status = CR.CPXXcompletelp(env, lp)
    check_status(env, status)
    return    


# Variables

def newcols(env, lp, obj, lb, ub, xctype, colname, enc=default_encoding):
    ccnt = max(len(obj), len(lb), len(ub), len(xctype), len(colname))
    obj = LAU.double_C_array(obj)
    lb  = LAU.double_C_array(lb)
    ub  = LAU.double_C_array(ub)
    status = CR.CPXXnewcols(env, lp, ccnt, obj.array,
                            lb.array, ub.array,
                            cpx_decode(xctype, enc),
                            [cpx_decode(x, enc) for x in colname])
    check_status(env, status)
    return

def addcols(env, lp, ccnt, nzcnt, obj, cmat, lb, ub, colname, enc=default_encoding):
    obj = LAU.double_C_array(obj)
    lb  = LAU.double_C_array(lb)
    ub  = LAU.double_C_array(ub)
    status = CR.CPXXaddcols(env, lp, ccnt, nzcnt, obj.array, cmat._mat,
                            lb.array, ub.array,
                            [cpx_decode(x, enc) for x in colname])
    check_status(env, status)
    return

def delcols(env, lp, begin, end):
    status = CR.CPXXdelcols(env, lp, begin, end)
    check_status(env, status)
    return [0]

def delsetcols(env, lp, delstat):
    delstat_array = LAU.int_list_to_array(delstat)
    status = CR.CPXXdelsetcols(env, lp, delstat_array)
    check_status(env, status)
    return LAU.int_array_to_list(delstat_array, len(delstat))

def chgbds(env, lp, indices, lu, bd):
    status = CR.CPXXchgbds(env, lp, len(indices), LAU.int_list_to_array(indices),
                          lu, LAU.double_list_to_array(bd))
    check_status(env, status)
    return

def chgcolname(env, lp, indices, newnames, enc=default_encoding):
    status = CR.CPXXchgcolname(env, lp, len(indices),
                               LAU.int_list_to_array(indices),
                               [cpx_decode(x, enc) for x in newnames])
    check_status(env, status)
    return

def chgctype(env, lp, indices, xctype):
    status = CR.CPXXchgctype(env, lp, len(indices),
                             LAU.int_list_to_array(indices),
                             cpx_decode(xctype, default_encoding))
    check_status(env, status)
    return

def getcolindex(env, lp, colname, enc=default_encoding):
    index = CR.intPtr()
    status = CR.CPXXgetcolindex(env, lp, cpx_decode_noop3(colname, enc), index)
    check_status(env, status)
    return index.value()

def getcolname(env, lp, begin, end, enc=default_encoding):
    inout_list = [0, begin, end]
    status = CR.CPXXgetcolname(env, lp, inout_list)
    if status == CR.CPXERR_NO_NAMES:
        return []
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list == [0]:
        return [cpx_encode_noop3(x, enc) for x in [""] * (end - begin + 1)]
    inout_list.extend([begin, end])
    status = CR.CPXXgetcolname(env, lp, inout_list)
    check_status(env, status)
    return [cpx_encode_noop3(x, enc) for x in inout_list]

def getctype(env, lp, begin, end):
    numcols = end - begin + 1
    ctype = ""
    for i in range(numcols):
        ctype = "".join([ctype, _getstrbfr()])
    status = CR.CPXXgetctype(env, lp, ctype, begin, end)
    if status == CR.CPXERR_NOT_MIP:
        return ""
    check_status(env, status)
    return ctype[:numcols]

def getlb(env, lp, begin, end):
    lblen = end - begin + 1
    lb    = _safeDoubleArray(lblen)
    status = CR.CPXXgetlb(env, lp, lb, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(lb, lblen)

def getub(env, lp, begin, end):
    ublen = end - begin + 1
    ub    = _safeDoubleArray(ublen)
    status = CR.CPXXgetub(env, lp, ub, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(ub, ublen)

def getcols(env, lp, begin, end):
    inout_list = [0, begin, end]
    status = CR.CPXXgetcols(env, lp, inout_list)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list == [0]:
        return ([0] * (end - begin + 1), [], [])
    inout_list.extend([begin, end])
    status = CR.CPXXgetcols(env, lp, inout_list)
    check_status(env, status)
    return tuple(inout_list)
    
def copyprotected(env, lp, indices):
    status = CR.CPXXcopyprotected(env, lp, len(indices), LAU.int_list_to_array(indices))
    check_status(env, status)
    return    

def getprotected(env, lp):
    count   = CR.intPtr()
    surplus = CR.intPtr()
    indices = LAU.int_list_to_array([])
    pspace  = 0
    status  = CR.CPXXgetprotected(env, lp, count, indices, pspace, surplus)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if surplus.value() == 0:
        return []
    pspace   = -surplus.value()
    indices = _safeIntArray(pspace)
    status  = CR.CPXXgetprotected(env, lp, count, indices, pspace, surplus)
    check_status(env, status)
    return LAU.int_array_to_list(indices, pspace)

def tightenbds(env, lp, indices, lu, bd):
    status = CR.CPXXtightenbds(env, lp, len(indices), LAU.int_list_to_array(indices),
                              lu, LAU.double_list_to_array(bd))
    check_status(env, status)
    return

# Linear Constraints

def newrows(env, lp, rhs, sense, rngval, rowname, enc=default_encoding):
    rcnt = max(len(rhs), len(sense), len(rngval), len(rowname))
    rhs = LAU.double_C_array(rhs)
    rng = LAU.double_C_array(rngval)
    status = CR.CPXXnewrows(env, lp, rcnt, rhs.array,
                            cpx_decode(sense, enc),
                            rng.array,
                            [cpx_decode(x, enc) for x in rowname])
    check_status(env, status)
    return

def addrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, rmat, colname, rowname,
            enc=default_encoding):
    rhs = LAU.double_C_array(rhs)
    status = CR.CPXXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs.array, 
                            cpx_decode(sense, enc), rmat._mat, colname,
                            [cpx_decode(x, enc) for x in rowname])
    check_status(env, status)
    return

def delrows(env, lp, begin, end):
    status = CR.CPXXdelrows(env, lp, begin, end)
    check_status(env, status)
    return [0]

def delsetrows(env, lp, delstat):
    delstat_array = LAU.int_list_to_array(delstat)
    status = CR.CPXXdelsetrows(env, lp, delstat_array)
    check_status(env, status)
    return LAU.int_array_to_list(delstat_array, len(delstat))

def chgrowname(env, lp, indices, newnames, enc=default_encoding):
    status = CR.CPXXchgrowname(env, lp, len(indices),
                               LAU.int_list_to_array(indices),
                               [cpx_decode(x, enc) for x in newnames])
    check_status(env, status)
    return

def chgcoef(env, lp, i, j, newvalue):
    status = CR.CPXXchgcoef(env, lp, i, j, newvalue)
    check_status(env, status)
    return

def chgcoeflist(env, lp, rowlist, collist, vallist):
    status = CR.CPXXchgcoeflist(env, lp, len(rowlist),
                               LAU.int_list_to_array(rowlist),
                               LAU.int_list_to_array(collist),
                               LAU.double_list_to_array(vallist))
    check_status(env, status)
    return

def chgrhs(env, lp, indices, values):
    status = CR.CPXXchgrhs(env, lp, len(indices),
                           LAU.int_list_to_array(indices),
                           LAU.double_list_to_array(values))
    check_status(env, status)
    return

def chgrngval(env, lp, indices, values):
    ind = LAU.int_C_array(indices)
    val = LAU.double_C_array(values)
    status = CR.CPXXchgrngval(env, lp, len(indices), ind.array, val.array)
    check_status(env, status)
    return

def chgsense(env, lp, indices, senses):
    status = CR.CPXXchgsense(env, lp, len(indices),
                             LAU.int_list_to_array(indices),
                             cpx_decode(senses, default_encoding))
    check_status(env, status)
    return

def getrhs(env, lp, begin, end):
    rhslen = end - begin + 1
    rhs = _safeDoubleArray(rhslen)
    status = CR.CPXXgetrhs(env, lp, rhs, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(rhs, rhslen)

def getsense(env, lp, begin, end):
    inout_list = [begin, end]
    status = CR.CPXXgetsense(env, lp, inout_list)
    check_status(env, status)
    return cpx_encode(inout_list[0][:end - begin + 1], default_encoding)

def getrngval(env, lp, begin, end):
    rngvallen = end - begin + 1
    rngval    = _safeDoubleArray(rngvallen)
    status = CR.CPXXgetrngval(env, lp, rngval, begin, end)
    if status == CR.CPXERR_NO_RNGVAL:
        return []
    check_status(env, status)
    return LAU.double_array_to_list(rngval, rngvallen)

def getrowname(env, lp, begin, end, enc=default_encoding):
    inout_list = [0, begin, end]
    status = CR.CPXXgetrowname(env, lp, inout_list)
    if status == CR.CPXERR_NO_NAMES:
        return []
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list == [0]:
        return [cpx_encode_noop3(x, enc) for x in [""] * (end - begin + 1)]
    inout_list.extend([begin, end])
    status = CR.CPXXgetrowname(env, lp, inout_list)
    check_status(env, status)
    return [cpx_encode_noop3(x, enc) for x in inout_list]

def getcoef(env, lp, i, j):
    coef = CR.doublePtr()
    status = CR.CPXXgetcoef(env, lp, i, j, coef)
    check_status(env, status)
    return coef.value()

def getrowindex(env, lp, rowname, enc=default_encoding):
    index = CR.intPtr()
    status = CR.CPXXgetrowindex(env, lp, cpx_decode_noop3(rowname, enc), index)
    check_status(env, status)
    return index.value()

def getrows(env, lp, begin, end):
    inout_list = [0, begin, end]
    status = CR.CPXXgetrows(env, lp, inout_list)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list == [0]:
        return ([0] * (end - begin + 1), [], [])
    inout_list.extend([begin, end])
    status = CR.CPXXgetrows(env, lp, inout_list)
    check_status(env, status)
    return tuple(inout_list)

def getnumnz(env, lp):
    return CR.CPXXgetnumnz(env, lp)

def addlazyconstraints(env, lp, rhs, sense, rmatbeg, rmatind, rmatval, names,
                       enc=default_encoding):
    status = CR.CPXXaddlazyconstraints(env, lp, len(rmatbeg), len(rmatind),
                                       LAU.double_list_to_array(rhs),
                                       cpx_decode(sense, enc),
                                       LAU.int_list_to_array(rmatbeg),
                                       LAU.int_list_to_array(rmatind),
                                       LAU.double_list_to_array(rmatval),
                                       [cpx_decode(x, enc) for x in names])
    check_status(env, status)
    return

def addusercuts(env, lp, rhs, sense, rmatbeg, rmatind, rmatval, names,
                enc=default_encoding):
    status = CR.CPXXaddusercuts(env, lp, len(rmatbeg), len(rmatind),
                                LAU.double_list_to_array(rhs),
                                cpx_decode(sense, enc),
                                LAU.int_list_to_array(rmatbeg),
                                LAU.int_list_to_array(rmatind),
                                LAU.double_list_to_array(rmatval),
                                [cpx_decode(x, enc) for x in names])
    check_status(env, status)
    return

def freelazyconstraints(env, lp):
    status = CR.CPXXfreelazyconstraints(env, lp)
    check_status(env, status)
    return

def freeusercuts(env, lp):
    status = CR.CPXXfreeusercuts(env, lp)
    check_status(env, status)
    return

# SOS

def copysos(env, lp, sostype, sosbeg, sosind, soswt, sosname, enc=default_encoding):
    status = CR.CPXXcopysos(env, lp, len(sosbeg), len(sosind), sostype,
                            LAU.int_list_to_array(sosbeg),
                            LAU.int_list_to_array(sosind),
                            LAU.double_list_to_array(soswt),
                            [cpx_decode(x, enc) for x in sosname])
    check_status(env, status)
    return

def addsos(env, lp, sostype, sosbeg, sosind, soswt, sosnames, enc=default_encoding):
    status = CR.CPXXaddsos(env, lp, len(sosbeg), len(sosind), sostype,
                           LAU.int_list_to_array(sosbeg),
                           LAU.int_list_to_array(sosind),
                           LAU.double_list_to_array(soswt),
                           [cpx_decode(x, enc) for x in sosnames])
    check_status(env, status)
    return

def delsetsos(env, lp, delset):
    status = CR.CPXXdelsetsos(env, lp, LAU.int_list_to_array(delset))
    check_status(env, status)
    return

def getsos_info(env, lp, begin, end):
    numsos  = (end - begin + 1)
    numnz   = CR.intPtr()
    sostype = ""
    for i in range(numsos):
        sostype = "".join([sostype, _getstrbfr()])
    sosbeg  = _safeIntArray(numsos)
    space   = 0
    sosind  = LAU.int_list_to_array([])
    soswt   = LAU.double_list_to_array([])
    surplus = CR.intPtr()
    status = CR.CPXXgetsos(env, lp, numnz, sostype, sosbeg, sosind, soswt,
                          space, surplus, begin, end)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    return (sostype[0:numsos], -surplus.value())
    
def getsos(env, lp, begin, end):
    numsos  = (end - begin + 1)
    numnz   = CR.intPtr()
    sostype = ""
    for i in range(numsos):
        sostype = "".join([sostype, _getstrbfr()])
    sosbeg  = _safeIntArray(numsos)
    space   = 0
    sosind  = LAU.int_list_to_array([])
    soswt   = LAU.double_list_to_array([])
    surplus = CR.intPtr()
    status = CR.CPXXgetsos(env, lp, numnz, sostype, sosbeg, sosind, soswt,
                          space, surplus, begin, end)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if surplus.value() == 0:
        return ("1" * numsos, [0] * numsos, [], [])
    space = -surplus.value()
    sosind  = _safeIntArray(space)
    soswt   = _safeDoubleArray(space)
    status = CR.CPXXgetsos(env, lp, numnz, sostype, sosbeg, sosind, soswt,
                          space, surplus, begin, end)
    check_status(env, status)
    return (sostype[0:numsos],
            LAU.int_array_to_list(sosbeg, numsos),
            LAU.int_array_to_list(sosind, space),
            LAU.double_array_to_list(soswt, space))
    
def getsosindex(env, lp, name, enc=default_encoding):
    index  = CR.intPtr()
    status = CR.CPXXgetsosindex(env, lp, cpx_decode_noop3(name, enc), index)
    check_status(env, status)
    return index.value()

def getsosname(env, lp, begin, end, enc=default_encoding):
    if end < begin:
        return []
    inout_list = [0, begin, end]
    status = CR.CPXXgetsosname(env, lp, inout_list)
    if status == CR.CPXERR_NO_NAMES:
        return []
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list == [0]:
        return [cpx_encode_noop3(x, enc) for x in [""] * (end - begin + 1)]
    inout_list.extend([begin, end])
    status = CR.CPXXgetsosname(env, lp, inout_list)
    check_status(env, status)
    return [cpx_encode_noop3(x, enc)for x in inout_list]

# Indicator Constraints

def addindconstr(env, lp, indvar, complemented, rhs, sense, linind, linval,
                 name, enc=default_encoding):
    status = CR.CPXXaddindconstr(env, lp, indvar, complemented, len(linind),
                                 rhs,
                                 cpx_decode(sense, enc),
                                 LAU.int_list_to_array(linind),
                                 LAU.double_list_to_array(linval),
                                 cpx_decode_noop3(name, enc))
    check_status(env, status)
    return

def getindconstr(env, lp, which):
    indvar = CR.intPtr()
    complemented = CR.intPtr()
    nzcnt = CR.intPtr()
    rhs = CR.doublePtr()
    sensebfr = cpx_decode(_getstrbfr(), default_encoding)
    space = 0
    linind = LAU.int_list_to_array([])
    linval = LAU.double_list_to_array([])
    surplus = CR.intPtr()
    status = CR.CPXXgetindconstr(env, lp, indvar, complemented, nzcnt, rhs,
                                 sensebfr,
                                 linind, linval, space, surplus, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if surplus.value() == 0:
        sensebfr = cpx_encode(sensebfr, default_encoding)
        return (indvar.value(), complemented.value(), rhs.value(),
                sensebfr[0], [], [])
    space = -surplus.value()
    linind = _safeIntArray(space)
    linval = _safeDoubleArray(space)
    status = CR.CPXXgetindconstr(env, lp, indvar, complemented, nzcnt, rhs,
                                 sensebfr,
                                 linind, linval, space, surplus, which)
    check_status(env, status)
    sensebfr = cpx_encode(sensebfr, default_encoding)
    return (indvar.value(), complemented.value(), rhs.value(), sensebfr[0],
            LAU.int_array_to_list(linind, space),
            LAU.double_array_to_list(linval, space))

def getindconstr_constant(env, lp, which):
    indvar = CR.intPtr()
    complemented = CR.intPtr()
    nzcnt = CR.intPtr()
    rhs = CR.doublePtr()
    sensebfr = cpx_decode(_getstrbfr(), default_encoding)
    space = 0
    linind = LAU.int_list_to_array([])
    linval = LAU.double_list_to_array([])
    surplus = CR.intPtr()
    status = CR.CPXXgetindconstr(env, lp, indvar, complemented, nzcnt, rhs,
                                 sensebfr,
                                 linind, linval, space, surplus, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    sensebfr = cpx_encode(sensebfr, default_encoding)
    return (indvar.value(), complemented.value(), - surplus.value(),
            rhs.value(), sensebfr[0])

def getindconstrindex(env, lp, name, enc=default_encoding):
    index = CR.intPtr()
    status = CR.CPXXgetindconstrindex(env, lp, cpx_decode_noop3(name, enc), index)
    check_status(env, status)
    return index.value()

def delindconstrs(env, lp, begin, end):
    status = CR.CPXXdelindconstrs(env, lp, begin, end)
    check_status(env, status)
    return
    
def getindconstrslack(env, lp, begin, end):
    slacklen = end - begin + 1
    slacks   = _safeDoubleArray(slacklen)
    status   = CR.CPXXgetindconstrslack(env, lp, slacks, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(slacks, slacklen)
    
def getindconstrname(env, lp, which, enc=default_encoding):
    inoutlist = [0]
    status = CR.CPXXgetindconstrname(env, lp, inoutlist, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inoutlist == [0]:
        return ""
    status = CR.CPXXgetindconstrname(env, lp, inoutlist, which)
    check_status(env, status)
    return cpx_encode_noop3(inoutlist[0], enc)

# Quadratic Constraints

def addqconstr(env, lp, rhs, sense, linind, linval, quadrow, quadcol, quadval,
               name, enc=default_encoding):
    status = CR.CPXXaddqconstr(env, lp, len(linind), len(quadrow), rhs,
                               cpx_decode(sense, enc),
                               LAU.int_list_to_array(linind),
                               LAU.double_list_to_array(linval),
                               LAU.int_list_to_array(quadrow),
                               LAU.int_list_to_array(quadcol),
                               LAU.double_list_to_array(quadval),
                               cpx_decode_noop3(name, enc))
    check_status(env, status)
    return

def getqconstr(env, lp, which):
    linnzcnt = CR.intPtr()
    quadnzcnt = CR.intPtr()
    rhs = CR.doublePtr()
    sense = _getstrbfr()
    inout_list = [0, 0]
    status = CR.CPXXgetqconstr(env, lp, linnzcnt, quadnzcnt, rhs,
                               cpx_decode(sense, default_encoding),
                               inout_list, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list == [0, 0]:
        return ()
    status = CR.CPXXgetqconstr(env, lp, linnzcnt, quadnzcnt, rhs,
                               cpx_decode(sense, default_encoding),
                               inout_list, which)
    check_status(env, status)
    out_list = [rhs.value(), sense[0]]
    out_list.extend(inout_list)
    return tuple(out_list)

def getqconstr_rhs_sense(env, lp, which):
    linnzcnt = CR.intPtr()
    quadnzcnt = CR.intPtr()
    rhs = CR.doublePtr()
    sense = _getstrbfr()
    status = CR.CPXXgetqconstr(env, lp, linnzcnt, quadnzcnt, rhs,
                               cpx_decode(sense, default_encoding),
                               [0, 0], which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    return (rhs.value(), sense[0])

def getqconstr_linnz(env, lp, which):
    linnzcnt = CR.intPtr()
    quadnzcnt = CR.intPtr()
    rhs = CR.doublePtr()
    sense = _getstrbfr()
    inout_list = [0, 0]
    status = CR.CPXXgetqconstr(env, lp, linnzcnt, quadnzcnt, rhs,
                               cpx_decode(sense, default_encoding),
                               inout_list, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    return inout_list[0]

def getqconstr_lin(env, lp, which):
    linnzcnt = CR.intPtr()
    quadnzcnt = CR.intPtr()
    rhs = CR.doublePtr()
    sense = _getstrbfr()
    inout_list = [0, 0]
    status = CR.CPXXgetqconstr(env, lp, linnzcnt, quadnzcnt, rhs,
                               cpx_decode(sense, default_encoding),
                               inout_list, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list[0] == 0:
        return ([], [])
    inout_list[1] = 0
    status = CR.CPXXgetqconstr(env, lp, linnzcnt, quadnzcnt, rhs,
                               cpx_decode(sense, default_encoding),
                               inout_list, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    return (inout_list[0], inout_list[1])

def getqconstr_quadnz(env, lp, which):
    linnzcnt = CR.intPtr()
    quadnzcnt = CR.intPtr()
    rhs = CR.doublePtr()
    sense = _getstrbfr()
    inout_list = [0, 0]
    status = CR.CPXXgetqconstr(env, lp, linnzcnt, quadnzcnt, rhs,
                               cpx_decode(sense, default_encoding),
                               inout_list, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    return inout_list[1]

def getqconstr_quad(env, lp, which):
    linnzcnt = CR.intPtr()
    quadnzcnt = CR.intPtr()
    rhs = CR.doublePtr()
    sense = _getstrbfr()
    inout_list = [0, 0]
    status = CR.CPXXgetqconstr(env, lp, linnzcnt, quadnzcnt, rhs,
                               cpx_decode(sense, default_encoding),
                               inout_list, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list[1] == 0:
        return ([], [])
    inout_list[0] = 0
    status = CR.CPXXgetqconstr(env, lp, linnzcnt, quadnzcnt, rhs,
                               cpx_decode(sense, default_encoding),
                               inout_list, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    return (inout_list[2], inout_list[3], inout_list[4])

def delqconstrs(env, lp, begin, end):
    status = CR.CPXXdelqconstrs(env, lp, begin, end)
    check_status(env, status)
    return [0]

def getqconstrindex(env, lp, name, enc=default_encoding):
    index = CR.intPtr()
    status = CR.CPXXgetqconstrindex(env, lp, cpx_decode_noop3(name, enc), index)
    check_status(env, status)
    return index.value()

def getqconstrslack(env, lp, begin, end):
    slacklen = end - begin + 1
    slacks   = _safeDoubleArray(slacklen)
    status   = CR.CPXXgetqconstrslack(env, lp, slacks, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(slacks, slacklen)

def getqconstrname(env, lp, which, enc=default_encoding):
    inoutlist = [0]
    status = CR.CPXXgetqconstrname(env, lp, inoutlist, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inoutlist == [0]:
        return ""
    status = CR.CPXXgetqconstrname(env, lp, inoutlist, which)
    check_status(env, status)
    return cpx_encode_noop3(inoutlist[0], enc)

    
# Objective

def copyobjname(env, lp, objname, enc=default_encoding):
    status = CR.CPXXcopyobjname(env, lp, objname)
    check_status(env, status)
    return

def chgobj(env, lp, indices, values):
    status = CR.CPXXchgobj(env, lp, len(indices), LAU.int_list_to_array(indices), LAU.double_list_to_array(values))
    check_status(env, status)
    return

def chgobjsen(env, lp, maxormin):
    status = CR.CPXXchgobjsen(env, lp, maxormin)
    check_status(env, status)
    return

def getobjsen(env, lp):
    return CR.CPXXgetobjsen(env, lp)

def getobj(env, lp, begin, end):
    objlen = end - begin + 1
    obj    = _safeDoubleArray(objlen)
    status = CR.CPXXgetobj(env, lp, obj, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(obj, objlen)

def getobjname(env, lp, enc=default_encoding):
    inoutlist = [0]
    status = CR.CPXXgetobjname(env, lp, inoutlist)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inoutlist == [0]:
        return ""
    status = CR.CPXXgetobjname(env, lp, inoutlist)
    check_status(env, status)
    return inoutlist[0]

def copyquad(env, lp, qmatbeg, qmatind, qmatval):
    qmatcnt = [qmatbeg[i+1] - qmatbeg[i] for i in range(len(qmatbeg) - 1)]
    qmatcnt.append(len(qmatind) - qmatbeg[-1])
    status = CR.CPXXcopyquad(env, lp, LAU.int_list_to_array(qmatbeg), LAU.int_list_to_array(qmatcnt),
                            LAU.int_list_to_array(qmatind), LAU.double_list_to_array(qmatval))
    check_status(env, status)
    return

def copyqpsep(env, lp, qsepvec):
    status = CR.CPXXcopyqpsep(env, lp, LAU.double_list_to_array(qsepvec))
    check_status(env, status)
    return    

def chgqpcoef(env, lp, row, col, value):
    status = CR.CPXXchgqpcoef(env, lp, row, col, value)
    check_status(env, status)
    return    
    
def getquad(env, lp, begin, end):
    nzcnt   = CR.intPtr()
    ncols   = end - begin + 1
    qmatbeg = _safeIntArray(ncols)
    qmatind = LAU.int_list_to_array([])
    qmatval = LAU.double_list_to_array([])
    space   = 0
    surplus = CR.intPtr()
    status = CR.CPXXgetquad(env, lp, nzcnt, qmatbeg, qmatind, qmatval, space, surplus, begin, end)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if surplus.value() == 0:
        return []
    space = -surplus.value()
    qmatind = _safeIntArray(space)
    qmatval = _safeDoubleArray(space)
    status = CR.CPXXgetquad(env, lp, nzcnt, qmatbeg, qmatind, qmatval,
                            space, surplus, begin, end)
    check_status(env, status)
    return (LAU.int_array_to_list(qmatbeg, ncols),
            LAU.int_array_to_list(qmatind, space),
            LAU.double_array_to_list(qmatval, space))

def getqpcoef(env, lp, row, col):
    val = CR.doublePtr()
    status = CR.CPXXgetqpcoef(env, lp, row, col, val)
    check_status(env, status)
    return val.value()

def getnumquad(env, lp):
    return CR.CPXXgetnumquad(env, lp)

def getnumqpnz(env, lp):
    return CR.CPXXgetnumqpnz(env, lp)



# Optimizing Problems

# Accessing LP results

def solution(env, lp):
    numcols  = CR.CPXXgetnumcols(env, lp)
    numrows  = CR.CPXXgetnumrows(env, lp)
    xlen     = numcols
    pilen    = numrows
    slacklen = numrows
    djlen    = numcols
    lpstat   = CR.intPtr()
    objval   = CR.doublePtr()
    x        = _safeDoubleArray(xlen)
    pi       = _safeDoubleArray(pilen)
    slack    = _safeDoubleArray(slacklen)
    dj       = _safeDoubleArray(djlen)
    # first try to get all the values
    status  = CR.CPXXsolution(env, lp, lpstat, objval, x, pi, slack, dj)
    if status == CR.CPXERR_NOT_FOR_MIP or status == CR.CPXERR_NO_DUAL_SOLN:
        # if pi and dj are not available, NULL them
        pi    = CR.cvar.CPX_NULL
        dj    = CR.cvar.CPX_NULL
        pilen = 0
        djlen = 0
        status = CR.CPXXsolution(env, lp, lpstat, objval, x, pi, slack, dj)
    if status == CR.CPXERR_NO_INT_X:
        # if x and slacks are not available, NULL them
        x        = CR.cvar.CPX_NULL
        slack    = CR.cvar.CPX_NULL
        xlen     = 0
        slacklen = 0
        status  = CR.CPXXsolution(env, lp, lpstat, objval, x, pi, slack, dj)
    check_status(env, status)        
    return (lpstat.value(), objval.value(),
            LAU.double_array_to_list(x, xlen),
            LAU.double_array_to_list(pi, pilen),
            LAU.double_array_to_list(slack, slacklen),
            LAU.double_array_to_list(dj, djlen))

def solninfo(env, lp):
    lpstat = CR.intPtr()
    stype  = CR.intPtr()
    pfeas  = CR.intPtr()
    dfeas  = CR.intPtr()
    status = CR.CPXXsolninfo(env, lp, lpstat, stype, pfeas, dfeas)
    check_status(env, status)        
    return (lpstat.value(), stype.value(), pfeas.value(), dfeas.value())

def getstat(env, lp):
    return CR.CPXXgetstat(env, lp)

def getmethod(env, lp):
    return CR.CPXXgetmethod(env, lp)

def getobjval(env, lp):
    objval = CR.doublePtr()
    status = CR.CPXXgetobjval(env, lp, objval)
    check_status(env, status)        
    return objval.value()

def getx(env, lp, begin, end):
    xlen = end - begin + 1
    x    = _safeDoubleArray(xlen)
    status = CR.CPXXgetx(env, lp, x, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(x, xlen)

def getnumcores(env):
    numcores = CR.intPtr()
    status = CR.CPXXgetnumcores(env, numcores)
    check_status(env, status)
    return numcores.value()

def getax(env, lp, begin, end):
    axlen = end - begin + 1
    ax    = _safeDoubleArray(axlen)
    status = CR.CPXXgetax(env, lp, ax, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(ax, axlen)

def getxqxax(env, lp, begin, end):
    qaxlen = end - begin + 1
    qax    = _safeDoubleArray(qaxlen)
    status = CR.CPXXgetxqxax(env, lp, qax, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(qax, qaxlen)

def getpi(env, lp, begin, end):
    pilen = end - begin + 1
    pi    = _safeDoubleArray(pilen)
    status = CR.CPXXgetpi(env, lp, pi, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(pi, pilen)

def getslack(env, lp, begin, end):
    slacklen = end - begin + 1
    slack    = _safeDoubleArray(slacklen)
    status = CR.CPXXgetslack(env, lp, slack, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(slack, slacklen)

def getdj(env, lp, begin, end):
    djlen = end - begin + 1
    dj    = _safeDoubleArray(djlen)
    status = CR.CPXXgetdj(env, lp, dj, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(dj, djlen)

def getqconstrdslack(env, lp, qind):
    inout_list = [0, qind]
    status = CR.CPXXgetqconstrdslack(env, lp, inout_list)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list == [0]:
        return ([], [])
    inout_list.extend([qind])
    status = CR.CPXXgetqconstrdslack(env, lp, inout_list)
    check_status(env, status)
    return tuple(inout_list)


# Infeasibility

def getrowinfeas(env, lp, x, begin, end):
    infeasoutlen = end - begin + 1
    infeasout    = _safeDoubleArray(infeasoutlen)
    status = CR.CPXXgetrowinfeas(env, lp, LAU.double_list_to_array(x),
                                 infeasout, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(infeasout, infeasoutlen)

def getcolinfeas(env, lp, x, begin, end):
    infeasoutlen = end - begin + 1
    infeasout    = _safeDoubleArray(infeasoutlen)
    status = CR.CPXXgetcolinfeas(env, lp, LAU.double_list_to_array(x),
                                 infeasout, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(infeasout, infeasoutlen)

def getqconstrinfeas(env, lp, x, begin, end):
    infeasoutlen = end - begin + 1
    infeasout    = _safeDoubleArray(infeasoutlen)
    status = CR.CPXXgetqconstrinfeas(env, lp, LAU.double_list_to_array(x),
                                     infeasout, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(infeasout, infeasoutlen)

def getindconstrinfeas(env, lp, x, begin, end):
    infeasoutlen = end - begin + 1
    infeasout    = _safeDoubleArray(infeasoutlen)
    status = CR.CPXXgetindconstrinfeas(env, lp, LAU.double_list_to_array(x),
                                       infeasout, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(infeasout, infeasoutlen)

def getsosinfeas(env, lp, x, begin, end):
    infeasoutlen = end - begin + 1
    infeasout    = _safeDoubleArray(infeasoutlen)
    status = CR.CPXXgetsosinfeas(env, lp, LAU.double_list_to_array(x),
                                 infeasout, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(infeasout, infeasoutlen)

# Basis

def getbase_c(env, lp):
    numcols = CR.CPXXgetnumcols(env, lp)
    cstat   = _safeIntArray(numcols)
    rstat   = LAU.int_list_to_array([])
    status = CR.CPXXgetbase(env, lp, cstat, rstat)
    check_status(env, status)
    return LAU.int_array_to_list(cstat, numcols)

def getbase_r(env, lp):
    numrows = CR.CPXXgetnumrows(env, lp)
    cstat   = LAU.int_list_to_array([])
    rstat   = _safeIntArray(numrows)
    status = CR.CPXXgetbase(env, lp, cstat, rstat)
    check_status(env, status)
    return LAU.int_array_to_list(cstat, numrows)

def getbase(env, lp):
    numcols = CR.CPXXgetnumcols(env, lp)
    numrows = CR.CPXXgetnumrows(env, lp)
    cstat   = _safeIntArray(numcols)
    rstat   = _safeIntArray(numrows)
    status = CR.CPXXgetbase(env, lp, cstat, rstat)
    check_status(env, status)
    return (LAU.int_array_to_list(cstat, numcols), LAU.int_array_to_list(rstat, numrows))

def getbhead(env, lp):
    headlen = CR.CPXXgetnumrows(env, lp)
    head    = _safeIntArray(headlen)
    x       = _safeDoubleArray(headlen)
    status  = CR.CPXXgetbhead(env, lp, head, x)
    check_status(env, status)
    return (LAU.int_array_to_list(head, headlen),
            LAU.double_array_to_list(x, headlen))

def mbasewrite(env, lp, filename):
    status = CR.CPXXmbasewrite(env, lp, filename)
    check_status(env, status)
    return

def getijrow(env, lp, i, row_or_column):
    row = CR.intPtr()
    if row_or_column == 'r' or row_or_column == 'R':
        status = CR.CPXXgetijrow(env, lp, i, -1, row)
    elif row_or_column == 'c' or row_or_column == 'C':
        status = CR.CPXXgetijrow(env, lp, -1, i, row)
    if status == CR.CPXERR_INDEX_NOT_BASIC:
        return -1
    else:
        check_status(env, status)
    return row.value()

def getpnorms(env, lp):
    numcols = CR.CPXXgetnumcols(env, lp)
    numrows = CR.CPXXgetnumrows(env, lp)
    cnorm   = _safeDoubleArray(numcols)
    rnorm   = _safeDoubleArray(numrows)
    length  = CR.intPtr()
    status  = CR.CPXXgetpnorms(env, lp, cnorm, rnorm, length)
    check_status(env, status)
    return (LAU.double_array_to_list(cnorm, length.value()),
            LAU.double_array_to_list(rnorm, numrows))

def getdnorms(env, lp):
    numrows = CR.CPXXgetnumrows(env, lp)
    norm    = _safeDoubleArray(numrows)
    head    = _safeIntArray(numrows)
    length  = CR.intPtr()
    status  = CR.CPXXgetdnorms(env, lp, norm, head, length)
    check_status(env, status)
    return (LAU.double_array_to_list(norm, length.value()),
            LAU.int_array_to_list(head, length.value()))

def getbasednorms(env, lp):
    numcols = CR.CPXXgetnumcols(env, lp)
    numrows = CR.CPXXgetnumrows(env, lp)
    cstat   = _safeIntArray(numcols)
    rstat   = _safeIntArray(numrows)
    dnorm   = _safeDoubleArray(numrows)
    status = CR.CPXXgetbasednorms(env, lp, cstat, rstat, dnorm)
    check_status(env, status)
    return (LAU.int_array_to_list(cstat, numcols),
            LAU.int_array_to_list(rstat, numrows),
            LAU.double_array_to_list(dnorm, numrows))

def getpsbcnt(env, lp):
    return CR.CPXXgetpsbcnt(env, lp)

def getdsbcnt(env, lp):
    return CR.CPXXgetdsbcnt(env, lp)





def getdblquality(env, lp, what):
    quality = CR.doublePtr()
    status = CR.CPXXgetdblquality(env, lp, quality, what)
    check_status(env, status)
    return quality.value()

def getintquality(env, lp, what):
    quality = CR.intPtr()
    status = CR.CPXXgetintquality(env, lp, quality, what)
    check_status(env, status)
    return quality.value()



# Sensitivity Analysis Results 

def boundsa_lower(env, lp, begin, end):
    listlen = end - begin + 1
    lblower = _safeDoubleArray(listlen)
    lbupper = _safeDoubleArray(listlen)
    ublower = LAU.double_list_to_array([])
    ubupper = LAU.double_list_to_array([])
    status = CR.CPXXboundsa(env, lp, begin, end, lblower, lbupper,
                            ublower, ubupper)
    check_status(env, status)
    return (LAU.double_array_to_list(lblower, listlen),
            LAU.double_array_to_list(lbupper, listlen))

def boundsa_upper(env, lp, begin, end):
    listlen = end - begin + 1
    lblower = LAU.double_list_to_array([])
    lbupper = LAU.double_list_to_array([])
    ublower = _safeDoubleArray(listlen)
    ubupper = _safeDoubleArray(listlen)
    status = CR.CPXXboundsa(env, lp, begin, end, lblower, lbupper,
                            ublower, ubupper)
    check_status(env, status)
    return (LAU.double_array_to_list(ublower, listlen),
            LAU.double_array_to_list(ubupper, listlen))

def boundsa(env, lp, begin, end):
    listlen = end - begin + 1
    lblower = _safeDoubleArray(listlen)
    lbupper = _safeDoubleArray(listlen)
    ublower = _safeDoubleArray(listlen)
    ubupper = _safeDoubleArray(listlen)
    status = CR.CPXXboundsa(env, lp, begin, end, lblower, lbupper,
                            ublower, ubupper)
    check_status(env, status)
    return (LAU.double_array_to_list(lblower, listlen),
            LAU.double_array_to_list(lbupper, listlen),
            LAU.double_array_to_list(ublower, listlen),
            LAU.double_array_to_list(ubupper, listlen))

def objsa(env, lp, begin, end):
    listlen = end - begin + 1
    lower   = _safeDoubleArray(listlen)
    upper   = _safeDoubleArray(listlen)
    status = CR.CPXXobjsa(env, lp, begin, end, lower, upper)
    check_status(env, status)
    return (LAU.double_array_to_list(lower, listlen),
            LAU.double_array_to_list(upper, listlen))

def rhssa(env, lp, begin, end):
    listlen = end - begin + 1
    lower   = _safeDoubleArray(listlen)
    upper   = _safeDoubleArray(listlen)
    status = CR.CPXXrhssa(env, lp, begin, end, lower, upper)
    check_status(env, status)
    return (LAU.double_array_to_list(lower, listlen),
            LAU.double_array_to_list(upper, listlen))


# Conflicts

def refinemipstartconflictext(env, lp, mipstartindex, grppref, grpbeg, grpind, grptype):
    grpcnt = len(grppref)
    concnt = len(grpind)
    sigint_swap()
    status = CR.CPXXrefinemipstartconflictext(env, lp, mipstartindex, grpcnt, concnt,
                                             LAU.double_list_to_array(grppref),
                                             LAU.int_list_to_array(grpbeg),
                                             LAU.int_list_to_array(grpind), grptype)
    sigint_swap()
    check_status(env, status)
    return

def refineconflict(env, lp):
    confnumrows = CR.intPtr()
    confnumcols = CR.intPtr()
    sigint_swap()
    status      = CR.CPXXrefineconflict(env, lp, confnumrows, confnumcols)
    sigint_swap()
    check_status(env, status)
    return (confnumrows.value(), confnumcols.value())

def refineconflictext(env, lp, grppref, grpbeg, grpind, grptype):
    grpcnt = len(grppref)
    concnt = len(grpind)
    sigint_swap()
    status = CR.CPXXrefineconflictext(env, lp, grpcnt, concnt,
                                     LAU.double_list_to_array(grppref), LAU.int_list_to_array(grpbeg),
                                     LAU.int_list_to_array(grpind), grptype)
    sigint_swap()
    check_status(env, status)
    return

def getconflict(env, lp, n_r, n_c):
    c_stat = CR.intPtr()
    n_rows = CR.intPtr()
    n_cols = CR.intPtr()
    rowind = _safeIntArray(n_r)
    rowstt = _safeIntArray(n_r)
    colind = _safeIntArray(n_c)
    colstt = _safeIntArray(n_c)
    status = CR.CPXXgetconflict(env, lp, c_stat, rowind, rowstt, n_rows, colind, colstt, n_cols)
    check_status(env, status)
    return (LAU.int_array_to_list(rowind, n_rows.value()),
            LAU.int_array_to_list(rowstt, n_rows.value()),
            LAU.int_array_to_list(colind, n_cols.value()),
            LAU.int_array_to_list(colstt, n_cols.value()))

def getconflictext(env, lp, begin, end):
    grpstatlen = end - begin + 1
    grpstat    = _safeIntArray(grpstatlen)
    status = CR.CPXXgetconflictext(env, lp, grpstat, begin, end)
    check_status(env, status)
    return LAU.int_array_to_list(grpstat, grpstatlen)

def clpwrite(env, lp, filename):
    status = CR.CPXXclpwrite(env, lp, filename)
    check_status(env, status)
    return

# Problem Modification Routines

# File Reading Routines

# File Writing Routines

def solwrite(env, lp, filename):
    status = CR.CPXXsolwrite(env, lp, filename)
    check_status(env, status)
    return

# Message Handling Routines

# Advanced LP Routines

def binvcol(env, lp, j):
    xlen   = CR.CPXXgetnumrows(env, lp)
    x      = _safeDoubleArray(xlen)
    status = CR.CPXXbinvcol(env, lp, j, x)
    check_status(env, status)
    return LAU.double_array_to_list(x, xlen)

def binvrow(env, lp, i):
    ylen   = CR.CPXXgetnumrows(env, lp)
    y      = _safeDoubleArray(ylen)
    status = CR.CPXXbinvrow(env, lp, i, y)
    check_status(env, status)
    return LAU.double_array_to_list(y, ylen)

def binvacol(env, lp, j):
    xlen   = CR.CPXXgetnumrows(env, lp)
    x      = _safeDoubleArray(xlen)
    status = CR.CPXXbinvacol(env, lp, j, x)
    check_status(env, status)
    return LAU.double_array_to_list(x, xlen)

def binvarow(env, lp, i):
    zlen   = CR.CPXXgetnumcols(env, lp)
    z      = _safeDoubleArray(zlen)
    status = CR.CPXXbinvarow(env, lp, i, z)
    check_status(env, status)
    return LAU.double_array_to_list(z, zlen)

def ftran(env, lp, x):
    x_array = LAU.double_list_to_array(x)
    status = CR.CPXXftran(env, lp, x_array)
    check_status(env, status)
    return LAU.double_array_to_list(x_array, len(x))

def btran(env, lp, y):
    y_array = LAU.double_list_to_array(y)
    status = CR.CPXXbtran(env, lp, y_array)
    check_status(env, status)
    return LAU.double_array_to_list(y_array, len(y))

def getobjoffset(env, lp):
    objoffset = CR.doublePtr()
    status = CR.CPXXgetobjoffset(env, lp, objoffset)
    check_status(env, status)
    return objoffset.value()

def copybasednorms(env, lp, cstat, rstat, dnorm):
    status = CR.CPXXcopybasednorms(env, lp,
                                  LAU.int_list_to_array(cstat), LAU.int_list_to_array(rstat),
                                  LAU.double_list_to_array(dnorm))
    check_status(env, status)
    return

def copydnorms(env, lp, norm, head):
    status = CR.CPXXcopydnorms(env, lp, LAU.double_list_to_array(norm),
                              LAU.int_list_to_array(head), len(norm))
    check_status(env, status)
    return

def killdnorms(lp):
    status = CR.CPXXkilldnorms(lp)
    check_status(env, status)
    return

def copypnorms(env, lp, cnorm, rnorm):
    status = CR.CPXXcopypnorms(env, lp, LAU.double_list_to_array(cnorm),
                              LAU.double_list_to_array(rnorm), len(cnorm))
    check_status(env, status)
    return

def killpnorms(lp):
    status = CR.CPXXkillpnorms(lp)
    check_status(env, status)
    return


# Advanced Solution functions

def getgrad(env, lp, j):
    numrows = CR.CPXXgetnumrows(env, lp)
    head    = _safeIntArray(numrows)
    y       = _safeDoubleArray(numrows)
    status = CR.CPXXgetgrad(env, lp, j, head, y)
    check_status(env, status)
    return (LAU.int_array_to_list(head, numrows),
            LAU.double_array_to_list(y, numrows))

def slackfromx(env, lp, x):
    numrows = CR.CPXXgetnumrows(env, lp)
    slack   = _safeDoubleArray(numrows)
    status  = CR.CPXXslackfromx(env, lp, LAU.double_list_to_array(x), slack)
    check_status(env, status)
    return (LAU.double_array_to_list(slack, numrows))

def qconstrslackfromx(env, lp, x):
    numqcon = CR.CPXXgetnumqconstrs(env, lp)
    slack   = _safeDoubleArray(numqcon)
    status  = CR.CPXXqconstrslackfromx(env, lp,
                                       LAU.double_list_to_array(x), slack)
    check_status(env, status)
    return (LAU.double_array_to_list(slack, numqcon))

def djfrompi(env, lp, pi):
    numcols = CR.CPXXgetnumcols(env, lp)
    dj      = _safeDoubleArray(numcols)
    status  = CR.CPXXdjfrompi(env, lp, LAU.double_list_to_array(pi), dj)
    check_status(env, status)
    return (LAU.double_array_to_list(dj, numcols))

def qpdjfrompi(env, lp, pi, x):
    numcols = CR.CPXXgetnumcols(env, lp)
    dj      = _safeDoubleArray(numcols)
    status  = CR.CPXXqpdjfrompi(env, lp, LAU.double_list_to_array(pi),
                               LAU.double_list_to_array(x), dj)
    check_status(env, status)
    return (LAU.double_array_to_list(dj, numcols))

def mdleave(env, lp, goodlist):
    goodlen   = len(goodlist)
    downratio = _safeDoubleArray(goodlen)
    upratio   = _safeDoubleArray(goodlen)
    status = CR.CPXXmdleave(env, lp, LAU.int_list_to_array(goodlist),
                            goodlen, downratio, upratio)
    check_status(env, status)
    return (LAU.double_array_to_list(downratio, goodlen),
            LAU.double_array_to_list(upratio, goodlen))

def qpindefcertificate(env, lp):
    certlen = CR.CPXXgetnumquad(env, lp)
    cert    = _safeDoubleArray(certlen)
    status  = CR.CPXXqpindefcertificate(env, lp, cert)
    check_status(env, status)
    return LAU.double_array_to_list(cert, certlen)

def dualfarkas(env, lp):
    ylen   = CR.CPXXgetnumrows(env, lp)
    y      = _safeDoubleArray(ylen)
    proof  = CR.doublePtr()
    status = CR.CPXXdualfarkas(env, lp, y, proof)
    check_status(env, status)
    return (LAU.double_array_to_list(y, ylen), proof.value())

def getijdiv(env, lp):
    idiv = CR.intPtr()
    jdiv = CR.intPtr()
    status = CR.CPXXgetijdiv(env, lp, idiv, jdiv)
    check_status(env, status)
    if idiv.value() != -1:
        return idiv.value() + getnumcols(env, lp)
    elif jdiv.value() != -1:
        return  jdiv.value()
    else: # problem is not diverging
        return -1
    
def getray(env, lp):
    zlen   = CR.CPXXgetnumcols(env, lp)
    z      = _safeDoubleArray(zlen)
    status = CR.CPXXgetray(env, lp, z)
    check_status(env, status)
    return LAU.double_array_to_list(z, zlen)



# Advanced Presolve Routines

def presolve(env, lp, method):
    status = CR.CPXXpresolve(env, lp, method)
    check_status(env, status)
    return

def freepresolve(env, lp):
    status  = CR.CPXXfreepresolve(env, lp)
    check_status(env, status)
    return

def getredlp(env, lp):
    redlp  = CR.CPXLPptrPtr()
    status = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    return redlp.value()

def crushx(env, lp, x):
    redlp  = CR.CPXLPptrPtr()
    status = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    if redlp.value() is None:
        raise CplexError("No presolved problem found")
    numcols = CR.CPXXgetnumcols(env, redlp.value())
    prex    = _safeDoubleArray(numcols)
    status  = CR.CPXXcrushx(env, lp, LAU.double_list_to_array(x), prex)
    check_status(env, status)
    return LAU.double_array_to_list(prex, numcols)

def uncrushx(env, lp, prex):
    numcols = CR.CPXXgetnumcols(env, lp)
    x       = _safeDoubleArray(numcols)
    status  = CR.CPXXuncrushx(env, lp, x, LAU.double_list_to_array(prex))
    check_status(env, status)
    return LAU.double_array_to_list(x, numcols)


def crushpi(env, lp, pi):
    redlp  = CR.CPXLPptrPtr()
    status = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    if redlp.value() is None:
        raise CplexError("No presolved problem found")
    numrows = CR.CPXXgetnumrows(env, redlp.value())
    prepi    = _safeDoubleArray(numrows)
    status  = CR.CPXXcrushpi(env, lp, LAU.double_list_to_array(pi), prepi)
    check_status(env, status)
    return LAU.double_array_to_list(prepi, numrows)

def uncrushpi(env, lp, prepi):
    numrows = CR.CPXXgetnumrows(env, lp)
    pi      = _safeDoubleArray(numrows)
    status  = CR.CPXXuncrushpi(env, lp, pi, LAU.double_list_to_array(prepi))
    check_status(env, status)
    return LAU.double_array_to_list(pi, numrows)

def qpuncrushpi(env, lp, prepi, x):
    numrows = CR.CPXXgetnumrows(env, lp)
    pi      = _safeDoubleArray(numrows)
    status  = CR.CPXXqpuncrushpi(env, lp, pi, LAU.double_list_to_array(prepi),
                                 LAU.double_list_to_array(x))
    check_status(env, status)
    return LAU.double_array_to_list(pi, numrows)

def crushform(env, lp, ind, val):
    plen    = CR.intPtr()
    poffset = CR.doublePtr()
    redlp   = CR.CPXLPptrPtr()
    status  = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    if redlp.value() is None:
        raise CplexError("No presolved problem found")
    numcols = CR.CPXXgetnumcols(env, redlp.value())
    pind    = _safeIntArray(numcols)
    pval    = _safeDoubleArray(numcols)
    status  = CR.CPXXcrushform(env, lp, len(ind),
                               LAU.int_list_to_array(ind),
                               LAU.double_list_to_array(val),
                               plen, poffset, pind, pval)
    check_status(env, status)
    return (poffset.value(), LAU.int_array_to_list(pind, plen.value()),
            LAU.double_array_to_list(pval, plen.value()))
    

def uncrushform(env, lp, pind, pval):
    length = CR.intPtr()
    offset = CR.doublePtr()
    maxlen = CR.CPXXgetnumcols(env, lp) + CR.CPXXgetnumrows(env, lp)
    ind    = _safeIntArray(maxlen)
    val    = _safeDoubleArray(maxlen)
    status = CR.CPXXuncrushform(env, lp, len(pind),
                                LAU.int_list_to_array(pind),
                                LAU.double_list_to_array(pval),
                                length, offset, ind, val)
    check_status(env, status)
    return (offset.value(), LAU.int_array_to_list(ind, length.value()),
            LAU.double_array_to_list(val, length.value()))
    
def getprestat_status(env, lp):
    redlp   = CR.CPXLPptrPtr()
    status  = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    if redlp.value() is None:
        raise CplexError("No presolved problem found")
    prestat = CR.intPtr()
    pcstat  = LAU.int_list_to_array([])
    prstat  = LAU.int_list_to_array([])
    ocstat  = LAU.int_list_to_array([])
    orstat  = LAU.int_list_to_array([])
    status  = CR.CPXXgetprestat(env, lp, prestat, pcstat, prstat, ocstat, orstat)
    check_status(env, status)
    return prestat.value()

def getprestat_r(env, lp):
    redlp   = CR.CPXLPptrPtr()
    status  = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    if redlp.value() is None:
        raise CplexError("No presolved problem found")
    nrows   = CR.CPXXgetnumrows(env, lp)
    prestat = CR.intPtr()
    pcstat  = LAU.int_list_to_array([])
    prstat  = _safeIntArray(nrows)
    ocstat  = LAU.int_list_to_array([])
    orstat  = LAU.int_list_to_array([])
    status  = CR.CPXXgetprestat(env, lp, prestat, pcstat, prstat, ocstat, orstat)
    check_status(env, status)
    return LAU.int_array_to_list(prstat, nrows)

def getprestat_c(env, lp):
    redlp   = CR.CPXLPptrPtr()
    status  = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    if redlp.value() is None:
        raise CplexError("No presolved problem found")
    ncols   = CR.CPXXgetnumcols(env, lp)
    prestat = CR.intPtr()
    pcstat  = _safeIntArray(ncols)
    prstat  = LAU.int_list_to_array([])
    ocstat  = LAU.int_list_to_array([])
    orstat  = LAU.int_list_to_array([])
    status  = CR.CPXXgetprestat(env, lp, prestat, pcstat, prstat, ocstat, orstat)
    check_status(env, status)
    return LAU.int_array_to_list(pcstat, ncols)

def getprestat_or(env, lp):
    redlp   = CR.CPXLPptrPtr()
    status  = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    if redlp.value() is None:
        raise CplexError("No presolved problem found")
    nprows  = CR.CPXXgetnumrows(env, redlp.value())
    prestat = CR.intPtr()
    pcstat  = LAU.int_list_to_array([])
    prstat  = LAU.int_list_to_array([])
    ocstat  = LAU.int_list_to_array([])
    orstat  = _safeIntArray(nprows)
    status  = CR.CPXXgetprestat(env, lp, prestat, pcstat, prstat, ocstat, orstat)
    check_status(env, status)
    return LAU.int_array_to_list(orstat, nprows)

def getprestat_oc(env, lp):
    redlp   = CR.CPXLPptrPtr()
    status  = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    if redlp.value() is None:
        raise CplexError("No presolved problem found")
    npcols  = CR.CPXXgetnumcols(env, redlp.value())
    prestat = CR.intPtr()
    pcstat  = LAU.int_list_to_array([])
    prstat  = LAU.int_list_to_array([])
    ocstat  = _safeIntArray(npcols)
    orstat  = LAU.int_list_to_array([])
    status  = CR.CPXXgetprestat(env, lp, prestat, pcstat, prstat, ocstat, orstat)
    check_status(env, status)
    return LAU.int_array_to_list(ocstat, npcols)

def getprestat(env, lp):
    redlp   = CR.CPXLPptrPtr()
    status  = CR.CPXXgetredlp(env, lp, redlp)
    check_status(env, status)
    if redlp.value() is None:
        raise CplexError("No presolved problem found")
    ncols   = CR.CPXXgetnumcols(env, lp)
    nrows   = CR.CPXXgetnumrows(env, lp)
    npcols  = CR.CPXXgetnumcols(env, redlp.value())
    nprows  = CR.CPXXgetnumrows(env, redlp.value())
    prestat = CR.intPtr()
    pcstat  = _safeIntArray(ncols)
    prstat  = _safeIntArray(nrows)
    ocstat  = _safeIntArray(npcols)
    orstat  = _safeIntArray(nprows)
    status  = CR.CPXXgetprestat(env, lp, prestat, pcstat, prstat, ocstat, orstat)
    check_status(env, status)
    return (prestat.value(), LAU.int_array_to_list(pcstat, ncols),
            LAU.int_array_to_list(prstat, nrows),
            LAU.int_array_to_list(ocstat, npcols),
            LAU.int_array_to_list(orstat, nprows))

def prechgobj(env, lp, ind, val):
    status = CR.CPXXprechgobj(env, lp, len(ind), LAU.int_list_to_array(ind),
                             LAU.double_list_to_array(val))
    check_status(env, status)
    return

def preaddrows(env, lp, rhs, sense, rmatbeg, rmatind, rmatval, names,
               enc=default_encoding):
    status = CR.CPXXpreaddrows(env, lp, len(rmatbeg), len(rmatind),
                               LAU.double_list_to_array(rhs),
                               cpx_decode(sense, enc),
                               LAU.int_list_to_array(rmatbeg),
                               LAU.int_list_to_array(rmatind),
                               LAU.double_list_to_array(rmatval),
                               [cpx_decode(x, enc) for x in names])
    check_status(env, status)
    return


# Copying Data

def chgmipstarts(env, lp, mipstartindices, beg, varindices, values, effortlevel):
    status = CR.CPXXchgmipstarts(env, lp, len(mipstartindices), LAU.int_list_to_array(mipstartindices),
                                len(varindices), LAU.int_list_to_array(beg),
                                LAU.int_list_to_array(varindices), LAU.double_list_to_array(values),
                                LAU.int_list_to_array(effortlevel))
    check_status(env, status)
    return
    




def addmipstarts(env, lp, beg, varindices, values, effortlevel, mipstartname, enc=default_encoding):
    status = CR.CPXXaddmipstarts(env, lp, len(beg), len(varindices),
                                 LAU.int_list_to_array(beg),
                                 LAU.int_list_to_array(varindices),
                                 LAU.double_list_to_array(values),
                                 LAU.int_list_to_array(effortlevel),
                                 [cpx_decode(x, enc) for x in mipstartname])
    check_status(env, status)
    return
    
def delmipstarts(env, lp, begin, end):
    status = CR.CPXXdelmipstarts(env, lp, begin, end)
    check_status(env, status)
    return [0]

def delsetmipstarts(env, lp, delstat):
    delstat_array = LAU.int_list_to_array(delstat)
    status = CR.CPXXdelsetmipstarts(env, lp, delstat_array)
    check_status(env, status)
    return LAU.int_array_to_list(delstat_array, len(delstat))
    
# Optimizing Problems

# Progress

def getitcnt(env, lp):
    return CR.CPXXgetitcnt(env, lp)

def getphase1cnt(env, lp):
    return CR.CPXXgetphase1cnt(env, lp)

def getsiftitcnt(env, lp):
    return CR.CPXXgetsiftitcnt(env, lp)

def getsiftphase1cnt(env, lp):
    return CR.CPXXgetsiftphase1cnt(env, lp)

def getbaritcnt(env, lp):
    return CR.CPXXgetbaritcnt(env, lp)

def getcrossppushcnt(env, lp):
    return CR.CPXXgetcrossppushcnt(env, lp)

def getcrosspexchcnt(env, lp):
    return CR.CPXXgetcrosspexchcnt(env, lp)

def getcrossdpushcnt(env, lp):
    return CR.CPXXgetcrossdpushcnt(env, lp)

def getcrossdexchcnt(env, lp):
    return CR.CPXXgetcrossdexchcnt(env, lp)

def getmipitcnt(env, lp):
    return CR.CPXXgetmipitcnt(env, lp)

def getnodecnt(env, lp):
    return CR.CPXXgetnodecnt(env, lp)

def getnodeleftcnt(env, lp):
    return CR.CPXXgetnodeleftcnt(env, lp)



# MIP Only solution interface

def getbestobjval(env, lp):
    objval = CR.doublePtr()
    status = CR.CPXXgetbestobjval(env, lp, objval)
    check_status(env, status)
    return objval.value()

def getcutoff(env, lp):
    cutoff = CR.doublePtr()
    status = CR.CPXXgetcutoff(env, lp, cutoff)
    check_status(env, status)
    return cutoff.value()
    
def getmiprelgap(env, lp):
    relgap = CR.doublePtr()
    status = CR.CPXXgetmiprelgap(env, lp, relgap)
    check_status(env, status)
    return relgap.value()

def getnumcuts(env, lp, cuttype):
    num    = CR.intPtr()
    status = CR.CPXXgetnumcuts(env, lp, cuttype, num)
    check_status(env, status)
    return num.value()

def getnodeint(env, lp):
    return CR.CPXXgetnodeint(env, lp)
    
# MIP Starts

def getmipstarts_size(env, lp, begin, end):
    beglen      = end - begin + 1
    beg         = LAU.int_list_to_array([])
    effortlevel = _safeIntArray(beglen)
    nzcnt       = CR.intPtr()
    surplus     = CR.intPtr()
    varindices  = LAU.int_list_to_array([])
    values      = LAU.double_list_to_array([])
    startspace  = 0
    status = CR.CPXXgetmipstarts(env, lp, nzcnt, beg, varindices, values,
                                effortlevel, startspace, surplus, begin, end)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    return -surplus.value()

def getmipstarts_effort(env, lp, begin, end):
    beglen      = end - begin + 1
    beg         = LAU.int_list_to_array([])
    effortlevel = _safeIntArray(beglen)
    nzcnt       = CR.intPtr()
    surplus     = CR.intPtr()
    varindices  = LAU.int_list_to_array([])
    values      = LAU.double_list_to_array([])
    startspace  = 0
    status = CR.CPXXgetmipstarts(env, lp, nzcnt, beg, varindices, values,
                                effortlevel, startspace, surplus, begin, end)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if surplus.value() == 0:
        return ([0] * (end - begin + 1), [], [], [0] * (end - begin + 1))
    startspace = -surplus.value()
    beg         = _safeIntArray(beglen)
    varindices  = _safeIntArray(startspace)
    values      = _safeDoubleArray(startspace)
    status = CR.CPXXgetmipstarts(env, lp, nzcnt, beg, varindices, values,
                                effortlevel, startspace, surplus, begin, end)
    check_status(env, status)
    return LAU.int_array_to_list(effortlevel, beglen)

def getmipstarts(env, lp, begin, end):
    beglen      = end - begin + 1
    beg         = LAU.int_list_to_array([])
    effortlevel = _safeIntArray(beglen)
    nzcnt       = CR.intPtr()
    surplus     = CR.intPtr()
    varindices  = LAU.int_list_to_array([])
    values      = LAU.double_list_to_array([])
    startspace  = 0
    status = CR.CPXXgetmipstarts(env, lp, nzcnt, beg, varindices, values,
                                effortlevel, startspace, surplus, begin, end)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if surplus.value() == 0:
        return ([0] * (end - begin + 1), [], [], [0] * (end - begin + 1))
    beg         = _safeIntArray(beglen)
    startspace  = -surplus.value()
    varindices  = _safeIntArray(startspace)
    values      = _safeDoubleArray(startspace)
    status = CR.CPXXgetmipstarts(env, lp, nzcnt, beg, varindices, values,
                                effortlevel, startspace, surplus, begin, end)
    check_status(env, status)
    return (LAU.int_array_to_list(beg, beglen),
            LAU.int_array_to_list(varindices, startspace),
            LAU.double_array_to_list(values, startspace),
            LAU.int_array_to_list(effortlevel, beglen))

def getmipstartname(env, lp, begin, end, enc=default_encoding):
    inout_list = [0, begin, end]
    status = CR.CPXXgetmipstartname(env, lp, inout_list)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inout_list == [0]:
        return [cpx_encode_noop3(x, enc) for x in [""] * (end - begin + 1)]
    inout_list.extend([begin, end])
    status = CR.CPXXgetmipstartname(env, lp, inout_list)
    check_status(env, status)
    return [cpx_encode_noop3(x, enc) for x in inout_list]

def getmipstartindex(env, lp, mipstartname, enc=default_encoding):
    index = CR.intPtr()
    status = CR.CPXXgetmipstartindex(env, lp,
                                     cpx_decode_noop3(mipstartname, enc),
                                     index)
    check_status(env, status)
    return index.value()

def readcopymipstarts(env, lp, filename):
    status = CR.CPXXreadcopymipstarts(env, lp, filename)
    check_status(env, status)
    return 

def writemipstarts(env, lp, filename, begin, end):
    status = CR.CPXXwritemipstarts(env, lp, filename, begin, end)
    check_status(env, status)
    return 
    
def getsubstat(env, lp):
    return CR.CPXXgetsubstat(env, lp)


def getsubmethod(env, lp):
    return CR.CPXXgetsubmethod(env, lp)







# for callback query methods 

def get_wherefrom(cbstruct):
    return CR.get_wherefrom(cbstruct)

cpxlong_callback_info = [
    CPX_CALLBACK_INFO_ITCOUNT_LONG,
    CPX_CALLBACK_INFO_CROSSOVER_PPUSH_LONG,
    CPX_CALLBACK_INFO_CROSSOVER_PEXCH_LONG,
    CPX_CALLBACK_INFO_CROSSOVER_DPUSH_LONG,
    CPX_CALLBACK_INFO_CROSSOVER_DEXCH_LONG,
    CPX_CALLBACK_INFO_PRESOLVE_AGGSUBST_LONG,
    CPX_CALLBACK_INFO_NODES_LEFT_LONG,
    CPX_CALLBACK_INFO_NODE_COUNT_LONG,
    CPX_CALLBACK_INFO_MIP_ITERATIONS_LONG,
    CPX_CALLBACK_INFO_PRESOLVE_COEFFS_LONG,
    ]

int_callback_info = [
    CPX_CALLBACK_INFO_PRIMAL_FEAS,
    CPX_CALLBACK_INFO_DUAL_FEAS,
    CPX_CALLBACK_INFO_CROSSOVER_SBCNT,
    CPX_CALLBACK_INFO_PRESOLVE_ROWSGONE,
    CPX_CALLBACK_INFO_PRESOLVE_COLSGONE,
    CPX_CALLBACK_INFO_MIP_FEAS,
    CPX_CALLBACK_INFO_PROBE_PHASE,
    CPX_CALLBACK_INFO_CLIQUE_COUNT,
    CPX_CALLBACK_INFO_COVER_COUNT,
    CPX_CALLBACK_INFO_DISJCUT_COUNT,
    CPX_CALLBACK_INFO_FLOWCOVER_COUNT,
    CPX_CALLBACK_INFO_FLOWPATH_COUNT,
    CPX_CALLBACK_INFO_FRACCUT_COUNT,
    CPX_CALLBACK_INFO_GUBCOVER_COUNT,
    CPX_CALLBACK_INFO_IMPLBD_COUNT,
    CPX_CALLBACK_INFO_MIRCUT_COUNT,
    CPX_CALLBACK_INFO_ZEROHALFCUT_COUNT,
    CPX_CALLBACK_INFO_MCFCUT_COUNT,
    CPX_CALLBACK_INFO_LANDPCUT_COUNT,
    CPX_CALLBACK_INFO_USERCUT_COUNT,
    CPX_CALLBACK_INFO_TABLECUT_COUNT,
    CPX_CALLBACK_INFO_SOLNPOOLCUT_COUNT,
    CPX_CALLBACK_INFO_MY_THREAD_NUM,
    CPX_CALLBACK_INFO_USER_THREADS
    ]

double_callback_info = [
    CPX_CALLBACK_INFO_ENDTIME,
    CPX_CALLBACK_INFO_ENDDETTIME,
    CPX_CALLBACK_INFO_PRIMAL_OBJ,
    CPX_CALLBACK_INFO_DUAL_OBJ,
    CPX_CALLBACK_INFO_PRIMAL_INFMEAS,
    CPX_CALLBACK_INFO_DUAL_INFMEAS,
    CPX_CALLBACK_INFO_BEST_INTEGER,
    CPX_CALLBACK_INFO_BEST_REMAINING,
    CPX_CALLBACK_INFO_FRACCUT_PROGRESS,
    CPX_CALLBACK_INFO_DISJCUT_PROGRESS,
    CPX_CALLBACK_INFO_FLOWMIR_PROGRESS,
    CPX_CALLBACK_INFO_CUTOFF,
    CPX_CALLBACK_INFO_PROBE_PROGRESS,
    CPX_CALLBACK_INFO_TUNING_PROGRESS,
    CPX_CALLBACK_INFO_MIP_REL_GAP,
    CPX_CALLBACK_INFO_STARTTIME,
    CPX_CALLBACK_INFO_STARTDETTIME,
    ]
    
LPptr_callback_info = [
    CPX_CALLBACK_INFO_USER_PROBLEM,
    ]

cpxlong_callback_node_info = [
    CPX_CALLBACK_INFO_NODE_SEQNUM_LONG,
    CPX_CALLBACK_INFO_NODE_NODENUM_LONG,
    CPX_CALLBACK_INFO_NODE_DEPTH_LONG,
    ]
    
int_callback_node_info = [
    CPX_CALLBACK_INFO_NODE_NIINF,
    CPX_CALLBACK_INFO_NODE_VAR,
    CPX_CALLBACK_INFO_NODE_SOS,
    CPX_CALLBACK_INFO_LAZY_SOURCE,
    ]

double_callback_node_info = [
    CPX_CALLBACK_INFO_NODE_SIINF,
    CPX_CALLBACK_INFO_NODE_ESTIMATE,
    CPX_CALLBACK_INFO_NODE_OBJVAL,
    ]

char_callback_node_info = [
    CPX_CALLBACK_INFO_NODE_TYPE
    ]

user_handle_callback_node_info = [
    CPX_CALLBACK_INFO_NODE_USERHANDLE
    ]

def getcallbackinfo(cbstruct, whichinfo):
    if whichinfo in int_callback_info:
        data = CR.intPtr()
    elif whichinfo in double_callback_info:
        data = CR.doublePtr()
    elif whichinfo in LPptr_callback_info:
        data = CR.CPXLPptrPtr()
    elif whichinfo in cpxlong_callback_info:
        data = CR.cpxlongPtr()
    else:
        raise CplexError("invalid value for whichinfo in _internal._procedural.getcallbackinfo")
    status = CR.CPXXgetcallbackinfo(cbstruct, whichinfo, data)
    check_status(None, status)
    return data.value()


def getcallbacknodelp(cbstruct):
    nodelp = CR.CPXLPptrPtr()
    status = CR.CPXXgetcallbacknodelp(cbstruct, nodelp)
    check_status(None, status)
    return nodelp.value()

def getcallbacklp(cbstruct):
    lp     = CR.CPXLPptrPtr()
    status = CR.CPXXgetcallbacklp(cbstruct, lp)
    check_status(None, status)
    return lp.value()


def gettime(env):
    time = CR.doublePtr()
    status = CR.CPXXgettime(env, time)
    check_status(env, status)
    return time.value()

def getdettime(env):
    time = CR.doublePtr()
    status = CR.CPXXgetdettime(env, time)
    check_status(env, status)
    return time.value()

def getcallbackincumbent(cbstruct, begin, end):
    xlen = end - begin + 1
    x    = _safeDoubleArray(xlen)
    status = CR.CPXXgetcallbackincumbent(cbstruct, x, begin, end)
    check_status(None, status)
    return LAU.double_array_to_list(x, xlen)

def getcallbackpseudocosts(cbstruct, begin, end):
    pclen = end - begin + 1
    uppc  = _safeDoubleArray(pclen)
    dnpc  = _safeDoubleArray(pclen)
    status = CR.CPXXgetcallbackpseudocosts(cbstruct, uppc, dnpc, begin, end)
    check_status(None, status)
    return (LAU.double_array_to_list(uppc, pclen),
            LAU.double_array_to_list(dnpc, pclen))

def getcallbacknodeintfeas(cbstruct, begin, end):
    feaslen = end - begin + 1
    feas    = _safeIntArray(feaslen)
    status = CR.CPXXgetcallbacknodeintfeas(cbstruct, feas, begin, end)
    check_status(None, status)
    return LAU.int_array_to_list(feas, feaslen)
    
def getcallbacknodelb(cbstruct, begin, end):
    lblen = end - begin + 1
    lb    = _safeDoubleArray(lblen)
    status = CR.CPXXgetcallbacknodelb(cbstruct, lb, begin, end)
    check_status(None, status)
    return LAU.double_array_to_list(lb, lblen)

def getcallbacknodeub(cbstruct, begin, end):
    ublen = end - begin + 1
    ub    = _safeDoubleArray(ublen)
    status = CR.CPXXgetcallbacknodeub(cbstruct, ub, begin, end)
    check_status(None, status)
    return LAU.double_array_to_list(ub, ublen)

def getcallbacknodeobjval(cbstruct):
    x = CR.doublePtr()
    status = CR.CPXXgetcallbacknodeobjval(cbstruct, x)
    check_status(None, status)
    return x.value()

def getcallbacknodex(cbstruct, begin, end):
    xlen = end - begin + 1
    x    = _safeDoubleArray(xlen)
    status = CR.CPXXgetcallbacknodex(cbstruct, x, begin, end)
    check_status(None, status)
    return LAU.double_array_to_list(x, xlen)


def getcallbacknodeinfo(cbstruct, node, which):
    if which in int_callback_node_info:
        data = CR.intPtr()
    elif which in cpxlong_callback_node_info:
        data = CR.cpxlongPtr()
    elif which in double_callback_node_info:
        data = CR.doublePtr()
    elif which in char_callback_node_info:
        data = _getstrbfr()
    elif which in user_handle_callback_node_info:
        data = []
    else:
        raise CplexError("invalid value for which in _internal._procedural.getcallbacknodeinfo")
    status = CR.CPXXgetcallbacknodeinfo(cbstruct, [node, which, data])
    check_status(None, status)
    if which in int_callback_node_info or which in double_callback_node_info or which in cpxlong_callback_node_info:
        return data.value()
    elif which in char_callback_node_info:
        return data[0]
    elif which in user_handle_callback_node_info:
        return data[0]

def callbacksetuserhandle(cbstruct, userhandle):
    data = []
    status = CR.CPXXcallbacksetuserhandle(cbstruct, [userhandle, data])
    check_status(None, status)
    return data[0]

def callbacksetnodeuserhandle(cbstruct, nodeindex, userhandle):
    data = []
    status = CR.CPXXcallbacksetnodeuserhandle(cbstruct, [nodeindex, userhandle, data])
    check_status(None, status)
    return data[0]

def getcallbackseqinfo(cbstruct, node, which):
    if which in int_callback_node_info:
        data = CR.intPtr()
    elif which in cpxlong_callback_node_info:
        data = CR.cpxlongPtr()
    elif which in double_callback_node_info:
        data = CR.doublePtr()
    elif which in char_callback_node_info:
        data = _getstrbfr()
    elif which in user_handle_callback_node_info:
        data = []
    else:
        raise CplexError("invalid value for which in _internal._procedural.getcallbackseqinfo")
    status = CR.CPXXgetcallbackseqinfo(cbstruct, [node, which, data])
    check_status(None, status)
    if which in int_callback_node_info or which in double_callback_node_info or which in cpxlong_callback_node_info:
        return data.value()
    elif which in char_callback_node_info:
        return data[0]
    elif which in user_handle_callback_node_info:
        return data[0]


int_sos_info = [
    CPX_CALLBACK_INFO_SOS_NUM,
    CPX_CALLBACK_INFO_SOS_SIZE,
    CPX_CALLBACK_INFO_SOS_IS_FEASIBLE,
    CPX_CALLBACK_INFO_SOS_MEMBER_INDEX,
    ]

double_sos_info = [
    CPX_CALLBACK_INFO_SOS_MEMBER_REFVAL,
    ]

char_sos_info = [
    CPX_CALLBACK_INFO_SOS_TYPE,
    ]


def getcallbacksosinfo(cbstruct, sosindex, member, which):
    if which in int_sos_info:
        data = CR.intPtr()
    elif which in double_sos_info:
        data = CR.doublePtr()
    elif which in char_sos_info:
        data = _getstrbfr()
    else:
        raise CplexError("invalid value for which in _internal._procedural.getcallbacksosinfo")
    status = CR.CPXXgetcallbacksosinfo(cbstruct, sosindex, member, which, data)
    check_status(None, status)
    if which in int_sos_info or which in double_sos_info:
        return data.value()
    elif which in char_sos_info:
        return data[0]

def cutcallbackadd(cbstruct, rhs, sense, ind, val, purgeable):
    status = CR.CPXXcutcallbackadd(cbstruct, len(ind), rhs,
                                   cpx_decode(sense, default_encoding),
                                   LAU.int_list_to_array(ind),
                                   LAU.double_list_to_array(val),
                                   purgeable)
    check_status(None, status)
    return

def cutcallbackaddlocal(cbstruct, rhs, sense, ind, val):
    status = CR.CPXXcutcallbackaddlocal(cbstruct, len(ind), rhs,
                                        cpx_decode(sense, default_encoding),
                                        LAU.int_list_to_array(ind),
                                        LAU.double_list_to_array(val))
    check_status(None, status)
    return

def branchcallbackbranchgeneral(cbstruct, ind, lu, bd, rhs, sense, matbeg,
                                matind, matval, nodeest, userhandle):
    seqnum = CR.cpxlongPtr()
    status = CR.CPXXbranchcallbackbranchgeneral(
        cbstruct, len(ind),
        LAU.int_list_to_array(ind),
        lu,
        LAU.double_list_to_array(bd),
        len(matbeg), len(matind),
        LAU.double_list_to_array(rhs),
        cpx_decode(sense, default_encoding),
        LAU.int_list_to_array(matbeg),
        LAU.int_list_to_array(matind),
        LAU.double_list_to_array(matval),
        nodeest, userhandle, seqnum)
    check_status(None, status)
    return seqnum.value()

def branchcallbackbranchasCPLEX(cbstruct, n, userhandle):
    seqnum = CR.cpxlongPtr()
    status = CR.CPXXbranchcallbackbranchasCPLEX(cbstruct, n, userhandle, seqnum)
    check_status(None, status)
    return seqnum.value()

def setpydel(env):
    status = CR.setpydel(env)
    check_status(env, status)
    return

def delpydel(env):
    status = CR.delpydel(env)
    check_status(env, status)
    return

def setpyterminate(env):
    CR.setpyterminate(env)

# Solution pool

def addsolnpooldivfilter(env, lp, lb, ub, ind, wts, ref, name,
                         enc=default_encoding):
    status = CR.CPXXaddsolnpooldivfilter(env, lp, lb, ub, len(ind),
                                         LAU.int_list_to_array(ind),
                                         LAU.double_list_to_array(wts),
                                         LAU.double_list_to_array(ref),
                                         cpx_decode_noop3(name, enc))
    check_status(env, status)
    return

def addsolnpoolrngfilter(env, lp, lb, ub, ind, val, name,
                         enc=default_encoding):
    status = CR.CPXXaddsolnpoolrngfilter(env, lp, lb, ub, len(ind),
                                         LAU.int_list_to_array(ind),
                                         LAU.double_list_to_array(val),
                                         cpx_decode_noop3(name, enc))
    check_status(env, status)
    return

def getsolnpooldivfilter_constant(env, lp, which):
    lb = CR.doublePtr()
    ub = CR.doublePtr()
    nzcnt = CR.intPtr()
    space = 0
    surplus = CR.intPtr()
    ind = LAU.int_list_to_array([])
    wts = LAU.double_list_to_array([])
    ref = LAU.double_list_to_array([])
    status = CR.CPXXgetsolnpooldivfilter(env, lp, lb, ub, nzcnt, ind, wts, ref,
                                         space, surplus, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    return (lb.value(), ub.value(), -surplus.value())

def getsolnpooldivfilter(env, lp, which):
    lb = CR.doublePtr()
    ub = CR.doublePtr()
    nzcnt = CR.intPtr()
    space = 0
    surplus = CR.intPtr()
    ind = LAU.int_list_to_array([])
    wts = LAU.double_list_to_array([])
    ref = LAU.double_list_to_array([])
    status = CR.CPXXgetsolnpooldivfilter(env, lp, lb, ub, nzcnt, ind, wts, ref,
                                        space, surplus, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    space = -surplus.value()
    ind = _safeIntArray(space)
    wts = _safeDoubleArray(space)
    ref = _safeDoubleArray(space)
    status = CR.CPXXgetsolnpooldivfilter(env, lp, lb, ub, nzcnt, ind, wts, ref,
                                        space, surplus, which)
    check_status(env, status)
    return (lb.value(),
            ub.value(),
            LAU.int_array_to_list(ind, space),
            LAU.double_array_to_list(wts, space),
            LAU.double_array_to_list(ref, space))

def getsolnpoolrngfilter_constant(env, lp, which):
    lb = CR.doublePtr()
    ub = CR.doublePtr()
    nzcnt = CR.intPtr()
    space = 0
    surplus = CR.intPtr()
    ind = LAU.int_list_to_array([])
    val = LAU.double_list_to_array([])
    status = CR.CPXXgetsolnpoolrngfilter(env, lp, lb, ub, nzcnt, ind, val,
                                        space, surplus, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    return (lb.value(), ub.value(), -surplus.value())

def getsolnpoolrngfilter(env, lp, which):
    lb = CR.doublePtr()
    ub = CR.doublePtr()
    nzcnt = CR.intPtr()
    space = 0
    surplus = CR.intPtr()
    ind = LAU.int_list_to_array([])
    val = LAU.double_list_to_array([])
    status = CR.CPXXgetsolnpoolrngfilter(env, lp, lb, ub, nzcnt, ind, val,
                                        space, surplus, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    space = -surplus.value()
    ind = _safeIntArray(space)
    val = _safeDoubleArray(space)
    status = CR.CPXXgetsolnpoolrngfilter(env, lp, lb, ub, nzcnt, ind, val,
                                        space, surplus, which)
    check_status(env, status)
    return (lb.value(), ub.value(), LAU.int_array_to_list(ind, space),
            LAU.double_array_to_list(val, space))

def delsolnpoolfilters(env, lp, begin, end):
    status = CR.CPXXdelsolnpoolfilters(env, lp, begin, end)
    check_status(env, status)
    return [0]

def delsetsolnpoolfilters(env, lp, delstat):
    delstat_array = LAU.int_list_to_array(delstat)
    status = CR.CPXXdelsetsolnpoolfilters(env, lp, delstat_array)
    check_status(env, status)
    return LAU.int_array_to_list(delstat_array, len(delstat))

def getsolnpoolfiltername(env, lp, which, enc=default_encoding):
    inoutlist = [0]
    status = CR.CPXXgetsolnpoolfiltername(env, lp, inoutlist, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inoutlist == [0]:
        return cpx_encode_noop3("", enc)
    status = CR.CPXXgetsolnpoolfiltername(env, lp, inoutlist, which)
    check_status(env, status)
    return cpx_encode_noop3(inoutlist[0], enc)

def getsolnpoolnumfilters(env, lp):
    return CR.CPXXgetsolnpoolnumfilters(env, lp)

def fltwrite(env, lp, filename):
    status = CR.CPXXfltwrite(env, lp, filename)
    check_status(env, status)
    return

def readcopysolnpoolfilters(env, lp, filename):
    status = CR.CPXXreadcopysolnpoolfilters(env, lp, filename)
    check_status(env, status)
    return

def getsolnpoolfilterindex(env, lp, colname, enc=default_encoding):
    index = CR.intPtr()
    status = CR.CPXXgetsolnpoolfilterindex(env, lp,
                                           cpx_decode_noop3(colname, enc),
                                           index)
    check_status(env, status)
    return index.value()

def getsolnpoolfiltertype(env, lp, index):
    type_   = CR.intPtr()
    status = CR.CPXXgetsolnpoolfiltertype(env, lp, type_, index)
    check_status(env, status)
    return type_.value()

def delsolnpoolsolns(env, lp, begin, end):
    status = CR.CPXXdelsolnpoolsolns(env, lp, begin, end)
    check_status(env, status)
    return [0]

def delsetsolnpoolsolns(env, lp, delstat):
    delstat_array = LAU.int_list_to_array(delstat)
    status = CR.CPXXdelsetsolnpoolsolns(env, lp, delstat_array)
    check_status(env, status)
    return LAU.int_array_to_list(delstat_array, len(delstat))

def getsolnpoolnumsolns(env, lp):
    return CR.CPXXgetsolnpoolnumsolns(env, lp)

def getsolnpoolnumreplaced(env, lp):
    return CR.CPXXgetsolnpoolnumreplaced(env, lp)

def getsolnpoolsolnindex(env, lp, colname, enc=default_encoding):
    index = CR.intPtr()
    status = CR.CPXXgetsolnpoolsolnindex(env, lp, cpx_decode(colname, enc),
                                         index)
    check_status(env, status)
    return index.value()

def getsolnpoolmeanobjval(env, lp):
    objval = CR.doublePtr()
    status = CR.CPXXgetsolnpoolmeanobjval(env, lp, objval)
    check_status(env, status)
    return objval.value()

def getsolnpoolsolnname(env, lp, which, enc=default_encoding):
    inoutlist = [0]
    status = CR.CPXXgetsolnpoolsolnname(env, lp, inoutlist, which)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    if inoutlist == [0]:
        return cpx_encode_noop3("", enc)
    status = CR.CPXXgetsolnpoolsolnname(env, lp, inoutlist, which)
    check_status(env, status)
    return cpx_encode_noop3(inoutlist[0], enc)

def solwritesolnpool(env, lp, soln, filename):
    status = CR.CPXXsolwritesolnpool(env, lp, soln, filename)
    check_status(env, status)
    return

def solwritesolnpoolall(env, lp, filename):
    status = CR.CPXXsolwritesolnpoolall(env, lp, filename)
    check_status(env, status)
    return

def getsolnpoolobjval(env, lp, soln):
    obj = CR.doublePtr()
    status = CR.CPXXgetsolnpoolobjval(env, lp, soln, obj)
    check_status(env, status)
    return obj.value()

def getsolnpoolx(env, lp, soln, begin, end):
    xlen = end - begin + 1
    x    = _safeDoubleArray(xlen)
    status = CR.CPXXgetsolnpoolx(env, lp, soln, x, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(x, xlen)

def getsolnpoolslack(env, lp, soln, begin, end):
    slacklen = end - begin + 1
    slack    = _safeDoubleArray(slacklen)
    status = CR.CPXXgetsolnpoolslack(env, lp, soln, slack, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(slack, slacklen)

def getsolnpoolqconstrslack(env, lp, soln, begin, end):
    qlen = end - begin + 1
    q    = _safeDoubleArray(qlen)
    status = CR.CPXXgetsolnpoolqconstrslack(env, lp, soln, q, begin, end)
    check_status(env, status)
    return LAU.double_array_to_list(q, qlen)

def getsolnpoolintquality(env, lp, soln, what):
    quality = CR.intPtr()
    status = CR.CPXXgetsolnpoolintquality(env, lp, soln, quality, what)
    check_status(env, status)
    return quality.value()

def getsolnpooldblquality(env, lp, soln, what):
    quality = CR.doublePtr()
    status = CR.CPXXgetsolnpooldblquality(env, lp, soln, quality, what)
    check_status(env, status)
    return quality.value()



# Initial data



def copystart(env, lp, cstat, rstat, cprim, rprim, cdual, rdual):
    status = CR.CPXXcopystart(env, lp,
                           LAU.int_list_to_array(cstat), LAU.int_list_to_array(rstat),
                           LAU.double_list_to_array(cprim), LAU.double_list_to_array(rprim),
                           LAU.double_list_to_array(cdual), LAU.double_list_to_array(rdual))
    check_status(env, status)
    return

def readcopybase(env, lp, filename):
    status = CR.CPXXreadcopybase(env, lp, filename)
    check_status(env, status)
    return

def getorder(env, lp):
    count = CR.intPtr()
    surplus = CR.intPtr()
    space = 0
    ind = LAU.int_list_to_array([])
    pri = LAU.int_list_to_array([])
    dir_ = LAU.int_list_to_array([])
    status = CR.CPXXgetorder(env, lp, count, ind, pri, dir_, space, surplus)
    if status != CR.CPXERR_NEGATIVE_SURPLUS:
        check_status(env, status)
    space = -surplus.value()
    ind = _safeIntArray(space)
    pri = _safeIntArray(space)
    dir_ = _safeIntArray(space)
    status = CR.CPXXgetorder(env, lp, count, ind, pri, dir_, space, surplus)
    check_status(env, status)
    return (LAU.int_array_to_list(ind, space), LAU.int_array_to_list(pri, space),
            LAU.int_array_to_list(dir_, space))

def copyorder(env, lp, indices, priority, direction):
    status = CR.CPXXcopyorder(env, lp, len(indices), LAU.int_list_to_array(indices),
                             LAU.int_list_to_array(priority), LAU.int_list_to_array(direction))
    check_status(env, status)
    return

def readcopyorder(env, lp, filename):
    status = CR.CPXXreadcopyorder(env, lp, filename)
    check_status(env, status)
    return

def ordwrite(env, lp, filename):
    status = CR.CPXXordwrite(env, lp, filename)
    check_status(env, status)
    return    

def readcopysol(env, lp, filename):
    status = CR.CPXXreadcopysol(env, lp, filename)
    check_status(env, status)
    return

# handle the lock for parallel callbacks

def initlock(env):
    return CR.init_callback_lock(env)

def finitlock(env, lock):
    CR.finit_callback_lock(env, lock)
    return


# get problem statistics

def getprobstats(env, lp):
    rows_p = CR.intPtr()
    cols_p = CR.intPtr()
    objcnt_p = CR.intPtr()
    rhscnt_p = CR.intPtr()
    nzcnt_p = CR.intPtr()
    ecnt_p = CR.intPtr()
    gcnt_p = CR.intPtr()
    lcnt_p = CR.intPtr()
    rngcnt_p = CR.intPtr()
    ncnt_p = CR.intPtr()
    fcnt_p = CR.intPtr()
    xcnt_p = CR.intPtr()
    bcnt_p = CR.intPtr()
    ocnt_p = CR.intPtr()
    bicnt_p = CR.intPtr()
    icnt_p = CR.intPtr()
    scnt_p = CR.intPtr()
    sicnt_p = CR.intPtr()
    qpcnt_p = CR.intPtr()
    qpnzcnt_p = CR.intPtr()
    nqconstr_p = CR.intPtr()
    qrhscnt_p = CR.intPtr()
    qlcnt_p = CR.intPtr()
    qgcnt_p = CR.intPtr()
    quadnzcnt_p = CR.intPtr()
    linnzcnt_p = CR.intPtr()
    nindconstr_p = CR.intPtr()
    indrhscnt_p = CR.intPtr()
    indnzcnt_p = CR.intPtr()
    indcompcnt_p = CR.intPtr()
    indlcnt_p = CR.intPtr()
    indecnt_p = CR.intPtr()
    indgcnt_p = CR.intPtr()
    maxcoef_p = CR.doublePtr()
    mincoef_p = CR.doublePtr()
    minrhs_p = CR.doublePtr()
    maxrhs_p = CR.doublePtr()
    minrng_p = CR.doublePtr()
    maxrng_p = CR.doublePtr()
    minobj_p = CR.doublePtr()
    maxobj_p = CR.doublePtr()
    minlb_p = CR.doublePtr()
    maxub_p = CR.doublePtr()
    minqcoef_p = CR.doublePtr()
    maxqcoef_p = CR.doublePtr()
    minqcq_p = CR.doublePtr()
    maxqcq_p = CR.doublePtr()
    minqcl_p = CR.doublePtr()
    maxqcl_p = CR.doublePtr()
    minqcr_p = CR.doublePtr()
    maxqcr_p = CR.doublePtr()
    minind_p = CR.doublePtr()
    maxind_p = CR.doublePtr()
    minindrhs_p = CR.doublePtr()
    maxindrhs_p = CR.doublePtr()
    minlazy_p = CR.doublePtr()
    maxlazy_p = CR.doublePtr()
    minlazyrhs_p = CR.doublePtr()
    maxlazyrhs_p = CR.doublePtr()
    minucut_p = CR.doublePtr()
    maxucut_p = CR.doublePtr()
    minucutrhs_p = CR.doublePtr()
    maxucutrhs_p = CR.doublePtr()
    nsos_p = CR.intPtr()
    nsos1_p = CR.intPtr()
    sos1nmem_p = CR.intPtr()
    sos1type_p = CR.intPtr()
    nsos2_p = CR.intPtr()
    sos2nmem_p = CR.intPtr()
    sos2type_p = CR.intPtr()
    lazyrhscnt_p = CR.intPtr() 
    lazygcnt_p = CR.intPtr() 
    lazylcnt_p = CR.intPtr() 
    lazyecnt_p = CR.intPtr()
    lazycnt_p = CR.intPtr()
    lazynzcnt_p = CR.intPtr()
    ucutrhscnt_p = CR.intPtr() 
    ucutgcnt_p = CR.intPtr() 
    ucutlcnt_p = CR.intPtr() 
    ucutecnt_p = CR.intPtr()
    ucutcnt_p = CR.intPtr()
    ucutnzcnt_p = CR.intPtr()
    status = CR.CPXEgetprobstats(env, lp,
                                 rows_p,
                                 cols_p,
                                 objcnt_p,
                                 rhscnt_p,
                                 nzcnt_p,
                                 ecnt_p,
                                 gcnt_p,
                                 lcnt_p,
                                 rngcnt_p,
                                 ncnt_p,
                                 fcnt_p,
                                 xcnt_p,
                                 bcnt_p,
                                 ocnt_p,
                                 bicnt_p,
                                 icnt_p,
                                 scnt_p,
                                 sicnt_p,
                                 qpcnt_p,
                                 qpnzcnt_p,
                                 nqconstr_p,
                                 qrhscnt_p,
                                 qlcnt_p,
                                 qgcnt_p,
                                 quadnzcnt_p,
                                 linnzcnt_p,
                                 nindconstr_p,
                                 indrhscnt_p,
                                 indnzcnt_p,
                                 indcompcnt_p,
                                 indlcnt_p,
                                 indecnt_p,
                                 indgcnt_p,
                                 maxcoef_p,
                                 mincoef_p,
                                 minrhs_p,
                                 maxrhs_p,
                                 minrng_p,
                                 maxrng_p,
                                 minobj_p,
                                 maxobj_p,
                                 minlb_p,
                                 maxub_p,
                                 minqcoef_p,
                                 maxqcoef_p,
                                 minqcq_p,
                                 maxqcq_p,
                                 minqcl_p,
                                 maxqcl_p,
                                 minqcr_p,
                                 maxqcr_p,
                                 minind_p,
                                 maxind_p,
                                 minindrhs_p,
                                 maxindrhs_p,
                                 minlazy_p,
                                 maxlazy_p,
                                 minlazyrhs_p,
                                 maxlazyrhs_p,
                                 minucut_p,
                                 maxucut_p,
                                 minucutrhs_p,
                                 maxucutrhs_p,
                                 nsos_p,
                                 nsos1_p,
                                 sos1nmem_p,
                                 sos1type_p,
                                 nsos2_p,
                                 sos2nmem_p,
                                 sos2type_p,
                                 lazyrhscnt_p, 
                                 lazygcnt_p, 
                                 lazylcnt_p, 
                                 lazyecnt_p,
                                 lazycnt_p,
                                 lazynzcnt_p,
                                 ucutrhscnt_p, 
                                 ucutgcnt_p, 
                                 ucutlcnt_p, 
                                 ucutecnt_p,
                                 ucutcnt_p,
                                 ucutnzcnt_p)
    check_status(env, status)
    return [rows_p.value(),         #0
            cols_p.value(),         #1
            objcnt_p.value(),       #2
            rhscnt_p.value(),       #3
            nzcnt_p.value(),        #4
            ecnt_p.value(),         #5
            gcnt_p.value(),         #6
            lcnt_p.value(),         #7
            rngcnt_p.value(),       #8
            ncnt_p.value(),         #9
            fcnt_p.value(),         #10
            xcnt_p.value(),         #11
            bcnt_p.value(),         #12
            ocnt_p.value(),         #13
            bicnt_p.value(),        #14
            icnt_p.value(),         #15
            scnt_p.value(),         #16
            sicnt_p.value(),        #17
            qpcnt_p.value(),        #18
            qpnzcnt_p.value(),      #19
            nqconstr_p.value(),     #20
            qrhscnt_p.value(),      #21
            qlcnt_p.value(),        #22
            qgcnt_p.value(),        #23
            quadnzcnt_p.value(),    #24
            linnzcnt_p.value(),     #25
            nindconstr_p.value(),   #26
            indrhscnt_p.value(),    #27
            indnzcnt_p.value(),     #28
            indcompcnt_p.value(),   #29
            indlcnt_p.value(),      #30
            indecnt_p.value(),      #31
            indgcnt_p.value(),      #32
            maxcoef_p.value(),      #33
            mincoef_p.value(),      #34
            minrhs_p.value(),       #35
            maxrhs_p.value(),       #36
            minrng_p.value(),       #37
            maxrng_p.value(),       #38
            minobj_p.value(),       #39
            maxobj_p.value(),       #40
            minlb_p.value(),        #41
            maxub_p.value(),        #42
            minqcoef_p.value(),     #43
            maxqcoef_p.value(),     #44
            minqcq_p.value(),       #45
            maxqcq_p.value(),       #46
            minqcl_p.value(),       #47
            maxqcl_p.value(),       #48
            minqcr_p.value(),       #49
            maxqcr_p.value(),       #50
            minind_p.value(),       #51
            maxind_p.value(),       #52
            minindrhs_p.value(),    #53
            maxindrhs_p.value(),    #54
            minlazy_p.value(),      #55
            maxlazy_p.value(),      #56
            minlazyrhs_p.value(),   #57
            maxlazyrhs_p.value(),   #58
            minucut_p.value(),      #59
            maxucut_p.value(),      #60
            minucutrhs_p.value(),   #61
            maxucutrhs_p.value(),   #62
            nsos_p.value(),         #63
            nsos1_p.value(),        #64
            sos1nmem_p.value(),     #65
            sos1type_p.value(),     #66
            nsos2_p.value(),        #67
            sos2nmem_p.value(),     #68
            sos2type_p.value(),     #69
            lazyrhscnt_p.value(),   #70
            lazygcnt_p.value(),     #71
            lazylcnt_p.value(),     #72
            lazyecnt_p.value(),     #73
            lazycnt_p.value(),      #74
            lazynzcnt_p.value(),    #75
            ucutrhscnt_p.value(),   #76
            ucutgcnt_p.value(),     #77
            ucutlcnt_p.value(),     #78
            ucutecnt_p.value(),     #79
            ucutcnt_p.value(),      #80
            ucutnzcnt_p.value(),]   #81

# get histogram of non-zero row/column counts

def gethist(env, lp, key):
    if key == 'r':
        space = CR.CPXXgetnumcols(env, lp) + 1
    else:
        key = 'c'
        space = CR.CPXXgetnumrows(env, lp) + 1
    hist   = _safeIntArray(space)
    status = CR.CPXEgethist(env, lp, cpx_decode(key, default_encoding), hist)
    check_status(env, status)
    return LAU.int_array_to_list(hist, space)

# get solution quality metrics

def getqualitymetrics(env, lp, soln):
    space  = 26
    data   = _safeDoubleArray(space)
    ispace = 10
    idata  = _safeIntArray(ispace)
    status = CR.CPXEgetqualitymetrics(env, lp, soln, data, idata)
    check_status(env, status)
    return [LAU.int_array_to_list(idata, ispace),
            LAU.double_array_to_list(data, space)]
