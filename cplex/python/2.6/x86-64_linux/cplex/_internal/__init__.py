# --------------------------------------------------------------------------
# File: __init__.py 
# ---------------------------------------------------------------------------
# Licensed Materials - Property of IBM
# 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
# Copyright IBM Corporation 2008, 2014. All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with
# IBM Corp.
# ------------------------------------------------------------------------
"""


"""

import os
import sys


from . import _aux_functions
from . import _list_array_utils
from . import _ostream
from . import _procedural
from . import _constants
from . import _matrices
from . import _parameter_classes
from . import _parameter_hierarchy
from . import _subinterfaces
from . import _pycplex
from . import _parameters_auto
from ..exceptions import CplexError

__all__ = ["Environment", "_aux_functions", "_list_array_utils", "_ostream", "_procedural", "_constants", "_matrices", "_parameter_classes", "_subinterfaces", "_pycplex", "_parameters_auto"]


class ProblemType:
    """Types of problems the Cplex object can encapsulate.

       For explanations of the problem types, see those topics in the
       CPLEX User's Manual in the topic titled Continuous Optimization
       for LP, QP, and QCP or the topic titled Discrete Optimization 
       for MILP, FIXEDMILP, NODELP, NODEQP, MIQCP, NODEQCP.

    """
    LP         = _constants.CPXPROB_LP
    MILP       = _constants.CPXPROB_MILP
    fixed_MILP = _constants.CPXPROB_FIXEDMILP
    node_LP    = _constants.CPXPROB_NODELP
    QP         = _constants.CPXPROB_QP
    MIQP       = _constants.CPXPROB_MIQP
    fixed_MIQP = _constants.CPXPROB_FIXEDMIQP
    node_QP    = _constants.CPXPROB_NODEQP
    QCP        = _constants.CPXPROB_QCP
    MIQCP      = _constants.CPXPROB_MIQCP
    node_QCP   = _constants.CPXPROB_NODEQCP
    def __getitem__(self, item):
        """Converts a constant to a string.

        >>> import cplex
        >>> c = cplex.Cplex()
        >>> c.problem_type.LP
        0
        >>> c.problem_type[0]
        'LP'

        """
        if item == _constants.CPXPROB_LP:
            return 'LP'
        if item == _constants.CPXPROB_MILP:
            return 'MILP'
        if item == _constants.CPXPROB_FIXEDMILP:
            return 'fixed_MILP'
        if item == _constants.CPXPROB_NODELP:
            return 'node_LP'
        if item == _constants.CPXPROB_QP:
            return 'QP'
        if item == _constants.CPXPROB_MIQP:
            return 'MIQP'
        if item == _constants.CPXPROB_FIXEDMIQP:
            return 'fixed_MIQP'
        if item == _constants.CPXPROB_NODEQP:
            return 'node_QP'
        if item == _constants.CPXPROB_QCP:
            return 'QCP'
        if item == _constants.CPXPROB_MIQCP:
            return 'MIQCP'
        if item == _constants.CPXPROB_NODEQCP:
            return 'node_QCP'

class Environment(object):

    """non-public


    """

    def __init__(self):
        """non-public"""
        self._e = _procedural.openCPLEX()
        self.parameters = _parameter_classes.RootParameterGroup(self, _parameter_hierarchy.root_members)
        _procedural.setpyterminate(self._e)
        _procedural.set_status_checker()
        self.__lock = _procedural.initlock(self._e)
        self._callback_exception = None
        self._callbacks = []
        self.set_results_stream(sys.stdout)
        self.set_warning_stream(sys.stderr)
        self.set_error_stream(sys.stderr)
        self.set_log_stream(sys.stdout)
        
    def __del__(self):
        """non-public"""
        if hasattr(self, "_Environment__results_stream"):
            self.__delete_stream(0, self.__results_stream)
        if hasattr(self, "_Environment__warning_stream"):
            self.__delete_stream(1, self.__warning_stream)
        if hasattr(self, "_Environment__error_stream"):
            self.__delete_stream(2, self.__error_stream)
        if hasattr(self, "_Environment__log_stream"):
            self.__delete_stream(3, self.__log_stream)
        if hasattr(self, "_Environment__lock") and hasattr(self, "_e"):
            _procedural.finitlock(self._e, self.__lock)
        if hasattr(self, "_e"):
            _procedural.closeCPLEX(self._e)


    def _needs_delete_callback(self, callback_instance):
        """non-public"""
        # If the user has registered any callback that may change
        # the user data at a node then we need to register the
        # delete callback.
        # The Control, Node, and Incumbent callbacks have the set_node_data
        # method (and all who inherit from these).
        return hasattr(callback_instance, "set_node_data")

    def register_callback(self, callback_class):
        """Registers a callback for use when solving.

        callback_class must be a proper subclass of one of the
        callback classes defined in the module callbacks.  It must
        override the __call__ method with a method that has signature
        __call__(self) -> None.  If callback_class is a subclass of
        more than one callback class, it will only be called when its
        first superclass is called.  register_callback returns the
        instance of callback_class registered for use.  Any previously
        registered callback of the same class will no longer be
        registered.

        """
        cb = callback_class(self)
        if cb._cb_type_string is None:
            raise CplexError(str(callback_class) +
                             " is not a subclass of a subclassable Callback class.")
        if hasattr(cb, "_unregister"):
            cb._cb_set_function(self._e, None)
        else:
            # Count the callbacks that are installed and require a
            # delete callback.
            num_delete = 0
            for c in self._callbacks:
                if self._needs_delete_callback(c):
                    num_delete = num_delete + 1
            old_cb = getattr(self, "_" + cb._cb_type_string + "_callback", None)
            if old_cb is not None:
                self._callbacks.remove(old_cb)
            setattr(self, "_" + cb._cb_type_string + "_callback", cb)
            if cb._cb_type_string == "MIP_info":
                cb._cb_set_function(self._e, self._MIP_info_callback)
            else:
                cb._cb_set_function(self._e, self)
            self._callbacks.append(cb)
            if self._needs_delete_callback(cb) and num_delete < 1:
                # We need a delete callback and did not have one
                # before -> install it.
                _procedural.setpydel(self._e)
        return cb

    def unregister_callback(self, callback_class):
        """Unregisters a callback.

        callback_class must be one of the callback classes defined in
        the module callback or a subclass of one of them.  This method 
        unregisters any previously registered callback of the same
        class.  If callback_class is a subclass of more than one
        callback class, this method unregisters only the callback of the
        same type as its first superclass.  unregister_callback
        returns the instance of callback_class just unregistered.

        """
        cb = callback_class(self)
        current_cb = getattr(self, "_" + cb._cb_type_string + "_callback", None)
        if current_cb is not None:
            # Count the number of installed callbacks that require
            # need a delete callback
            num_delete = 0
            for c in self._callbacks:
                if self._needs_delete_callback(c):
                    num_delete = num_delete + 1
            if self._needs_delete_callback(current_cb) and num_delete < 2:
                # We are about to remove the last callback that requires
                # a delete callback.
                _procedural.delpydel(self._e)
            self._callbacks.remove(current_cb)
            class do_nothing(callback_class):
                def __init__(self, env):
                    callback_class.__init__(self, env)
                    self._unregister = True
                def __call__(self):
                    return
            self.register_callback(do_nothing)
        return current_cb

    def __add_stream(self, which_channel, stream):
        """non-public"""
        channel = _procedural.getchannels(self._e)[which_channel]
        _procedural.addfuncdest(self._e, channel, stream)

    def __delete_stream(self, which_channel, stream):
        """non-public"""
        channel = _procedural.getchannels(self._e)[which_channel]
        _procedural.delfuncdest(self._e, channel, stream)
        del stream

    def set_results_stream(self, results_file, fn=None):
        """Specifies where results will be printed.

        The first argument must be either a file-like object (that is, an
        object with a write method and a flush method) or the name of
        a file to be written to.  Use None as the first argument to
        suppress output.

        The second optional argument is a function that takes a string
        as input and returns a string.  If specified, strings sent to
        this stream will be processed by this function before being
        written.

        Returns the stream to which results will be written.  To write
        to this stream, use the write() method of this object.

        """
        which = 0
        if hasattr(self, "_Environment__results_stream"):
            self.__delete_stream(which, self.__results_stream)
        self.__results_stream = _ostream.OutputStream(results_file, self, fn)
        self.__add_stream(which, self.__results_stream)
        return self.__results_stream
    
    def set_warning_stream(self, warning_file, fn=None):
        """Specifies where warnings will be printed.

        The first argument must be either a file-like object (that is, an
        object with a write method and a flush method) or the name of
        a file to be written to.  Use None as the first argument to
        suppress output.

        The second optional argument is a function that takes a string
        as input and returns a string.  If specified, strings sent to
        this stream will be processed by this function before being
        written.

        Returns the stream to which warnings will be written.  To write
        to this stream, use the write() method of this object.

        """
        which = 1
        if hasattr(self, "_Environment__warning_stream"):
            self.__delete_stream(which, self.__warning_stream)
        self.__warning_stream = _ostream.OutputStream(warning_file, self, fn)
        self.__add_stream(which, self.__warning_stream)
        return self.__warning_stream

    def set_error_stream(self, error_file, fn=None):
        """Specifies where errors will be printed.

        The first argument must be either a file-like object (that is, an
        object with a write method and a flush method) or the name of
        a file to be written to.  Use None as the first argument to
        suppress output.

        The second optional argument is a function that takes a string
        as input and returns a string.  If specified, strings sent to
        this stream will be processed by this function before being
        written.

        Returns the stream to which errors will be written.  To write
        to this stream, use the write() method of this object.

        """
        which = 2
        if hasattr(self, "_Environment__error_stream"):
            self.__delete_stream(which, self.__error_stream)
        self.__error_stream = _ostream.OutputStream(error_file, self, fn)
        self.__error_stream._error_string = None
        self.__add_stream(which, self.__error_stream)
        return self.__error_stream

    def set_log_stream(self, log_file, fn=None):
        """Specifies where the log will be printed.

        The first argument must be either a file-like object (that is, an
        object with a write method and a flush method) or the name of
        a file to be written to.  Use None as the first argument to
        suppress output.

        The second optional argument is a function that takes a string
        as input and returns a string.  If specified, strings sent to
        this stream will be processed by this function before being
        written.

        Returns the stream to which the log will be written.  To write
        to this stream, use this object's write() method.

        """
        which = 3
        if hasattr(self, "_Environment__log_stream"):
            self.__delete_stream(which, self.__log_stream)
        self.__log_stream = _ostream.OutputStream(log_file, self, fn)
        self.__add_stream(which, self.__log_stream)
        return self.__log_stream

    def get_version(self):
        """Returns a string specifying the version of CPLEX."""
        return _procedural.version(self._e)

    def get_versionnumber(self):
        """Returns an integer specifying the version of CPLEX.

        The version of CPLEX is in the format vvrrmmff, where vv is
        the version, rr is the release, mm is the modification, and ff
        is the fixpack number. For example, for CPLEX version 12.5.0.1
        the returned value is 12050001.
        """
        return _procedural.versionnumber(self._e)

    def get_num_cores(self):
        """Returns the number of cores on this machine."""
        return _procedural.getnumcores(self._e)

    def get_time(self):
        """Returns a timestamp in CPU or wallclock seconds from CPLEX."""
        return _procedural.gettime(self._e)

    def get_dettime(self):
        """Returns the current deterministic time in ticks."""
        return _procedural.getdettime(self._e)
