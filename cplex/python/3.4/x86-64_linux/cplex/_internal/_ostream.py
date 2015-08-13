# --------------------------------------------------------------------------
# File: _ostream.py 
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

import weakref
from ._procedural import check_status
from ..exceptions import CplexError, CplexSolverError
from .. import six

class OutputStream(object):

    """Class to parse and write strings to a file object.

    Methods:
    __init__(self, outputfile, fn = None)
    __del__(self)
    write(self)
    flush(self)
    
    """
    
    def __init__(self, outputfile, env, fn = None):
        """OutputStream constructor.

        outputfile must provide methods write(self, str) and
        flush(self).

        If fn is specified, it must be a fuction with signature
        fn(str) -> str.

        """
        self._env = weakref.proxy(env)
        self._fn = fn
        self._is_valid = False
        if isinstance(outputfile, six.string_types):
            self._file = open(outputfile, "w")
        else:
            self._file = outputfile
        if self._file is not None:
            if not hasattr(self._file, "write"):
                raise CplexError("Output object must have write method")
            elif not callable(self._file.write):
                raise CplexError("Output object must have write method")
            if not hasattr(self._file, "flush"):
                raise CplexError("Output object must have flush method")
            elif not callable(self._file.flush):
                raise CplexError("Output object must have flush method")
        self._is_valid = True

    def __del__(self):
        """OutputStream destructor."""
        # If something bad happened in the constructor, then don't flush.
        if self._is_valid:
            self.flush()

    def _write_wrap(self, str_):
        try:
            self._terminate = 0
            self.write(str_)
            self.flush()
            if hasattr(self, "_error_string"):
                msg = self._error_string
                if msg is not None:
                    if not msg.startswith("CPLEX Error  1006"):
                        self._error_string = None
                        raise CplexError("ERROR", msg)
        except Exception as exc:
            self._env._callback_exception = exc
            check_status._pyenv = self._env
            self._terminate = 1
            

    def write(self, str_):
        """Parses and writes a string.

        If self._fn is not None, self._fn(str_) is passed to
        self._file.write.  Otherwise, str_ is passed to self._file.write

        """
        if self._file is None:
            return
        if str_ is None:
            str_ = ""
        if self._fn is None:
            self._file.write(str_)
        else:
            self._file.write(self._fn(str_))

    def flush(self):
        """Flushes the buffer."""
        if self._file is not None:
            self._file.flush()

        
