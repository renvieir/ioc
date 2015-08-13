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


"""Exceptions raised by the CPLEX Python API.

   For documentation of CPLEX error codes, see the group
   optim.cplex.errorcodes in the reference 
   manual of the CPLEX Callable Library, and the topic 
   Interpreting Error Codes in the Overview of the APIs.

"""


from . import error_codes


__all__ = ["CplexError", "CplexSolverError", "error_codes"]


class CplexError(Exception):

    """Class for exceptions raised by the CPLEX Python API."""

    pass


class CplexSolverError(CplexError):

    """Class for errors returned by the Callable Library functions.

    self.args[0] : A string describing the error.

    self.args[1] : The address of the environment that raised the error.

    self.args[2] : The integer status code of the error.

    """

    def __str__(self):
        return self.args[0]
    

