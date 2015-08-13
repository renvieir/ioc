# --------------------------------------------------------------------------
# File: _matrices.py 
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

from ..exceptions import CplexError
from . import _procedural
from .. import six

class SparsePair(object):

    """A class for storing sparse vector data.

    An instance of this class has two attributes, ind and val.  ind
    specifies the indices and val specifies the values.  ind and val
    must be sequences of the same length.  In general, ind may contain
    any identifier; for example, when a SparsePair object is passed to
    Cplex.linear_constraints.add, its ind attribute may be a list
    containing both variable names and variable indices.
    
    """
    
    def __init__(self, ind = [], val = []):
        """Constructor for SparsePair.

        Takes two arguments, ind and val; ind specifies the indices that
        the SparsePair refers to, and val specifies the float values 
        associated with those indices; ind and val must have the same
        length.

        """
        self.ind = ind
        self.val = val
        if not self.isvalid():
            raise CplexError("Inconsistent input data to SparsePair")

    def __repr__(self):
        """Representation method of SparsePair."""
        return "".join(["SparsePair(ind = ",
                        repr(self.ind),
                        ", val = ",
                        repr(self.val), ")"])

    def isvalid(self):
        """Tests that ind and val have the same length."""
        return len(self.ind) == len(self.val)


class _C_HBMatrix(object):

    """non-public


    """

    def __init__(self, lolmat, env_lp, r_c, enc):
        """non-public"""
        self._mat = _procedural.Pylolmat_to_CHBmat(lolmat, env_lp, r_c, enc)
        if self._mat is None:
            raise MemoryError
        elif len(self._mat) == 2:
            raise CplexError(" %d: Invalid name -- '%s'\n" % tuple(self._mat))
        elif len(self._mat) == 1:
            raise TypeError(" invalid matrix input type -- ", self._mat[0])
        
    def __del__(self):
        """non-public"""
        if hasattr(self, "_mat") and self._mat is not None:
            if len(self._mat) == 4:
                _procedural.free_CHBmat(self._mat)
    
    def _get_nnz(self):
        """non-public"""
        return self._mat[3]


class _HBMatrix(object):

    """non-public


    """

    def __init__(self, matrix = None):
        """non-public"""
        self.matbeg = []
        self.matind = []
        self.matval = []
        if matrix is not None:
            for vector in matrix:
                if isinstance(vector, SparsePair):
                    v0 = vector.ind
                    v1 = vector.val
                else:
                    v0 = vector[0]
                    v1 = vector[1]
                if len(v0) != len(v1):
                    raise CplexError("Inconsistent input data to _HBMatrix")
                self.matbeg.append(len(self.matind))
                self.matind.extend(v0)
                self.matval.extend(v1)
            
    def __len__(self):
        """non-public"""
        return len(self.matbeg)

    def __getitem__(self, key):
        """non-public"""
        if isinstance(key, six.integer_types):
            if key < 0:
                key += len(self)
            begin = self.matbeg[key]
            if key == len(self) - 1:
                end = len(self.matind)
            else:
                end = self.matbeg[key + 1]
            return SparsePair(self.matind[begin:end], self.matval[begin:end])
        elif isinstance(key, slice):
            start, stop, step = key.start, key.stop, key.step
            if start is None:
                start = 0
            if stop is None or stop > len(self):
                stop = len(self)
            if step is None:
                step = 1
            return [self[i] for i in range(start, stop, step)]
        else:
            raise TypeError


class SparseTriple(object):

    """A class for storing sparse matrix data.

    An instance of this class has three attributes, ind1, ind2, and
    val.  ind1 and ind2 specify the indices and val specifies the
    values.  ind1, ind2, and val must be sequences of the same length.
    In general, ind1 and ind2 may contain any identifier; for example, when a
    SparseTriple object is passed to Cplex.quadratic_constraints.add,
    its ind1 attribute may be a list containing both variable names
    and variable indices.

    """

    def __init__(self, ind1 = [], ind2 = [], val = []):
        """Constructor for SparseTriple.

        Takes three arguments, ind1, ind2 and val, specifying the
        indices that the SparseTriple refers to and the float values
        associated with those indices, respectively.  ind1, ind2, and
        val must all have the same length.

        """
        self.ind1 = ind1
        self.ind2 = ind2
        self.val  = val
        if not self.isvalid():
            raise CplexError("Inconsistent input data to SparseTriple")

    def __repr__(self):
        """Representation method of SparseTriple."""
        return "".join(["SparseTriple(ind1 = ",
                        repr(self.ind1),
                        ", ind2 = ",
                        repr(self.ind2),
                        ", val = ",
                        repr(self.val), ")"])

    def isvalid(self):
        """Tests that ind1, ind2, and val have the same length."""
        return (len(self.ind1) == len(self.ind2) and
                len(self.ind1) == len(self.val))
