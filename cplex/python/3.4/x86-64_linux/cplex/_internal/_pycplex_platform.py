# ------------------------------------------------------------------------------
# Licensed Materials - Property of IBM
# 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
# Copyright IBM Corporation 2008, 2014. All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with
# IBM Corp.
# ------------------------------------------------------------------------------
"""
"""

import platform

from sys import version_info

ERROR_STRING = "CPLEX 12.6.1.0 is not compatible with this version of Python."

if platform.system() in ('Darwin', 'Linux', 'AIX', 'Windows', 'Microsoft'):
    if version_info < (3, 4, 0):
        raise Exception(ERROR_STRING)
    elif version_info < (3, 5, 0):
        from cplex._internal.py34_cplex1261 import *
    else:
        raise Exception(ERROR_STRING)
else:
    raise Exception("The CPLEX Python API is not supported on this platform.")
