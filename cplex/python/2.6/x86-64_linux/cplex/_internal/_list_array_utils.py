# --------------------------------------------------------------------------
# File: _list_array_utils.py 
# ---------------------------------------------------------------------------
# Licensed Materials - Property of IBM
# 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
# Copyright IBM Corporation 2008, 2014. All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with
# IBM Corp.
# ------------------------------------------------------------------------

from . import _pycplex as CPX

# int_list_to_C_array    = CPX.int_list_to_C_array
# double_list_to_C_array = CPX.double_list_to_C_array

def int_list_to_array(input):
    length = len(input)
    if length == 0:
        return CPX.cvar.CPX_NULL
    output = CPX.intArray(length)
    for i in range(length):
        output[i] = input[i]
    return output

def int_list_to_array_trunc_int32(input):
    int32_min = -2147483648
    int32_max =  2147483647
    length = len(input)
    if length == 0:
        return CPX.cvar.CPX_NULL
    output = CPX.intArray(length)
    for i in range(length):
        if input[i] > int32_max:
            output[i] = int32_max
        elif input[i] < int32_min:
            output[i] = int32_min
        else:
            output[i] = input[i]
    return output

def double_list_to_array(input):
    length = len(input)
    if length == 0:
        return CPX.cvar.CPX_NULL
    output = CPX.doubleArray(length)
    for i in range(length):
        output[i] = input[i]
    return output

def int_array_to_list(input, length):
    output = []
    for i in range(length):
        output.append(input[i])
    return output

def double_array_to_list(input, length):
    output = []
    for i in range(length):
        output.append(input[i])
    return output


class int_C_array(object):

    def __init__(self, list_):
        self.array = CPX.int_list_to_C_array(list_)

    def __del__(self):
        CPX.free_int_C_array(self.array)


class double_C_array(object):

    def __init__(self, list_):
        self.array = CPX.double_list_to_C_array(list_)
        if self.array == "error":
            del self.array
            raise TypeError("Non-float value in float input list")

    def __del__(self):
        if hasattr(self, "array"):
            CPX.free_double_C_array(self.array)


