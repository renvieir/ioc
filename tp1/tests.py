#!/usr/bin/python

import sys
from nearest_insertion_tsp import *
from two_opt_heuristic import *

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print 'exemplo de uso: python tests.py start_node improvement'
    else:
        node = int(sys.argv[-2])
        '''
        cost_matrix = [[None,9,7,None,15,10],
                       [9,None,2,None,8,12],
                       [7,2,None,6,None,None],
                       [None,None,6,None,3,None],
                       [15,8,None,3,None,24],
                       [10,12,None,None,24,None]]
        '''
        cost_matrix = [[None,6,100,100,8,5],
                       [6,None,9,100,100,20],
                       [100,9,None,7,5,20],
                       [100,100,7,None,3,10],
                       [8,100,5,3,None,12],
                       [5,20,20,10,12,None]]

        path = nearest(cost_matrix,node,[node])
        print 'path = ', path
        if path != None:
            print 'cost = ', cost_calc(cost_matrix,path)

    improvement = int(sys.argv[-1])
    print two_opt_heuristic(cost_matrix, path,improvement)