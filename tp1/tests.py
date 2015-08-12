#!/usr/bin/python
'''
    Renato Pacheco Vieira 2150283
    Referencias:
        http://www2.isye.gatech.edu/~mgoetsch/cali/VEHICLE/TSP/TSP009__.HTM
        https://en.wikipedia.org/wiki/2-opt
'''

import sys
from nearest_insertion_tsp import *
from two_opt_heuristic import *

costs = []
costs.append([[None,9,7,None,15,10],
               [9,None,2,None,8,12],
               [7,2,None,6,None,None],
               [None,None,6,None,3,None],
               [15,8,None,3,None,24],
               [10,12,None,None,24,None]])
costs.append([[None,6,100,100,8,5],
               [6,None,9,100,100,20],
               [100,9,None,7,5,20],
               [100,100,7,None,3,10],
               [8,100,5,3,None,12],
               [5,20,20,10,12,None]])

if __name__ == "__main__":
    print '             example of use: python tests.py <matrix> <start_node> <improvement>'
    print ''
    if len(sys.argv) != 4:
        print len(sys.argv),'arguments received.'
        print '4 arguments expected.'
    else:

        m = int(sys.argv[-3])

        if not m < len(costs):
            print '<matrix> must be less then ',len(costs)

        else:
            node = int(sys.argv[-2])

            if not node < len(costs[m]):
                print '<start_node> must be less then ',len(costs[m])

            else:


                # Escolher abaixo entre 0 e 1 para mudar a matriz usada
                matrix = costs[1]

                print 'Nearest Insertion Solution'
                # path = nearest_insertion(cost_matrix,node,[node])
                path = iterative_nearest_insertion(matrix,node)
                print 'path = ', path
                if path != None:
                    print 'cost = ', cost_calc(matrix,path)
                else:
                    print 'No Solution'

                print ''
                improvement = int(sys.argv[-1])
                path = two_opt_heuristic(matrix, path,improvement)
                print 'Improved Solution after ',improvement,'trys'
                print 'path = ', path
                print 'cost = ', cost_calc(matrix,path)
                print ''