#!/usr/bin/python

import sys

def nearest(cost_matrix, node=0, path=[]):
    if len(path) == len(cost_matrix):
        if path[0] in cost_matrix[node]:
            return path
        else:
            return None
    else:
        new_path = None
        s = cost_matrix[node]
        neighboors = sorted(range(len(s)), key=lambda k: s[k])
        for v in neighboors:
            vertex = s[v]
            if vertex != None and v not in path:
                new_path = nearest(cost_matrix,v,path+[v])
                if new_path!= None:
                    return new_path
        return None

def cost_calc(cost_matrix, path):
    cost = 0
    for i in range(len(path)-1):
        cost += cost_matrix[path[i]][path[i+1]]
    return cost

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print 'exemplo de uso: python nearest_insertion.py 0'
    else:
        node = int(sys.argv[-1])
        cost_matrix = [[None,9,7,None,15,10],[9,None,2,None,8,12],[7,2,None,6,None,None],[None,None,6,None,3,None],[15,8,None,3,None,24],[10,12,None,None,24,None],]
        path = nearest(cost_matrix,node)
        print 'path = ', path
        if path != None:
            print 'cost = ', cost_calc(cost_matrix,path)
