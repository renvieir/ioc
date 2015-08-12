'''
    Renato Pacheco Vieira 2150283
    Referencias:
        http://www2.isye.gatech.edu/~mgoetsch/cali/VEHICLE/TSP/TSP009__.HTM
        https://en.wikipedia.org/wiki/2-opt
'''

def cost_calc(cost_matrix, path):
    cost = 0
    for i in range(len(path)-1):
        value = cost_matrix[path[i]][path[i+1]]
        cost += value
    value = cost_matrix[path[i]][path[i+1]]
    cost += value
    return cost


def iterative_nearest_insertion(cost_matrix, node):
    unvisited = range(len(cost_matrix))
    unvisited.remove(node)
    last = node

    tour = [node]

    while unvisited != []:
        next = nearest(cost_matrix, last, unvisited)
        tour.append(next)
        unvisited.remove(next)
        last = next
    return tour


def nearest(cost_matrix, node, unvisited):
    neighboors = cost_matrix[node]
    neighboors_index = sorted(range(len(neighboors)), key=lambda k: neighboors[k])

    for index in neighboors_index:
        vertex_value = neighboors[index]
        if vertex_value != None and index in unvisited:
            return index
    return None