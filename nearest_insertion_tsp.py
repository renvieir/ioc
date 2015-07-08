def nearest(cost_matrix, node=0, path=[]):
    neighboors = cost_matrix[node]
    neighboors_index = sorted(range(len(neighboors)), key=lambda k: neighboors[k])

    if len(path) == len(cost_matrix):
        if path[0] in neighboors_index:
            return path
        else:
            return None
    else:
        for index in neighboors_index:
            vertex_value = neighboors[index]
            if vertex_value != None and index not in path:
                new_path = nearest(cost_matrix,index,path+[index])
                if new_path != None:
                    return new_path
        return None

def cost_calc(cost_matrix, path):
    cost = 0
    for i in range(len(path)-1):
        value = cost_matrix[path[i]][path[i+1]]
        print cost,'+',value,'=',cost+value
        cost += value
    value = cost_matrix[path[i]][path[i+1]]
    cost += value
    print cost,'+',value,'=',cost+value
    return cost