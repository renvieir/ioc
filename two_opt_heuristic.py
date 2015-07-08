def _cost_calc(cost_matrix, path):
    cost = 0
    for i in range(len(path)-1):
        value = cost_matrix[path[i]][path[i+1]]
        # Se tiver um caminho que nao existe entao custo eh infinito
        if not value:
            return None
        cost += value
    value = cost_matrix[path[i]][path[i+1]]
    cost += value
    return cost

def two_opt_swap(route, i, k):
    new_route = route[:i-1]
    r = route[i-1:k]
    r = r[::-1]
    new_route += r
    new_route += route[k:]

    return new_route

def two_opt_heuristic(cost_matrix, route, max):
    improvement = 0
    while improvement < max:
        best_distance = _cost_calc(cost_matrix, route)
        for i, k in enumerate(range(1,len(route))):
            new_route = two_opt_swap(route,i,k)
            new_distance = _cost_calc(cost_matrix,new_route)

            if new_distance != None and new_distance < best_distance:
                route = new_route
                best_distance = new_distance
                improvement = 0

        improvement+=1
    return route
