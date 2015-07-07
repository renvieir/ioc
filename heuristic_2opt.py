# 2optSwap(route, i, k) {
#       1. take route[1] to route[i-1] and add them in order to new_route
#       2. take route[i] to route[k] and add them in reverse order to new_route
#       3. take route[k+1] to end and add them in order to new_route
#       return new_route;
#   }

# repeat until no improvement is made {
#       start_again:
#       best_distance = calculateTotalDistance(existing_route)
#       for (i = 0; i < number of nodes eligible to be swapped - 1; i++) {
#           for (k = i + 1; k < number of nodes eligible to be swapped; k++) {
#               new_route = 2optSwap(existing_route, i, k)
#               new_distance = calculateTotalDistance(new_route)
#               if (new_distance < best_distance) {
#                   existing_route = new_route
#                   goto start_again
#               }
#           }
#       }
#   }
