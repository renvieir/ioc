from utils import *

def test_default(argv):
    """Local search for the Travelling Saleman Problem: sample usage."""
    #
    # test the functions:
    #

    # random.seed(1)    # uncomment for having always the same behavior

    if len(argv) == 1:
        # create a graph with several cities' coordinates
        coord = [(4,0),(5,6),(8,3),(4,4),(4,1),(4,10),(4,7),(6,8),(8,1)]
        n, D = mk_matrix(coord, distL2) # create the distance matrix
        instance = "toy problem"
    else:
        instance = argv[1]
        n, coord, D = read_tsplib(instance)     # create the distance matrix
        # n, coord, D = read_tsplib('INSTANCES/TSP/eil51.tsp')  # create the distance matrix

    # function for printing best found solution when it is found
    from time import clock
    init = clock()
    def report_sol(obj, s=""):
        print "cpu:%g\tobj:%g\ttour:%s" % \
              (clock(), obj, s)


    print "*** travelling salesman problem ***"
    print

    # random construction
    print "random construction + local search:"
    tour = randtour(n)     # create a random tour
    z = length(tour, D)     # calculate its length
    print "random:", tour, z, '  -->  ',
    z = localsearch(tour, z, D)      # local search starting from the random tour
    print tour, z
    print

    # greedy construction
    print "greedy construction with nearest neighbor + local search:"
    for i in range(n):
        tour = nearest_neighbor(n, i, D)     # create a greedy tour, visiting city 'i' first
        z = length(tour, D)
        print "nneigh:", tour, z, '  -->  ',
        z = localsearch(tour, z, D)
        print tour, z
    print

    # multi-start local search
    print "random start local search:"
    niter = 100
    tour,z = multistart_localsearch(niter, n, D, report_sol)
    assert z == length(tour, D)
    print "best found solution (%d iterations): z = %g" % (niter, z)
    print tour

if __name__ == "__main__":
    import sys
    test_default(sys.argv)

