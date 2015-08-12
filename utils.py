def read_tsplib(filename):
    "basic function for reading a TSP problem on the TSPLIB format"
    "NOTE: only works for explicit values in upper row matrix format"
    f = open(filename, 'r');

    line = f.readline()
    while line.find("EDGE_WEIGHT_SECTION") == -1:
        line = f.readline()

    matrix=[]
    while 1:
        line = f.readline()
        if line.find("EOF") != -1: break

        values=line.split(' ')



        matrix.append(values)
    return

