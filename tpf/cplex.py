import math


class TSPReader:
    def __init__(self, filename):
        self.info = {}
        self.filename = filename
        self.dist = None

    def parse(self):
        self.__read_file()
        self.__calculate_distance_matrix()

    def __read_file(self):
        "basic function for reading a TSP problem on the TSPLIB format"
        "NOTE: only works for 2D euclidean, manhattan distances or Geographical distance"
        f = open(self.filename, 'r')

        line = f.readline()
        while line.find("EDGE_WEIGHT_TYPE") == -1:
            line = f.readline()

        if line.find("EUC_2D") != -1:
            self.dist = self.__euclidean_distance
        elif line.find("MAN_2D") != -1:
            self.dist = self.__manhattan_distance
        elif line.find("GEO") != -1:
            self.dist = self.__geographical_distance
        else:
            print "cannot deal with non-euclidean or non-manhattan distances"
            raise Exception

        while line.find("NODE_COORD_SECTION") == -1:
            line = f.readline()

        xy_positions = []
        while 1:
            line = f.readline()
            if line.find("EOF") != -1: break
            (i, x, y) = line.split()
            x = float(x)
            y = float(y)
            xy_positions.append((x, y))

        self.info['coordinates'] = xy_positions

    def __calculate_distance_matrix(self):
        """Compute a distance matrix for a set of points.

        Uses function 'dist' to calculate distance between
        any two points.  Parameters:
        -coord -- list of tuples with coordinates of all points, [(x1,y1),...,(xn,yn)]
        -dist -- distance function
        """
        coord = self.info['coordinates']
        n = len(coord)
        D = {}  # dictionary to hold n times n matrix
        for i in range(n - 1):
            for j in range(i + 1, n):
                (x1, y1) = coord[i]
                (x2, y2) = coord[j]
                D[i, j] = self.dist((x1, y1), (x2, y2))
                D[j, i] = D[i, j]
        # return n,D
        self.info['distance_matrix'] = D

    def __euclidean_distance(self, (x1, y1), (x2, y2)):
        """Compute the L2-norm (Euclidean) distance between two points.

        The two points are located on coordinates (x1,y1) and (x2,y2),
        sent as parameters"""
        xdiff = x2 - x1
        ydiff = y2 - y1
        return int(math.sqrt(xdiff * xdiff + ydiff * ydiff))

    def __manhattan_distance(self, (x1, y1), (x2, y2)):
        """Compute the L1-norm (Manhattan) distance between two points.

        The two points are located on coordinates (x1,y1) and (x2,y2),
        sent as parameters"""
        return int(abs(x2 - x1) + abs(y2 - y1))

    def __geographical_distance(self, (x1, y1), (x2, y2)):
        """Compute the Geographical distance between latitude, longitude points.

        The two points are located on coordinates (x1,y1) and (x2,y2),
        sent as parameters"""
        def radian_coordinate(position):
            pi = 3.141592
            deg = int(position)
            min = position - deg
            radian_position = pi * (float(deg) + (5.0 * min)) / 180.0
            return radian_position

        # The radius of Earth in kilometers
        rrr = 6378.388

        q1 = math.cos(radian_coordinate(y1) - radian_coordinate(y2))
        q2 = math.cos(radian_coordinate(x1) - radian_coordinate(x2))
        q3 = math.cos(radian_coordinate(x1) + radian_coordinate(x2))

        geo_distance = int(rrr * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0)

        return geo_distance


if __name__ == '__main__':
    print 'starting cplex solving'

    # tsplib_file = 'data/a280.tsp'
    # tsplib_file = 'data/berlin52.tsp'
    tsplib_file = 'data/burma14.tsp'

    reader = TSPReader(tsplib_file)
    reader.parse()

    for items in reader.info.items():
        print '\n', items

        # TODO: DONE! 1. read tsp data
        # TODO: 2. parse data to cplex format
        # TODO: 3. solve first time
        # TODO: 4. do branch and bound while not integer
        # TODO: 5. do branch and cut while not integer
