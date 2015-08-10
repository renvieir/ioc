def initial_random_tour():
    pass


def escolha_x(param, param1):   # Step 6
                                # (a) if t2i is joined to t1 , the resulting configuration is a tour, T', and
                                # (b) x i != y s for all s < i
    pass


def escolha_y(param, param1):
    pass


def escolha_t(param):
    pass


def lin_kernighan(G):
    T = []
    goto = 0
    while True:
        T = initial_random_tour()       # Step 1

        while True:                     # Step 2
            i = 0
            t1 = escolha_t(1)

            while True:                 # Step 3
                x_1 = escolha_x(1,2)

                while True:             # Step 4
                    y_1 = escolha_y(2,3)

                    while True:         # Step 5
                        i = i+1

                        while True:     # Step 6
                            escolha_x(t[2*i -1],t[2*i])

                            if T_linha < T:
                                goto = 2
                                break                   # Goto Step 2

                            while True: # Step 7
                                y = escolha_y(t[2*i],t[2*i+1])

                                if y:           # if such y exists goto Step 5
                                    goto = 5
                                    break
                                if untried_y(2):    # Step 8
                                    i = 2
                                    continue            # Goto Step 7

                                if untried_x(2):    # Step 9
                                    i=2
                                    goto = 6
                                    break

                                if untried_y(1):    # Step 10
                                    i = 1
                                    goto = 4
                                    break

                                if untried_x(1):    # Step 11
                                    i = 1
                                    goto = 3
                                    break

                                if untried_t(1):    # Step 12
                                    goto = 2
                                    break

                                if stop:            # Step 13
                                    return T

                                goto=1
                                break

                            if goto <= 5:
                                break
                        if goto <= 4:
                            break
                    if goto <= 3:
                        break
                if goto <= 2:
                    break
            if goto <= 1:
                break

        T = []                  # 1. Generate a random starting tour T
        i = 0

        best_gain = 0           # 2. set G*=0 [G* is the best improvement made so far]
        t1 = None               # choose any node t1
        x1 = None               # and let x1 be one of the edges of T adjacents to t1.
        i = 1                   # Let i = 1

        y1 = (t2, t3)           # 3. From the other endpoint of x1 (t2) choose y1=(t2,t3) with g1 > 0.
        if y1:                  # If no such y1 exists go to step 6(d)

            i = i +1                            # 4. Let i = i + 1
            x[i] = (t2*i -1, t2*i)              # (a) Choose xi
            y[i] = choose_y(i)                  # (b) Choose yi, some available link at endpoint t2i shared with xi subject to (c), (d) and (e)

            # (c) x[i] cannot be a link previously joined (a y[j], j<i)
            # and y[i] cannot be a link previously borken

            # (d) G[i] = sum_from_i_equal_1_to_j(g_i > 0) [Gain criterion]

            # (e) in order to ensure that feasability criterion of (a) can be satisfied at y+1,
            # the y_i chosen mus permit the breaking os and x_i+1

            # (f) before U_i is constructed, whe check is closing up by joining

    return T