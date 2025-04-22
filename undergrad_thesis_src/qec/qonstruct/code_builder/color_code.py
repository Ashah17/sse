"""
    author: Suhas Vittal
    date:   1 December 2023

    Color code tanner graphs.
"""

from qonstruct.code_builder.base import *

import networkx as nx

def make_hexagonal(d: int, both_at_once=False) -> nx.Graph:
    """
        This function will return a tanner graph for a hexagonal color
        code (weight-6 RGB plaquettes).
    """
    gr = tanner_init()

    # Some geometric properties of the code:
    offset, side_len = 2, (3*d - 1)//2
    # In our initial pass, we will go through the structure, and identify
    # each qubit.
    loc_map, check_locs, obs = {}, [], []
    n = 0
    for r in range(side_len):
        row_off = offset
        for c in range(r+1):
            if row_off == 0:
                check_locs.append((r, c))
            else:
                loc_map[(r, c)] = n
                add_data_qubit(gr, n)
                if c == 0:
                    # This is the left edge of the triangle, which is a
                    # X and Z operator.
                    obs.append(n)
                n += 1
            row_off = (row_off+1) % 3
        offset = (offset+1) % 3
    add_observable(gr, obs, 'x')
    add_observable(gr, obs, 'z')
    # Now, create all the checks and complete the Tanner graph.
    get_loc = lambda _r, _c: loc_map[(_r, _c)] if (_r, _c) in loc_map else None

    plaq = 0
    for (i, j) in check_locs:
        a = get_loc(i-1, j-1)
        b = get_loc(i-1, j)
        c = get_loc(i, j+1)
        d = get_loc(i+1, j+1)
        e = get_loc(i+1, j)
        f = get_loc(i, j-1)

        add_plaquette(gr, plaq, [a, b, c, d, e, f], i%3)
        for stabilizer in ['x', 'z']:
            # Create check vertex.
            if both_at_once:
                if stabilizer == 'z':
                    schedule_order = [b, c, d, a, f, e, None]
                else:
                    schedule_order = [None, b, a, f, c, d, e]
            else:
                schedule_order = [b, c, e, d, a, f]
            add_check(gr, n, stabilizer, [a, b, c, d, e, f],
                    color=i%3, plaquette=plaq, schedule_order=schedule_order)
            set_plaquette(gr, n, plaq)
            n += 1
        plaq += 1
    return gr

def make_hycc_d4(both_at_once=False) -> nx.Graph:
    gr = tanner_init()

    red_operators = [
        [0, 1, 2, 3, 4, 5],
        [6, 7, 8, 9, 10, 11],
        [12, 13, 14, 15, 16, 17],
        [18, 19, 20, 21, 22, 23]
    ]

    blue_operators = [
        [2, 4, 8, 10, 21, 19, 13, 15],
        [18, 20, 14, 12, 3, 5, 9, 11]
    ]

    green_operators = [
        [0, 2, 6, 8, 23, 21, 15, 17],
        [20, 22, 14, 16, 3, 1, 7, 9]
    ]

    z_obs_list = [
        [11, 10, 4, 5],
        [4, 5, 12, 13],
        [7, 6, 0, 1],
        [0, 1, 16, 17],
        [9, 11, 18, 20],
        [8, 10, 4, 2],
        [22, 20, 9, 7],
        [6, 8, 2, 0]
    ]

    x_obs_list = [
        [11, 10, 4, 5],
        [4, 5, 12, 13],
        [7, 6, 0, 1],
        [0, 1, 16, 17],
        [9, 11, 18, 20],
        [8, 10, 4, 2],
        [22, 20, 9, 7],
        [6, 8, 2, 0]
    ]
    # Add all data qubits to the graph.
    for i in range(24):
        add_data_qubit(gr, i)
    # Add observables:
    for i in range(8):
        add_observable(gr, x_obs_list[i], 'x')
        add_observable(gr, z_obs_list[i], 'z')
    # Create checks:
    n, plaq = 0, 0
    for (c, check_list) in enumerate([red_operators, green_operators, blue_operators]):
        for support in check_list:
            add_plaquette(gr, plaq, support, c)
            for stabilizer in ['x', 'z']:
                add_check(gr, n, stabilizer, support, color=c, plaquette=plaq, schedule_order=[])
                set_plaquette(gr, n, plaq)
                n += 1
            plaq += 1
    # Now, we must determine the schedules. We will do so algorithmically.
    from qonstruct.scheduling import compute_syndrome_extraction_schedule
    compute_syndrome_extraction_schedule(gr)
    return gr
