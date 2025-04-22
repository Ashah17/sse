"""
    author: Suhas Vittal
    date:   1 December 2023

    Surface code tanner graphs.
"""

from qonstruct.code_builder.base import *

import networkx as nx

def make_rotated(d: int) -> nx.Graph:
    """
        This function will return a tanner graph for a rotated surface
        code.
    """
    gr = tanner_init()
    # Add data qubits and observables.
    n = 0
    for _ in range(d*d):
        add_data_qubit(gr, n)
        n += 1
    add_observable(gr, [x for x in range(d)], 'x')
    add_observable(gr, [(d-1)+x*d for x in range(d)], 'z')
    # Add checks. Do the boundary checks first, then the bulk.
    # Z checks:
    for r in range(0, d-1, 2):
        q1, q2 = r, r+1
        add_check(gr, n, 'z', [q1, q2], schedule_order=[q2, q1, None, None])
        q1, q2 = d*(d-1) + r+1, d*(d-1) + r+2
        add_check(gr, n+1, 'z', [q1, q2], schedule_order=[None, None, q2, q1])
        n += 2
    for c in range(0, d-1, 2):
        q1, q2 = c*d + d-1, (c+1)*d + d-1
        add_check(gr, n, 'x', [q1, q2], schedule_order=[None, None, q2, q1])
        q1, q2 = (c+1)*d, (c+2)*d
        add_check(gr, n+1, 'x', [q1, q2], schedule_order=[q2, q1, None, None])
        n += 2
    # Now do bulk:
    for r in range(0, d-1):
        for c in range(0, d-1):
            x, y, z, w = r + c*d, r + (c+1)*d, r+1 + c*d, r+1 + (c+1)*d
            if (r+c) % 2 == 1:
                add_check(gr, n, 'z', [x, y, z, w], schedule_order=[w, y, z, x])
            else:
                add_check(gr, n, 'x', [x, y, z, w], schedule_order=[w, z, y, x])
            n += 1
    return gr
