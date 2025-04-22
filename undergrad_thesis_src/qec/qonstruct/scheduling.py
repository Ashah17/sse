"""
    author: Suhas Vittal
    date:   12 January 2024

    This file contains code for algorithmically computing syndrome extraction schedules.
    Requires cplex.
"""

from qonstruct.code_builder.base import *

import networkx as nx

import docplex.mp
from docplex.mp.model import Model

def compute_syndrome_extraction_schedule(tanner_graph: nx.Graph):
    """
        Modifies the tanner_graph by setting the scheduling order for
        each check vertex. This is done greedily with repeated calls
        to CPLEX.
    """
    M = 100000

    checks = tanner_graph.graph['checks']['all']
    max_check_weight = max(tanner_graph.degree(x) for x in checks)
    for (i, ch) in enumerate(checks):
        support = get_support(tanner_graph, ch)
        # Develop the program. The objective of the program is to minimize the depth of the circuit.
        prog = Model(name='ext')

        max_of_all = prog.integer_var(name='max_of_all') # This is the depth of the circuit.

        qu_var_map = {}
        # Create a variable for each data qubit in the support.
        for q in support:
            v = prog.integer_var(name=str(q))
            qu_var_map[q] = v
            # Add constraints for q.
            prog.add_constraint(v >= 1)
            prog.add_constraint(v <= 2*max_check_weight)
            prog.add_constraint(max_of_all >= v)
        # Add uniqueness constraints.
        for (ii, q1) in enumerate(support):
            v1 = qu_var_map[q1]
            for jj in range(ii+1, len(support)):
                q2 = support[jj]
                v2 = qu_var_map[q2]
                # Add the constraint.
                prog.add_constraint(v1 != v2)
        # Add commutativity and timing constraints.
        for j in range(0, i):
            other_ch = checks[j]
            schedule_order = tanner_graph.nodes[other_ch]['schedule_order']

            ind_sum_array = []
            for (t, q) in enumerate(schedule_order):
                if q is None or q not in support:
                    continue
                v = qu_var_map[q]
                prog.add_constraint(v != t) 
                # If ch and other_ch are different types, add commutativity constraints.
                if tanner_graph.nodes[ch]['node_type'] == tanner_graph.nodes[other_ch]['node_type']:
                    continue
                y = prog.binary_var()
                prog.add_constraint(v - t <= M*y)
                prog.add_constraint(t - v <= M*(1-y))
                ind_sum_array.append(y)
            # Require that the ind_sum_array is even.
            y = prog.integer_var()
            prog.add_constraint(y >= 0)
            prog.add_constraint(prog.sum(ind_sum_array) == 2*y)
        prog.minimize(max_of_all)
        soln = prog.solve()
        depth = int(round(soln.get_var_value(max_of_all)))

        schedule_order = [None for _ in range(depth)]
        for q in support:
            v = qu_var_map[q]
            t = int(round(soln.get_var_value(v)))
            # As t is 1-indexed, the actual time is t-1.
            schedule_order[t-1] = q
        tanner_graph.nodes[ch]['schedule_order'] = schedule_order

