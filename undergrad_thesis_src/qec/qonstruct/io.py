"""
    author: Suhas Vittal
    date:   2 January 2024
"""

from qonstruct.code_builder.base import *

import networkx as nx

def read_tanner_graph_file(input_file: str) -> nx.Graph:
    """
        This reads a Tanner graph from an input file. Each line should be
        as follows:
            <O?><X/Z><number>,<data-qubit>,<data-qubit>,...
        If a line begins with an 'O', it is considered an observable declaration.
        Otherwise, it is a check declaration. Example (of not a real code):

            Z0,0,1,2
            Z1,2,3,4
            X0,1,2,3
            X2,0,2,4
            OZ0,0,1,3
            OX0,0,2,3
    """
    reader = open(input_file, 'r')
    lines = reader.readlines()
    reader.close()

    gr = tanner_init()

    n = 0
    check_list_map = {'x': [], 'z': []}
    for ln in lines:
        line_data = ln.split(',')
        ch_s = line_data[0]
        # Get type:
        support = [int(x) for x in line_data[1:]]
        for x in support:
            if not gr.has_node(x):
                add_data_qubit(gr, x)
                n = max(n, x)
        if ch_s[0] == 'O':
            # This is an observable.
            node_type = 'x' if ch_s[1] == 'X' or ch_s[1] == 'x' else 'z'
            add_observable(gr, support, node_type)
        else:
            node_type = 'x' if ch_s[0] == 'X' or ch_s[0] == 'x' else 'z'
            check_list_map[node_type].append(support)
    n = n+1
    for s in ['x', 'z']:
        for support in check_list_map[s]:
            add_check(gr, n, s, support)
            n += 1

    return gr

def write_tanner_graph_file(tanner_graph: nx.Graph, output_file: str) -> None:
    """
        This code will write a Tanner graph to an output file.
    """
    writer = open(output_file, 'w')
    # Map data qubit nodes to indices.
    x_ctr, z_ctr = 0, 0
    for s in tanner_graph.graph['checks']['all']:
        if get_check_property(tanner_graph, s, 'node_type') == 'x':
            writer.write('X%d' % x_ctr)
            x_ctr += 1
        else:
            writer.write('Z%d' % z_ctr)
            z_ctr += 1
        for d in tanner_graph.neighbors(s):
            writer.write(',%d' % tanner_graph.nodes[d]['label'])
        writer.write('\n')
    # Write logical observables of the code.
    x_obs_list, z_obs_list = tanner_graph.graph['obs_list']['x'], tanner_graph.graph['obs_list']['z']
    for (i, obs) in enumerate(x_obs_list):
        writer.write('OX%d' % i)
        for d in obs:
            writer.write(',%d' % d)
        writer.write('\n')
    for (i, obs) in enumerate(z_obs_list):
        writer.write('OZ%d' % i)
        for d in obs:
            writer.write(',%d' % d)
        writer.write('\n')
    writer.close()

