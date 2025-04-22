"""
    author: Suhas Vittal
    date:   21 January 2024
"""

from qonstruct.code_builder.base import *

import networkx as nx
import numpy as np

def make_support_graph(tanner_graph: nx.Graph, op_type: str) -> nx.Graph:
    gr = nx.Graph()
    for x in tanner_graph.nodes():
        if tanner_graph.nodes[x]['node_type'] == op_type:
            gr.add_node(x)
    for x in gr.nodes():
        for y in gr.nodes():
            if x == y:
                continue
            if len(list(nx.common_neighbors(tanner_graph, x, y))) > 0:
                gr.add_edge(x, y)
    return gr

def compute_code_distance(tanner_graph: nx.Graph, op_type: str) -> int:
    """
        Computes the code distance of a code given the Tanner graph.
    """
    gr = make_support_graph(tanner_graph, op_type)
    d = np.inf

    for obs in tanner_graph.graph['obs_list'][op_type]:
        _obs = [get_data_qubit_node(tanner_graph, x) for x in obs]
        crossing_checks = set()
        for q in _obs:
            for v in tanner_graph.neighbors(q):
                if tanner_graph.nodes[v]['node_type'] == op_type:
                    crossing_checks.add(v)

        gr_dup = nx.Graph()
        # Duplicate every node in gr.
        n = max(x for x in tanner_graph.nodes())+1
        for x in gr.nodes():
            gr_dup.add_node(x)
            gr_dup.add_node(x+n)
        # Add edges as well.
        for (x, y) in gr.edges():
            if x in crossing_checks or y in crossing_checks:
                gr_dup.add_edge(x, y+n)
                gr_dup.add_edge(x+n, y)
                print('<crossing> adding edges (%d, %d) and (%d, %d)' % (x, y+n, x+n, y))
            else:
                gr_dup.add_edge(x, y)
                gr_dup.add_edge(x+n, y+n)
                print('adding edges (%d, %d) and (%d, %d)' % (x, y, x+n, y+n))
        # Compute shortest path distance between all v and v+n.
        paths = dict(nx.all_pairs_shortest_path(gr_dup))
        for v in crossing_checks:
            if v+n in paths[v]:
                print('%d to %d:' % (v, v+n), paths[v][v+n])
            else:
                print('%d to %d: not found' % (v, v+n))
        _d = min(len(paths[v][v+n]) if v+n in paths[v] else np.inf for v in crossing_checks)
        d = min(d, _d)
        break
    return d
