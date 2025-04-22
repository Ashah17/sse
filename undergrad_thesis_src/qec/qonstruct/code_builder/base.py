"""
    author: Suhas Vittal
    date:   17 January 2024

    Utility functions for constructing codes.
"""

import networkx as nx

from collections import defaultdict

def tanner_init() -> nx.Graph:
    gr = nx.Graph()
    # Initialize structures for the graph.
    gr.graph['data_qubits'] = []
    gr.graph['checks'] = {'all': [], 'x': [], 'z': []}
    gr.graph['obs_list'] = {'x': [], 'z': []}
    # These are specific to color codes:
    gr.graph['plaquettes'] = []
    gr.graph['plaquette_support_map'] = {}
    gr.graph['plaquette_color_map'] = {}
    gr.graph['plaquette_check_map'] = defaultdict(list)
    return gr

def add_data_qubit(gr: nx.Graph, q: int, **kwargs) -> None:
    gr.add_node(q, node_type='data')
    for (k, v) in kwargs.items():
        gr.nodes[q][k] = v
    gr.graph['data_qubits'].append(q)

def add_check(gr: nx.Graph, check: int, check_type: str, support: list[int], **kwargs) -> None:
    gr.add_node(check, node_type=check_type, schedule_order=[])
    for (k, v) in kwargs.items():
        gr.nodes[check][k] = v
    gr.graph['checks']['all'].append(check)
    gr.graph['checks'][check_type].append(check)
    # Add edges to data qubits:
    for q in support:
        if not gr.has_node(q):
            add_data_qubit(gr, q)
        gr.add_edge(check, q)

def add_to_support(gr: nx.Graph, check: int, q: int) -> None:
    n = gr.graph['_check_to_index'][check]
    x = gr.graph['_data_to_index'][q]
    gr.add_edge(n, x)

def get_support(gr: nx.Graph, check: int) -> list[int]:
    return [x for x in gr.neighbors(check)]

def add_observable(gr: nx.Graph, observable: list[int], obs_type: str) -> None:
    gr.graph['obs_list'][obs_type].append(observable)

# Specific Code Functions:

def add_plaquette(gr: nx.Graph, plaquette: int, support: list[int], color: int) -> None:
    gr.graph['plaquettes'].append(plaquette)
    gr.graph['plaquette_support_map'][plaquette] = support
    gr.graph['plaquette_color_map'][plaquette] = color

def set_plaquette(gr: nx.Graph, check: int, plaquette: int) -> None:
    if plaquette not in gr.graph['plaquettes']:
        add_plaquette(gr, plaquette)
    gr.graph['plaquette_check_map'][plaquette].append(check)

def get_plaquette(gr: nx.Graph, check: int) -> int:
    return gr.nodes[check]['plaquette']

