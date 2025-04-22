"""
    author: Suhas Vittal
    date:   1 December 2023

    Pre-written generation code.
"""

from qonstruct.qes.manager import QesManager

from qonstruct.codes import qecolor

def qes_hexagonal_color_code(
        output_file: str,
        distance: int,
        rounds: int,
        memory='z',
        use_flags=True,
        both_at_once=False
):
    code = qecolor.make_hexagonal_tanner_graph(distance, both_at_once=both_at_once)
    mgr = QesManager(code, use_plaquettes_instead_of_checks=(not both_at_once))
    mgr.memory = memory

    if use_flags:
        plaquettes = code.graph['plaquettes']
        for plq in plaquettes:
            support = code.graph['plaquette_support_map'][plq]
            for i in range(0, len(support), 2):
                q1 = support[i]
                q2 = support[i+1]
                if q1 is None or q2 is None:
                    continue
                mgr.add_flags_to(q1, q2, plq)
    # Write file.
    mgr.fopen(output_file)
    mgr.write_memory_experiment(rounds)
    mgr.fclose()
