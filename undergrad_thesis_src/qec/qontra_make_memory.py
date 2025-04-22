"""
    author: Suhas Vittal
    date:   21 January 2024
"""

import qonstruct.io
import qonstruct.scheduling

from qonstruct.code_builder.surface_code import make_rotated

from qonstruct.parsing.cmd import *
from qonstruct.qes.manager import QesManager

from sys import argv

arg_list = parse(argv)

code_file = try_get_string(arg_list, 'tanner')
output_file = try_get_string(arg_list, 'qes')
rounds = try_get_int(arg_list, 'rounds')
is_memory_x = 'x' in arg_list

if code_file == 'rsc':
    d = try_get_int(arg_list, 'd')
    gr = make_rotated(d)
else:
    gr = qonstruct.io.read_tanner_graph_file(code_file)
    qonstruct.scheduling.compute_syndrome_extraction_schedule(gr)

mgr = QesManager(gr)
mgr.memory = 'x' if is_memory_x else 'z'

mgr.fopen(output_file)
mgr.write_memory_experiment(rounds)
mgr.fclose()
