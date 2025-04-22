"""
    author: Suhas Vittal
    date:   29 January 2024
"""

def concat_args(arr):
    return ','.join(str(x) for x in arr)

def write_memory(fout, gr, cycles):
    data_qubits = gr.graph['data_qubits']
    all_checks = gr.graph['checks']['all']
    mem_checks = gr.graph['checks']['z']

    n_meas_per_cycle = len(all_checks) * rounds_per_cycle
    n_events_per_cycle = len(mem_checks) * rounds_per_cycle

    data_to_meas_index = { d : len(all_checks)+i for (i,d) in enumerate(data_qubits) }
    check_to_meas_index = { ch : i for (i,ch) in enumerate(all_checks) }

    fout.write('reset %s;\n' % concat_args([*data_qubits, *all_checks]))
    fout.write('call \"ext0\";\n')
    fout.write('mshift %d;\n' % (n_meas_per_cycle - len(all_checks)))
    fout.write('eoffset %d;\n' % n_events_per_cycle)
    for _ in range(cycles-1):
        fout.write('call \"ext1\";\n')
        fout.write('mshift %d;\n' % n_meas_per_cycle)
        fout.write('eoffset %d;\n' % n_events_per_cycle)

    fout.write('measure %s;\n' % concat_args(data_qubits))
    for (i, ch) in enumerate(mem_checks):
        fout.write('event %d,%d,%s;\n'
                % (i, check_to_meas_index[ch], concat_args([data_to_meas_index[x] for x in gr.neighbors(ch)])))
    for (i, obs) in enumerate(gr.graph['obs_list']['z']):
        fout.write('obs %d,%s;\n' % (i, concat_args([data_to_meas_index[x] for x in obs])))

def write_ext(fout, gr, is_first_round=False, eoff=0, moff=0):
    data_qubits = gr.graph['data_qubits']
    all_checks = gr.graph['checks']['all']
    x_checks = gr.graph['checks']['x']
    mem_checks = gr.graph['checks']['z']

    check_to_meas_index = { ch : i for (i,ch) in enumerate(all_checks) }

    fout.write('@annotation timing_error\n'
               'h %s;\n' % concat_args(x_checks))
    k = 0
    while True:
        cnot_args = []
        for ch in all_checks:
            if k >= len(gr.nodes[ch]['schedule_order']):
                continue
            q = gr.nodes[ch]['schedule_order'][k]
            if q is None:
                continue
            if gr.nodes[ch]['node_type'] == 'z':
                cnot_args.extend([q, ch])
            else:
                cnot_args.extend([ch, q])
        if len(cnot_args) == 0:
            break
        fout.write('cx %s;\n' % concat_args(cnot_args))
        k += 1
    fout.write('h %s;\n' % concat_args(x_checks))
    fout.write('measure %s;\n' % concat_args(all_checks))
    fout.write('reset %s;\n' % concat_args(all_checks))
    for (i, ch) in enumerate(mem_checks):
        e_args = [check_to_meas_index[ch]+moff]
        if not is_first_round:
            e_args.append(check_to_meas_index[ch]+len(all_checks)+moff)
        fout.write('event %d,%s;\n' % (i+eoff, concat_args(e_args)))

if __name__ == '__main__':
    import qonstruct.io
    import qonstruct.scheduling
    from qonstruct.code_builder.surface_code import make_rotated

    from qonstruct.parsing.cmd import *

    from sys import argv

    import os

    arg_list = parse(argv[3:])

    code_file = argv[1]
    output_folder = argv[2]
    rounds_per_cycle = try_get_int(arg_list, 'r')
    cycles = try_get_int(arg_list, 'cycles')
    # Check if the user is requesting a prebuilt code (i.e. rotated surface code (rsc))
    if code_file == 'rsc':
        d = try_get_int(arg_list, 'd')
        gr = make_rotated(d)
    else:
        gr = qonstruct.io.read_tanner_graph_file(code_file)
        qonstruct.scheduling.compute_syndrome_extraction_schedule(gr)
    main_file = '%s/main.qes' % output_folder
    ext0_file = '%s/ext0.qes' % output_folder
    ext1_file = '%s/ext1.qes' % output_folder

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    with open(main_file, 'w') as fout:
        write_memory(fout, gr, cycles)
    with open(ext0_file, 'w') as fout:
        moff, eoff = 0, 0
        for r in range(rounds_per_cycle):
            write_ext(fout, gr, is_first_round=(r==0), moff=moff, eoff=eoff)
            if r > 0:
                moff += len(gr.graph['checks']['all'])
            eoff += len(gr.graph['checks']['z'])
    with open(ext1_file, 'w') as fout:
        moff, eoff = 0, 0
        for r in range(rounds_per_cycle):
            write_ext(fout, gr, is_first_round=False, moff=moff, eoff=eoff)
            if 'sse' in arg_list:
                if r >= rounds_per_cycle//2:
                    fout.write('return_if \"$reoz\";\n')
            moff += len(gr.graph['checks']['all'])
            eoff += len(gr.graph['checks']['z'])
    # Write .ini file.
    ini_file = '%s/config.ini' % output_folder
    total_events = (rounds_per_cycle*cycles + 1) * len(gr.graph['checks']['z'])
    total_obs = len(gr.graph['obs_list']['z'])
    rreoz_events = (rounds_per_cycle//2) * len(gr.graph['checks']['z'])
    with open(ini_file, 'w') as fout:
        fout.write('ext0 = ext0.qes\n')
        fout.write('ext1 = ext1.qes\n')
        fout.write('[Config]\n')
        fout.write('record_events_until = %d\n' % total_events)
        fout.write('record_obs_until = %d\n' % total_obs)
        fout.write('[SSE]\n')
        fout.write('$reoz_track_last_n_events = %d\n' % rreoz_events)

