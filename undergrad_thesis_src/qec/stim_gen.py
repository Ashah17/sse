import stim
import sinter


def gen_rsc(d, p):
    circ = stim.Circuit.generated(
        'surface_code:rotated_memory_z',
        rounds=10,
        distance=d,
        after_clifford_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p
    )

    yield sinter.Task(
        circuit = circ
    )

    with open('rsc_d3_circ.tex_z_exp', 'w') as f:
        f.write(str(circ)) #write to stim file


def main():
    d = 3
    p = 0.001
    sample = sinter.collect(
        num_workers = 1,
        max_shots=1_000_000,
        max_errors=1000,
        tasks=gen_rsc(d, p),
        decoders=['pymatching'],
    )

    print(sinter.CSV_HEADER)
    for s in sample:
        print(s.to_csv_line())

if __name__ == '__main__':
    main()

# with open('rsc_d3_circ.tex', 'w') as f:
#     f.write(str(circ)) #write to stim file

#now decode with sinter
