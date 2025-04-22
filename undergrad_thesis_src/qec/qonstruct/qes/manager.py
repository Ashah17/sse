"""
    author: Suhas Vittal
    date:   1 December 2023
"""

import networkx as nx

from collections import defaultdict

def concat(arr: list[int], delimiter=','):
    return delimiter.join(str(x) for x in arr)

class QesManager:
    def __init__(self, tanner_graph: nx.Graph, use_plaquettes_instead_of_checks=False):
        self.code = tanner_graph
        if use_plaquettes_instead_of_checks:
            self.n = len(tanner_graph.graph['data']) + len(tanner_graph.graph['plaquettes'])
        else:
            self.n = tanner_graph.number_of_nodes()

        # Public parameters:
        self.memory = 'z'

        # Simulation structures
        self.meas_ctr_map = {'curr': 0}
        self.event_ctr = 0
        self.obs_ctr = 0

        # Flag data structures:
        #   - flag_assignment_map[x] for some parity qubit x contains
        #       { data : flag } and an entry 'all' which contains all flags
        #       used in the syndrome extraction of x.
        #   - flag_edge_map[flag] corresponds to the parity qubits
        #       affected by the flag. The data contained in the entry is
        #       as follows:
        #           base: the owning stabilizer
        #           edges: [a, b] lists of parity qubits that may have a hook error.
        self.flag_qubits = []
        self.flag_assignment_map = {}
        self.flag_edge_map = {}

        # Parameters that should be set via a function:
        self._mode_use_plaquettes_instead_of_checks = use_plaquettes_instead_of_checks

    def add_flags_to(self, q1: int, q2: int, check: int) -> int:
        """
            Adds flags to qubits q1 and q2 during the measurement of the specified check.
            Returns the flag qubit id.
        """
        if check not in self.flag_assignment_map:
            self.flag_assignment_map[check] = {'all': []}
        # Compute flag edges.
        self.flag_edge_map[self.n] = {'base': check, 'edges': []}
        affected_checks = set()
        _tmp = [(q1, q2), (q2, q1)]
        for (qa, qb) in _tmp:
            for x in self.code.neighbors(qa):
                if not self.code.has_edge(x, qb):
                    if self._mode_use_plaquettes_instead_of_checks:
                        # We only need one stabilizer.
                        affected_checks.add(self.code.nodes[x]['plaquette'])
                    else:
                        stabilizer = 'x' if self.code.nodes[check]['node_type'] == 'z' else 'z'
                        if self.code.nodes[x]['node_type'] != stabilizer:
                            affected_checks.add(x)
        affected_checks = list(affected_checks)
        if len(affected_checks) > 2 or len(affected_checks) == 0:
            del self.flag_edge_map[self.n]
            return -1
        if len(affected_checks) == 2:
            c1, c2 = affected_checks
            self.flag_edge_map[self.n]['edges'].append((c1, c2))
        else:
            c1 = affected_checks[0]
            # Other is the boundary.
            self.flag_edge_map[self.n]['edges'].append((c1, -1))
        # If everything checks out, update flag assignment map.
        self.flag_assignment_map[check]['all'].append(self.n)
        self.flag_assignment_map[check][q1] = self.n
        self.flag_assignment_map[check][q2] = self.n
        self.flag_qubits.append(self.n)
        self.n += 1

    def get_detection_events(self, rounds: int, offset=0) -> tuple[list[int], dict[int, list[int]]]:
        """
            Returns a list of qubits which will have detection events in order,
                and a mapping from qubits to detection events.

            If offset is nonzero, then any detection events less than offset are ignored.
        """
        detection_event_order = []
        # These are the parity qubits.
        if self._mode_use_plaquettes_instead_of_checks:
            # The memory type changes the order. For color codes, we do
            # X then Z.
            if self.memory == 'x':
                detection_event_order.extend(self.code.graph['plaquettes'])
                detection_event_order.extend(self.flag_qubits)
            else:
                detection_event_order.extend(self.flag_qubits)
                detection_event_order.extend(self.code.graph['plaquettes'])
        else:
            detection_event_order.extend(self.code.graph['checks'][self.memory])
            detection_event_order.extend(self.flag_qubits)
        detection_events = []
        detection_event_map = defaultdict(list)
        for r in range(rounds):
            for (i, x) in enumerate(detection_event_order):
                k = r*len(detection_event_order) + i
                if k < offset:
                    continue
                detection_events.append(x)
                detection_event_map[x].append(k)
        # We also have epilogue events.
        base_k = rounds*len(detection_event_order)
        if self._mode_use_plaquettes_instead_of_checks:
            checks = self.code.graph['plaquettes']
        else:
            checks = self.code.graph['checks'][self.memory]
        for (i, x) in enumerate(checks):
            detection_events.append(x)
            detection_event_map[x].append(base_k + i)
        return detection_events, detection_event_map

    def fopen(self, filename: str):
        self._writer = open(filename, 'w')

    def fclose(self):
        self._writer.close()

    def write_memory_experiment(self, rounds: int):
        data_qubits = self.code.graph['data_qubits']
        if self._mode_use_plaquettes_instead_of_checks:
            plaquettes = self.code.graph['plaquettes']
        else:
            mem_checks = self.code.graph['checks'][self.memory]
        detection_events, detection_event_map = self.get_detection_events(rounds)
        #
        # PROLOGUE
        #
        self.big_comment('PROLOGUE')
        self.reset(data_qubits)
        if self.memory == 'x':
            self.h(data_qubits)
        #
        # BODY
        #
        self.big_comment('BODY')
        for r in range(rounds):
            prev_meas_ctr_map = self.meas_ctr_map.copy()
            if self._mode_use_plaquettes_instead_of_checks:
                # This will be an X then Z stabilizer syndrome extraction round. Each
                # plaquette will correspond to two stabilizers.
                for s in ['x', 'z']:
                    self.annotation('inject_timing_error')
                    self._plaquette_syndrome_extraction(s)
                    # Create detection events here.
                    if self.memory != s:
                        # These are flag detection events.
                        for fq in self.flag_qubits:
                            event_no = detection_event_map[fq][r]
                            m = self.meas_ctr_map[fq]
                            # Get next detection event.
                            self.annotation('flag')
                            self.event(event_no, m)
                    else:
                        # These are stabilizer detection events.
                        for pq in plaquettes:
                            event_no = detection_event_map[pq][r]
                            if r == 0:
                                m = self.meas_ctr_map[pq]
                                self.property('color', self.code.graph['plaquette_color_map'][pq])
                                self.event(event_no, m)
                            else:
                                m1, m2 = self.meas_ctr_map[pq], prev_meas_ctr_map[pq]
                                self.property('color', self.code.graph['plaquette_color_map'][pq])
                                self.event(event_no, m1, m2)
            else:
                self.annotation('inject_timing_error')
                self._standard_syndrome_extraction()
                # Create detection events.
                for pq in mem_checks:
                    event_no = detection_event_map[pq][r]
                    if r == 0:
                        m = self.meas_ctr_map[pq]
                        self.event(event_no, m)
                    else:
                        m1, m2 = self.meas_ctr_map[pq], prev_meas_ctr_map[pq]
                        self.event(event_no, m1, m2)
        #
        # EPILOGUE
        #
        if self.memory == 'x':
            self.h(data_qubits)
        self.measure(data_qubits)
        # Create detection events and observables.
        if self._mode_use_plaquettes_instead_of_checks:
            for pq in plaquettes:
                event_no = detection_event_map[pq][rounds]
                support = self.code.graph['plaquette_support_map'][pq]
                meas_array = [self.meas_ctr_map[x] for x in support if x is not None]
                meas_array.append(self.meas_ctr_map[pq])
                self.property('color', self.code.graph['plaquette_color_map'][pq])
                self.event(event_no, *meas_array)
        else:
            for pq in mem_checks:
                event_no = detection_event_map[pq][rounds]
                support = self.code.nodes[pq]['schedule_order']
                meas_array = [self.meas_ctr_map[x] for x in support if x is not None]
                meas_array.append(self.meas_ctr_map[pq])
                self.event(event_no, *meas_array)
        # Write observables.
        obs_list = self.code.graph['obs_list'][self.memory]
        for obs in obs_list:
            meas_array = [self.meas_ctr_map[x] for x in obs]
            self.auto_obs(*meas_array)

    def h(self, operands: list[int]):
        self._op('h', operands)

    def cx(self, operands: list[int]):
        self._op('cx', operands)

    def measure(self, operands: list[int], update_meas_ctr_map=True):
        self._op('measure', operands)
        k = self.meas_ctr_map['curr']
        for (i, x) in enumerate(operands):
            if update_meas_ctr_map:
                self.meas_ctr_map[x] = k+i
        self.meas_ctr_map['curr'] = k + len(operands)

    def reset(self, operands: list[int]):
        self._op('reset', operands)
    
    def event(self, event_no: int, *meas_operands: list[int]):
        self._op('event', [event_no, *meas_operands])

    def obs(self, obs_no: int, *meas_operands: list[int]):
        self._op('obs', [obs_no, *meas_operands])
    
    def auto_event(self, *meas_operands: list[int]):
        self._op('event', [self.event_ctr, *meas_operands])
        self.event_ctr += 1

    def auto_obs(self, *meas_operands: list[int]):
        self._op('obs', [self.obs_ctr, *meas_operands])
        self.obs_ctr += 1

    def _op(self, opname: str, operands: list[int]):
        if len(operands) == 0:
            return
        self._writer.write('%s %s;\n' % (opname, concat(operands)))

    def annotation(self, name: str):
        self._writer.write('@annotation %s\n' % name)

    def property(self, name: str, value: int|float):
        self._writer.write('@property %s %s\n' % (name, str(value)))

    def comment(self, text: str):
        self._writer.write('# %s\n' % text)

    def big_comment(self, text: str):
        self._writer.write('#\n# %s \n#\n' % text)

    def _standard_syndrome_extraction(self):
        checks = self.code.graph['checks']['all']
        x_checks = self.code.graph['checks']['x']

        self.reset(self.flag_qubits)
        self.reset(checks)

        self.h(x_checks)

        self.comment('FLAG SETUP')
        self._flag_hadamards(checks)
        self._flag_cnots(checks)

        self.comment('DATA CNOTS')
        self._data_cnots(checks)

        self.comment('FLAG TEARDOWN')
        self._flag_cnots(checks)
        self._flag_hadamards(checks)

        self.h(x_checks)
        self.measure(self.flag_qubits)
        self.measure(checks)

    def _plaquette_syndrome_extraction(self, stabilizer: str):
        plaquettes = self.code.graph['plaquettes']

        self.reset(self.flag_qubits)
        self.reset(plaquettes)

        if stabilizer == 'x':
            self.h(plaquettes)

        self.comment('FLAG SETUP (%s)' % stabilizer)
        self._flag_hadamards(plaquettes, stabilizer)
        self._flag_cnots(plaquettes, stabilizer)

        self.comment('DATA CNOTS (%s)' % stabilizer)
        self._data_cnots(plaquettes, stabilizer)

        self.comment('FLAG TEARDOWN (%s)' % stabilizer)
        self._flag_cnots(plaquettes, stabilizer)
        self._flag_hadamards(plaquettes, stabilizer)

        if stabilizer == 'x':
            self.h(plaquettes)
        self.measure(self.flag_qubits, update_meas_ctr_map=(self.memory != stabilizer))
        self.measure(plaquettes, update_meas_ctr_map=(self.memory == stabilizer))

    def _flag_hadamards(self, checks: list[int], plaquette_stabilizer=''):
        if self._mode_use_plaquettes_instead_of_checks and plaquette_stabilizer == 'x':
            return # Nothing to be done.
        depth = 0
        while True: 
            h_list = []
            for pq in checks:
                if pq not in self.flag_assignment_map:
                    continue
                flag_list = self.flag_assignment_map[pq]['all']
                if depth >= len(flag_list):
                    continue
                fq = flag_list[depth]
                if self._mode_use_plaquettes_instead_of_checks:
                    s = 'z'
                else:
                    s = self.code.nodes[pq]['node_type']
                if s == 'z':
                    h_list.append(fq)
            if len(h_list) == 0:
                break
            self.h(h_list)
            depth += 1

    def _flag_cnots(self, checks: list[int], plaquette_stabilizer=''):
        # Note: plaquette_stabilizer is unused if _mode_use_plaquette_instead_of_checks is not set.
        depth = 0
        while True:
            cx_list = []
            for pq in checks:
                if pq not in self.flag_assignment_map:
                    continue
                flag_list = self.flag_assignment_map[pq]['all']
                if depth >= len(flag_list):
                    continue
                fq = flag_list[depth]
                if self._mode_use_plaquettes_instead_of_checks:
                    s = plaquette_stabilizer
                else:
                    s = self.code.nodes[pq]['node_type']
                if s == 'x':
                    cx_list.extend([pq, fq])
                else:
                    cx_list.extend([fq, pq])
            if len(cx_list) == 0:
                break
            self.cx(cx_list)
            depth += 1

    def _data_cnots(self, checks: list[int], plaquette_stabilizer=''):
        depth = 0
        while True:
            cx_list = []
            for pq in checks:
                if self._mode_use_plaquettes_instead_of_checks:
                    support = self.code.graph['plaquette_support_map'][pq]
                else:
                    support = self.code.nodes[pq]['schedule_order']
                if depth >= len(support):
                    continue
                dq = support[depth]
                if dq is None:
                    continue
                cx_list.extend(self._get_cnot_targets(dq, pq, plaquette_stabilizer))
            if len(cx_list) == 0:
                break
            self.cx(cx_list)
            depth += 1

    def _get_cnot_targets(self, dq: int, pq: int, plaquette_stabilizer='') -> tuple[int, int]:
        if self._mode_use_plaquettes_instead_of_checks:
            s = plaquette_stabilizer
        else:
            s = self.code.nodes[pq]['node_type']
        if pq not in self.flag_assignment_map or dq not in self.flag_assignment_map[pq]:
            fq = pq
        else:
            fq = self.flag_assignment_map[pq][dq]
        if s == 'x':
            return fq, dq
        else:
            return dq, fq

