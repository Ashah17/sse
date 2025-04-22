/*
 *  author: Suhas Vittal
 *  date:   23 January 2024
 * */

#ifndef QONTRA_FULL_SYSTEM_SIM_h
#define QONTRA_FULL_SYSTEM_SIM_h

#include "qontra/experiments/histogram.h"
#include "qontra/ext/qes.h"
#include "qontra/sim/base/state_sim.h"
#include "qontra/tables.h"

#include <deque>
#include <map>
#include <vector>
#include <tuple>

namespace qontra {

size_t  get_register_index(std::string);

struct program_status_t {
    program_status_t(uint64_t shots)
        :pc(0),
        pc_specific_map(),
        branch_and_wait_reach_map(),
        return_if_waiting_trials(shots)
    {}
        
    // PC structures:
    uint64_t                        pc;
    std::map<uint64_t, uint64_t>    pc_specific_map;
    // Branch tracking structures:
    std::map<uint64_t, stim::simd_bits<SIMD_WIDTH>>
        branch_and_wait_reach_map;
    stim::simd_bits<SIMD_WIDTH>
        return_if_waiting_trials;
};

class FullSystemSimulator {
public:
    FullSystemSimulator();

    // run_program steps through the program in batches. The output is a shot histogram
    // of all possible outputs. The keys of the histogram are the measured observables in each trial,
    // and the values are the occurrences of each key.
    template <class SIM>
    histogram_t<uint64_t> run_program(const qes::Program<>&, uint64_t shots);

    void load_subroutine(std::string name, const qes::Program<>&);

    struct {
        // Simulation configuration:
        uint64_t    distance;
        ErrorTable  errors;
        TimeTable   timing;
        std::string stim_output_file;
        std::string syndrome_output_folder;
        std::string data_output_file;
    
        // Normally, the simulator will track the max written event. However, to force the simulator
        // to write events that may not be written, these two parameters can be set.
        int64_t     record_events_until = -1;
        int64_t     record_obs_until = -1;

        // Microarchitecture:
        size_t  n_registers = 64;

        struct {
            size_t  rreoz_track_last_n_events = 0;
        } sse;
    } config;

    statistic_t<uint64_t> logical_errors;

    // qes::Program<>      create_program(int d, double measurement_percentage); //ADDED THIS

    qes::Program<>      create_program(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs); //ADDED THIS

    qes::Program<>      create_program_masked_selective(int d, double p, int policy, int exp_type, int numRounds); //2 params: double p masked, int for policy code, int for mem_exp type
    //pass in numRounds too

    qes::Program<>      create_program_stability_exp(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs, int stabilityOrControl);

    // std::tuple<qes::Program<>, int numEvents, int numMeasurements, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, int roundNum>
    // return type ^^

    std::tuple<qes::Program<>, qes::Program<>, int, int, std::map<int, int>, std::map<int, int>, int>      stability_syndrome_prologue(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs, int numEvents, int numMeasurements, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, int roundNum, qes::Program<> program);
    std::tuple<qes::Program<>, qes::Program<>, int, int, std::map<int, int>, std::map<int, int>, int>      stability_syndrome_extraction_round(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs, int numEvents, int numMeasurements, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, int roundNum, qes::Program<> program);
    std::tuple<qes::Program<>, qes::Program<>, int, int, std::map<int, int>, std::map<int, int>, int>      stability_syndrome_epilogue(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs, int numEvents, int numMeasurements, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, int roundNum, qes::Program<> program);

    std::map<int, std::vector<int>>  create_p_to_d_map(int d); //helper added too

    std::map<int, std::vector<int>> create_p_to_d_map_for_stability_exp(int d); //for stab exp

    std::map<int, double>      read_dets(std::string path); //new method to read dets files

    void    run_batch(const qes::Program<>&, uint64_t shots_in_batch); //making public

    stim::simd_bit_table<SIMD_WIDTH>    syndrome_table; //making public

    void    execute_routine(const qes::Program<>& program); //making public


    std::tuple<std::set<int>, std::set<int>, std::set<int>, std::set<int>, qes::Program<>, int> create_masked_prologue(int d, qes::Program<> program, int instruction_number);
    std::tuple<qes::Program<>, int, int, int, std::map<int, int>, std::map<int, int>, std::vector<int>, qes::Program<>, int>    create_masked_se_round(int d, qes::Program<> program, int numEvents, int numMeasurements, int numRounds, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, std::set<int> stabsThisRound, std::set<int> allStabs, std::set<int> innerStabs, std::set<int> xStabs, std::set<int> zStabs, std::vector<int> firstRound, int instruction_number, int disconnectRound);

    qes::Program<> create_masked_epilogue(int d, qes::Program<> program, int numEvents, int numMeasurements, int round, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, std::set<int> selective_stabs, std::set<int> allStab, std::set<int> innerStabilizers, std::set<int> xStabilizers, std::set<int> zStabilizers, std::vector<int> firstRound, int instruction_number);


    //making these public

    stim::simd_bit_table<SIMD_WIDTH>    observable_table;
    stim::simd_bit_table<SIMD_WIDTH>    observable_table_cpy;
private:
    
    void    write_stats(uint64_t batchno);

    void    read_next_instruction(const qes::Program<>& from, program_status_t&);

    fp_t    do_gate(const qes::Instruction<>&, program_status_t&, int64_t trial=-1);
    fp_t    do_measurement(const qes::Instruction<>&, program_status_t&, int64_t trial=-1);

    void    create_event_or_obs(const qes::Instruction<>&, program_status_t&);

    void    recalibrate_timing(void);

    void    inject_timing_error(void);
    void    inject_idling_error_positive(std::vector<uint64_t> on_qubits, int64_t trial=-1);
    void    inject_idling_error_negative(std::vector<uint64_t> not_on_qubits, int64_t trial=-1);

    void    snapshot(void);
    void    rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH>);

    stim::simd_bits_range_ref<SIMD_WIDTH> get_register(std::string);

    // Architectural structures:
    //
    // The register_file is implemented as bitvectors (so each register is
    // actually a 1-bit value). This is done because general arithmetic (adds,
    // multiplies, etc.) are generally hard to implement efficiently, whereas
    // bitwise arithmetic (which we expect is the common case) is efficient to
    // implement.
    stim::simd_bit_table<SIMD_WIDTH>    register_file;

    // Timing structures:
    //
    // elapsed_time is the time elapsed to the current point of execution.
    //
    // shot_time_delta_map is to accomodate events that may have taken less or
    // more time to execute certain instructions. As this is not a common issue,
    // we maintain a table of trials which have this delta. 
    //
    // Both elapsed_time and shot_time_delta_map are reset upon injecting a
    // timing error.
    fp_t                            elapsed_time;
    std::map<uint64_t, fp_t>        shot_time_delta_map;

    // Tracking structures:
    //
    // syndrome_table and observable_table are the outputs written to the trace
    // file.
    //
    // meas_ctr is the number of measurements that have been recorded.
    //
    // meas_offset is used to record measurements at "meas_ctr + meas_offset".
    // This is useful when measurements are repetitive and need not have a long
    // lifetime.
    //
    // event_offset is akin to meas_offset, but for events.
    
    

    int64_t     meas_ctr;
    int64_t     meas_offset;
    int64_t     event_offset;

    uint64_t    max_event_written_to;
    uint64_t    max_obs_written_to;

    // Stim structures:
    DetailedStimCircuit     sample_circuit;
    bool                    is_recording_stim_instructions;

    // Other simulation structures:
    uint64_t                n_qubits;
    uint64_t                current_shots;
    uptr<StateSimulator>    base_sim;
    histogram_t<uint64_t>   shot_histogram;

    std::map<std::string, qes::Program<>>   subroutine_map;

    // Table copies for snapshots and rollbacks:
    stim::simd_bit_table<SIMD_WIDTH>    register_file_cpy;
    stim::simd_bit_table<SIMD_WIDTH>    syndrome_table_cpy;
    //moved obs table public

    // Policy specific stuff below:
    struct {
        std::deque<uint64_t>   rreoz_events;
    } sse;
};

}   // qontra

#include "full_system_sim.inl"

#endif  // QONTRA_FULL_SYSTEM_SIM_h
