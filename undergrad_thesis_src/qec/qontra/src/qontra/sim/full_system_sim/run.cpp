/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include "qontra/sim/full_system_sim.h"

#include <fstream>
#include <iostream>

#include <stdio.h>

namespace qontra {

FullSystemSimulator::FullSystemSimulator()
    :register_file(1, 1),
    elapsed_time(0),
    shot_time_delta_map(),
    syndrome_table(1, 1),
    observable_table(1, 1),
    meas_ctr(0),
    meas_offset(0),
    event_offset(0),
    max_event_written_to(0),
    max_obs_written_to(0),
    sample_circuit(),
    is_recording_stim_instructions(false),
    n_qubits(1),
    current_shots(0),
    base_sim(nullptr),
    shot_histogram(),
    subroutine_map(),
    register_file_cpy(1, 1),
    syndrome_table_cpy(1, 1),
    observable_table_cpy(1, 1),
    logical_errors(MPI_SUM, 1)
{}

// void
// FullSystemSimulator::run_batch_orig(const qes::Program<>& program, uint64_t shots) {
//     base_sim->reset_sim();
//     base_sim->shots = shots;
//     current_shots = shots;
//     // Reset all structures.
//     register_file.clear();
    
//     elapsed_time = 0.0;
//     shot_time_delta_map.clear();

//     syndrome_table.clear();
//     observable_table.clear();
    
//     meas_ctr = 0;
//     meas_offset = 0;
//     event_offset = 0;
//     max_event_written_to = 0;
//     max_obs_written_to = 0;

//     execute_routine(program);
// }

// void
// FullSystemSimulator::write_stats(uint64_t batchno) {
//     // Write stim file.

//     std::string file_name = config.stim_output_file + "/" + config.stim_output_file + "_shot_" + std::to_string(batchno) + ".stim";

//     if (is_recording_stim_instructions) {
//         std::ofstream stim_fout(file_name);
//         stim_fout << sample_circuit.str() << std::endl; //printing without noise so clean
//     }

//     // Write syndromes and observables to a file.
//     //
//     // First, we must combine the two tables into one table before writing.
//     const uint64_t n_det = config.record_events_until >= 0 ? config.record_events_until : max_event_written_to,
//                     n_obs = config.record_obs_until >= 0 ? config.record_obs_until : max_obs_written_to;

//     stim::simd_bits<SIMD_WIDTH> ref(n_det+n_obs);
//     stim::simd_bit_table<SIMD_WIDTH> output_trace = syndrome_table.concat_major(observable_table, n_det, n_obs);
//     // Check popcnts.
//     std::string filename = config.syndrome_output_folder + "/" + get_batch_filename(batchno);
//     FILE* trace_fout = fopen(filename.c_str(), "w");
//     stim::write_table_data(trace_fout,
//                             current_shots,
//                             n_det+n_obs,
//                             ref,
//                             output_trace,
//                             stim::SampleFormat::SAMPLE_FORMAT_DETS,
//                             'D',
//                             'D',
//                             0);
//     fclose(trace_fout);
//     // Write other data to output_file here:

//     // Update histogram.
//     auto obs_tr = observable_table.transposed();
//     for (uint64_t t = 0; t < current_shots; t++) {
//         const size_t n_obs_u64 = (max_obs_written_to >> 6) + 1;
//         std::vector<uint64_t> obs_vec(obs_tr[t].u64, obs_tr[t].u64 + n_obs_u64);
//         shot_histogram[obs_vec]++;
//     }
// }

void
FullSystemSimulator::run_batch(const qes::Program<>& program, uint64_t shots) {
    base_sim->reset_sim();
    base_sim->shots = shots;
    current_shots = shots;
    // Reset all structures.
    register_file.clear();
    
    elapsed_time = 0.0;
    shot_time_delta_map.clear();

    // for (int i = 0; i < G_RECORD_SPACE_SIZE; i++) {
    //     if (syndrome_table[i].popcnt() != 0) {
    //         std::cout<<"WE GOT A 1 BROOOOO"<<std::endl;
    //     }
    // }

    // for (int i = 0; i < G_RECORD_SPACE_SIZE; i++) {
    //     if (syndrome_table[i].popcnt() != 0) {
    //         std::cout<<"WE GOT A 1"<<std::endl;
    //     }
    // }

    syndrome_table.clear();
    observable_table.clear();
    
    meas_ctr = 0;
    meas_offset = 0;
    event_offset = 0;
    max_event_written_to = 0;
    max_obs_written_to = 0;

    execute_routine(program);
}

void
FullSystemSimulator::write_stats(uint64_t batchno) {
    // Write stim file.
    if (is_recording_stim_instructions) {
        std::ofstream stim_fout(config.stim_output_file);
        stim_fout << sample_circuit.str() << std::endl;
    }
    // Write syndromes and observables to a file.
    //
    // First, we must combine the two tables into one table before writing.
    const uint64_t n_det = config.record_events_until >= 0 ? config.record_events_until : max_event_written_to,
                    n_obs = config.record_obs_until >= 0 ? config.record_obs_until : max_obs_written_to;

    stim::simd_bits<SIMD_WIDTH> ref(n_det+n_obs);
    stim::simd_bit_table<SIMD_WIDTH> output_trace = syndrome_table.concat_major(observable_table, n_det, n_obs);

    // std::ofstream outFile("actual_flipped_frames.txt", std::ios::app);

    // outFile << "AND THE OBS IS " << observable_table.str() << std::endl;

    // outFile.close();

    // Check popcnts.
    std::string filename = config.syndrome_output_folder + "/" + get_batch_filename(batchno);
    FILE* trace_fout = fopen(filename.c_str(), "w");
    stim::write_table_data(trace_fout,
                            current_shots,
                            n_det+n_obs,
                            ref,
                            output_trace,
                            stim::SampleFormat::SAMPLE_FORMAT_DETS,
                            'L',
                            'D',
                            0);
    fclose(trace_fout);
    // Write other data to output_file here:


    // Update histogram.
    auto obs_tr = observable_table.transposed();


    for (uint64_t t = 0; t < current_shots; t++) {
        const size_t n_obs_u64 = (max_obs_written_to >> 6) + 1;
        std::vector<uint64_t> obs_vec(obs_tr[t].u64, obs_tr[t].u64 + n_obs_u64);
        shot_histogram[obs_vec]++;
    }
}

// void
// FullSystemSimulator::write_stats(uint64_t batchno) {
//     // Write stim file for each individual shot to a specific folder
    
//     std::string file_name = config.stim_output_file + "/" + config.stim_output_file + "_shot_" + std::to_string(batchno) + ".stim";

//     if (is_recording_stim_instructions) {
//         //creating a new stim file each time for write_stats based on batchno
//         //std::ofstream stim_fout(config.stim_output_file);
//         std::ofstream stim_fout(file_name);
//         stim_fout << sample_circuit.str() << std::endl;
//     }
//     // Write syndromes and observables to a file.
//     //
//     // First, we must combine the two tables into one table before writing.
//     const uint64_t n_det = config.record_events_until >= 0 ? config.record_events_until : max_event_written_to,
//                     n_obs = config.record_obs_until >= 0 ? config.record_obs_until : max_obs_written_to;

//     stim::simd_bits<SIMD_WIDTH> ref(n_det+n_obs);
//     stim::simd_bit_table<SIMD_WIDTH> output_trace = syndrome_table.concat_major(observable_table, n_det, n_obs);
//     // Check popcnts.
//     std::string filename = config.syndrome_output_folder + "/" + get_batch_filename(batchno);
//     FILE* trace_fout = fopen(filename.c_str(), "w");
//     stim::write_table_data(trace_fout,
//                             current_shots,
//                             n_det+n_obs,
//                             ref,
//                             output_trace,
//                             stim::SampleFormat::SAMPLE_FORMAT_DETS,
//                             'D',
//                             'D',
//                             0);
//     fclose(trace_fout);
//     // Write other data to output_file here:

//     // Update histogram.
//     auto obs_tr = observable_table.transposed();
//     for (uint64_t t = 0; t < current_shots; t++) {
//         const size_t n_obs_u64 = (max_obs_written_to >> 6) + 1;
//         std::vector<uint64_t> obs_vec(obs_tr[t].u64, obs_tr[t].u64 + n_obs_u64);
//         shot_histogram[obs_vec]++;
//     }
// }

}   // qontra
