/*
 *  author: Suhas Vittal
 *  date:   11 January 2024
 * */

#include "qontra/experiments.h"
#include "qontra/experiments/stats.h"
#include "qontra/ext/qes.h"
#include <vector>

#include "qontra/graph/decoding_graph.h" //for decoding graph

#include "qontra/decoder/mwpm.h"

static int i = 0;

namespace qontra {

inline DetailedStimCircuit
make_circuit(std::string qes_file, fp_t p) {
    return make_circuit(qes::from_file(qes_file), p);
}

inline DetailedStimCircuit
make_circuit(const qes::Program<> program, fp_t p) {
    const size_t n = get_number_of_qubits(program);

    tables::ErrorAndTiming et;
    et.e_g1q *= 0.1;
    et.e_idle *= 0.1;
    et = et * (1000*p);

    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);
    return DetailedStimCircuit::from_qes(program, errors, timing);
}

inline memory_result_t
memory_experiment(Decoder* dec, memory_config_t config) {
    return memory_experiment(dec, config, [] (shot_payload_t x) {}, [] (Decoder::result_t x) {});
}

static int shot_count = 0; 
template <class PROLOGUE, class EPILOGUE> memory_result_t
memory_experiment(Decoder* dec, memory_config_t config, PROLOGUE p_cb, EPILOGUE e_cb) {
    const DetailedStimCircuit circuit = dec->get_circuit();
    const size_t n_obs = circuit.count_observables();

    statistic_t<uint64_t> logical_errors(MPI_SUM, n_obs+1);
    statistic_t<uint64_t> hw_sum(MPI_SUM), hw_sqr_sum(MPI_SUM), hw_max(MPI_MAX);
    statistic_t<fp_t> t_sum(MPI_SUM), t_sqr_sum(MPI_SUM), t_max(MPI_MAX);

    int count = 0;
    auto dec_cb = [&] (shot_payload_t payload)
    {
        
        // std::cout<<payload.syndrome<<std::endl;

        stim::simd_bits<SIMD_WIDTH> syndrome(payload.syndrome),
                                    obs(payload.observables);
        p_cb({syndrome, obs});
        const size_t hw = syndrome.popcnt();

        auto size = syndrome.num_bits_padded();
        int index = 0;

        // for (int i = 0; i < size; i++) {
        //     if (hw != 1) {
        //         std::cout<< i << "_____" << hw << std::endl;
        //     }
        // }

        // std::cout<<hw<<std::endl;
        // std::cout<<"hw"<<std::endl;

        // Update HW statistics and skip the trial if the HW is too small
        
        // and filtering is enabled.
        // std::cout<<count<<std::endl;
        // // // std::cout<<hw<<std::endl;
        // count+=1;

        // std::cout<<syndrome<<std::endl;

        hw_sum += hw;
        hw_sqr_sum += sqr(hw);
        hw_max.scalar_replace_if_better_extrema(hw);

        
        if (G_FILTER_OUT_SYNDROMES && hw <= G_FILTERING_HAMMING_WEIGHT) {
            return;
        }
        // Decode syndrome

        // std::ofstream outFile("flipped_stabilizer_FA.txt", std::ios::app);

        // auto size = syndrome.num_bits_padded();

        // for (int i = 0; i < size; i++) {
        //     if (syndrome[i] == 1) {
        //         outFile << "Actual flip at: " << i << std::endl;
        //     }
        // }

        // std::cout << "decoding error shot " << shot_count++ << std::endl;
        auto res = dec->decode_error(syndrome);

        std::ofstream outputFile("slw_correction_vectors.txt", std::ios::app);

        outputFile << res.corr.str() << std::endl; //correction string

        outputFile.close();

        // auto ea = res.error_assignments;

        // for (auto a : ea) {
        //     std::cout << "hi" << std::endl;
        // }

        // for (auto assign : ea) {
        //     auto x = std::get<0>(assign);
        //     auto y = std::get<1>(assign);
        //     stim::simd_bits<SIMD_WIDTH> local = std::get<2>(assign);

        //     std::cout<<x<<" WITH " << y << " AND LOCAL IS" << std::endl;

        //     std::cout<<"innerPass"<<std::endl;
        // }
    

        // std::cout<<"OBS AMOUNT" << circuit.count_observables() <<std::endl;

        // std::cout<<res.corr.str()<<std::endl;

        // std::vector<uint64_t> detectors = get_nonzero_detectors(syndrome);
        
        // graph::DecodingGraph dec_graph = circuit->decoding_graph;

        // uint n_vertices = detectors.size();

        // for (size_t i = 0; i < n_vertices; i++) {
        //     uint64_t di = detectors[i];
        //     uint64_t dj = detectors[j];

        //     auto vi = decoding_graph.get_vertex(di);
        //     auto vj = decoding_graph.get_vertex(dj);
        //     auto error_data = decoding_graph.get_error_chain_data(vi, vj);

        //     for (auto f : error_data.frame_changes) {
        //         std::cout<<f<<std::endl;
        //     }
        // }

        // std::cout<<"hi"<<std::endl;



        //do work here with the dec and the res too

        //decoding exp now!

        //stim circuit stored in circuit

        //decode_error made a res too

        // std::vector<uint64_t> flipped_dets = dec->get_nonzero_detectors(syndrome);

        

        // for (auto det : flipped_dets) {
        //     outFile<<det<<std::endl;
        // }

        

        //go thru syndrome to see orig

        


        // outFile<<res.corr<<std::endl;

        // auto error_matches = res.error_assignments;

        // for (auto pair : error_matches) {
        //     auto first = std::get<0>(pair);
        //     auto second = std::get<1>(pair);

        //     outFile << "first " << first << "second " << second << std::endl;
        // }

        // size_t total = res.corr.num_bits_padded();

        // for (int i = 0; i < total; i++) {
        //     if (res.corr[i] == 1) {
        //         outFile << "flipped at " << i << std::endl;
        //     }
        // }

        // outFile<<"pass?"<<std::endl;

        // outFile.close();

        // qontra::graph::DecodingGraph decGraph = qontra::graph::to_decoding_graph(circuit, qontra::graph::DecodingGraph::Mode::NORMAL);

        logical_errors[0] += (bool) (payload.observables != res.corr);
        t_sum += res.exec_time;
        t_sqr_sum += sqr(res.exec_time);
        t_max.scalar_replace_if_better_extrema(res.exec_time);

        for (size_t i = 0; i < n_obs; i++) {
            logical_errors[i+1] += (bool)(res.corr[i] != obs[i]);
        }
        e_cb(res);
    };

    uint64_t shots = config.shots;

    // std::cout<<shots<<std::endl;

    if (shots == 0) {
        shots = read_syndrome_trace(config.trace_folder, circuit, dec_cb);
    } else {
        generate_syndromes(circuit, shots, dec_cb);
    }

    // Collect results across all processors.
    logical_errors.reduce();

    hw_sum.reduce();
    hw_sqr_sum.reduce();
    hw_max.reduce();

    t_sum.reduce();
    t_sqr_sum.reduce();
    t_max.reduce();
    // Compute means and variances/std. deviations.

    statistic_t<fp_t> logical_error_rate = logical_errors.get_mean(shots);
    statistic_t<fp_t> hw_mean = hw_sum.get_mean(shots);
    statistic_t<fp_t> t_mean = t_sum.get_mean(shots);

    statistic_t<fp_t> hw_std = hw_sqr_sum.get_std(hw_mean, shots);
    statistic_t<fp_t> t_std = t_sqr_sum.get_std(t_mean, shots);

    std::vector<fp_t> logical_error_rate_by_obs = logical_error_rate.slice(1, n_obs+1).vec();

    memory_result_t res = {
        logical_error_rate.at(),
        hw_mean.at(),
        hw_std.at(),
        hw_max.at(),
        t_mean.at(),
        t_std.at(),
        t_max.at(),
        // Additional statistics:
        logical_error_rate_by_obs
    };
    return res;
}

//new code for shot by shot (sbs) for memory experiment with shot by shot decoding: above commented out

inline std::vector<memory_result_t>
memory_experiment_sbs(Decoder* dec, memory_config_t config) {
    return memory_experiment_sbs(dec, config, [] (shot_payload_t x) {}, [] (Decoder::result_t x) {});
}

template <class PROLOGUE, class EPILOGUE> std::vector<memory_result_t>
memory_experiment_sbs(Decoder* dec, memory_config_t config, PROLOGUE p_cb, EPILOGUE e_cb) {

    const DetailedStimCircuit circuit = dec->get_circuit();
    const size_t n_obs = circuit.count_observables();

    statistic_t<uint64_t> logical_errors(MPI_SUM, n_obs+1);
    statistic_t<uint64_t> hw_sum(MPI_SUM), hw_sqr_sum(MPI_SUM), hw_max(MPI_MAX);
    statistic_t<fp_t> t_sum(MPI_SUM), t_sqr_sum(MPI_SUM), t_max(MPI_MAX);

    std::vector<memory_result_t> results;

    // int ones = 0;
    // int zeroes = 0;

    auto dec_cb = [&] (shot_payload_t payload)
    {
        stim::simd_bits<SIMD_WIDTH> syndrome(payload.syndrome),
                                    obs(payload.observables);

        p_cb({syndrome, obs});
        const size_t hw = syndrome.popcnt();
        // Update HW statistics and skip the trial if the HW is too small
        // and filtering is enabled.
        hw_sum += hw;
        hw_sqr_sum += sqr(hw);
        hw_max.scalar_replace_if_better_extrema(hw);
        
        if (G_FILTER_OUT_SYNDROMES && hw <= G_FILTERING_HAMMING_WEIGHT) {
            return;
        }
        // Decode syndrome

        //all of these changed from += to just =
        
        auto res = dec->decode_error(syndrome); 
        logical_errors[0] = (bool) (payload.observables != res.corr);
        t_sum.fill(res.exec_time);
        t_sqr_sum.fill(res.exec_time);
        t_max.scalar_replace_if_better_extrema(res.exec_time);

        for (size_t i = 0; i < n_obs; i++) {
            logical_errors[i+1] = (bool)(res.corr[i] != obs[i]);
        }

        e_cb(res);

        //res creation: tweaked as it's only one shot

        // Collect results across all processors.

        logical_errors.reduce();
        hw_sum.reduce();
        hw_sqr_sum.reduce();
        hw_max.reduce();
        t_sum.reduce();
        t_sqr_sum.reduce();
        t_max.reduce();

        uint64_t shots = 1; //this line should take care of the statistics getting tweaked
        //because mean will be itself and therefore std dev will just be 0

        // Compute means and variances/std. deviations.
        statistic_t<fp_t> logical_error_rate = logical_errors.get_mean(shots);
        statistic_t<fp_t> hw_mean = hw_sum.get_mean(shots);
        statistic_t<fp_t> t_mean = t_sum.get_mean(shots);
        statistic_t<fp_t> hw_std = hw_sqr_sum.get_std(hw_mean, shots);
        statistic_t<fp_t> t_std = t_sqr_sum.get_std(t_mean, shots);

        std::vector<fp_t> logical_error_rate_by_obs = logical_error_rate.slice(1, n_obs+1).vec();
    
        // if (logical_errors[0] == 1) {
        //     ones++;
        // }
        // if (logical_errors[1] == 1) {
        //     ones++;
        // }
        // if (logical_errors[0] == 0) {
        //     zeroes++;
        // }
        // if (logical_errors[1] == 0) {
        //     zeroes++;
        // }
        // std::cout<<"ONES"<<std::endl;

        // std::cout<<ones<<std::endl;

        // std::cout<<"ZEROES"<<std::endl;
        // std::cout<<zeroes<<std::endl;

        memory_result_t indivResult = {
            logical_error_rate.at(),
            hw_mean.at(),
            hw_std.at(),
            hw_max.at(),
            t_mean.at(),
            t_std.at(),
            t_max.at(),
            // Additional statistics:
            logical_error_rate_by_obs
        };

        results.push_back(indivResult); //add the result to the results array
    };

    uint64_t shots = config.shots;

    // std::cout << shots << std::endl;

    int syndromeIndex = 0;

    if (shots == 0) {
        while (syndromeIndex < 3073) {
            shots = read_syndrome_trace_individually(config.trace_folder, circuit, dec_cb, syndromeIndex);
            syndromeIndex++;
        }
    } else {
        generate_syndromes(circuit, shots, dec_cb);
    }

    //moved the res creation from below into the callback function so we can create
    //an individual result object for each call to each syndrome

    // std::cout << "second" << std::endl;
    // std::cout << shots << std::endl;

    // Collect results across all processors.
    // logical_errors.reduce();

    // hw_sum.reduce();
    // hw_sqr_sum.reduce();
    // hw_max.reduce();

    // t_sum.reduce();
    // t_sqr_sum.reduce();
    // t_max.reduce();
    // // Compute means and variances/std. deviations.
    // statistic_t<fp_t> logical_error_rate = logical_errors.get_mean(shots);
    // statistic_t<fp_t> hw_mean = hw_sum.get_mean(shots);
    // statistic_t<fp_t> t_mean = t_sum.get_mean(shots);

    // statistic_t<fp_t> hw_std = hw_sqr_sum.get_std(hw_mean, shots);
    // statistic_t<fp_t> t_std = t_sqr_sum.get_std(t_mean, shots);

    // std::vector<fp_t> logical_error_rate_by_obs = logical_error_rate.slice(1, n_obs+1).vec();

    // memory_result_t res = {
    //     logical_error_rate.at(),
    //     hw_mean.at(),
    //     hw_std.at(),
    //     hw_max.at(),
    //     t_mean.at(),
    //     t_std.at(),
    //     t_max.at(),
    //     // Additional statistics:
    //     logical_error_rate_by_obs
    // };

    return results; //return array of results instead
}

}   // qontra
