// /*
//  *  author: Suhas Vittal
//  *  date:   30 January 2024
//  * */

// #include <qontra/decoder/pymatching.h>
// #include <qontra/experiments.h>
// #include <qontra/experiments/memory.h>
// #include <qontra/ext/stim.h>

// #include <vtils/filesystem.h>

// #include <iostream>

// #include <mpi.h>
// #include <stdio.h>

// #include <vector>

// using namespace qontra;
// using namespace vtils;

// int main(int argc, char* argv[]) {
//     std::string help = "usage: ./sse_memory <trace-folder> <stim-file> <output-file>";
//     if (argc < 3) {
//         std::cerr << help << std::endl;
//         return 1;
//     }
//     std::string trace_folder(argv[1]),
//                 stim_file(argv[2]),
//                 output_file(argv[3]);

//     MPI_Init(NULL, NULL);
//     int world_rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//     configure_optimal_batch_size();

//     FILE* fin = fopen(stim_file.c_str(), "r");
//     stim::Circuit error_model = stim::Circuit::from_file(fin);
//     fclose(fin);

//     PyMatching dec(error_model);

//     G_FILTER_OUT_SYNDROMES = false;
    
//     memory_config_t config = { 0, trace_folder };

//     std::vector<memory_result_t> results = memory_experiment_sbs(&dec, config); //changed to the sbs method, so we store in a vector instead

//     int ones = 0;
//     int zeroes = 0;

//     for (int i = 0; i < results.size(); i++) {
//         if (results.at(i).logical_error_rate == 0) {
//             zeroes++;
//         } else {
//             ones++;
//         }
//     }

//     std::cout<<"0:"<<std::endl;
//     std::cout<<zeroes<<std::endl;
//     std::cout<<"1:"<<std::endl;
//     std::cout<<ones<<std::endl;
//     double mean = ones/3073.00;

//     std::cout<<"LER:"<<std::endl;
//     std::cout<<mean<<std::endl;

//     //ITS EITHER 0 or 1 WHICH MAKES ENSE

//     // Write to output file.
//     bool write_header = false;
//     if (!file_exists(output_file)) {
//         safe_create_directory(get_parent_directory(output_file));
//         write_header = true;
//     }

//     if (world_rank == 0) {
//         std::ofstream fout(output_file, std::ios::app);
//         if (write_header) {
//             fout << "Trace Folder,Logical Error Rate" << std::endl;
//         }
//         // fout << get_basename(trace_folder) << ","
//         //     << res.logical_error_rate << std::endl;
//     }

//     MPI_Finalize();
//     return 0;
// }

/*
 *  author: Suhas Vittal
 *  date:   30 January 2024
 * */

#include <qontra/decoder/pymatching.h>
#include <qontra/experiments.h>
#include <qontra/experiments/memory.h>
#include <qontra/ext/stim.h>

#include <vtils/filesystem.h>

#include <iostream>

#include <mpi.h>
#include <stdio.h>

using namespace qontra;
using namespace vtils;

int main(int argc, char* argv[]) {
    std::string help = "usage: ./sse_memory <trace-folder> <stim-file> <output-file>";
    if (argc < 3) {
        std::cerr << help << std::endl;
        return 1;
    }
    std::string trace_folder(argv[1]),
                stim_file(argv[2]),
                output_file(argv[3]);

    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    configure_optimal_batch_size();

    FILE* fin = fopen(stim_file.c_str(), "r");
    stim::Circuit error_model = stim::Circuit::from_file(fin);
    fclose(fin);

    PyMatching dec(error_model);

    G_FILTER_OUT_SYNDROMES = false;
    
    memory_config_t config = { 0, trace_folder };
    memory_result_t res = memory_experiment(&dec, config);

    // Write to output file.
    bool write_header = false;
    if (!file_exists(output_file)) {
        safe_create_directory(get_parent_directory(output_file));
        write_header = true;
    }

    if (world_rank == 0) {
        std::ofstream fout(output_file, std::ios::app);
        if (write_header) {
            fout << "Trace Folder,Logical Error Rate" << std::endl;
        }
        fout << get_basename(trace_folder) << ","
            << res.logical_error_rate << std::endl;
    }

    MPI_Finalize();
    return 0;
}
