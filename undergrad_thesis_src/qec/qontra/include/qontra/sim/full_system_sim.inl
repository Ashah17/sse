/*
 *  author: Suhas Vittal
 *  date:   23 January 2024
 * */

#include "qontra/experiments.h"

#include <vtils/filesystem.h>

#include <qontra/decoder/pymatching.h>
#include <qontra/decoder/mwpm.h>
#include <qontra/experiments.h>
#include <qontra/experiments/memory.h>
#include <qontra/ext/stim.h>
#include <vtils/cmd_parse.h>

#include <fstream>
#include <iostream>
#include <cctype>

#include <string>
#include <regex>

#include <tuple>
#include <array>

#include <mpi.h>

#include <cstdlib>

#include <algorithm>

using namespace vtils;


namespace qontra {

inline size_t
get_register_index(std::string r) {
    const size_t REGISTER_SPECIAL_START = 32;
    // Check if the register is special, otherwise it is straightforward to get the index.
    if (r == "$rbrk") {
        return REGISTER_SPECIAL_START + 0;
    } else if (r == "$reoz") {
        return REGISTER_SPECIAL_START + 1;
    } else {
        int id = std::stoi(r.substr(2));
        return static_cast<size_t>(id);
    }
}

template <class SIM> histogram_t<uint64_t>
FullSystemSimulator::run_program(const qes::Program<>& program, uint64_t shots) {
    int world_rank = 0, world_size = 1;
    if (G_USE_MPI) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }
    // Execute program simulation in batches:
    const uint64_t local_shots = shots / world_size + static_cast<int>(world_rank==0)*(shots % world_size);

    // Set up all structures.
    is_recording_stim_instructions = true;
    n_qubits = get_number_of_qubits(program);
    base_sim = uptr<SIM>(new SIM(n_qubits, G_SHOTS_PER_BATCH));

    base_sim->set_seed(G_BASE_SEED + world_rank);

    register_file = stim::simd_bit_table<SIMD_WIDTH>(config.n_registers, G_SHOTS_PER_BATCH);
    syndrome_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);
    observable_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);

    shot_histogram.clear();

    // Create the stats files if they do not exist:
    using namespace vtils;
    if (!file_exists(config.syndrome_output_folder)) safe_create_directory(config.syndrome_output_folder);
    if (!file_exists(config.stim_output_file)) {
        safe_create_directory(get_parent_directory(config.stim_output_file.c_str()));
    }
    if (!file_exists(config.data_output_file)) {
        safe_create_directory(get_parent_directory(config.data_output_file.c_str()));
    }

    // uint64_t shots_remaining = local_shots;
    
    uint64_t shots_remaining = 1;
    uint64_t batchno = world_rank;

    // std::cout<<"shots IS " << shots_remaining<<std::endl;
    

    // int count = 0;

    while (shots_remaining) {
        const uint64_t shots_this_batch = shots_remaining < G_SHOTS_PER_BATCH ? shots_remaining : G_SHOTS_PER_BATCH;
        run_batch(program, shots_this_batch);
        write_stats(batchno);

        // std::cout<<batchno<<std::endl;

        //see syndrome table

        // if (syndrome_table[0].popcnt() > 0) {
        //     count++;
        // }

        // std::cout<<syndrome_table[0].str()<<std::endl;

        is_recording_stim_instructions = false;

        // shots_remaining -= shots_this_batch;
        shots_remaining -= 1;
        batchno += world_size;
    }
    // std::cout<<"pop count is "<< count << std::endl;


    // int count = 0;
    // int amount = 0;

    // for (int i = 0; i < G_RECORD_SPACE_SIZE; i++) {
    //     amount+=1;
    //     if (syndrome_table[i].popcnt() > 0) {
    //         count++;
    //         std::cout<<"i is "<< i<< "syndrome is" << syndrome_table[i].popcnt()<<std::endl;
    //     }
    // }
    // // std::cout<<amount<<std::endl;

    // std::cout<<count<<std::endl;


    shot_histogram = histogram_reduce(shot_histogram);
    return shot_histogram;

    
}

// template <class SIM> histogram_t<uint64_t>
// FullSystemSimulator::run_program(const qes::Program<>& program, uint64_t shots) {
//     int world_rank = 0, world_size = 1;
//     if (G_USE_MPI) {
//         MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//         MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//     }

//     // Execute program simulation in batches:

//     const uint64_t local_shots = shots / world_size + static_cast<int>(world_rank==0)*(shots % world_size);

//     // const uint64_t local_shots = 1

//     // Set up all structures.
//     is_recording_stim_instructions = true;
//     n_qubits = get_number_of_qubits(program);
//     base_sim = uptr<SIM>(new SIM(n_qubits, G_SHOTS_PER_BATCH));

//     // std::cout<< "WUBITS" << std::endl;
//     // std::cout<< n_qubits << std::endl;

//     std::cout<<"GOT wat above"<<std::endl;

//     base_sim->set_seed(G_BASE_SEED + world_rank);

//     std::cout<<"GOT way above"<<std::endl;

//     register_file = stim::simd_bit_table<SIMD_WIDTH>(config.n_registers, G_SHOTS_PER_BATCH);
//     syndrome_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);
//     observable_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);

//     shot_histogram.clear();

//     // Create the stats files if they do not exist:
//     using namespace vtils;
//     if (!file_exists(config.syndrome_output_folder)) safe_create_directory(config.syndrome_output_folder);
//     if (!file_exists(config.stim_output_file)) {
//         safe_create_directory(get_parent_directory(config.stim_output_file.c_str()));
//     }
//     if (!file_exists(config.data_output_file)) {
//         safe_create_directory(get_parent_directory(config.data_output_file.c_str()));
//     }

//     uint64_t shots_remaining = local_shots;
//     uint64_t batchno = world_rank;
//     std::cout<<"GOT above"<<std::endl;

//     while (shots_remaining) {
//         std::cout<<"GOT into loop??????"<<std::endl;
//         const uint64_t shots_this_batch = shots_remaining < G_SHOTS_PER_BATCH ? shots_remaining : G_SHOTS_PER_BATCH;
//         std::cout<<shots_remaining<<std::endl;
        
//         // const uint64_t shots_this_batch = 1; //1 shot per batch 
//         run_batch(program, shots_this_batch);

//         //wrote stim file for circuit, and will do so for each batch
//         //also wrote to the syndrome_table, so can decode

//         // stim::simd_bits_range_ref<SIMD_WIDTH> syndrome(syndrome_table, 1);

//         // std::string stimFile = config.stim_output_file; //stim file name: will be same every time

//         // std::string qesFile = "qes_circuit";

//         // std::string conv_comm = "./converter " + stimFile + " " + qesFile + " --p 1e-3";
//         // can't use a converter: it's only qes to stim, not the other way around

//         // DetailedStimCircuit circuit = sample_circuit; //sample_circuit is what is written to stim file

//         // // dec.decode_error(syndrome_table.to); //decode syndrome from syndrome table

//         // std::cout<<shots_remaining<<std::endl;
//         // std::cout<<"Made decoder"<<std::endl;

//         // auto syndrome = syndrome_table.transposed()[0]; //tranpose and take 0th - will give you shot data

//         // auto nonzero = get_nonzero_detectors_(syndrome, 124);

//         // if (syndrome.popcnt() > -5) {
//         //     //to save time we can skip hw 0 syndromes
//         //     std::cout<<"GOT into HERE"<<std::endl;
//         //     PyMatching dec(circuit); //make decoder for the circuit
//         //     auto res = dec.decode_error(syndrome);
//         //     if (res.corr[0] != observable_table.transposed()[0][0]) {
                
//         //     }
//         //     std::cout<<"detectors:" <<std::endl;

//         //     for (auto x : nonzero) {
//         //         std::cout << " " << x;
//         //     }

//         //     std::cout << std::endl;

//         //     std::cout<<batchno<<std::endl;

//         write_stats(batchno);

//         //     logical_errors += 1;
//         //     // logical_errors += (bool)(res.corr[0] != observable_table.transposed()[0][0]); //transpose and take 0th index again for obs
//         //     // if (res.corr[0] != observable_table.transposed()[0][0]) {
//         //     //     std::cout<<"logical error"<<std::endl;
//         //     // }
//         // }

//         //stop here and see if decoder still throws error

//         is_recording_stim_instructions = false; //leave this true so that a circuit is written for every batch, not just first

//         shots_remaining -= shots_this_batch;
//         batchno += world_size;

//         // std::ofstream ofs(stimFile, std::ofstream::out | std::ofstream::trunc); //truncation mode will clear file 
//         // ofs.close();
//         // std::ofstream ofs2(qesFile, std::ofstream::out | std::ofstream::trunc); //also clear qes file: need new ofs object tho
//         // ofs2.close();
//         std::cout<<"GOT HERE"<<std::endl;
//     }

//     // logical_errors.reduce();

//     // std::cout<<"LERs"<<std::endl;

//     // double ler = logical_errors.get_mean(shots).at(); //mean over all shots

//     // std::cout<<ler<<std::endl;

//     shot_histogram = histogram_reduce(shot_histogram);
//     return shot_histogram;
// }

inline void
FullSystemSimulator::load_subroutine(std::string name, const qes::Program<>& program) {
    subroutine_map[name] = program;
}

inline void
FullSystemSimulator::execute_routine(const qes::Program<>& program) {
    program_status_t status(current_shots);

    // std::cout<<"HIII"<<std::endl;
    // std::cout<<program.size()<<std::endl;

    while (status.pc < program.size()) {
        const qes::Instruction<>& instruction = program.at(status.pc);

        // std::cout<<instruction.get_name()<<std::endl;
        // std::cout<< qes:Instruction::print_instr(instruction, True) << std::endl;
        // Check if the instruction is requesting the execution of the
        // microcode.
        if (instruction.get_name() == "call") {
            // Switch into the subroutine.
            std::string subroutine_name = instruction.get<std::string>(0);
            execute_routine(subroutine_map.at(subroutine_name));
            status.pc++;
        } else {
            // std::cout<<"hitttinggg"<<std::endl;
            read_next_instruction(program, status);
            // if (status.return_if_waiting_trials[0]) {
            //     std::cout<<"breaking"<<std::endl;
            //     break; //
            // }
        }

    }
}

//**ADDING THIS**//

//map syntax

// std::map<char, char> my_map = {
//     {'A', '1'},
//     {...}
// }

// my_map["a"] = 101;

//loop start to end with .begin() and .end(), .first for key, .second for value
//.size() for end, .erase("key") to remove key
//.insert({key, element}) to add a pair


//so let's do the same process as we did with the python generation, but now in c++
//insight is to have partial measurements, so we can improve the code by adding a second parameter


//after this function is written, next step is to make sure integration is there with rest of methods
//main integration is with read_next_instruction method: ensure that the program (instruction vector) is built correctly


//STEPS:

//EXTO0:
//1) hadamard the x paritys
//2) create a qubit checks map: parity to the data qubits it's checking


//3) z parity boundary: first (n - 1)/2 are z boundary
//4) even ones: A is (zParity[counter] - n**2 + 1, zParity[counter] - n**2), C is (zParity[counter] - n + 1, zParity[counter] - n)
//5) odd ones: B is (zParity[counter] - n**2 + 1, zParity[counter] - n**2), D is (zParity[counter] - n + 1, zParity[counter] - n)


//6) x parity boundary: next (n - 1)/2 are x boundary
//7) even ones: C is (2*n - 1 + adder, n - 1 + adder), D is (2*n + adder, n + adder), adder += 2n
//8) odd ones: A is (2*n - 1 + adder, n - 1 + adder), B is (2*n + adder, n + adder), adder += 2n

//9) regular parity: for inner, do all in one loop


//creating a function to build the map

inline std::map<int, std::vector<int>> 
FullSystemSimulator::create_p_to_d_map(int d) {
    // int total = pow(d, 2) + pow (d, 2) - 1;
    int zStart = pow(d, 2); //take next (d-1) for z boundary
    int xStart = pow(d, 2) + (d-1); //take next (d-1) for x bound
    int innerStart = pow(d, 2) + 2*(d-1);

    // std::cout<< "INNER SATTR" << std::endl;
    // std::cout<< innerStart << std::endl;

    std::map<int, std::vector<int>> map = std::map<int, std::vector<int>>(); //use for parity to data for checks

    //zBoundary: first d-1

    // for (int i = zStart; i < zStart + (d-1); i++) {
    //     //new if statement
    //     std::vector<int> withData = std::vector<int>();
    //     //odds are on left
    //     if (i%2 == 1) {
    //         withData.push_back((i - pow(d, 2) + 1));
    //         withData.push_back((i - pow(d, 2)));
    //     } else {
    //         withData.push_back((i - d + 1));
    //         withData.push_back((i - d));
    //     }
        
    //     map[i] = withData;
    // }

    //xBoundary: next d-1

    // int adder = 0;

    // for (int i = xStart; i < xStart + (d-1); i++) {
    //     std::vector<int> withData = std::vector<int>();
    //     //odds are on bottom
    //     if (i%2 == 1) {
    //         withData.push_back((2*d - 1 + adder));
    //         withData.push_back((d - 1 + adder));
    //     } else {
    //         withData.push_back((2*d + adder));
    //         withData.push_back((d + adder));
    //         adder += 2*d; //the gap logic, even will always happen second
    //     }
    //     map[i] = withData;
    // }

    //NOW THE Z NEEDS ADDER LOGIC N X HAS PLAIN LOGIC

    int adder = 0;

    for (int i = zStart; i < zStart + (d-1); i++) {
        std::vector<int> withData = std::vector<int>();
        //odds are on left
        if (i%2 == 1) {
            withData.push_back((adder));
            withData.push_back((d + adder));
            adder += 2*d; //the gap logic, even will always happen second
        } else {
            withData.push_back((adder - 1));
            withData.push_back((d - 1 + adder));
        }
        map[i] = withData;
    }

    //cool x gud

    adder = 1; //use this but as a "subtracter" and adder

    for (int i = xStart; i < xStart + (d-1); i++) {
        //new if statement
        std::vector<int> withData = std::vector<int>();
        //odds are on top
        if (i%2 == 1) {
            withData.push_back(d - adder);
            withData.push_back(d - adder - 1);
        } else {
            withData.push_back((d-1)*(d) + adder);
            withData.push_back((d-1)*(d) + adder - 1);
            adder+=2;
        }
        
        map[i] = withData;
    }

    //cool z gud


    //INNER can stay the same - just make the order BADC for X and BDAC for Z

    //inner now, boundaries are fixed

    for (int i = 0; i < d - 1; i++) {
        //for every row

        //didn't need if either: it was redunant in old code, just using ranges fixes it

        //switching order to DCBA for x, DBCA for z
        
        //the DBCA and DCBA order works out: alternates b/w x, z row by row

        int start = 0; //starting of first square
        int end = d+1; //end of first square

        for (int j = 0; j < (d-1)/2; j++) {
            std::vector<int> withData = std::vector<int>();

            if (i%2 == 0) {
                //even row, order is x then z

                withData.push_back(int(start + 1 + (i*d))); //B
                withData.push_back(int(start + (i*d))); //A
                withData.push_back(int(start + end + (i*d))); //D
                withData.push_back(int(start + end - 1 + (i*d))); //C
            

                //+(i*d) for the row shift
                
                map[innerStart] = withData;

                innerStart++; //next one
                withData.clear(); //to make z

                start++; //shift square 1

                withData.push_back(int(start + 1 + (i*d))); //B

                withData.push_back(int(start + end + (i*d))); //D

                withData.push_back(int(start + (i*d))); //A
                
                withData.push_back(int(start + end - 1 + (i*d))); //C
                
                
                map[innerStart] = withData;

                start++; //shift square

                // withData.push_back(int(start2 + adder + i));
                // withData.push_back(int(start2 + adder - 1 + i));
                // withData.push_back(int(start2 + adder - d + i));
                // withData.push_back(int(start2 + adder - d - 1 + i));

                
                // map[innerStart] = withData;

                // innerStart++; //next one
                // withData.clear(); //to make z

                // withData.push_back(int(start1 + adder + i));
                // withData.push_back(int(start1 + adder - d + i));
                // withData.push_back(int(start1 + adder - 1 + i));
                // withData.push_back(int(start1 + adder - d - 1 + i));

                
                
                // map[innerStart] = withData;
                // adder2 += 2;

                innerStart++; //next
            } else {
                withData.push_back(int(start + 1 + (i*d))); //B
                withData.push_back(int(start + end + (i*d))); //D
                withData.push_back(int(start + (i*d))); //A
                
                withData.push_back(int(start + end - 1 + (i*d))); //C
                
                
                
                
                map[innerStart] = withData;

                innerStart++; //next one
                withData.clear(); //to make x
                start++;

                withData.push_back(int(start + 1 + (i*d))); //B
                withData.push_back(int(start + (i*d))); //A

                withData.push_back(int(start + end + (i*d))); //D
                withData.push_back(int(start + end - 1 + (i*d))); //C
                
                
                
                
                
                map[innerStart] = withData;
                // adder2 += 2;
                start++;

                innerStart++; //next


                // withData.push_back(int(start2 + adder + i));
                // withData.push_back(int(start2 + adder - d + i));
                // withData.push_back(int(start2 + adder - 1 + i));
                // withData.push_back(int(start2 + adder - d - 1 + i));

                // map[innerStart] = withData;

                // innerStart++; //next one
                // withData.clear(); //to make z

                // withData.push_back(int(start1 + adder + i));
                // withData.push_back(int(start1 + adder - 1 + i));
                // withData.push_back(int(start1 + adder - d + i));
                // withData.push_back(int(start1 + adder - d - 1 + i));

                
                // map[innerStart] = withData;
            }
    
            adder += 2*d; //gap for next cell in row
        }
    }

    //all parity done

    std::ofstream outFile("themap.txt");

    for (auto parity : map) {
        outFile << "KEY: " << parity.first << "--> " << std::endl;
        for (auto data : parity.second) {
            outFile << data << std::endl;
        }
    }

    outFile.close();


    return map; //keys are all parity qubits, values are the data qubits they need to cnot with

    //when doing within function, z should be data to parity, x should be parity to data (order for control, target)

}

inline qes::Program<>
FullSystemSimulator::create_program(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs) {

    std::vector<int> selectors;
    //program is a vector of instructions: simple end goal
    
    std::vector<qes::Instruction<>> program; //append to this
    qes::Instruction<> curr; //build this instruction every time: don't need the default types
    std::vector<int> operands; //use this for every instruction, clear after using

    //reset all qubits

    int total_qubits = pow(d, 2) + pow(d, 2) - 1;

    for (int i = 0; i < total_qubits; i++) {
        operands.push_back(i);
    }

    std::set<int> zStabilizers = std::set<int>();

    std::set<int> xStabilizers = std::set<int>();
    std::set<int> innerStabilizers = std::set<int>();
    std::vector<int> firstRound; //save first round

    std::map<int, std::vector<int>> p_to_d = create_p_to_d_map(d); //use for parity to data for checks


    for (int dRound = 0; dRound < 1; dRound++) {
        // d outer rounds: keep it at 1 now, should be d+1? But causing seg faults

        int numEvents = 0;
        int newNumEvents = 0; //var to reset the numevents for ext1
        // std::vector<int> lastRoundChecks; //don't use this anymore

        std::map<int, int> latestMeasurement; //map for <stabilizer:latest measurement number>
        std::map<int, int> previousMeasurement; //map for <stabilizer:previous measurement number>



        curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);
        program.push_back(curr);
        //std::cout<<curr<<std::endl;
        operands.clear(); //reset all

        //go through map: the key is the parity, the value is the data qubits

        //10 inner passes fixed: first 1 is ext0, next 9 are ext1 - only difference is the events

        for (int j = pow(d, 2); j < pow(d, 2) + (d - 1); j++) {
            zStabilizers.insert(j);
        }


        for (int j = pow(d, 2) + (d-1); j < pow(d, 2) + 2*(d-1); j++) {
            xStabilizers.insert(j);
        }

        // for (auto i : zStabilizers) {
        //     // if (innerStabilizers.contains(i)) {
        //         std::cout<<i<<std::endl;
        //     // }  
        // }

        // std::cout<<"break"<<std::endl;

        // for (auto i : xStabilizers) {
        //     // if (innerStabilizers.contains(i)) {
        //         std::cout<<i<<std::endl;
        //     // }  
        // }

        int rowC = 0; //keep track of this: divide by d and cast to int to get the row num

        for (int j = pow(d, 2) + 2*(d-1); j < total_qubits; j++) {
            int rowNum = (int)(rowC/(d-1)); //changed logic

            innerStabilizers.insert(j);

            //pattern: if (rowNum + j)%2 == 0, then it is a z stabilizer, else it is x

            if ((rowNum + j)%2 == 0) {
                //z stabilizer
                zStabilizers.insert(j);
            } else {
                xStabilizers.insert(j);
            }

            rowC++;
        }
        
        // for (auto i : zStabilizers) {
        //     if (innerStabilizers.contains(i)) {
        //         std::cout<<i<<std::endl;
        //     }  
        // }
        

        // std::cout<<"break2"<<std::endl;

        // for (auto i : xStabilizers) {
        //     if (innerStabilizers.contains(i)) {
        //         std::cout<<i<<std::endl;
        //     }
        // }

        // for (auto x : zStabilizers) {
        //     std::cout<<x<<std::endl;
        // }




        

        //zStabs populated

        int numMeasurements = 0;

        int selectiveOrControl = 1; //0 means control, 1 means selective

        std::random_device rand; //rand
        std::mt19937 gen(rand()); //seed
        std::uniform_real_distribution<> dist(0.0,1.0); //0 to 1 dist
        
        for (int round = 0; round < 10*d; round++) {

            //set a flag here by checking the round: check every 3 rounds -> round%d == d-1
            
            selectiveOrControl = 1;

            if ((round%d) == (d-1)) {
                //check the dict

                // if (stab_freqs.find(round) != stab_freqs.end()) {
                //     //if key in: as default behav of c++ map is to add default value
                //     double numFlipped = stab_freqs[round];

                //     double freq = numFlipped/4096; //turn into prob 

                //     std::ofstream outFile5("theRands.txt", std::ios::app);

                double randomNum = dist(gen); //create random num from dist

                    // outFile5 <<"RAND : " << randomNum<< " round : " << round <<std::endl;

                if (randomNum > 0.5) {//then it is error, so do control
                    selectiveOrControl = 1; //set to 0 for next 3 rounds till recheck
                    // std::cout<<"HII"<<std::endl;
                } else {
                    selectiveOrControl = 1; //else 1
                    // std::cout<<"buhhhh"<<std::endl;
                }
                    // selectors.push_back(selectiveOrControl);
                // }

                // std::cout<<round<<std::endl;
                // std::cout<<randomNum<<std::endl;
                // std::cout<<freq<<std::endl;
                // std::cout<<selectiveOrControl<<std::endl;
                // std::cout<<"pass"<std::endl;
            }

            // if (round%(1*d) == 0) {
            //     //set to 0 every 2d rounds (to measure alls)
            //     selectiveOrControl = 0;
            // }

            int selected = stabilizers_each_round[round].size(); //num stabilizers to measure: key = round

            //the 10 syndrome extraction passes

            //h the x stabilizers

            for (int xStab : xStabilizers) {
                operands.push_back(xStab);
            }

            curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);
            curr.put("timing_error"); //with the timing error
            program.push_back(curr);
            //std::cout<<curr<<std::endl;
            operands.clear(); //reset all

            //CNOT EVENTS

            //adding selective cnot logic

            int i = 0;
            

            //CHANGED J2 LOGIC
            int j1 = pow(d, 2);
            int j2 = (pow(d, 2) + (d/2) * 4) - 1;
            

            // selectiveOrControl = 1;

            while (i < 4) {        

                //x
                // if (i==2) {
                //     //reset the j2 and j1 to get the other pairing for the last 2 lines
                //     j1 = pow(d, 2);
                //     j2 = (pow(d, 2) + (d/2) * 4) - 1;
                //     second = 1; //
                // }

                if (i==2) {
                    j1++;
                    j2--;
                }

                std::vector<int> values;

                if (round==0 || selectiveOrControl == 0) {
                    for (int x = 0; x < (d/2); x++) {
                        operands.push_back(j2);
                        values = p_to_d.at(j2); //the values
                        operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                        j2-=2; //do d/2 each time
                    }
                    // if (i!=1) {
                    j2+=2*(d/2);
                    // }
                } else {
                    if (stabilizers_each_round[round].contains(j2)) {
                        for (int x = 0; x < (d/2); x++) {
                        operands.push_back(j2);
                        values = p_to_d.at(j2); //the values
                        operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                        j2-=2; //do d/2 each time
                    }
                    // if (i!=1) {
                    j2+=2*(d/2);
                    }
                }

                //z

                if (round==0 || selectiveOrControl == 0) {
                    //1st round is normal cnots
                    //ALSO if selectiveOrControl is 0 then do all
                    for (int x = 0; x < (d/2); x++) {
                        values = p_to_d.at(j1); //the values
                        operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                        operands.push_back(j1);
                        j1+=2;
                    }
                    // if (i!=1) {
                    j1-=(2*(d/2));
                    // }
                } else {
                    if (stabilizers_each_round[round].contains(j1)) {
                        for (int x = 0; x < (d/2); x++) {
                        values = p_to_d.at(j1); //the values
                        operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                        operands.push_back(j1);
                        j1+=2;
                    }
                    // if (i!=1) {
                    j1-=(2*(d/2));
                    }
                }
                
                //inner below

                for (int j = pow(d, 2) + 2*(d-1); j < total_qubits; j++) {
                    if (round==0 || selectiveOrControl == 0) {
                        if (zStabilizers.contains(j)) {
                            operands.push_back(p_to_d.at(j)[i%4]);
                            operands.push_back(j);
                        } else {
                            operands.push_back(j);
                            operands.push_back(p_to_d.at(j)[i%4]);
                        }
                    } else {
                        if (stabilizers_each_round[round].contains(j)) {
                            if (zStabilizers.contains(j)) {
                                operands.push_back(p_to_d.at(j)[i%4]);
                                operands.push_back(j);
                            } else {
                                operands.push_back(j);
                                operands.push_back(p_to_d.at(j)[i%4]);
                            }
                        } 
                    }
                    
                }

                //operands line done

                // for (auto op : operands) {
                //     std::cout<<op<<std::endl;
                // }

                // std::cout<<"pass"<<std::endl;

                curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);
                program.push_back(curr); //instruction made and added to program
                //std::cout<<curr<<std::endl;

                i++;

                operands.clear(); //clear for next loop
            }

            //the cnots are done

            //h the x stabs again

            for (int xStab : xStabilizers) {
                operands.push_back(xStab);
            }

            curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);
            program.push_back(curr);
            //std::cout<<curr<<std::endl;
            operands.clear(); //reset all

            //measure stabilizers and reset

            //only do below if rounds < numRounds

            if (round==0 || selectiveOrControl == 0) {
                //measure all stabilizers in first round

                //OR IF CONTROL BIT IS OFF

                previousMeasurement = latestMeasurement; //set it before changing latest

                for (int i = pow(d, 2); i < total_qubits; i++) {
                    operands.push_back(i);
                    latestMeasurement[i] = numMeasurements; //set the latest measurement to this num of measurement
                    numMeasurements++;
                }

                //add to the map
            } else {
                //or else dif measurement patterns
                previousMeasurement = latestMeasurement; //set it before changing latest
                // for (int i = pow(d, 2); i < total_qubits; i++) {
                //     operands.push_back(i);
                //     latestMeasurement[i] = numMeasurements;
                //     numMeasurements++;
                // }

                for (int i : stabilizers_each_round[round]) {
                    operands.push_back(i);
                    latestMeasurement[i] = numMeasurements;
                    numMeasurements++;
                }
            }
            
            // for (auto pair : latestMeasurement) {
            //     std::cout<< pair.first << "-->" << pair.second << std::endl;
            // }
            

            curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);
            program.push_back(curr); //added measurements
            //std::cout<<curr<<std::endl;
            operands.clear(); //for resets

            for (int i = pow(d, 2); i < total_qubits; i++) {
                operands.push_back(i);
            }

            curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);
            program.push_back(curr); //added resets: all stabilizers
            //std::cout<<curr<<std::endl;
            operands.clear();
            

            //events now

            if (round == 0) {
                //initial ext0 events

                //all events only for z stabilizer measurements
                for (int i : zStabilizers) {
                    operands.push_back(numEvents);
                    operands.push_back(latestMeasurement[i]); //push the latest z stabilizer measurement
                    // std::cout << i << "--> " << latestMeasurement[i] << std::endl;
                    curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
                    program.push_back(curr);
                    //std::cout<<curr<<std::endl;
                    operands.clear();
                    numEvents++;

                    //save these for epilogue events
                    firstRound.push_back(latestMeasurement[i]);
                }

            } else {
                //rest of rounds
                std::set<int> allStabs = std::set<int>();
                allStabs.insert(xStabilizers.begin(), xStabilizers.end());
                allStabs.insert(zStabilizers.begin(), zStabilizers.end());

                // for (auto c : allStabs) {
                //     std::cout<<c<<std::endl;
                // }

                // std::cout<<"pass"<<std::endl;
                for (int i : allStabs) {
                    operands.push_back(numEvents);
                    operands.push_back(previousMeasurement[i]); //u want the previous measurement too
                    operands.push_back(latestMeasurement[i]);
                    
                    curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
                    program.push_back(curr);

                    operands.clear();
                    numEvents++;
                }
            }
        }

        //everything else done: epilogue now

        //measure data qubits

        std::vector<int> observable;

        for (int i = 0; i < pow(d, 2); i++) {
            operands.push_back(i);

            if (i < d) {
                observable.push_back(numMeasurements); //for final obs
            }
            // if (numMeasurements%d == (d-1)) {
                
            // }
            numMeasurements++;
        }

        // for (int i = pow(d,2)-d; i < pow(d,2);i++) {
        //     observable.push_back(numMeasurements); //for final obs
        //     numMeasurements++;
        // }

        curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);
        program.push_back(curr);
        //std::cout<<curr<<std::endl;
        operands.clear();

        //final events

        int index = 0;

        //x memory exp so we care about z stabilizers 

        std::vector<int> save;

        int adder = 1;

        bool allDoneOnce = false;

        //ok so subtract the previous d^2 data q measurements then just add from map each time

        numMeasurements -= pow(d,2); 

        for (int z : zStabilizers) {
            
            if (z < pow(d, 2) + 2*(d-1)) {
                //that means it's a boundary z, so only 2
                // operands.push_back(index); //remove the eoffset
                // operands.push_back(firstRound[index]); //remove the mshift
                operands.push_back(numEvents);
                operands.push_back(latestMeasurement[z]); //moved this from below push back num events

                
                std::vector<int> values = p_to_d[z];
                
                // operands.push_back(values[1] + numMeasurements - 1);
                // operands.push_back(values[0] + numMeasurements - 1);

                operands.push_back(numMeasurements + values[1]);
                operands.push_back(numMeasurements + values[0]);
                // operands.push_back(numMeasurements - 1 - values[2]);
                // operands.push_back(numMeasurements - 1 - values[0]);


            } else {
                //4 values
                // operands.push_back(index);
                // operands.push_back(firstRound[index]); //removed shifts
                operands.push_back(numEvents);
                operands.push_back(latestMeasurement[z]);//moved
                

                std::vector<int> values = p_to_d.at(z); //the values

                operands.push_back(numMeasurements + values[0]);
                operands.push_back(numMeasurements + values[2]);
                operands.push_back(numMeasurements + values[1]);
                operands.push_back(numMeasurements + values[3]);

            }
            // std::cout << "LINE" << std::endl;
            curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
            program.push_back(curr);
            //std::cout<<curr<<std::endl;
            
            index++;
            numEvents++;

            operands.clear();
        }

        //events done

        //obs

        operands.push_back(0);

        //making the bottom logic more generic: it's the newest data qubit measures

        // operands.push_back(242);
        // operands.push_back(245);
        // operands.push_back(248);

        // for (int obs : observable) {
        //     operands.push_back(obs);
        // }

        for (int i = 0; i < d; i++) {
            operands.push_back(numMeasurements + i);
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("obs", operands);
        program.push_back(curr);
        ////std::cout<<curr<<std::endl;
        operands.clear();

        // operands.push_back(1);

        // for (int i = d; i < 2*d; i++) {
        //     operands.push_back(numMeasurements + i);
        // }

        // curr = qes::Instruction<qes::any_t, qes::any_t>("obs", operands);
        // program.push_back(curr);
        // ////std::cout<<curr<<std::endl;
        // operands.clear();

        // operands.push_back(2);

        // for (int i = 2*d; i < 3*d; i++) {
        //     operands.push_back(numMeasurements + i);
        // }

        // curr = qes::Instruction<qes::any_t, qes::any_t>("obs", operands);
        // program.push_back(curr);
        // ////std::cout<<curr<<std::endl;
        // operands.clear();

        // for (int i = d+1; i < 2*d; i++) {
        //     operands.push_back(numMeasurements + i);
        // }

        // for (int i =  2*d; i < 3*d; i++) {
        //     operands.push_back(numMeasurements + i);
        // }

        

        //tried with multiple obs too

        

        // int obsCount = 0;

        // for (int obs : observable) {
        //     operands.push_back(obsCount++);
        //     operands.push_back(obs);
        //     curr = qes::Instruction<qes::any_t, qes::any_t>("obs", operands);
        //     program.push_back(curr);
        //     ////std::cout<<curr<<std::endl;
        //     operands.clear();
        // }

        //above was interesting multiple obs
    }
    
    //done now go back and inject timing errors

    //DONE (except annotation)

    // for (int i = 0; i < SIMD_WIDTH; i++) {
    //     std::cout<<syndrome_table[i].str()<<std::endl;
    // }

    // std::cout<<syndrome_table.data<<std::endl;

    if (stab_freqs.size() != 0) {
        std::ofstream outFile("instructions.txt"); //saving to file to analyze

        // for (const auto& pair : stab_freqs) {
        //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
        // }

        for (auto instr : program) {
            outFile << instr << std::endl;
        }

        // for (auto sel : selectors) {
        //     outFile <<sel<<std::endl;
        // }

        outFile.close();
    } else {
        std::ofstream outFile2("final_project_surface_code_old_code.txt"); //saving to file to analyze

        // for (const auto& pair : stab_freqs) {
        //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
        // }

        for (auto instr : program) {
            outFile2 << instr << std::endl;
        }

        // for (auto sel : selectors) {
        //     outFile <<sel<<std::endl;
        // }
 
        outFile2.close();
    }
    

    // int count = 0;

    // for (int i = 0; i < G_RECORD_SPACE_SIZE; i++) {
    //     if (syndrome_table[i].popcnt() > 0) {
    //         count++;
    //     }
    //     // std::cout<<syndrome_table[i].popcnt()<<std::endl;
    // }

    // std::cout<<count<<std::endl;

    

    // for (int i = 1; i < 33; i++) {
    //     std::string regnum = "$r";
    //     regnum.append(std::to_string(i));
    //     std::cout<< regnum <<std::endl;
    //     stim::simd_bits_range_ref<SIMD_WIDTH> r = get_register(regnum);
    //     std::cout<< r <<std::endl;
    // }


    return program;
}

// -- code below is to generate stability exp, first map then circuit --

// inline std::map<int, std::vector<int>>
// FullSystemSimulator::create_p_to_d_map_for_stability_exp(int d) {
//     // int total = pow(d, 2) + pow (d, 2) - 1;
//     // int zStart = pow(d, 2); //take next (d-1) for z boundary
//     // int xStart = pow(d, 2) + (d-1); //take next (d-1) for x bound
//     // int innerStart = pow(d, 2) + 2*(d-1);

//     int zStart = pow(d,2); //all outer r z
//     int innerStart = pow(d,2)+2*d; //so inner starts with x after 2d outer z qubits


//     std::map<int, std::vector<int>> map = std::map<int, std::vector<int>>(); //use for parity to data for checks

//     //NOW THE Z NEEDS ADDER LOGIC N X HAS PLAIN LOGIC

//     // int adder = 0;

//     // for (int i = zStart; i < zStart + (d-1); i++) {
//     //     std::vector<int> withData = std::vector<int>();
//     //     //odds are on left
//     //     if (i%2 == 1) {
//     //         withData.push_back((adder));
//     //         withData.push_back((d + adder));
//     //         adder += 2*d; //the gap logic, even will always happen second
//     //     } else {
//     //         withData.push_back((adder - 1));
//     //         withData.push_back((d - 1 + adder));
//     //     }
//     //     map[i] = withData;
//     // }



//     //***DOING FOR A D2 CODE*** NOT SCALABLE **

//     //the outer z logic
//     int adder = 0;

//     for (int i = zStart; i < zStart + 2*d; i++) {
//         std::vector<int> withData = std::vector<int>();

//         //first half r on sides

//         if (i < zStart + d) {
//             if (i%2==0) {
//                 //left side
//                 withData.push_back((adder));
//                 withData.push_back((adder+d));
//             } else {
//                 //right side
//                 withData.push_back((adder+1));
//                 withData.push_back((adder+d+1));
//                 //now bouncing back to left so jump
//             }
//         } else {
//             //next half r on top/bottoms
//             if (i%2==0) {
//                 //top
//                 withData.push_back((adder));
//                 withData.push_back((adder+1));
//             } else {
//                 withData.push_back((d));
//                 withData.push_back((d+1));
//             }
//         }

//         map[i] = withData;
//         withData.clear();
//     }

//     //outer z gud

//     // adder = 1; //use this but as a "subtracter" and adder

//     // for (int i = xStart; i < xStart + (d-1); i++) {
//     //     //new if statement
//     //     std::vector<int> withData = std::vector<int>();
//     //     //odds are on top
//     //     if (i%2 == 1) {
//     //         withData.push_back(d - adder);
//     //         withData.push_back(d - adder - 1);
//     //     } else {
//     //         withData.push_back((d-1)*(d) + adder);
//     //         withData.push_back((d-1)*(d) + adder - 1);
//     //         adder+=2;
//     //     }
        
//     //     map[i] = withData;
//     // }

//     // //cool z gud


//     // //INNER can stay the same - just make the order BADC for X and BDAC for Z

//     // //inner now, boundaries are fixed

//     // for (int i = 0; i < d - 1; i++) {
//     //     //for every row

//     //     //didn't need if either: it was redunant in old code, just using ranges fixes it

//     //     //switching order to DCBA for x, DBCA for z
        
//     //     //the DBCA and DCBA order works out: alternates b/w x, z row by row

//     //     int start = 0; //starting of first square
//     //     int end = d+1; //end of first square

//     //     for (int j = 0; j < (d-1)/2; j++) {
//     //         std::vector<int> withData = std::vector<int>();

//     //         if (i%2 == 0) {
//     //             //even row, order is x then z

//     //             withData.push_back(int(start + 1 + (i*d))); //B
//     //             withData.push_back(int(start + (i*d))); //A
//     //             withData.push_back(int(start + end + (i*d))); //D
//     //             withData.push_back(int(start + end - 1 + (i*d))); //C
            

//     //             //+(i*d) for the row shift
                
//     //             map[innerStart] = withData;

//     //             innerStart++; //next one
//     //             withData.clear(); //to make z

//     //             start++; //shift square 1

//     //             withData.push_back(int(start + 1 + (i*d))); //B

//     //             withData.push_back(int(start + end + (i*d))); //D

//     //             withData.push_back(int(start + (i*d))); //A
                
//     //             withData.push_back(int(start + end - 1 + (i*d))); //C
                
                
//     //             map[innerStart] = withData;

//     //             start++; //shift square

//     //             innerStart++; //next
//     //         } else {
//     //             withData.push_back(int(start + 1 + (i*d))); //B
//     //             withData.push_back(int(start + end + (i*d))); //D
//     //             withData.push_back(int(start + (i*d))); //A
                
//     //             withData.push_back(int(start + end - 1 + (i*d))); //C
                
                
                
                
//     //             map[innerStart] = withData;

//     //             innerStart++; //next one
//     //             withData.clear(); //to make x
//     //             start++;

//     //             withData.push_back(int(start + 1 + (i*d))); //B
//     //             withData.push_back(int(start + (i*d))); //A

//     //             withData.push_back(int(start + end + (i*d))); //D
//     //             withData.push_back(int(start + end - 1 + (i*d))); //C
    
                
                
//     //             map[innerStart] = withData;
//     //             // adder2 += 2;
//     //             start++;

//     //             innerStart++; //next

//     //         }
    
//     //         adder += 2*d; //gap for next cell in row
//     //     }
//     // }

//     //inner is just the 1 x

//     std::vector<int> withData = std::vector<int>();

//     for (int i = 0; i < 4; i++) {
//         withData.push_back(i);
//     }

//     map[8] = withData;
    

//     //all parity done

//     std::ofstream outFile("themap_stability.txt");

//     for (auto parity : map) {
//         outFile << "KEY: " << parity.first << "--> " << std::endl;
//         for (auto data : parity.second) {
//             outFile << data << std::endl;
//         }
//     }

//     outFile.close();


//     return map; //keys are all parity qubits, values are the data qubits they need to cnot with

//     //when doing within function, z should be data to parity, x should be parity to data (order for control, target)

// }

// inline qes::Program<>
// FullSystemSimulator::create_program_stability_exp(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs) {
//     std::map<int, std::vector<int>> p_to_d = create_p_to_d_map_for_stability_exp(d); //use for parity to data for checks

//     std::vector<int> selectors;
//     //program is a vector of instructions: simple end goal
    
//     std::vector<qes::Instruction<>> program; //append to this
//     qes::Instruction<> curr; //build this instruction every time: don't need the default types
//     std::vector<int> operands; //use this for every instruction, clear after using

//     //reset all qubits

//     int totalQubits = 2*pow(d,2)+1; //2d^2 qubits, just add 1 cuz doing < logic

//     std::set<int> zStabilizers = std::set<int>();

//     std::set<int> xStabilizers = std::set<int>();
//     std::set<int> allStabilizers = std::set<int>();
//     std::vector<int> firstRound; //save first round

//     int numEvents = 0;
//     int numMeasurements = 0;

//     std::map<int, int> latestMeasurement; //map for <stabilizer:latest measurement number>
//     std::map<int, int> previousMeasurement; //map for <stabilizer:previous measurement number>

//     for (int i = 0; i < totalQubits; i++) {
//         operands.push_back(i); //all qubits
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);

//     program.push_back(curr); //reset all

//     operands.clear();

//     for (int i = 4; i < 9; i ++) {
//         allStabilizers.insert(i);
//     }

//     xStabilizers.insert(8);

//     for (int i = 4; i < 8; i ++) {
//         zStabilizers.insert(i);
//     }

//     //stabilizer sets populated

//     //h all stabs

//     // for (auto i : allStabilizers) {
//     //     operands.push_back(i);
//     // }

//     for (int i = 1; i < totalQubits; i++) {
//         operands.push_back(i);
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//     program.push_back(curr);

//     operands.clear();

//     //cnots

//     //1

//     operands.push_back(8);
//     operands.push_back(0);
//     operands.push_back(2);
//     operands.push_back(7);
//     operands.push_back(1);
//     operands.push_back(5);

//     curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//     program.push_back(curr);

//     operands.clear();

//     for (int i = 0; i < 2; i++) {
//         operands.push_back(i);
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//     program.push_back(curr);

//     operands.clear();

//     //2

//     operands.push_back(0);
//     operands.push_back(4);
//     operands.push_back(8);
//     operands.push_back(2);
//     operands.push_back(3);
//     operands.push_back(7);

//     curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//     program.push_back(curr);

//     operands.clear();

//     for (int i = 1; i < 3; i++) {
//         operands.push_back(i);
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//     program.push_back(curr);

//     operands.clear();

//     //3

//     operands.push_back(0);
//     operands.push_back(6);
//     operands.push_back(8);
//     operands.push_back(1);
//     operands.push_back(3);
//     operands.push_back(5);

//     curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//     program.push_back(curr);

//     operands.clear();

//     for (int i = 2; i < 4; i++) {
//         operands.push_back(i);
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//     program.push_back(curr);

//     operands.clear();

//     //4

//     operands.push_back(2);
//     operands.push_back(4);
//     operands.push_back(1);
//     operands.push_back(6);
//     operands.push_back(8);
//     operands.push_back(3);

//     curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//     program.push_back(curr);

//     operands.clear();

//     for (int i = 3; i < totalQubits; i++) {
//         operands.push_back(i);
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//     program.push_back(curr);

//     operands.clear();

//     //measure all stabs

//     previousMeasurement = latestMeasurement; //set it before changing latest

//     for (auto i : allStabilizers) {
//         operands.push_back(i);
//         latestMeasurement[i] = numMeasurements;
//         numMeasurements++;
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);

//     program.push_back(curr);
//     operands.clear();

//     //event for x stab

//     operands.push_back(numEvents);
//     operands.push_back(5);
//     curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
//     program.push_back(curr);
//     operands.clear();
//     numEvents++;

//     //NO RESETTING

//     //NOW SYNDROME EXTRACT, do 10 rounds

//     for (int i = 0; i < 10; i++) {
//         //reset all stabs

//         for (auto i : allStabilizers) {
//             operands.push_back(i);
//         }

//         curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);

//         program.push_back(curr);

//         operands.clear();

//         //h all stabs

//         // for (auto i : allStabilizers) {
//         //     operands.push_back(i);
//         // }
//         // for (int i = 0; i < totalQubits; i++) {
//         //     operands.push_back(i);
//         // }
//         operands.push_back(0);
//         operands.push_back(4);
//         operands.push_back(5);
//         operands.push_back(6);
//         operands.push_back(7);
//         operands.push_back(8);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//         //cnots

//         //1

//         operands.push_back(8);
//         operands.push_back(0);
//         operands.push_back(2);
//         operands.push_back(7);
//         operands.push_back(1);
//         operands.push_back(5);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//         program.push_back(curr);

//         operands.clear();

//         operands.push_back(0);
//         operands.push_back(1);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//         //2

//         operands.push_back(0);
//         operands.push_back(4);
//         operands.push_back(8);
//         operands.push_back(2);
//         operands.push_back(3);
//         operands.push_back(7);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//         program.push_back(curr);

//         operands.clear();

//         operands.push_back(1);
//         operands.push_back(2);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//         //3

//         operands.push_back(0);
//         operands.push_back(6);
//         operands.push_back(8);
//         operands.push_back(1);
//         operands.push_back(3);
//         operands.push_back(5);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//         program.push_back(curr);

//         operands.clear();

//         operands.push_back(2);
//         operands.push_back(3);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//         //4

//         operands.push_back(2);
//         operands.push_back(4);
//         operands.push_back(1);
//         operands.push_back(6);
//         operands.push_back(8);
//         operands.push_back(3);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//         program.push_back(curr);

//         operands.clear();

//         operands.push_back(3);
//         operands.push_back(4);
//         operands.push_back(5);
//         operands.push_back(6);
//         operands.push_back(7);
//         operands.push_back(8);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//         //cnots done

//         //measure all stabs

//         previousMeasurement = latestMeasurement; //set it before changing latest

//         for (auto i : allStabilizers) {
//             operands.push_back(i);
//             latestMeasurement[i] = numMeasurements;
//             numMeasurements++;
//         }

//         curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);

//         program.push_back(curr);
//         operands.clear();

//         //events for all stabs

//         for (auto i : allStabilizers) {
//             operands.push_back(numEvents);
//             operands.push_back(previousMeasurement[i]); //got set cuz of prologue and then curr round
//             operands.push_back(latestMeasurement[i]);

//             curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
//             program.push_back(curr);
//             operands.clear();
//             numEvents++;
//         }
//     }

//     //shud be done, now epilogue, same as one more round
//     //but measure all at end

//     for (auto i : allStabilizers) {
//         operands.push_back(i);
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);

//     program.push_back(curr);

//     operands.clear();

//     //h all stabs

//     // for (auto i : allStabilizers) {
//     //     operands.push_back(i);
//     // }

//     operands.push_back(0);
//         operands.push_back(4);
//         operands.push_back(5);
//         operands.push_back(6);
//         operands.push_back(7);
//         operands.push_back(8);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//         //cnots

//         //1

//         operands.push_back(8);
//         operands.push_back(0);
//         operands.push_back(2);
//         operands.push_back(7);
//         operands.push_back(1);
//         operands.push_back(5);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//         program.push_back(curr);

//         operands.clear();

//         operands.push_back(0);
//         operands.push_back(1);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//         //2

//         operands.push_back(0);
//         operands.push_back(4);
//         operands.push_back(8);
//         operands.push_back(2);
//         operands.push_back(3);
//         operands.push_back(7);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//         program.push_back(curr);

//         operands.clear();

//         operands.push_back(1);
//         operands.push_back(2);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//         //3

//         operands.push_back(0);
//         operands.push_back(6);
//         operands.push_back(8);
//         operands.push_back(1);
//         operands.push_back(3);
//         operands.push_back(5);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//         program.push_back(curr);

//         operands.clear();

//         operands.push_back(2);
//         operands.push_back(3);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//         //4

//         operands.push_back(2);
//         operands.push_back(4);
//         operands.push_back(1);
//         operands.push_back(6);
//         operands.push_back(8);
//         operands.push_back(3);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

//         program.push_back(curr);

//         operands.clear();

//         operands.push_back(0);
//         operands.push_back(1);
//         operands.push_back(2);
//         operands.push_back(4);
//         operands.push_back(5);
//         operands.push_back(6);
//         operands.push_back(7);
//         operands.push_back(8);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

//         program.push_back(curr);

//         operands.clear();

//     //cnots done

//     //measure ALL QUBITS - this the dif for epilogue

//     previousMeasurement = latestMeasurement; //set it before changing latest

//     for (int i = 0; i < totalQubits; i++) {
//         operands.push_back(i);
//         latestMeasurement[i] = numMeasurements;
//         numMeasurements++;
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);

//     program.push_back(curr);
//     operands.clear();

//     //events for all stabs, logic shud be same


//     for (auto i : allStabilizers) {
//         operands.push_back(numEvents);
//         operands.push_back(previousMeasurement[i]); //got set cuz of prologue and then curr round
//         operands.push_back(latestMeasurement[i]);

//         curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
//         program.push_back(curr);
//         operands.clear();
//         numEvents++;
//     }
    

//     //event for the final x stabilizer (with data qubits measure)

//     operands.push_back(numEvents);
//     operands.push_back(latestMeasurement[6]); //latestMeasurement for x stabilizer (6)
    
//     for (int i = 0; i < 4; i++) {
//         //latest measurements for data qubits
//         operands.push_back(latestMeasurement[i]);
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
//     program.push_back(curr);
//     operands.clear();
//     numEvents++;

//     //obs is z stabilizer measures

//     operands.push_back(0);

//     for (auto i : zStabilizers) {
//         operands.push_back(latestMeasurement[i]);
//     }

//     curr = qes::Instruction<qes::any_t, qes::any_t>("obs", operands);
//     program.push_back(curr);
//     operands.clear();

//     std::ofstream outFile2("instructions_stability_experiment.txt"); //saving to file to analyze

//         // for (const auto& pair : stab_freqs) {
//         //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
//         // }

//     for (auto instr : program) {
//         outFile2 << instr << std::endl;
//     }

//     // for (auto sel : selectors) {
//     //     outFile <<sel<<std::endl;
//     // }

//     outFile2.close();

//     //done

//     return program;
// }



// inline std::map<int, double> FullSystemSimulator::read_dets(std::string pathfake){
//     //takes in path to read the dets files

//     std::map<int, double> stab_freqs;

//     //define regex pattern to match "D##" and extract info

//     for (int i = 0; i < 4096; i++) {
//         // std::string path = "/Users/aryan/Desktop/quantumResearch/qec/qontra/build/temp_trace/";
//         std::string path = "/Users/aryan/Desktop/quantumResearch/qec/qontra/build/temp_trace/";
//         path.append("batch_");
//         path.append(std::to_string(i));
//         path.append(".dets");

//         // std::cout<<path<<std::endl;
        
//         std::ifstream file(path); //open with param path to file
        
//         std::string line;

//         std::set<int> rounds_flipped = std::set<int>(); //keep set of rounds flipped, clear per shot

//         while (std::getline(file, line)) {
//             //go thru each line

//             std::vector<std::string> terms; //to add to from curr line

//             std::istringstream stream(line);
//             std::string word;

//             while (stream >> word) {
//                 terms.push_back(word);
//             }

//             //now terms has the word shot n then flipped stabilizers so remove shot

//             terms.erase(terms.begin());

//             for (int j = 0; j < terms.size(); j++) {
//                 std::string curr = terms[j];
                
//                 curr.erase(0,1); //remove the D

//                 // std::cout<<curr<<std::endl;

//                 //now add curr to dict

//                 int currNum = std::stoi(curr);

//                 currNum /= 4; //int divide by 4 to make key the round num

//                 ///add currNum to set: no dupes

//                 rounds_flipped.insert(currNum); //insert into set: no dupes

//                 //below is code for # of stabs flipped

//                 // auto key = stab_freqs.find(currNum); //find the key, returns ref to iterator

//                 // int val = 0; //default val

//                 // if (key != stab_freqs.end()) {
//                 //     //if not then keeps the default 0 count, else get the value
//                 //     val = key->second; //deref the iterator for val
//                 // }

//                 // val += 1;

//                 // stab_freqs[currNum] = val; //set the val
//             } 
//         }

//         //below is code for binary IF stab is flipped

//         for (auto round : rounds_flipped) {
//             auto key = stab_freqs.find(round); //find the key, returns ref to iterator

//             int val = 0; //default val

//             if (key != stab_freqs.end()) {
//                 //if not then keeps the default 0 count, else get the value
//                 val = key->second; //deref the iterator for val
//             }

//             val += 1;

//             stab_freqs[round] = val; //set the val
//         }

//         rounds_flipped.clear();

//     }


//     //reading properly, now feed into map

//     return stab_freqs;
// }

inline std::map<int, std::vector<int>>
FullSystemSimulator::create_p_to_d_map_for_stability_exp(int d) {

    //hard coding map for d4 stability exp :/

    std::set<int> zOuter = {16,17,18,22,26,30,31,32};
    std::set<int> x = {19,21,24,27,29};
    std::set<int> zInner = {20,23,25,28};

    std::map<int, std::vector<int>> map = std::map<int, std::vector<int>>(); //use for parity to data for checks

    std::vector<int> withData = std::vector<int>();

    //16
    withData.push_back(0);
    withData.push_back(1);
    map[16] = withData;
    withData.clear();

    //17
    withData.push_back(2);
    withData.push_back(3);
    map[17] = withData;
    withData.clear();

    //18
    withData.push_back(0);
    withData.push_back(4);
    map[18] = withData;
    withData.clear();

    //19
    withData.push_back(4);
    withData.push_back(0);
    withData.push_back(5);
    withData.push_back(1);
    map[19] = withData;
    withData.clear();

    //20
    withData.push_back(5);
    withData.push_back(6);
    withData.push_back(1);
    withData.push_back(2);
    map[20] = withData;
    withData.clear();

    //21
    withData.push_back(6);
    withData.push_back(2);
    withData.push_back(7);
    withData.push_back(3);
    map[21] = withData;
    withData.clear();

    //22
    withData.push_back(3);
    withData.push_back(7);
    map[22] = withData;
    withData.clear();

    //23
    withData.push_back(8);
    withData.push_back(9);
    withData.push_back(4);
    withData.push_back(5);
    map[23] = withData;
    withData.clear();

    //24
    withData.push_back(9);
    withData.push_back(5);
    withData.push_back(10);
    withData.push_back(6);
    map[24] = withData;
    withData.clear();

    //25
    withData.push_back(10);
    withData.push_back(11);
    withData.push_back(6);
    withData.push_back(7);
    map[25] = withData;
    withData.clear();

    //26
    withData.push_back(8);
    withData.push_back(12);
    map[26] = withData;
    withData.clear();

    //27
    withData.push_back(12);
    withData.push_back(8);
    withData.push_back(13);
    withData.push_back(9);
    map[27] = withData;
    withData.clear();

    //28
    withData.push_back(13);
    withData.push_back(14);
    withData.push_back(9);
    withData.push_back(10);
    map[28] = withData;
    withData.clear();

    //29
    withData.push_back(14);
    withData.push_back(10);
    withData.push_back(15);
    withData.push_back(11);
    map[29] = withData;
    withData.clear();

    //30
    withData.push_back(11);
    withData.push_back(15);
    map[30] = withData;
    withData.clear();

    //31
    withData.push_back(12);
    withData.push_back(13);
    map[31] = withData;
    withData.clear();

    //32
    withData.push_back(14);
    withData.push_back(15);
    map[32] = withData;
    withData.clear();

    //done for d4

    std::ofstream outFile("stability_d4_map.txt");

    for (auto parity : map) {
        outFile << "KEY: " << parity.first << "--> " << std::endl;
        for (auto data : parity.second) {
            outFile << data << std::endl;
        }
    }

    outFile.close();


    return map; //keys are all parity qubits, values are the data qubits they need to cnot with

    //when doing within function, z should be data to parity, x should be parity to data (order for control, target)
}

inline qes::Program<>
FullSystemSimulator::create_program_stability_exp(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs, int selectiveOrControl) {
    std::map<int, std::vector<int>> p_to_d = create_p_to_d_map_for_stability_exp(d); //use for parity to data for checks

    //program is a vector of instructions: simple end goal
    
    std::vector<qes::Instruction<>> program; //append to this
    qes::Instruction<> curr; //build this instruction every time: don't need the default types
    std::vector<int> operands; //use this for every instruction, clear after using

    //reset all qubits

    int totalQubits = 2*pow(d,2)+1; //2d^2 qubits, just add 1 cuz doing < logic

    std::set<int> zOuter = {16,17,18,22,26,30,31,32};
    std::set<int> x = {19,21,24,27,29}; //only inner
    std::set<int> zInner = {20,23,25,28};
    std::set<int> zStabilizers = zOuter;
    zStabilizers.insert(zInner.begin(), zInner.end()); //all zstabs
    std::set<int> allStabilizers = zStabilizers;
    allStabilizers.insert(x.begin(), x.end()); //all stabs

    //stabilizer sets populated

    //using sets for h gates as they're weird - 5 sets

    std::set<int> h1 = {1,3,4,6,7,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    std::set<int> h2 = {0,1,2,3,5,6,8,9,10,11};
    std::set<int> h3 = {1,3,4,11,12,14};
    std::set<int> h4 = {4,5,6,7,9,10,12,13,14,15};
    std::set<int> h5 = {5,7,10,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};

    //syndrome extract and epilogue have different h1, epilogue has different h5

    std::set<int> h1dif = {0,2,5,8,10,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    std::set<int> h5dif = {0,1,2,3,4,6,8,9,11,12,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};

    int numEvents = 0;
    int numMeasurements = 0;

    std::map<int, int> latestMeasurement; //map for <stabilizer:latest measurement number>
    std::map<int, int> previousMeasurement; //map for <stabilizer:previous measurement number>

    //reset all

    for (int i = 0; i < totalQubits; i++) {
        operands.push_back(i); //all qubits
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);

    program.push_back(curr); 

    operands.clear();
    

    //h1

    for (auto q : h1) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //CNOT EVENTS - prologue so do all

    std::map<int, int> cnot_check_indices = std::map<int, int>(); //init all to 0

    for (int i = 16; i < 33; i++) {
        cnot_check_indices[i] = 0;
    }

    //cuz map is correct order, so go thru the indices
    
    //cnot_line1
    for (int i = 19; i < 33; i++) {
        if (i != 26) {
            //skip 26
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);

    operands.clear();

    //h2

    for (auto q : h2) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //cnot_line2
    for (int i = 16; i < 31; i++) {
        if (i != 18 && i != 26) {
            //skip 18, 26
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);

    operands.clear();

    //h3

    for (auto q : h3) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //cnot_line3
    for (int i = 18; i < 33; i++) {
        if (i != 22 && i != 30) {
            //skip 22, 30
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);

    operands.clear();

    //h4

    for (auto q : h4) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //cnot_line4
    for (int i = 16; i < 30; i++) {
        if (i != 22) {
            //skip 22
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);

    operands.clear();

    //h5

    for (auto q : h5) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //measure all stabs

    previousMeasurement = latestMeasurement; //set it before changing latest

    for (auto i : allStabilizers) {
        operands.push_back(i);
        latestMeasurement[i] = numMeasurements;
        numMeasurements++;
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);

    program.push_back(curr);
    operands.clear();

    //event for x stabs

    for (auto xQ : x) {
        operands.push_back(numEvents);
        operands.push_back(latestMeasurement[xQ]);
        curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
        program.push_back(curr);
        operands.clear();
        numEvents++;
    }
    

    //NO RESETTING

    //NOW SYNDROME EXTRACT, do 10 rounds

    for (int round = 0; round < 10; round++) {
        //reset all stabs

        for (auto i : allStabilizers) {
            operands.push_back(i);
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);

        program.push_back(curr);

        operands.clear();

        //h1 (DIF ONE)

        for (auto q : h1dif) {
            operands.push_back(q);
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

        program.push_back(curr);

        operands.clear();

        //CNOT EVENTS - CHECK IF IN THE SELECTIVE SET CUZ SYNDROME EXTRACT

        std::map<int, int> cnot_check_indices = std::map<int, int>(); //init all to 0

        for (int i = 16; i < 33; i++) {
            cnot_check_indices[i] = 0;
        }

        //cuz map is correct order, so j go thru the indices
        
        //cnot_line1
        for (int i = 19; i < 33; i++) {
            if (i != 26) {
                //skip 26
                //do the correct index of all from map
                //SELECTIVE LOGIC
                if (stabilizers_each_round[round].contains(i)) {
                    int index = cnot_check_indices[i];
                    //if z then value first then i
                    if (zStabilizers.contains(i)) {
                        operands.push_back(p_to_d[i][index]);
                        operands.push_back(i);
                    } else {
                        operands.push_back(i);
                        operands.push_back(p_to_d[i][index]);
                    }
                    cnot_check_indices[i]++; //increment for next time
                }
            }
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

        program.push_back(curr);

        operands.clear();

        //h2

        for (auto q : h2) {
            operands.push_back(q);
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

        program.push_back(curr);

        operands.clear();

        //cnot_line2
        for (int i = 16; i < 31; i++) {
            if (i != 18 && i != 26) {
                //skip 18, 26
                //do the correct index of all from map
                //SELECTIVE LOGIC
                if (stabilizers_each_round[round].contains(i)) {
                    int index = cnot_check_indices[i];
                    //if z then value first then i
                    if (zStabilizers.contains(i)) {
                        operands.push_back(p_to_d[i][index]);
                        operands.push_back(i);
                    } else {
                        operands.push_back(i);
                        operands.push_back(p_to_d[i][index]);
                    }
                    cnot_check_indices[i]++; //increment for next time
                }
            }
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

        program.push_back(curr);

        operands.clear();

        //h3

        for (auto q : h3) {
            operands.push_back(q);
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

        program.push_back(curr);

        operands.clear();

        //cnot_line3
        for (int i = 18; i < 33; i++) {
            if (i != 22 && i != 30) {
                //skip 22, 30
                //do the correct index of all from map
                //SELECTIVE LOGIC
                if (stabilizers_each_round[round].contains(i)) {
                    int index = cnot_check_indices[i];
                    //if z then value first then i
                    if (zStabilizers.contains(i)) {
                        operands.push_back(p_to_d[i][index]);
                        operands.push_back(i);
                    } else {
                        operands.push_back(i);
                        operands.push_back(p_to_d[i][index]);
                    }
                    cnot_check_indices[i]++; //increment for next time
                }
            }
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

        program.push_back(curr);

        operands.clear();

        //h4

        for (auto q : h4) {
            operands.push_back(q);
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

        program.push_back(curr);

        operands.clear();

        //cnot_line4
        for (int i = 16; i < 30; i++) {
            if (i != 22) {
                //skip 22
                //do the correct index of all from map
                //SELECTIVE LOGIC
                if (stabilizers_each_round[round].contains(i)) {
                    int index = cnot_check_indices[i];
                    //if z then value first then i
                    if (zStabilizers.contains(i)) {
                        operands.push_back(p_to_d[i][index]);
                        operands.push_back(i);
                    } else {
                        operands.push_back(i);
                        operands.push_back(p_to_d[i][index]);
                    }
                    cnot_check_indices[i]++; //increment for next time
                }
            }
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

        program.push_back(curr);

        operands.clear();

        //h5

        for (auto q : h5) {
            operands.push_back(q);
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

        program.push_back(curr);

        operands.clear();

        //measure all stabs ONLY IF IN SELECTIVE ROUND

        previousMeasurement = latestMeasurement; //set it before changing latest

        for (auto i : allStabilizers) {
            if (stabilizers_each_round[round].contains(i)) {
                operands.push_back(i);
                latestMeasurement[i] = numMeasurements;
                numMeasurements++;
            }
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);

        program.push_back(curr);
        operands.clear();

        //events for all stabs - logic stays same

        for (auto i : allStabilizers) {
            operands.push_back(numEvents);
            operands.push_back(previousMeasurement[i]); //got set cuz of prologue and then curr round
            operands.push_back(latestMeasurement[i]);

            curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
            program.push_back(curr);
            operands.clear();
            numEvents++;
        }
    }

    //syndrome extract shud be done
    //now epilogue, same as one more round
    //but measure all at end

    //NOT SELECTIVE - JUST LIKE PROLOGUE

    for (auto i : allStabilizers) {
        operands.push_back(i);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);

    program.push_back(curr);

    operands.clear();

    //h1 (dif one)

    for (auto q : h1dif) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //CNOT EVENTS - epilogue so do all

    // std::map<int, int> cnot_check_indices = std::map<int, int>(); //init all to 0

    for (int i = 16; i < 33; i++) {
        cnot_check_indices[i] = 0;
    }

    //cuz map is correct order, so go thru the indices
    
    //cnot_line1
    for (int i = 19; i < 33; i++) {
        if (i != 26) {
            //skip 26
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);

    operands.clear();

    //h2

    for (auto q : h2) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //cnot_line2
    for (int i = 16; i < 31; i++) {
        if (i != 18 && i != 26) {
            //skip 18, 26
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);

    operands.clear();

    //h3

    for (auto q : h3) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //cnot_line3
    for (int i = 18; i < 33; i++) {
        if (i != 22 && i != 30) {
            //skip 22, 30
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);

    operands.clear();

    //h4

    for (auto q : h4) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //cnot_line4
    for (int i = 16; i < 30; i++) {
        if (i != 22) {
            //skip 22
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);

    operands.clear();

    //h5 (DIF FOR EPILOGUE)

    for (auto q : h5dif) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);

    operands.clear();

    //measure ALL QUBITS - this the dif for epilogue

    previousMeasurement = latestMeasurement; //set it before changing latest

    for (int i = 0; i < totalQubits; i++) {
        operands.push_back(i);
        latestMeasurement[i] = numMeasurements;
        numMeasurements++;
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);

    program.push_back(curr);
    operands.clear();

    //events for all stabs, logic shud be same

    for (auto i : allStabilizers) {
        operands.push_back(numEvents);
        operands.push_back(previousMeasurement[i]); //got set cuz of prologue and then curr round
        operands.push_back(latestMeasurement[i]);

        curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
        program.push_back(curr);
        operands.clear();
        numEvents++;
    }
    

    //events for x stabilizers and their data qubits (only inners)

    for (auto xQ : x) {
        operands.push_back(numEvents);
        operands.push_back(latestMeasurement[xQ]);
        for (auto dataWith : p_to_d[xQ]) {
            operands.push_back(latestMeasurement[dataWith]); //those got stored in map too!
        }
        curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
        program.push_back(curr);
        operands.clear();
        numEvents++;
    }

    //obs is z stabilizer measures now

    operands.push_back(0);

    for (auto i : zStabilizers) {
        operands.push_back(latestMeasurement[i]);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("obs", operands);
    program.push_back(curr);
    operands.clear();

    std::ofstream outFile2("d4_stability_experiment.txt"); //saving to file to analyze

        // for (const auto& pair : stab_freqs) {
        //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
        // }

    for (auto instr : program) {
        outFile2 << instr << std::endl;
    }

    // for (auto sel : selectors) {
    //     outFile <<sel<<std::endl;
    // }

    outFile2.close();

    //done

    return program;
}



inline stim::simd_bits_range_ref<SIMD_WIDTH>
FullSystemSimulator::get_register(std::string r) {
    // Handle certain registers differently. If it is a special register, we may need to
    // update the value here.
    if (r == "$reoz") {
        // Set $reoz to the inverse of the bitwise OR of the 
        // last config.sse.rreoz_track_n_det events.
        size_t k = get_register_index("$reoz");
        register_file[k].clear();
        for (uint64_t e : sse.rreoz_events) {
            register_file[k] |= syndrome_table[e];
            // std::cout<<syndrome_table[e][0];
        }
        register_file[k].invert_bits();
        // std::cout<<" " << register_file[k][0] << std::endl;
    }
    return register_file[get_register_index(r)];
}

inline void
FullSystemSimulator::recalibrate_timing() {
    // Compute min time when taking account delta.
    fp_t min_time_delta = 0.0;
    for (const auto& p : shot_time_delta_map) {
        min_time_delta = std::min(min_time_delta, p.second);
    }
    if (min_time_delta >= 0.0) return;
    // Otherwise, change everything up.
    std::map<uint64_t, fp_t> new_time_delta_map;
    for (uint64_t t = 0; t < current_shots; t++) {
        if (shot_time_delta_map.count(t)) {
            if (min_time_delta != shot_time_delta_map.at(t)) {
                new_time_delta_map[t] = shot_time_delta_map.at(t) - min_time_delta;
            }
        } else {
            new_time_delta_map[t] = -min_time_delta;
        }
    }
    elapsed_time += min_time_delta;
    shot_time_delta_map = std::move(new_time_delta_map);
}

inline void
FullSystemSimulator::snapshot() {
    base_sim->snapshot();

    register_file_cpy = register_file;
    syndrome_table_cpy = syndrome_table;
    observable_table_cpy = observable_table;
}

inline void
FullSystemSimulator::rollback_where(stim::simd_bits_range_ref<SIMD_WIDTH> pred) {
    base_sim->rollback_where(pred);
    copy_where(register_file_cpy, register_file, pred);
    copy_where(syndrome_table_cpy, syndrome_table, pred);
    copy_where(observable_table_cpy, observable_table, pred);
}

// template <class SIM> histogram_t<uint64_t>
// FullSystemSimulator::run_program(const qes::Program<>& program, uint64_t shots) {
//     int world_rank = 0, world_size = 1;
//     if (G_USE_MPI) {
//         MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//         MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//     }
//     // Execute program simulation in batches:

//     // const uint64_t local_shots = shots / world_size + static_cast<int>(world_rank==0)*(shots % world_size);

//     const uint64_t local_shots = 1;

//     G_SHOTS_PER_BATCH = 1; //set this to 1, as we have 1 shot per batch

//     // Set up all structures.
//     is_recording_stim_instructions = true;
//     n_qubits = get_number_of_qubits(program);
//     base_sim = uptr<SIM>(new SIM(n_qubits, G_SHOTS_PER_BATCH));

//     base_sim->set_seed(G_BASE_SEED + world_rank);

//     register_file = stim::simd_bit_table<SIMD_WIDTH>(config.n_registers, G_SHOTS_PER_BATCH);
//     syndrome_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);
//     observable_table = stim::simd_bit_table<SIMD_WIDTH>(G_RECORD_SPACE_SIZE, G_SHOTS_PER_BATCH);

//     shot_histogram.clear();

//     // Create the stats files if they do not exist:
//     using namespace vtils;
//     if (!file_exists(config.syndrome_output_folder)) safe_create_directory(config.syndrome_output_folder);
//     if (!file_exists(config.stim_output_file)) {
//         safe_create_directory(get_parent_directory(config.stim_output_file.c_str()));
//     }
//     if (!file_exists(config.data_output_file)) {
//         safe_create_directory(get_parent_directory(config.data_output_file.c_str()));
//     }

//     // uint64_t shots_remaining = local_shots;
//     // uint64_t batchno = world_rank;

//     uint64_t shot_num = 0; //loop thru shots rather than batches
//     uint64_t    total_shots = 4096; //hardcoded values
//     fp_t        p = 1e-3;
//     std::string data_file = config.data_output_file;

//     //we make circuit from the program passed in once we execute it - not qes file

//     while (shot_num < 4096) {
//         const uint64_t shots_this_batch = 1; //1 shot per batch

//         run_batch(program, shots_this_batch); //run program with 1 shot per batch
//         write_stats(shot_num); //write syndrome to file so that it can be read and decoded

//         DetailedStimCircuit circuit = make_circuit(program, p); //create circuit
//         // std::cout<<shot_num<<std::endl;

//         qontra::MWPMDecoder dec(circuit);

//         auto res = dec.decode_error(syndrome_table[0]);

//         //use the decoder_error method which takes in syndrome
        
//         //write the logical error rate from each res to a file
//         //for now just print

//         shot_num++;
//         syndrome_table.clear(); //clear everytime
//     }
        
//     shot_histogram = histogram_reduce(shot_histogram);
//     return shot_histogram;
// }



// -- creating incremental stability functions here -- //

inline std::tuple<qes::Program<>, qes::Program<>, int, int, std::map<int, int>, std::map<int, int>, int>
FullSystemSimulator::stability_syndrome_prologue(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs, int numEvents, int numMeasurements, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, int roundNum, qes::Program<> program){
    
    std::map<int, std::vector<int>> p_to_d = create_p_to_d_map_for_stability_exp(d); //use for parity to data for checks

    //program is a vector of instructions: simple end goal
    
    // std::vector<qes::Instruction<>> program; //append to this
    qes::Instruction<> curr; //build this instruction every time: don't need the default types
    std::vector<int> operands; //use this for every instruction, clear after using
    std::vector<qes::Instruction<>> curr_program; //CURR PROGRAM

    //reset all qubits

    int totalQubits = 2*pow(d,2)+1; //2d^2 qubits, just add 1 cuz doing < logic

    std::set<int> zOuter = {16,17,18,22,26,30,31,32};
    std::set<int> x = {19,21,24,27,29}; //only inner
    std::set<int> zInner = {20,23,25,28};
    std::set<int> zStabilizers = zOuter;
    zStabilizers.insert(zInner.begin(), zInner.end()); //all zstabs
    std::set<int> allStabilizers = zStabilizers;
    allStabilizers.insert(x.begin(), x.end()); //all stabs

    //stabilizer sets populated

    //using sets for h gates as they're weird - 5 sets

    std::set<int> h1 = {1,3,4,6,7,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    std::set<int> h2 = {0,1,2,3,5,6,8,9,10,11};
    std::set<int> h3 = {1,3,4,11,12,14};
    std::set<int> h4 = {4,5,6,7,9,10,12,13,14,15};
    std::set<int> h5 = {5,7,10,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};

    //syndrome extract and epilogue have different h1, epilogue has different h5

    std::set<int> h1dif = {0,2,5,8,10,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    std::set<int> h5dif = {0,1,2,3,4,6,8,9,11,12,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};


    //using these from params
    // int numEvents = 0;
    // int numMeasurements = 0;

    // std::map<int, int> latestMeasurement; //map for <stabilizer:latest measurement number>
    // std::map<int, int> previousMeasurement; //map for <stabilizer:previous measurement number>

    //reset all

    for (int i = 0; i < totalQubits; i++) {
        operands.push_back(i); //all qubits
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();
    

    //h1

    for (auto q : h1) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //CNOT EVENTS - prologue so do all

    std::map<int, int> cnot_check_indices = std::map<int, int>(); //init all to 0

    for (int i = 16; i < 33; i++) {
        cnot_check_indices[i] = 0;
    }

    //cuz map is correct order, so go thru the indices
    
    //cnot_line1
    for (int i = 19; i < 33; i++) {
        if (i != 26) {
            //skip 26
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h2

    for (auto q : h2) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //cnot_line2
    for (int i = 16; i < 31; i++) {
        if (i != 18 && i != 26) {
            //skip 18, 26
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h3

    for (auto q : h3) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //cnot_line3
    for (int i = 18; i < 33; i++) {
        if (i != 22 && i != 30) {
            //skip 22, 30
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h4

    for (auto q : h4) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //cnot_line4
    for (int i = 16; i < 30; i++) {
        if (i != 22) {
            //skip 22
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h5

    for (auto q : h5) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //measure all stabs

    previousMeasurement = latestMeasurement; //set it before changing latest

    for (auto i : allStabilizers) {
        operands.push_back(i);
        latestMeasurement[i] = numMeasurements;
        numMeasurements++;
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);

    program.push_back(curr);
    curr_program.push_back(curr);
    operands.clear();

    //event for x stabs

    for (auto xQ : x) {
        operands.push_back(numEvents);
        operands.push_back(latestMeasurement[xQ]);
        curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
        program.push_back(curr);
        curr_program.push_back(curr);
        operands.clear();
        numEvents++;
    }

    //wrap info to pass as a tuple

    return std::make_tuple(program, curr_program, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum); //prologue instructions + info

    //unpack with square brackets and auto within qontrasim
}

// -- syndrome extraction round -- //

inline std::tuple<qes::Program<>, qes::Program<>, int, int, std::map<int, int>, std::map<int, int>, int>
FullSystemSimulator::stability_syndrome_extraction_round(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs, int numEvents, int numMeasurements, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, int roundNum, qes::Program<> program){
    std::map<int, std::vector<int>> p_to_d = create_p_to_d_map_for_stability_exp(d); //use for parity to data for checks
    std::cout<<roundNum<<std::endl;
    //program is a vector of instructions: simple end goal
    
    std::vector<qes::Instruction<>> curr_program; //append to this
    qes::Instruction<> curr; //build this instruction every time: don't need the default types
    std::vector<int> operands; //use this for every instruction, clear after using

    //reset all qubits

    int totalQubits = 2*pow(d,2)+1; //2d^2 qubits, just add 1 cuz doing < logic

    std::set<int> zOuter = {16,17,18,22,26,30,31,32};
    std::set<int> x = {19,21,24,27,29}; //only inner
    std::set<int> zInner = {20,23,25,28};
    std::set<int> zStabilizers = zOuter;
    zStabilizers.insert(zInner.begin(), zInner.end()); //all zstabs
    std::set<int> allStabilizers = zStabilizers;
    allStabilizers.insert(x.begin(), x.end()); //all stabs

    //stabilizer sets populated

    //using sets for h gates as they're weird - 5 sets

    std::set<int> h1 = {1,3,4,6,7,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    std::set<int> h2 = {0,1,2,3,5,6,8,9,10,11};
    std::set<int> h3 = {1,3,4,11,12,14};
    std::set<int> h4 = {4,5,6,7,9,10,12,13,14,15};
    std::set<int> h5 = {5,7,10,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};

    //syndrome extract and epilogue have different h1, epilogue has different h5

    std::set<int> h1dif = {0,2,5,8,10,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    std::set<int> h5dif = {0,1,2,3,4,6,8,9,11,12,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};


    //using these from params
    // int numEvents = 0;
    // int numMeasurements = 0;

    // std::map<int, int> latestMeasurement; //map for <stabilizer:latest measurement number>
    // std::map<int, int> previousMeasurement; //map for <stabilizer:previous measurement number>

    for (auto i : allStabilizers) {
        operands.push_back(i);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h1 (DIF ONE)

    for (auto q : h1dif) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //CNOT EVENTS - CHECK IF IN THE SELECTIVE SET CUZ SYNDROME EXTRACT

    std::map<int, int> cnot_check_indices = std::map<int, int>(); //init all to 0

    for (int i = 16; i < 33; i++) {
        cnot_check_indices[i] = 0;
    }

    //cuz map is correct order, so j go thru the indices
    
    //cnot_line1
    for (int i = 19; i < 33; i++) {
        if (i != 26) {
            //skip 26
            //do the correct index of all from map
            //SELECTIVE LOGIC
            if (stabilizers_each_round[roundNum].contains(i)) {
                int index = cnot_check_indices[i];
                //if z then value first then i
                if (zStabilizers.contains(i)) {
                    operands.push_back(p_to_d[i][index]);
                    operands.push_back(i);
                } else {
                    operands.push_back(i);
                    operands.push_back(p_to_d[i][index]);
                }
                cnot_check_indices[i]++; //increment for next time
            }
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h2

    for (auto q : h2) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //cnot_line2
    for (int i = 16; i < 31; i++) {
        if (i != 18 && i != 26) {
            //skip 18, 26
            //do the correct index of all from map
            //SELECTIVE LOGIC
            if (stabilizers_each_round[roundNum].contains(i)) {
                int index = cnot_check_indices[i];
                //if z then value first then i
                if (zStabilizers.contains(i)) {
                    operands.push_back(p_to_d[i][index]);
                    operands.push_back(i);
                } else {
                    operands.push_back(i);
                    operands.push_back(p_to_d[i][index]);
                }
                cnot_check_indices[i]++; //increment for next time
            }
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h3

    for (auto q : h3) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //cnot_line3
    for (int i = 18; i < 33; i++) {
        if (i != 22 && i != 30) {
            //skip 22, 30
            //do the correct index of all from map
            //SELECTIVE LOGIC
            if (stabilizers_each_round[roundNum].contains(i)) {
                int index = cnot_check_indices[i];
                //if z then value first then i
                if (zStabilizers.contains(i)) {
                    operands.push_back(p_to_d[i][index]);
                    operands.push_back(i);
                } else {
                    operands.push_back(i);
                    operands.push_back(p_to_d[i][index]);
                }
                cnot_check_indices[i]++; //increment for next time
            }
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h4

    for (auto q : h4) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //cnot_line4
    for (int i = 16; i < 30; i++) {
        if (i != 22) {
            //skip 22
            //do the correct index of all from map
            //SELECTIVE LOGIC
            if (stabilizers_each_round[roundNum].contains(i)) {
                int index = cnot_check_indices[i];
                //if z then value first then i
                if (zStabilizers.contains(i)) {
                    operands.push_back(p_to_d[i][index]);
                    operands.push_back(i);
                } else {
                    operands.push_back(i);
                    operands.push_back(p_to_d[i][index]);
                }
                cnot_check_indices[i]++; //increment for next time
            }
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h5

    for (auto q : h5) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //measure all stabs ONLY IF IN SELECTIVE ROUND

    previousMeasurement = latestMeasurement; //set it before changing latest

    for (auto i : allStabilizers) {
        if (stabilizers_each_round[roundNum].contains(i)) {
            operands.push_back(i);
            latestMeasurement[i] = numMeasurements;
            numMeasurements++;
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);

    program.push_back(curr);
    curr_program.push_back(curr);
    operands.clear();

    //events for all stabs - logic stays same

    for (auto i : allStabilizers) {
        operands.push_back(numEvents);
        operands.push_back(previousMeasurement[i]); //got set cuz of prologue and then curr round
        operands.push_back(latestMeasurement[i]);

        curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
        program.push_back(curr);
        curr_program.push_back(curr);
        operands.clear();
        numEvents++;
    }

    return std::make_tuple(program, curr_program, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum++);

}

// -- epilogue -- //

inline std::tuple<qes::Program<>, qes::Program<>, int, int, std::map<int, int>, std::map<int, int>, int>
FullSystemSimulator::stability_syndrome_epilogue(int d, std::map<int, std::set<int>> stabilizers_each_round, std::map<int, double> stab_freqs, int numEvents, int numMeasurements, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, int roundNum, qes::Program<> program){
    std::map<int, std::vector<int>> p_to_d = create_p_to_d_map_for_stability_exp(d); //use for parity to data for checks

    //program is a vector of instructions: simple end goal
    
    std::vector<qes::Instruction<>> curr_program; //append to this
    qes::Instruction<> curr; //build this instruction every time: don't need the default types
    std::vector<int> operands; //use this for every instruction, clear after using

    //reset all qubits

    int totalQubits = 2*pow(d,2)+1; //2d^2 qubits, just add 1 cuz doing < logic

    std::set<int> zOuter = {16,17,18,22,26,30,31,32};
    std::set<int> x = {19,21,24,27,29}; //only inner
    std::set<int> zInner = {20,23,25,28};
    std::set<int> zStabilizers = zOuter;
    zStabilizers.insert(zInner.begin(), zInner.end()); //all zstabs
    std::set<int> allStabilizers = zStabilizers;
    allStabilizers.insert(x.begin(), x.end()); //all stabs

    //stabilizer sets populated

    //using sets for h gates as they're weird - 5 sets

    std::set<int> h1 = {1,3,4,6,7,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    std::set<int> h2 = {0,1,2,3,5,6,8,9,10,11};
    std::set<int> h3 = {1,3,4,11,12,14};
    std::set<int> h4 = {4,5,6,7,9,10,12,13,14,15};
    std::set<int> h5 = {5,7,10,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};

    //syndrome extract and epilogue have different h1, epilogue has different h5

    std::set<int> h1dif = {0,2,5,8,10,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    std::set<int> h5dif = {0,1,2,3,4,6,8,9,11,12,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};


    //using these from params
    // int numEvents = 0;
    // int numMeasurements = 0;

    // std::map<int, int> latestMeasurement; //map for <stabilizer:latest measurement number>
    // std::map<int, int> previousMeasurement; //map for <stabilizer:previous measurement number>

    for (auto i : allStabilizers) {
        operands.push_back(i);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h1 (dif one)

    for (auto q : h1dif) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //CNOT EVENTS - epilogue so do all

    std::map<int, int> cnot_check_indices = std::map<int, int>(); //init all to 0

    for (int i = 16; i < 33; i++) {
        cnot_check_indices[i] = 0;
    }

    //cuz map is correct order, so go thru the indices
    
    //cnot_line1
    for (int i = 19; i < 33; i++) {
        if (i != 26) {
            //skip 26
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h2

    for (auto q : h2) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //cnot_line2
    for (int i = 16; i < 31; i++) {
        if (i != 18 && i != 26) {
            //skip 18, 26
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h3

    for (auto q : h3) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //cnot_line3
    for (int i = 18; i < 33; i++) {
        if (i != 22 && i != 30) {
            //skip 22, 30
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h4

    for (auto q : h4) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //cnot_line4
    for (int i = 16; i < 30; i++) {
        if (i != 22) {
            //skip 22
            //do the correct index of all from map
            int index = cnot_check_indices[i];
            //if z then value first then i
            if (zStabilizers.contains(i)) {
                operands.push_back(p_to_d[i][index]);
                operands.push_back(i);
            } else {
                operands.push_back(i);
                operands.push_back(p_to_d[i][index]);
            }
            cnot_check_indices[i]++; //increment for next time
        }
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //h5 (DIF FOR EPILOGUE)

    for (auto q : h5dif) {
        operands.push_back(q);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);

    program.push_back(curr);
    curr_program.push_back(curr);

    operands.clear();

    //measure ALL QUBITS - this the dif for epilogue

    previousMeasurement = latestMeasurement; //set it before changing latest

    for (int i = 0; i < totalQubits; i++) {
        operands.push_back(i);
        latestMeasurement[i] = numMeasurements;
        numMeasurements++;
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);

    program.push_back(curr);
    curr_program.push_back(curr);
    operands.clear();

    //events for all stabs, logic shud be same

    for (auto i : allStabilizers) {
        operands.push_back(numEvents);
        operands.push_back(previousMeasurement[i]); //got set cuz of prologue and then curr round
        operands.push_back(latestMeasurement[i]);

        curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
        program.push_back(curr);
        curr_program.push_back(curr);
        operands.clear();
        numEvents++;
    }
    

    //events for x stabilizers and their data qubits (only inners)

    for (auto xQ : x) {
        operands.push_back(numEvents);
        operands.push_back(latestMeasurement[xQ]);
        for (auto dataWith : p_to_d[xQ]) {
            operands.push_back(latestMeasurement[dataWith]); //those got stored in map too!
        }
        curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
        program.push_back(curr);
        curr_program.push_back(curr);
        operands.clear();
        numEvents++;
    }

    //obs is z stabilizer measures now

    operands.push_back(0);

    for (auto i : zStabilizers) {
        operands.push_back(latestMeasurement[i]);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("obs", operands);
    program.push_back(curr);
    curr_program.push_back(curr);
    operands.clear();

    return std::make_tuple(program, curr_program, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum);

    // std::ofstream outFile2("d4_stability_experiment.txt"); //saving to file to analyze

    //     // for (const auto& pair : stab_freqs) {
    //     //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
    //     // }

    // for (auto instr : program) {
    //     outFile2 << instr << std::endl;
    // }

    // // for (auto sel : selectors) {
    // //     outFile <<sel<<std::endl;
    // // }

    // outFile2.close();

    //done
}



// -- beginning code for final project, simple/iterative masking schedule on a surface code


inline qes::Program<>
FullSystemSimulator::create_program_masked_selective(int d, double p, int policy, int mem, int numRounds) {

    //same circuit gen, except choose a random mask based on double p
    //policy chosen by int policy

    //-- num rounds is passed in now, but set epilogue flag based on that

    bool skipEpilogue = (numRounds != 10*d); //if 10*d then reg exp, or else we doing epilogue skip to get the obs thing

    std::vector<int> selectors;
    //program is a vector of instructions: simple end goal
    
    std::vector<qes::Instruction<>> program; //append to this
    qes::Instruction<> curr; //build this instruction every time: don't need the default types
    std::vector<int> operands; //use this for every instruction, clear after using

    //reset all qubits

    int total_qubits = pow(d, 2) + pow(d, 2) - 1;

    for (int i = 0; i < total_qubits; i++) {
        operands.push_back(i);
    }

    std::set<int> zStabilizers = std::set<int>();

    std::set<int> xStabilizers = std::set<int>();
    std::set<int> innerStabilizers = std::set<int>();
    std::vector<int> firstRound; //save first round

    std::set<int> allStab;

    std::map<int, std::vector<int>> p_to_d = create_p_to_d_map(d); //use for parity to data for checks

    int numEvents = 0;
    int newNumEvents = 0; //var to reset the numevents for ext1
    // std::vector<int> lastRoundChecks; //don't use this anymore

    std::map<int, int> latestMeasurement; //map for <stabilizer:latest measurement number>
    std::map<int, int> previousMeasurement; //map for <stabilizer:previous measurement number>
    


    curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);
    program.push_back(curr);
    //std::cout<<curr<<std::endl;
    operands.clear(); //reset all

    //go through map: the key is the parity, the value is the data qubits

    //10 inner passes fixed: first 1 is ext0, next 9 are ext1 - only difference is the events

    for (int j = pow(d, 2); j < pow(d, 2) + (d - 1); j++) {
        zStabilizers.insert(j);
        allStab.insert(j);
    }

    for (int j = pow(d, 2) + (d-1); j < pow(d, 2) + 2*(d-1); j++) {
        xStabilizers.insert(j);
        allStab.insert(j);
    }

    int rowC = 0; //keep track of this: divide by d and cast to int to get the row num

    for (int j = pow(d, 2) + 2*(d-1); j < total_qubits; j++) {
        int rowNum = (int)(rowC/(d-1)); //changed logic

        innerStabilizers.insert(j);

        //pattern: if (rowNum + j)%2 == 0, then it is a z stabilizer, else it is x

        if ((rowNum + j)%2 == 0) {
            //z stabilizer
            zStabilizers.insert(j);
        } else {
            xStabilizers.insert(j);
        }

        allStab.insert(j);

        rowC++;
    }

    std::set<int> outerStabs; //need for new selective logic

    auto it = std::set_difference(allStab.begin(), allStab.end(), 
                    innerStabilizers.begin(), innerStabilizers.end(), 
                    std::inserter(outerStabs, outerStabs.begin())); //adding with the inserter - cuz the dest needs to be alloc otherwise

    // for (auto i : outerStabs) {
    //     std::cout<< i < std::endl;
    // }

    // std::cout << "OuterStabs: ";
    // for (const int &val : outerStabs) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    int numMeasurements = 0;

    int selectiveOrControl = 1; //0 means control, 1 means selective


    std::random_device rand; //rand
    std::mt19937 gen(rand()); //seed
    // std::uniform_int_distribution<> dist(pow(d,2), 2*(pow(d,2)-1)); //d^2 to 2(d^2-1) inclusive

    std::vector<int> weights;
    for (int i = pow(d,2); i <= 2*(pow(d,2)-1); i++) {
        if (!outerStabs.contains(i)) {
            weights.push_back(6); // 6x as much weight to inner
        } else {
            weights.push_back(1); // Regular weight
        }
    }

    std::discrete_distribution<> dist(weights.begin(), weights.end());
    std::uniform_real_distribution<> dist2(0.25, 0.75);

        //1-p cuz we wanna mask p, so keep 1-p

    std::set<int> selective_stabs;

    //lets fix these "random" stabs for experimental data

    // selective_stabs.insert(9);
    // selective_stabs.insert(10);c
    // selective_stabs.insert(11);
    // selective_stabs.insert(12);
    // selective_stabs.insert(13);
    // selective_stabs.insert(14);
    // selective_stabs.insert(15);
    // selective_stabs.insert(16);


    std::map<int, std::set<int>> stabs_per_round; //map of round num to stabilizers to measure

    // if (policy == 0) {
    //     //simple
    //     for (int i = 0; i < numRounds; i++) {
    //         //new random per round
    //         selective_stabs.clear();
    //         while (selective_stabs.size() < num_stabs) {
    //             //needa do less than not != cuz it might go over by like 1 but chill
    //             //unique, so will get w out repeats

    //             int random_stab = dist(gen) + pow(d,2);
    //             //do the outer thing tho (if even outer add the -1, if odd outer add the +1)
    //             if (outerStabs.contains(random_stab)) {
    //                 if (random_stab%2==0) {
    //                     selective_stabs.insert(random_stab-1);
    //                 } else {
    //                     selective_stabs.insert(random_stab+1);
    //                 }
    //             }

    //             selective_stabs.insert(random_stab);

    //         }

    //         // selective_stabs = {25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48};

    //         // selective_stabs = {9,10,11,12, 13,14,15,16};

    //         // if (i%2 == 0) {
    //         //     selective_stabs = {9,10,11,12};
    //         // } else {
    //         //     selective_stabs = {13,14,15,16};
    //         // }

    //         // selective_stabs = {49,50,51,52,53,54,62,64,66,67,69,71,74,76,78,79,81,83,86,88,90,91,93,95,55,56,57,58,59,60,61,63,65,68,70,72,73,75,77,80,82,84,85,87,89,92,94,96, 55,56,57,58,59,60,61,63,65,68,70,72,73,75,77,80,82,84,85,87,89,92,94,96,49,50,51,52,53,54,62,64,66,67,69,71,74,76,78,79,81,83,86,88,90,91,93,95};

    //         stabs_per_round[i] = selective_stabs;
    //     }
        
    // } else {
    //     //iterative scheduling
    //     for (int i = 0; i < numRounds; i++) {
    //         if (i%10 == 0) {
    //             //then we'll unmask some %
    //             //see what round we at first
    //             int timestep = i/10;
    //             double portion_mask_to_keep = pow(10, -timestep); //keep this percent, like 10%, 1%, etc. (same as removing 1 - x)
    //             //bruh pretty much ends up being entire mask everytime
    //             //currently masking total - selective_stabs.size()
    //             int currently_masking = (allStab.size() - selective_stabs.size());
    //             int num_to_keep = int(portion_mask_to_keep*currently_masking);

    //             std::set<int> current_mask;

    //             std::set_difference(allStab.begin(), allStab.end(), selective_stabs.begin(), selective_stabs.end(), 
    //                                 std::inserter(current_mask, current_mask.begin())); //put dif into current mask

    //             //now only keep first num from current_mask

    //             std::set<int> new_mask;
    //             auto it = current_mask.begin(); //iterator thru mask
    //             for (int i = 0; i < num_to_keep; ++i) {
    //                 new_mask.insert(*it); //element at iterator
    //                 ++it; //move iterator
    //             }

    //             //get rid of new_mask from allstabs

    //             std::set<int> new_selective_stabs;

    //             std::set_difference(allStab.begin(), allStab.end(), new_mask.begin(), new_mask.end(), 
    //                                 std::inserter(new_selective_stabs, new_selective_stabs.begin())); //put dif into current mask


    //             //now set this round to new stabs
    //             stabs_per_round[i] = new_selective_stabs;

    //             // for (auto i : new_selective_stabs) {
    //             //     std::cout<<i<<std::endl;
    //             // }


    //         } else {
    //             //or else reg random
    //             selective_stabs.clear();
    //             while (selective_stabs.size() != num_stabs) {
    //             //unique, so will get w out repeats
    //                 int random_stab = dist(gen);
    //                 selective_stabs.insert(random_stab);
    //             }
    //             stabs_per_round[i] = selective_stabs;

    //             // if (i%2 == 0) {
    //             //     selective_stabs = {9,10,11,12};
    //             // } else {
    //             //     selective_stabs = {13,14,15,16};
    //             // }
                
    //             // stabs_per_round[i] = selective_stabs;
    //         }
    //     }
    // }

    // while (selective_stabs.size() > num_stabs) {
    //     //mighta gone over lets keep it to num_stabs so pop inner from end
    //     selective_stabs.erase(--selective_stabs.end());
    // }

    // for (const auto& pair : stabs_per_round) {
    //     std::cout << "Round " << pair.first << ": ";
    //     for (int value : pair.second) {
    //         std::cout << value << " ";
    //     }
    //     std::cout << std::endl;
    // }

     

    //now take stabs from stabs_per_round[i]


    //use 1 for even rounds, 2 for odd rounds


    //CHANGED FROM NUMROUNDS TO 10*D

    for (int round = 0; round < 10*d; round++) {
        //do the plus 1 cuz extra round at end with all
        //set a random mask before rounds - same per round (based on algo 2)
        //based on p value

        p = dist2(gen);
        int num_stabs = int((pow(d,2)-1) * (1-p)); //make it int to get num of total stabs based on p

        selective_stabs.clear();
        while (selective_stabs.size() < num_stabs) {
            //needa do less than not != cuz it might go over by like 1 but chill
            //unique, so will get w out repeats

            int random_stab = dist(gen) + pow(d,2);
            //do the outer thing tho (if even outer add the -1, if odd outer add the +1)
            if (outerStabs.contains(random_stab)) {
                if (random_stab%2==0) {
                    selective_stabs.insert(random_stab-1);
                } else {
                    selective_stabs.insert(random_stab+1);
                }
            }

            selective_stabs.insert(random_stab);

        }

        stabs_per_round[round] = selective_stabs;

        // if (round%2 == 0) {
        //     selective_stabs = selective_stabs1;
        // } else {
        //     selective_stabs = selective_stabs2;
        // }

        selective_stabs = stabs_per_round[round]; //set to this from map

        // for (auto i : selective_stabs) {
        //     std::cout<<i<<std::endl;
        // }

        if (1==1) {
            //simple policy - only measure selective_stabs for each round
            //so we applying the p mask
            for (int xStab : xStabilizers) {
                operands.push_back(xStab);
            }

            curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);
            curr.put("timing_error"); //with the timing error
            program.push_back(curr);
            //std::cout<<curr<<std::endl;
            operands.clear(); //reset all

            //CNOT EVENTS

            //adding selective cnot logic

            int i = 0;

            //CHANGED J2 LOGIC
            int j1 = pow(d, 2);
            int j2 = (pow(d, 2) + (d/2) * 4) - 1;
            
            // for (auto i : selective_stabs) {
            //     std::cout<<i<<std::endl;
            // }


            std::set<int> cnotsDone; //bruh dum solution but



            while (i < 4) {
                // cnotsDone.clear(); //per line
                if (i==2) {
                    j1++;
                    j2--;
                }

                std::vector<int> values;

                if (round == numRounds-1 || round == 0) {
                    //final round do ALL
                    for (int x = 0; x < (d/2); x++) {
                        operands.push_back(j2);
                        values = p_to_d.at(j2); //the values
                        operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                        j2-=2; //do d/2 each time
                    }
                    j2+=2*(d/2);
                } else {
                    loop2:
                    if (selective_stabs.contains(j2) && xStabilizers.contains(j2)) {
                        // std::cout<<j2<<std::endl;
                        for (int x = 0; x < (d/2); x++) {
                            // std::cout<<j2<<std::endl;
                            if (outerStabs.contains(j2) && selective_stabs.contains(j2) && p_to_d.contains(j2) && xStabilizers.contains(j2)) {
                                cnotsDone.insert(j2);
                                operands.push_back(j2);
                                values = p_to_d.at(j2); //the values
                                operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                            }
                            j2-=2; //do d/2 each time
                        }

                        j2+=2*(d/2);
                    } else {
                        j2-=2; //try it again
                        if (j2 >= *outerStabs.begin() && j2 <= *--outerStabs.end()) {
                            //sorted sets - deref the value at begin() for smallest, go one before the end
                            goto loop2;
                        }
                    }
                }

                //z

                if (round == numRounds - 1 || round == 0) {
                    //final round do all
                    for (int x = 0; x < (d/2); x++) {
                        values = p_to_d.at(j1); //the values
                        operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                        operands.push_back(j1);
                        j1+=2;
                    }
                    // if (i!=1) {
                    j1-=(2*(d/2));
                    // }
                } else {
                    //SELECTIVE NEEDA CHANGE LOGIC TF KINDA LOGIC THIS??

                    //take union of selective_stabs and Z and outer

                    //nvm

                    // std::set<int> selected_z;

                    // auto it1 = std::set_union(selective_stabs.begin(), selective_stabs.end(), zStabilizers.begin(), zStabilizers.end(), selected_z.begin());

                    // std::set<int> selected_z_outer;

                    // auto it2 = std::set_union(selected_z.begin(), selected_z.end(), outerStabs.begin(), outerStabs.end(), selected_z_outer.begin());

                    // for (auto stab : selected_z_outer) {
                    //     values = p_to_d.at(stab); //the values
                    //     operands.push_back(values[i%2]); //by
                    //     operands.push_back(stab);

                    // }
                    loop:
                    if (selective_stabs.contains(j1) && zStabilizers.contains(j1)) {
                        for (int x = 0; x < (d/2); x++) {
                            // std::cout<<j1<<std::endl;
                            if (outerStabs.contains(j1) && selective_stabs.contains(j1) && p_to_d.contains(j1) && zStabilizers.contains(j1)) {
                                cnotsDone.insert(j1);
                                values = p_to_d.at(j1); //the values
                                operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                                operands.push_back(j1);
                            }
                            j1+=2;
                        }
                        j1-=(2*(d/2));
                    } else {
                        j1+=2; //try it again
                        if (j1 >= *outerStabs.begin() && j1 <= *--outerStabs.end()) {
                            goto loop;
                        }
                    }
                }
                
                //inner below

                for (int j = pow(d, 2) + 2*(d-1); j < total_qubits; j++) {
                    if (round == numRounds - 1 || round == 0) {
                        if (zStabilizers.contains(j)) {
                            operands.push_back(p_to_d.at(j)[i%4]);
                            operands.push_back(j);
                        } else {
                            operands.push_back(j);
                            operands.push_back(p_to_d.at(j)[i%4]);
                        }
                    } else {
                        if (selective_stabs.contains(j)) {
                            if (zStabilizers.contains(j)) {
                                operands.push_back(p_to_d.at(j)[i%4]);
                                operands.push_back(j);
                            } else {
                                operands.push_back(j);
                                operands.push_back(p_to_d.at(j)[i%4]);
                            }
                        } 
                    }
                    
                }

                curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);
                program.push_back(curr); //instruction made and added to program
                //std::cout<<curr<<std::endl;

                i++;

                operands.clear(); //clear for next loop
            }

            //the cnots are done

            //h the x stabs again

            for (int xStab : xStabilizers) {
                operands.push_back(xStab);
            }

            curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);
            program.push_back(curr);
            //std::cout<<curr<<std::endl;
            operands.clear(); //reset all

            //measure stabilizers and reset

            if (1==1) {
            if (round == numRounds-1 || round == 0) {
                //measure all stabilizers in FINAL round

                previousMeasurement = latestMeasurement; //set it before changing latest

                for (int i = pow(d, 2); i < total_qubits; i++) {
                    operands.push_back(i);
                    latestMeasurement[i] = numMeasurements; //set the latest measurement to this num of measurement
                    numMeasurements++;
                }

                //add to the map
            } else {
                //or else dif measurement patterns
                previousMeasurement = latestMeasurement; //set it before changing latest
                // for (int i = pow(d, 2); i < total_qubits; i++) {
                //     operands.push_back(i);
                //     latestMeasurement[i] = numMeasurements;
                //     numMeasurements++;
                // }

                for (int i : selective_stabs) {
                    operands.push_back(i);
                    latestMeasurement[i] = numMeasurements;
                    numMeasurements++;
                }
            }
            
            // for (auto pair : latestMeasurement) {
            //     std::cout<< pair.first << "-->" << pair.second << std::endl;
            // }
            

            curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);
            program.push_back(curr); //added measurements
            //std::cout<<curr<<std::endl;
            operands.clear(); //for resets

            for (int i = pow(d, 2); i < total_qubits; i++) {
                operands.push_back(i);
            }

            curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);
            program.push_back(curr); //added resets: all stabilizers
            //std::cout<<curr<<std::endl;
            operands.clear();

            //events now

            if (round == 0) {
                //initial ext0 events

                std::set<int> allStabs = std::set<int>();
                allStabs.insert(xStabilizers.begin(), xStabilizers.end());
                allStabs.insert(zStabilizers.begin(), zStabilizers.end());

                //all events only for z stabilizer measurements
                for (int i : zStabilizers) {
                    operands.push_back(numEvents);
                    operands.push_back(latestMeasurement[i]); //push the latest z stabilizer measurement
                    // std::cout << i << "--> " << latestMeasurement[i] << std::endl;
                    curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
                    program.push_back(curr);
                    //std::cout<<curr<<std::endl;
                    operands.clear();
                    numEvents++;

                    //save these for epilogue events
                    firstRound.push_back(latestMeasurement[i]);
                }

            } else {
                //rest of rounds
                std::set<int> allStabs = std::set<int>();
                allStabs.insert(xStabilizers.begin(), xStabilizers.end());
                allStabs.insert(zStabilizers.begin(), zStabilizers.end());

                // for (auto c : allStabs) {
                //     std::cout<<c<<std::endl;
                // }

                // std::cout<<"pass"<<std::endl;
                for (int i : allStabs) {
                    operands.push_back(numEvents);
                    operands.push_back(previousMeasurement[i]); //u want the previous measurement too
                    operands.push_back(latestMeasurement[i]);
                    
                    curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
                    program.push_back(curr);

                    operands.clear();
                    numEvents++;
                }
            }
            }
        }

    }

    //everything else done: epilogue now

    // -- ONLY DOING IF EPILOGUE FLAG IS SET

    //measure data qubits

    std::vector<int> observable;

    for (int i = 0; i < pow(d, 2); i++) {
        operands.push_back(i);

        if (i < d) {
            observable.push_back(numMeasurements); //for final obs
        }
        //obs is j any row cuz x mem
        numMeasurements++;
    }

    // for (int i = pow(d,2)-d; i < pow(d,2);i++) {
    //     observable.push_back(numMeasurements); //for final obs
    //     numMeasurements++;
    // }

    curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);
    program.push_back(curr);
    //std::cout<<curr<<std::endl;
    operands.clear();

    //final events


    //-- ONLY DO THESE IF SKIP EPILOGUE IS OFF
    
    std::cout << numRounds << std::endl;
    std::cout << skipEpilogue << std::endl;

    if (1==1) {
        int index = 0;

        //x memory exp so we care about z stabilizers

        //ok so subtract the previous d^2 data q measurements then just add from map each time

        numMeasurements -= pow(d,2); 

        for (int z : zStabilizers) {
            
            if (z < pow(d, 2) + 2*(d-1)) {
                //that means it's a boundary z, so only 2
                // operands.push_back(index); //remove the eoffset
                // operands.push_back(firstRound[index]); //remove the mshift
                operands.push_back(numEvents);
                operands.push_back(latestMeasurement[z]); //moved this from below push back num events

                
                std::vector<int> values = p_to_d[z];
                
                // operands.push_back(values[1] + numMeasurements - 1);
                // operands.push_back(values[0] + numMeasurements - 1);

                operands.push_back(numMeasurements + values[1]);
                operands.push_back(numMeasurements + values[0]);
                // operands.push_back(numMeasurements - 1 - values[2]);
                // operands.push_back(numMeasurements - 1 - values[0]);


            } else {
                //4 values
                // operands.push_back(index);
                // operands.push_back(firstRound[index]); //removed shifts
                operands.push_back(numEvents);
                operands.push_back(latestMeasurement[z]);//moved
                

                std::vector<int> values = p_to_d.at(z); //the values

                operands.push_back(numMeasurements + values[0]);
                operands.push_back(numMeasurements + values[2]);
                operands.push_back(numMeasurements + values[1]);
                operands.push_back(numMeasurements + values[3]);

            }
            // std::cout << "LINE" << std::endl;
            curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
            program.push_back(curr);
            //std::cout<<curr<<std::endl;
            
            index++;
            numEvents++;

            operands.clear();
        }

    }

    //obs

    operands.push_back(0);

    operands.insert(operands.end(), observable.begin(), observable.end()); //right?

    // for (int i = 0; i < d; i++) {
    //     operands.push_back(numMeasurements + i);
    // }

    curr = qes::Instruction<qes::any_t, qes::any_t>("obs", operands);
    program.push_back(curr);

    operands.clear();


    if (p == 0) {
        //control
        std::ofstream outFile("final_project_surface_code_control.txt"); //saving to file to analyze

        for (auto instr : program) {
            outFile << instr << std::endl;
        }
        
        outFile.close();

    } else {
        //any non zero p

        std::ofstream outFile2("error_after_x_rounds.txt"); //saving to file to analyze

        for (auto instr : program) {
            outFile2 << instr << std::endl;
        }
        
        outFile2.close();
    }


    return program;
}





// spring 25, adaptive syndrome extraction

inline std::tuple<std::set<int>, std::set<int>, std::set<int>, std::set<int>, qes::Program<>, int>
FullSystemSimulator::create_masked_prologue(int d, qes::Program<> program, int instruction_number) {
    
    std::vector<int> selectors;

    //appending to the program passed in

    // std::vector<qes::Instruction<>> curr_program;
    qes::Instruction<> curr;
    std::vector<int> operands;

    std::map<int, std::vector<int>> p_to_d = create_p_to_d_map(d); //use for parity to data for checks

    int total_qubits = pow(d, 2) + pow(d, 2) - 1;

    for (int i = 0; i < total_qubits; i++) {
        operands.push_back(i);
    }

    std::set<int> zStabilizers = std::set<int>();
    std::set<int> xStabilizers = std::set<int>();
    std::set<int> innerStabilizers = std::set<int>();
    std::set<int> allStab = std::set<int>();

    //latest/prev measurements passed in

    curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);
    program.push_back(curr);
    instruction_number++;
    operands.clear(); //reset all

    for (int j = pow(d, 2); j < pow(d, 2) + (d - 1); j++) {
        zStabilizers.insert(j);
        allStab.insert(j);
    }

    for (int j = pow(d, 2) + (d-1); j < pow(d, 2) + 2*(d-1); j++) {
        xStabilizers.insert(j);
        allStab.insert(j);
    }

    int rowC = 0; //keep track of this: divide by d and cast to int to get the row num

    for (int j = pow(d, 2) + 2*(d-1); j < total_qubits; j++) {
        int rowNum = (int)(rowC/(d-1)); //changed logic

        innerStabilizers.insert(j);

        //pattern: if (rowNum + j)%2 == 0, then it is a z stabilizer, else it is x

        if ((rowNum + j)%2 == 0) {
            //z stabilizer
            zStabilizers.insert(j);
        } else {
            xStabilizers.insert(j);
        }

        allStab.insert(j);

        rowC++;
    }


    std::set<int> outerStabs; //need for new selective logic

    auto it = std::set_difference(allStab.begin(), allStab.end(), 
                    innerStabilizers.begin(), innerStabilizers.end(), 
                    std::inserter(outerStabs, outerStabs.begin())); //adding with the inserter - cuz the dest needs to be alloc otherwise

    std::ofstream outFile("dynamic_sse.txt");

    for (auto instr : program) {
        outFile << instr << std::endl;
    }

    outFile.close();

    return std::make_tuple(allStab, innerStabilizers, zStabilizers, xStabilizers, program, instruction_number); //return instr_number too
}

inline std::tuple<qes::Program<>, int, int, int, std::map<int, int>, std::map<int, int>, std::vector<int>, qes::Program<>, int>
FullSystemSimulator::create_masked_se_round(int d, qes::Program<> program, int numEvents, int numMeasurements, int round, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, std::set<int> selective_stabs, std::set<int> allStab, std::set<int> innerStabilizers, std::set<int> xStabilizers, std::set<int> zStabilizers, std::vector<int> firstRound, int instruction_number, int disconnectRound) {

    //using the stabsThisRound PARAM as the stabs to measure (handle that externally)

    qes::Program<> currProgram; //add to this, return this at end to cuz needa decode indiv rounds
    qes::Instruction<> curr;
    std::vector<int> operands;
    int total_qubits = pow(d, 2) + pow(d, 2) - 1;
    int numRounds = 10*d; //shud match the external caller

    std::map<int, std::vector<int>> p_to_d = create_p_to_d_map(d); //use for parity to data for checks

    // for (const auto& pair : previousMeasurement) {
    //     std::cout << "KeyPPP: " << pair.first << ", Value: " << pair.second << '\n';
    // }
    // for (const auto& pair : latestMeasurement) {
    //     std::cout << "KeyLLL: " << pair.first << ", Value: " << pair.second << '\n';
    // }

    std::set<int> outerStabs; //need for new selective logic

    auto it = std::set_difference(allStab.begin(), allStab.end(), 
            innerStabilizers.begin(), innerStabilizers.end(), 
            std::inserter(outerStabs, outerStabs.begin())); //adding with the inserter - cuz the dest needs to be alloc otherwise

    for (int xStab : xStabilizers) {
        operands.push_back(xStab); //should be zStabs for x mem exp
    }
    
    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);
    curr.put("timing_error"); //with the timing error
    program.push_back(curr);
    currProgram.push_back(curr);
    instruction_number++;
    operands.clear(); //reset all

    int i = 0;

    int j1 = pow(d, 2);
    int j2 = (pow(d, 2) + (d/2) * 4) - 1;
    
    std::set<int> cnotsDone; 

    // for (auto i : selective_stabs) {
    //     std::cout << i << std::endl;
    // }


    while (i < 4) {
        // cnotsDone.clear(); //per line
        if (i==2) {
            j1++;
            j2--;
        }

        std::vector<int> values;

        if (round == numRounds-1 || round == 0) {
            //final round and first round do ALL
            for (int x = 0; x < (d/2); x++) {
                operands.push_back(j2);
                values = p_to_d.at(j2); //the values
                operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                j2-=2; //do d/2 each time
            }
            j2+=2*(d/2);
        } else {
            loop2:
            if (selective_stabs.contains(j2) && xStabilizers.contains(j2)) {
                // std::cout<<j2<<std::endl;
                for (int x = 0; x < (d/2); x++) {
                    // std::cout<<j2<<std::endl;
                    if (outerStabs.contains(j2) && selective_stabs.contains(j2) && p_to_d.contains(j2) && xStabilizers.contains(j2)) {
                        cnotsDone.insert(j2);
                        operands.push_back(j2);
                        values = p_to_d.at(j2); //the values
                        operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                    }
                    j2-=2; //do d/2 each time
                }

                j2+=2*(d/2);
            } else {
                j2-=2; //try it again
                if (j2 >= *outerStabs.begin() && j2 <= *--outerStabs.end()) {
                    //sorted sets - deref the value at begin() for smallest, go one before the end
                    goto loop2;
                }
            }
        }

        //z

        if (round == numRounds - 1 || round == 0) {
            //final round do all
            for (int x = 0; x < (d/2); x++) {
                values = p_to_d.at(j1); //the values
                operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                operands.push_back(j1);
                j1+=2;
            }
            // if (i!=1) {
            j1-=(2*(d/2));
            // }
        } else {
            loop:
            if (selective_stabs.contains(j1) && zStabilizers.contains(j1)) {
                for (int x = 0; x < (d/2); x++) {
                    // std::cout<<j1<<std::endl;
                    if (outerStabs.contains(j1) && selective_stabs.contains(j1) && p_to_d.contains(j1) && zStabilizers.contains(j1)) {
                        cnotsDone.insert(j1);
                        values = p_to_d.at(j1); //the values
                        operands.push_back(values[i%2]); //i%2 so if 0th and 2nd row add first connection, 1st and 3rd row add second connection
                        operands.push_back(j1);
                    }
                    j1+=2;
                }
                j1-=(2*(d/2));
            } else {
                j1+=2; //try it again
                if (j1 >= *outerStabs.begin() && j1 <= *--outerStabs.end()) {
                    goto loop;
                }
            }
        }
        
        //inner below

        for (int j = pow(d, 2) + 2*(d-1); j < total_qubits; j++) {
            if (round == numRounds - 1 || round == 0) {
                if (zStabilizers.contains(j)) {
                    operands.push_back(p_to_d.at(j)[i%4]);
                    operands.push_back(j);
                } else {
                    operands.push_back(j);
                    operands.push_back(p_to_d.at(j)[i%4]);
                }
            } else {
                if (selective_stabs.contains(j)) {
                    if (zStabilizers.contains(j)) {
                        operands.push_back(p_to_d.at(j)[i%4]);
                        operands.push_back(j);
                    } else {
                        operands.push_back(j);
                        operands.push_back(p_to_d.at(j)[i%4]);
                    }
                }
            }
            
        }

        curr = qes::Instruction<qes::any_t, qes::any_t>("cx", operands);
        program.push_back(curr); //instruction made and added to program
        currProgram.push_back(curr);
        instruction_number++;
        //std::cout<<curr<<std::endl;

        i++;

        operands.clear(); //clear for next loop
    }

    //the cnots are done

    //h the x stabs again

    for (int xStab : xStabilizers) {
        operands.push_back(xStab);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("h", operands);
    program.push_back(curr);
    currProgram.push_back(curr);
    instruction_number++;
    //std::cout<<curr<<std::endl;
    operands.clear(); //reset all

    //measure stabilizers and reset

    if (round == numRounds-1 || round == 0) {
        //measure all stabilizers in FINAL round

        previousMeasurement = latestMeasurement; //set it before changing latest

        for (int i = pow(d, 2); i < total_qubits; i++) {
            operands.push_back(i);
            latestMeasurement[i] = numMeasurements; //set the latest measurement to this num of measurement
            numMeasurements++;
        }

        //add to the map
    } else {
        //or else dif measurement patterns
        previousMeasurement = latestMeasurement; //set it before changing latest

        for (int i : selective_stabs) {
            operands.push_back(i);
            latestMeasurement[i] = numMeasurements;
            numMeasurements++;
        }
    }
    

    curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);
    program.push_back(curr); //added measurements
    currProgram.push_back(curr);
    instruction_number++;
    //std::cout<<curr<<std::endl;
    operands.clear(); //for resets

    for (int i = pow(d, 2); i < total_qubits; i++) {
        operands.push_back(i);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("reset", operands);
    program.push_back(curr); //added resets: all stabilizers
    currProgram.push_back(curr);
    instruction_number++;
    //std::cout<<curr<<std::endl;
    operands.clear();

    //events now

    if (round == 0) {
        //initial ext0 events

        std::set<int> allStabs = std::set<int>();
        allStabs.insert(xStabilizers.begin(), xStabilizers.end());
        allStabs.insert(zStabilizers.begin(), zStabilizers.end());

        //all events only for z stabilizer measurements

        for (int i : zStabilizers) {
            operands.push_back(numEvents);
            operands.push_back(latestMeasurement[i]); //push the latest z stabilizer measurement
            // std::cout << i << "--> " << latestMeasurement[i] << std::endl;
            curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
            program.push_back(curr);
            currProgram.push_back(curr);
            instruction_number++;
            //std::cout<<curr<<std::endl;
            operands.clear();
            numEvents++;

            //save these for epilogue events
            firstRound.push_back(latestMeasurement[i]);
        }

    } else {
        //rest of rounds
        std::set<int> allStabs = std::set<int>();
        allStabs.insert(xStabilizers.begin(), xStabilizers.end());
        allStabs.insert(zStabilizers.begin(), zStabilizers.end());

        // int disconnect = 3;



        if (round%disconnectRound==0) {
            //just dont do anything
            //disconnect (no prev)
            // for (int i : zStabilizers) {
            //     operands.push_back(numEvents);
            //     operands.push_back(latestMeasurement[i]);
            //     // operands.push_back(latestMeasurement[i]); //with itself wont do anything right?

            //     curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
            //     program.push_back(curr);
            //     currProgram.push_back(curr);
            //     instruction_number++;

            //     operands.clear();
            //     numEvents++;
            // }

        } else {
            for (int i : allStabs) {
                operands.push_back(numEvents);
                operands.push_back(previousMeasurement[i]); //u want the previous measurement too
                operands.push_back(latestMeasurement[i]);
                
                curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
                program.push_back(curr);
                currProgram.push_back(curr);
                instruction_number++;

                operands.clear();
                numEvents++;
            }
        }
        
    }

    round++;
    

    std::ofstream outFile("dynamic_sse_rounds.txt");

    for (auto instr : program) {
        outFile <<instr << std::endl;
    }

    outFile.close();

    return std::make_tuple(program, numEvents, numMeasurements, round, previousMeasurement, latestMeasurement, firstRound, currProgram, instruction_number); //instr number too
    
}




inline qes::Program<>
FullSystemSimulator::create_masked_epilogue(int d, qes::Program<> program, int numEvents, int numMeasurements, int round, std::map<int, int> latestMeasurement, std::map<int, int> previousMeasurement, std::set<int> selective_stabs, std::set<int> allStab, std::set<int> innerStabilizers, std::set<int> xStabilizers, std::set<int> zStabilizers, std::vector<int> firstRound, int instruction_number) {

    std::vector<int> observable;

    std::map<int, std::vector<int>> p_to_d = create_p_to_d_map(d); //use for parity to data for checks

    qes::Instruction<> curr; //build this instruction every time: don't need the default types
    std::vector<int> operands; //use this for every instruction, clear after using

    for (int i = 0; i < pow(d, 2); i++) {
        operands.push_back(i);

        if (i < d) {
            observable.push_back(numMeasurements); //for final obs
        }
        // if (numMeasurements%d == (d-1)) {
            
        // }
        numMeasurements++;
    }

    // for (int i = pow(d,2)-d; i < pow(d,2);i++) {
    //     observable.push_back(numMeasurements); //for final obs
    //     numMeasurements++;
    // }

    curr = qes::Instruction<qes::any_t, qes::any_t>("measure", operands);
    program.push_back(curr);
    instruction_number++;
    //std::cout<<curr<<std::endl;
    operands.clear();

    //final events

    int index = 0;

    //x memory exp so we care about z stabilizers

    //ok so subtract the previous d^2 data q measurements then just add from map each time

    numMeasurements -= pow(d,2); 

    for (int z : zStabilizers) {
        
        if (z < pow(d, 2) + 2*(d-1)) {
            //that means it's a boundary z, so only 2
            // operands.push_back(index); //remove the eoffset
            // operands.push_back(firstRound[index]); //remove the mshift
            operands.push_back(numEvents);
            operands.push_back(latestMeasurement[z]); //moved this from below push back num events

            
            std::vector<int> values = p_to_d[z];
            
            // operands.push_back(values[1] + numMeasurements - 1);
            // operands.push_back(values[0] + numMeasurements - 1);

            operands.push_back(numMeasurements + values[1]);
            operands.push_back(numMeasurements + values[0]);
            // operands.push_back(numMeasurements - 1 - values[2]);
            // operands.push_back(numMeasurements - 1 - values[0]);


        } else {
            //4 values
            // operands.push_back(index);
            // operands.push_back(firstRound[index]); //removed shifts
            operands.push_back(numEvents);
            operands.push_back(latestMeasurement[z]);//moved
            

            std::vector<int> values = p_to_d.at(z); //the values

            operands.push_back(numMeasurements + values[0]);
            operands.push_back(numMeasurements + values[2]);
            operands.push_back(numMeasurements + values[1]);
            operands.push_back(numMeasurements + values[3]);

        }
        // std::cout << "LINE" << std::endl;
        curr = qes::Instruction<qes::any_t, qes::any_t>("event", operands);
        program.push_back(curr);
        instruction_number++;
        //std::cout<<curr<<std::endl;
        
        index++;
        numEvents++;

        operands.clear();
    }

    //events done

    //obs

    operands.push_back(0);

    for (int i = 0; i < d; i++) {
        operands.push_back(numMeasurements + i);
    }

    curr = qes::Instruction<qes::any_t, qes::any_t>("obs", operands);
    program.push_back(curr);
    instruction_number++;


    operands.clear();

    return program;
}

}   // qontra
