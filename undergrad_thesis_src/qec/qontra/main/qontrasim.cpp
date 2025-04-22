// /*
//  *  author: Suhas Vittal
//  *  date:   28 January 2024
//  * */

// #include <qontra/ext/qes.h>
// #include <qontra/sim/base/clifford_sim.h>
// #include <qontra/sim/base/frame_sim.h>
// #include <qontra/sim/full_system_sim.h>

// #include <vtils/cmd_parse.h>
// #include <vtils/filesystem.h>
// #include <vtils/ini_parse.h>
// #include <vtils/timer.h>

// #include <tuple>

// #include <filesystem>

// using namespace qontra;
// using namespace vtils;

// int main(int argc, char* argv[]) {
//     MPI_Init(NULL, NULL);
//     int world_rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//     std::string help = "usage: ./qontrasim <qes-file> "
//                         "--shots <shots> "
//                         "--trace-out <trace-folder> "
//                         "--stim-out <stim-file> "
//                         "--data-out <data-file>\n"
//                         "Optional:\n"
//                         "\t--config <ini-file>";

//     CmdParser pp(argc, argv, 1);
//     pp.help = help;
//     if (pp.option_set("h")) {
//         std::cerr << help << std::endl;
//         return 1;
//     }
//     configure_optimal_batch_size();
//     G_SHOTS_PER_BATCH = 1;

//     // Initialize simulator:
//     FullSystemSimulator sim;

//     // Get command line inputs.
//     std::string qes_file(argv[1]);

//     uint64_t shots;

//     std::string ini_file;
//     pp.get("shots", shots, true);
//     pp.get("trace-out", sim.config.syndrome_output_folder, true);
//     pp.get("stim-out", sim.config.stim_output_file, true);
//     pp.get("data-out", sim.config.data_output_file, true);

//     fp_t p = 0.0;
//     pp.get("p", p);

//     if (pp.get("config", ini_file)) {
//         std::string parent_dir = get_parent_directory(ini_file);

//         IniParser ini(ini_file);
//         const auto& ini_map = ini.get_ini_map();
//         // First load in the subroutines.
//         for (const auto& p : ini_map.at("__ANON__")) {
//             qes::Program<> subroutine = qes::from_file(parent_dir + "/" + p.second);
//             sim.load_subroutine(p.first, subroutine);
//         }
//         // Set simulation config if possible.
//         ini.get("Config", "record_events_until", sim.config.record_events_until);
//         ini.get("Config", "record_obs_until", sim.config.record_obs_until);
//         ini.get("SSE", "$reoz_track_last_n_events", sim.config.sse.rreoz_track_last_n_events);
//     }
//     // qes::Program<> main_program = qes::from_file(qes_file);

//     // std::map<int, std::set<int>> stabilizers_each_round; //populate stabilizers for each round in map

//     // -- d5 stabilizer layouts -- //

//     // int numRounds = 10*5; //10d rounds

//     // std::set<int> stabilizers;

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%2==0) {
//     //         stabilizers = {25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48};
//     //     } else {
//     //         stabilizers = {25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48};
//     //     }
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

//     // int numRounds = 10*5; //10d rounds

//     // std::set<int> stabilizers;

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%4==0) {
//     //         stabilizers = {25,26,27,28};
//     //     } else if (i%4==1) {
//     //         stabilizers = {29, 30, 31, 32};
//     //     } else if (i%4==2) {
//     //         stabilizers = {34,36,37,39,42,44,45,47};
//     //     } else if (i%4==3) {
//     //         stabilizers = {33, 35, 38, 40, 41, 43, 46, 48};
//     //     }
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%5==0) {
//     //         stabilizers = {25,26,27,28,33}; //need to have even amount of outer
//     //     } else if (i%5==1) {
//     //         stabilizers = {29,30,31,32,34}; //need to have even amount of outer
//     //     } else if (i%5==2) {
//     //         stabilizers = {35,36,37,38,39};
//     //     } else if (i%5==3) {
//     //         stabilizers = {40,41,42,43,44};
//     //     } else if (i%5==4) {
//     //         stabilizers = {45,46,47,48};
//     //     }
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

// //     for (int i = 0; i < numRounds; i++) {
// //     if (i % 20 == 1) {
// //         stabilizers = {25, 26, 27, 28}; // need to have even amount of outer
// //     } else if (i % 20 == 2) {
// //         stabilizers= {29, 30, 31, 32};
// //     } else if (i % 20 == 3) {
// //         stabilizers= {33};
// //     } else if (i % 20 == 4) {
// //         stabilizers= {34};
// //     } else if (i % 20 == 5) {
// //         stabilizers= {35};
// //     } else if (i % 20 == 6) {
// //         stabilizers= {36};
// //     } else if (i % 20 == 7) {
// //         stabilizers = {37}; // Incremented condition and value
// //     } else if (i % 20 == 8) {
// //         stabilizers = {38}; // Incremented condition and value
// //     } else if (i % 20 == 9) {
// //         stabilizers = {39}; // Incremented condition and value
// //     } else if (i % 20 == 10) {
// //         stabilizers = {40}; // Incremented condition and value
// //     } else if (i % 20 == 11) {
// //         stabilizers = {41}; // Incremented condition and value
// //     } else if (i % 20 == 12) {
// //         stabilizers = {42}; // Incremented condition and value
// //     } else if (i % 20 == 13) {
// //         stabilizers = {43}; // Incremented condition and value
// //     } else if (i % 20 == 14) {
// //         stabilizers = {44}; // Incremented condition and value
// //     } else if (i % 20 == 15) {
// //         stabilizers = {45}; // Incremented condition and value
// //     } else if (i % 20 == 16) {
// //         stabilizers = {46}; // Incremented condition and value
// //     } else if (i % 20 == 17) {
// //         stabilizers = {47}; // Incremented condition and value
// //     } else if (i % 20 == 18) {
// //         stabilizers = {48}; // Incremented condition and value
// //     }

// //     stabilizers_each_round[i] = stabilizers;
// // }


//     // -- d7 stabilizer layouts -- //

//     // int numRounds = 70; //10d rounds

//     // std::set<int> stabilizers;

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%2==0) {
//     //         stabilizers = {49,50,51,52,53,54,62,64,66,67,69,71,74,76,78,79,81,83,86,88,90,91,93,95,55,56,57,58,59,60,61,63,65,68,70,72,73,75,77,80,82,84,85,87,89,92,94,96};
//     //     } else {
//     //         stabilizers = {55,56,57,58,59,60,61,63,65,68,70,72,73,75,77,80,82,84,85,87,89,92,94,96,49,50,51,52,53,54,62,64,66,67,69,71,74,76,78,79,81,83,86,88,90,91,93,95};
//     //     }
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

//     // -- d4 stability exp layouts -- //

//     // int numRounds = 10; //10d rounds

//     // std::set<int> stabilizers;

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%10==0) {
//     //         stabilizers = {16,17,18,22,26,30,31,32};
//     //     }else if (i%10==1) {
//     //         stabilizers = {19};
//     //     } else if (i%10==2) {
//     //         stabilizers = {20};
//     //     } else if (i%10==3) {
//     //         stabilizers = {21};
//     //     } else if (i%10==4) {
//     //         stabilizers = {23};
//     //     } else if (i%10==5) {
//     //         stabilizers = {24};
//     //     } else if (i%10==6) {
//     //         stabilizers = {25};
//     //     } else if (i%10==7) {
//     //         stabilizers = {27};
//     //     } else if (i%10==8) {
//     //         stabilizers = {28};
//     //     } else if (i%10==9) {
//     //         stabilizers = {29};
//     //     }
//     //     stabilizers_each_round[i] = stabilizers;
//     // }
    

//     // for (int j = 16; j < 33; j++) {
//     //     stabilizers.insert(j);
//     // }

//     // for (int i = 0; i < numRounds; i++) {
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

//     // -- d3 stabilizer layouts -- //

//     // int numRounds = 30; //10d rounds

//     // std::set<int> stabilizers;

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%2==0) {
//     //         stabilizers = {9, 10, 11, 12, 13, 14, 15, 16};
//     //     } else {
//     //         stabilizers = {9, 10, 11, 12, 13, 14, 15, 16};
//     //     }
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%2==0) {
//     //         stabilizers = {9, 10, 15, 14, 11,12};
//     //     } else {
//     //         stabilizers = {13, 16};
//     //     }
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%3==0) {
//     //         stabilizers = {9,10,14,15};
//     //     } else if (i%3 == 1) {
//     //         stabilizers = {11,12};
//     //     } else {
//     //         stabilizers = {13,16};
//     //     }
        
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%6==0) {
//     //         stabilizers = {9,10};
//     //     } else if (i%6 == 1) {
//     //         stabilizers = {11,12};
//     //     } else if (i%6==2) {
//     //         stabilizers = {13};
//     //     } else if (i%6==3) {
//     //         stabilizers = {14};
//     //     } else if (i%6==4) {
//     //         stabilizers = {15};
//     //     }else if (i%6==5) {
//     //         stabilizers = {16};
//     //     }
        
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

//     // std::map<int, double> frequency_stabilizers; //empty map, pass into the control round

//     // // //clear dets from prior experiment for safety

//     // --program creations function calls below-- //

//     // qes::Program<> main_program = sim.create_program(2, stabilizers_each_round, frequency_stabilizers);

//     // qes::Program<> main_program = sim.create_program_stability_exp(4, stabilizers_each_round, frequency_stabilizers, 0);

//     //**RUNNING 2 TIMES NOW. ONCE FOR CODE ABOVE, ONCE FOR CODE BELOW AFTER CODE ABOVE HAS MADE CONTROL DETS FILES**

//     // std::map<int, double> frequency_stabilizers;

//     // std::string path = "/Users/aryan/Desktop/quantumResearch/qec/qontra/build/temp_trace_control/";

//     // frequency_stabilizers = sim.read_dets(path); //now populate from dets files
    
//     // // for (const auto& pair : frequency_stabilizers) {
//     // //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
//     // // }

//     // int numRounds = 10*3; //10d rounds

//     // std::map<int, std::set<int>> stabilizers_each_round; //populate stabilizers for each round in map
//     // std::set<int> stabilizers;

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%2==0) {
//     //         stabilizers = {9,10,14,15};
//     //     } else {
//     //         stabilizers = {11,12,13,16};
//     //     }
//     //     stabilizers_each_round[i] = stabilizers;
//     // }

//     // qes::Program<> main_program = sim.create_program(3, stabilizers_each_round, frequency_stabilizers); 

//     // continue with retroactive

//     // for (const auto& pair : frequency_stabilizers) {
//     //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
//     // }
    
//     // Setup error model.



//     // -- incremental stability code generation -- //

//     // calling run_program individually per shot on the incremental program

//     std::map<int, double> frequency_stabilizers; //empty - shud get rid of this later
//     std::map<int, std::set<int>> stabilizers_each_round; //populate stabilizers for each round in map

//     int numRounds = 30; //fix to 30 rounds

//     // std::map<int, std::set<int>> selective_map;
//     // std::map<int, std::set<int>> regular_map;

//     // //create the maps above loop

//     // std::set<int> stabilizers_prior;

//     // for (int j = 16; j < 33; j++) {
//     //     stabilizers_prior.insert(j);
//     // }

//     // for (int i = 0; i < numRounds; i++) {
//     //     regular_map[i] = stabilizers_prior;
//     // }

//     // stabilizers_prior.clear();

//     // //now selective based on policy

//     // for (int i = 0; i < numRounds; i++) {
//     //     if (i%2==0) {
//     //         stabilizers_prior = {16,17,18,22,26,30,31,32};
//     //     } else {
//     //         stabilizers_prior = {19,20,21,23,24,25,27,28,29};
//     //     }
//     //     selective_map[i] = stabilizers_prior; //each round is this stabilizers
//     // }

//     //set the selective policy here - the "surrounding" stabilizers

//     std::map<int, std::set<int>> surrounding_stabs;

//     std::map<int, std::vector<int>> p_to_d_to_change = sim.create_p_to_d_map_for_stability_exp(4);

//     std::map<int, std::set<int>> p_to_d;

//     for (int i = 16; i <= 32; i++) {
//         std::set<int> curr;
//         for (auto q : p_to_d_to_change[i]) {
//             curr.insert(q);
//         }
//         p_to_d[i] = curr;
//     }



    // for (int i = 16; i <= 32; i++) {
    //     std::set<int> curr_data = p_to_d[i]; //check where else has this

    //     std::set<int> surrounding;

    //     for (int j = 16; j <= 32; j++) {
    //         if (i != j) {
    //             //use set intersection - see if any data qubits in common
//                 std::set<int> intersect; //common elements in here

//                 std::set_intersection(
//                     curr_data.begin(), curr_data.end(), p_to_d[j].begin(), p_to_d[j].end(), std::inserter(intersect, intersect.begin())
//                 );
                
//                 //if any common data qubit, add j
//                 if (!intersect.empty()) {
//                     surrounding.insert(j);
//                 }
//             }
//         }

//         //now added for all data q of curr
//         surrounding_stabs[i] = surrounding;
//     }

//     std::ofstream outFile("surrounding_stabs.txt");

//     for (auto parity : surrounding_stabs) {
//         outFile << "KEY: " << parity.first << "--> " << std::endl;
//         for (auto data : parity.second) {
//             outFile << data << std::endl;
//         }
//     }

//     outFile.close();

//     for (int shot = 0; shot < 1; shot++) {
//         // std::cout<<"HI"<<std::endl;
//         // std::cout<<"HI"<<std::endl;
//         //init params
//         std::cout<<shot<<std::endl;
//         qes::Program<> main_program; //add to this incrementally: it's just a vector of instructions so
//         qes::Program<> curr_program; //doesn't matter for prologue
//         int numEvents = 0;
//         int numMeasurements = 0;
//         std::map<int,int> latestMeasurement;
//         std::map<int,int> previousMeasurement;
//         int roundNum = 0;


//         //first init with prologue
//         //unpack tuple with auto and brackets

//         //params: d, stabilizers_each_round, stab_freqs, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum, program

//         std::tie(main_program, curr_program, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum) = sim.stability_syndrome_prologue(
//             4, stabilizers_each_round, frequency_stabilizers, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum, main_program);

//         //run

//         //setup before running
//         tables::ErrorAndTiming et;
//         et *= p*1000;

//         const uint64_t n = get_number_of_qubits(curr_program); //num qubits same
//         tables::populate(n, sim.config.errors, sim.config.timing, et);
//         sim.run_program<FrameSimulator>(curr_program, 1); //1 shot

//         //access the measures now
//         //put into a map: qubit to measurement

//         std::map<int, int> q_to_m;
//         std::map<int, int> previous_q_to_m;

//         //all stabs so qubits 16 to 32

//         int index = 0;

//         for (int i = 16; i <= 32; i++) {
//             q_to_m[i] = qontra::FrameSimulator::measures[index++]; //the respective measurement in vector (-16 cuz starts at 0)
//             // std::cout<<i << " is " << std::cout << ""
//         }

//         // for (auto i : qontra::FrameSimulator::measures) {
//         //     std::cout<<i<<std::endl;
//         // }

//         // for (auto i : qontra::FrameSimulator::measures) {
//         //     std::cout<<i<<std::endl;
//         // }
        

//         //k so it runs now nice, syndrome table checking is insig here tho

//         // std::ofstream outFile("currStabilityProgram.txt"); //saving to file to analyze

//         //     // for (const auto& pair : stab_freqs) {
//         //     //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
//         //     // }

//         // for (auto instr : curr_program) {
//         //     outFile << instr << std::endl;
//         // }

//         // // for (auto sel : selectors) {
//         // //     outFile <<sel<<std::endl;
//         // // }
//         // outFile.close();

//         // write_stats(batchno); //batchno starts at 1, up 1 each time
//         // batchno++;

//         int roundsSinceError = 1000; //set to high num - saying that first round is error free
//         int error = 0; //init to no error

//         int eventsToCheck = 5;
        
//         std::set<int> flipped_stabs; //keep flipped_stabs from round in here


//         std::ofstream outFile5;
//         outFile5.open("round_by_round_measured_stabs.txt", std::ios::app);

//         for (int round = 0; round < numRounds; round++) {
//             //call syndrome extract individually
//             //measure after each call

//             //if error: stick control for next 3
//             //else continue with selective

//             std::map<int, std::set<int>> stabilizers_each_round; //populate stabilizers for each round in map
//             std::set<int> stabilizers;

//             qes::Program<> curr_program; //for each round, run this to find errors

//             //error check IF roundsSinceError >= d

//             //control round all stabilizers:
            
//             if (error) {
//                 for (int j = 16; j < 33; j++) {
//                     stabilizers.insert(j);
//                 }

//                 for (int i = 0; i < numRounds; i++) {
//                     stabilizers_each_round[i] = stabilizers;
//                 }

//                 std::tie(main_program, curr_program, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum) = sim.stability_syndrome_extraction_round(
//             4, stabilizers_each_round, frequency_stabilizers, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum, main_program);

//             } else {
//                 //can do selective

//                 // if (round%2==0) {
//                 //     stabilizers = {16,17,18,22,26,30,31,32};
//                 // } else {
//                 //     stabilizers = {19,20,21,23,24,25,27,28,29};
//                 // }
//                 // stabilizers_each_round[round] = stabilizers; //each round is this stabilizers

//                 // if (round%3==0) {
//                 //     stabilizers = {16,17,18,22,26,30,31,32};
//                 // } else if (round%3==1) {
//                 //     stabilizers = {20,23,25,28};
//                 // } else {
//                 //     stabilizers = {19,21,24,27,29};
//                 // }

//                 // if (round%4==0) {
//                 //     stabilizers = {16,17,18,19};
//                 // } else if (round%4==1) {
//                 //     stabilizers = {20,21,22,23};
//                 // } else if (round%4==2) {
//                 //     stabilizers = {24,25,26,27};
//                 // } else {
//                 //     stabilizers = {28,29,30,31,32};
//                 // }
//                 // stabilizers_each_round[round] = stabilizers; //each round is this stabilizers


//                 for (int i = 0; i < numRounds; i++) {
                    
//                     if (i%2==0) {
//                         stabilizers = {16,17,18,22,26,30,31,32, };
//                     } else {
//                         stabilizers = {19,20,21,23,24,25,27,28,29};
//                     }
                    
//                 }

//                 // if (round%10==0) {
//                 //     stabilizers = {16,17,18,22,26,30,31,32};
//                 // }else if (round%10==1) {
//                 //     stabilizers = {19};
//                 // } else if (round%10==2) {
//                 //     stabilizers = {20};
//                 // } else if (round%10==3) {
//                 //     stabilizers = {21};
//                 // } else if (round%10==4) {
//                 //     stabilizers = {23};
//                 // } else if (round%10==5) {
//                 //     stabilizers = {24};
//                 // } else if (round%10==6) {
//                 //     stabilizers = {25};
//                 // } else if (round%10==7) {
//                 //     stabilizers = {27};
//                 // } else if (round%10==8) {
//                 //     stabilizers = {28};
//                 // } else if (round%10==9) {
//                 //     stabilizers = {29};
//                 // }
//                 stabilizers_each_round[round] = stabilizers;

//                 // stabilizers_each_round = selective_map;

//                 // for (auto key : stabilizers_each_round) {
//                 //     std::cout<<key<<stabilizers_each_round[key]<<std::endl;
//                 // }

//                 std::tie(main_program, curr_program, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum) = sim.stability_syndrome_extraction_round(
//             4, stabilizers_each_round, frequency_stabilizers, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum, main_program);
//             }



//             // -- NEW ADPATIVE LOGIC -- 

//             //changing the above logic - new adaptive so skip the round if no error

//             //first syndrome extract just set to doing measure ALL

//             // if (round % 5 == 0) {
//             //     stabilizers = {16,17,18,22,26,30,31,32, 19,20,21,23,24,25,27,28,29};

//             //     stabilizers_each_round[round] = stabilizers; //each round is this stabilizers

//             //     std::tie(main_program, curr_program, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum) = sim.stability_syndrome_extraction_round(
//             //             4, stabilizers_each_round, frequency_stabilizers, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum, main_program);
//             //     goto hi;
//             // }

//             // if (!error) {
//             //     goto hi; //skip that
//             //     //keep the last rounds stabs and measure
//             // } else {
//             //     //if error then do ADAPTIVE SELECTIVE:

//             //     stabilizers.clear();

//             //     std::set<int> selective_stabs;

//             //     outFile5 << "----ROUND " << round << "----" << std::endl;

//             //     for (auto q : flipped_stabs) {
//             //         outFile5 << "THE FLIPPED STAB WAS: " << q << std::endl;
//             //         for (auto i : surrounding_stabs[q]) {
//             //             stabilizers.insert(i); //no dupes cuz set
//             //             outFile5 << "AND SO WE MEASURED: " << i << std::endl;
//             //         }
//             //     }


//             //     // if (round%2==0) {
//             //     //     stabilizers = {16,17,18,22,26,30,31,32};
//             //     // } else {
//             //     //     stabilizers = {19,20,21,23,24,25,27,28,29};
//             //     // }

//             //     // stabilizers = {16,17,18,22,26,30,31,32, 19,20,21,23,24,25,27,28,29};

//             //     stabilizers_each_round[round] = stabilizers; //each round is the stabilizers set

//             //     std::tie(main_program, curr_program, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum) = sim.stability_syndrome_extraction_round(
//             //                 4, stabilizers_each_round, frequency_stabilizers, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum, main_program);
//             // }`

//             //"RUN" the curr_program - only need measurement table loaded (done with run_batch and write_stats)
            
//             // std::cout<<"for syndrome round" << round << std::endl;

//             // for (int i = 0; i < G_RECORD_SPACE_SIZE; i++) {
//             //     if (sim.syndrome_table[i].popcnt() == 0) {
//             //         std::cout<<"WE GOT A 1 BROOOOO"<<std::endl;
//             //     }
//             // }
//             hi:

//             sim.run_program<FrameSimulator>(curr_program, 1); //1 shot run program again

//             // for (auto i : qontra::FrameSimulator::measures) {
//             //     std::cout<<i<<std::endl;
//             // }

//             //update the measure map per round, set previous to curr too

//             // for (auto i : qontra::FrameSimulator::measures) {
//             //     std::cout<<i<<std::endl;
//             // }

        
//             previous_q_to_m = q_to_m; //can't just set this BRUHHHH

//             int index = 0;

//             // std::cout<<"18 IS " << qontra::FrameSimulator::measures[2]<<std::endl;

//             // for (auto i : qontra::FrameSimulator::measures) {
//             //     std::cout << " WUT DO U MENA ITS " << i << std::endl;
//             // }

//             for (auto i : stabilizers) {
//                 q_to_m[i] = qontra::FrameSimulator::measures[index++];
//                 std::cout<< i << " was " << previous_q_to_m[i] << std::endl;
//                 std::cout<< i << " is " << q_to_m[i] << std::endl;
//             }

//             //check to set for error or not at the END of round makes more sense

//             // for (auto event : stabilizers) {
//             //     //check the stabs we measured this round to see if any r dif from their prev measures
//             //     if (q_to_m[event] != previous_q_to_m[event]) {
//             //         //THEN MEASURED ERROR!!!
//             //         error = 1; //say that there is an error
//             //         roundsSinceError = 0; //new error this round
//             //         break; //don't care if multiple errors - just checking presence is enough
//             //     } else {
//             //         if (roundsSinceError >= 3) {
//             //             error = 0; //or else no error
//             //         }
//             //     }
//             // }

//             //new adaptive logic - if error, let's skip the measurements and cnots altog
//             //so check every round

//             // flipped_stabs.clear(); //want for this round so
 
//             for (auto event : stabilizers) {
//                 if (q_to_m[event] != previous_q_to_m[event]) {
//                     //add that event to the flipped stabs set
//                     flipped_stabs.insert(event);
//                     error = 1;
//                     roundsSinceError = 0;
//                     break;
//                 } else {
//                     error = 0; //set error to 0 asap, not the roundssinceerror
//                 }
//             }
            

//             std::cout<<"we had an error " << error << std::endl;
            
//             roundsSinceError++;

//             roundNum++;


//             // std::ofstream outFile2("ASPmanualEvents.txt"); //saving to file to analyze
            
//             // outFile2 << "this is the previous: " << std::endl;
//             // for (auto i : previous_q_to_m) {
//             //     outFile2<<"qubit was "<< i << "and the measurement was " << previous_q_to_m[i]<<std::endl;
//             // }
//             // outFile2 << "this is the new one: " << std::endl;
//             // for (auto i : q_to_m) {
//             //     outFile2<<"qubit was "<< i << "and the measurement was " << q_to_m[i]<<std::endl;
//             // } 

//             // outFile2.close();
            

//             // write_stats(batchno); //batchno starts at 1, up 1 each time
//             //WILL WRITE EVENTS TO SYNDROME TABLE **

//             //see the syndrome table -> to file

//             // for (const auto& pair : stab_freqs) {
//             //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
//             // }

//             // for (int i = 0; i < G_RECORD_SPACE_SIZE; i++) {
//             //     // if (syndrome_table[i].popcnt() > 0) {
//             //     //     count++;
//             //     // }
//             //     outFile2 << sim.syndrome_table[i].popcnt()<<std::endl;
//             // }
            
//             // for (auto sel : selectors) {
//             //     outFile <<sel<<std::endl;
//             // }

//             // batchno++;
//         }
//         outFile5.close();

//         // std::cout<<qontra::FullSystemSimulator::syndrome_table[35].popcnt()<<std::endl;

//         //epilogue now

//         // std::cout<<"Wait broke out" << std::endl;

//         std::tie(main_program, curr_program, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum) = sim.stability_syndrome_epilogue(
//                 4, stabilizers_each_round, frequency_stabilizers, numEvents, numMeasurements, latestMeasurement, previousMeasurement, roundNum, main_program);

//         //now care about main_program - the incrementally built

//         // std::filesystem::create_directory("adaptive_programs_per_shot"); //create folder for each shot cuz circuit dif per shot

//         // std::string fileName = "adaptive_programs_per_shot/adaptiveStabilityProgram_shot_" + std::to_string(shot) + ".txt"; //file name
        
//         std::ofstream outFile("new_approach_adaptiveStabilityProgram.txt"); //saving to file to analyze
//         // std::ofstream outFile(fileName);

//         // for (const auto& pair : stab_freqs) {
//         //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
//         // }

//         for (auto instr : main_program) {
//             outFile << instr << std::endl;
//         }

//         // for (auto sel : selectors) {
//         //     outFile <<sel<<std::endl;
//         // }
//         outFile.close();

//         //and now finally running main_program below
        
//         //ALL OF BELOW CODE ALSO GOES WITHIN LOOP - DO SHOT BY SHOT CUZ ALL SHOT CIRCS DIF


//         // error model is with the control simulator: constant error model
//         // tables::ErrorAndTiming et2;
//         // et2 *= p*1000;

//         // const uint64_t n2 = get_number_of_qubits(main_program);
//         // tables::populate(n2, sim.config.errors, sim.config.timing, et2);

//         // //below code is for the control code -> error table for that so
//         // //not useful info below

//         // histogram_t<uint64_t> shots_hist = sim.run_program<FrameSimulator>(main_program, 1); //1 shot run
//         // histogram_t<double> norm_hist = histogram_normalize(shots_hist);
        

//         // if (world_rank == 0) {
//         //     std::cout << "Shots Histogram ------------------------\n"
//         //             << shots_hist << "\n"
//         //             << "Probability Histogram ------------------\n"
//         //             << norm_hist << std::endl;
//         // }
//     }

//     MPI_Finalize();
//     return 0; //post loop
    
// }


// -- resetting code back to regular qontra sim - no incremental generation

/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

// #include <qontra/ext/qes.h>
// #include <qontra/sim/base/clifford_sim.h>
// #include <qontra/sim/base/frame_sim.h>
// #include <qontra/sim/full_system_sim.h>

// #include <vtils/cmd_parse.h>
// #include <vtils/filesystem.h>
// #include <vtils/ini_parse.h>
// #include <vtils/timer.h>

// #include <cstdlib>

// using namespace qontra;
// using namespace vtils;


// int main(int argc, char* argv[]) {
//     MPI_Init(NULL, NULL);
//     int world_rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//     std::string help = "usage: ./qontrasim <qes-file> "
//                         "--shots <shots> "
//                         "--trace-out <trace-folder> "
//                         "--stim-out <stim-file> "
//                         "--data-out <data-file>\n"
                        
//                         "Optional:\n"
//                         "--rounds <numRounds> " //adding a --rounds param (default j 10d)
//                         "\t--config <ini-file>";


//     CmdParser pp(argc, argv, 1);
//     pp.help = help;
//     if (pp.option_set("h")) {
//         std::cerr << help << std::endl;
//         return 1;
//     }
//     configure_optimal_batch_size();

//     // Initialize simulator:
//     FullSystemSimulator sim;

//     // Get command line inputs.
//     std::string qes_file(argv[1]);

//     uint64_t shots;
//     std::string ini_file;
//     pp.get("shots", shots, true);
//     pp.get("trace-out", sim.config.syndrome_output_folder, true);
//     pp.get("stim-out", sim.config.stim_output_file, true);
//     pp.get("data-out", sim.config.data_output_file, true);


//     fp_t p = 0.0;
//     pp.get("p", p);
//     if (pp.get("config", ini_file)) {
//         std::string parent_dir = get_parent_directory(ini_file);

//         IniParser ini(ini_file);
//         const auto& ini_map = ini.get_ini_map();
//         // First load in the subroutines.
//         for (const auto& p : ini_map.at("__ANON__")) {
//             qes::Program<> subroutine = qes::from_file(parent_dir + "/" + p.second);
//             sim.load_subroutine(p.first, subroutine);
//         }
//         // Set simulation config if possible.
//         ini.get("Config", "record_events_until", sim.config.record_events_until);
//         ini.get("Config", "record_obs_until", sim.config.record_obs_until);
//         ini.get("SSE", "$reoz_track_last_n_events", sim.config.sse.rreoz_track_last_n_events);
//     }
//     // qes::Program<> main_program = qes::from_file(qes_file);

//     std::map<int, std::set<int>> stabilizers_each_round;

//     std::set<int> stabilizers;


//     int d = 3;
//     double err = 0.75;
//     int numRounds = 10*d;
//     if (pp.option_set("rounds")) {
//         //then set
//         pp.get("rounds", numRounds, true);
//     }
    
//     for (int i = 0; i < numRounds; i++) {
//         // if (i == numRounds-1) {
//         // }
//         // else {
//         //     stabilizers = {10, 11, 12, 14};
//         // }
//         stabilizers = {9,10,11,12,13,14,15,16}; //all in last round trying

//         // stabilizers = {25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48};

        
//         stabilizers_each_round[i] = stabilizers;
//     }

//     std::map<int, double> frequency_stabilizers; //empty map, pass into the control round

//     // qes::Program<> main_program = sim.create_program(5, stabilizers_each_round, frequency_stabilizers);

//     //-- **get the numRounds from command line so we can bash this


//     qes::Program<> main_program = sim.create_program_masked_selective(d, err, 0, 0, numRounds); //new param - 0 for z mem, 1 for x mem

//     //pass in numRounds too now

//     // Setup error model.
//     tables::ErrorAndTiming et;
//     et *= p*1000;

//     //print the et

//     const uint64_t n = get_number_of_qubits(main_program);
//     tables::populate(n, sim.config.errors, sim.config.timing, et);

//     qontra::tables::printErrorAndTiming(et); //its correct

//     histogram_t<uint64_t> shots_hist = sim.run_program<FrameSimulator>(main_program, shots);
//     histogram_t<double> norm_hist = histogram_normalize(shots_hist);

//     if (world_rank == 0) {
//         std::cout << "Shots Histogram ------------------------\n"
//                 << shots_hist << "\n"
//                 << "Probability Histogram ------------------\n"
//                 << norm_hist << std::endl;
//     }
//     MPI_Finalize();
//     return 0;
// }


// ------- INCREMENTAL GENERATION BELOW --------


#include <qontra/ext/qes.h>
#include <qontra/sim/base/clifford_sim.h>
#include <qontra/sim/base/frame_sim.h>
#include <qontra/sim/full_system_sim.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>
#include <vtils/ini_parse.h>
#include <vtils/timer.h>

//trying cache
#include <queue>
#include <deque>


using namespace qontra;
using namespace vtils;

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::string help = "usage: ./qontrasim <qes-file> "
                        "--shots <shots> "
                        "--trace-out <trace-folder> "
                        "--stim-out <stim-file> "
                        "--data-out <data-file>\n"
                        "Optional:\n"
                        "\t--config <ini-file>"
                        "--disconnect <disconnect-round>";
    CmdParser pp(argc, argv, 1);
    pp.help = help;
    if (pp.option_set("h")) {
        std::cerr << help << std::endl;
        return 1;
    }
    configure_optimal_batch_size();

    // Initialize simulator:
    FullSystemSimulator sim;

    // Get command line inputs.
    std::string qes_file(argv[1]);

    uint64_t shots;
    std::string ini_file;
    pp.get("shots", shots, true);
    pp.get("trace-out", sim.config.syndrome_output_folder, true);
    pp.get("stim-out", sim.config.stim_output_file, true);
    pp.get("data-out", sim.config.data_output_file, true);
    

    fp_t p = 0.0;
    pp.get("p", p);
    if (pp.get("config", ini_file)) {
        std::string parent_dir = get_parent_directory(ini_file);

        IniParser ini(ini_file);
        const auto& ini_map = ini.get_ini_map();
        // First load in the subroutines.
        for (const auto& p : ini_map.at("__ANON__")) {
            qes::Program<> subroutine = qes::from_file(parent_dir + "/" + p.second);
            sim.load_subroutine(p.first, subroutine);
        }
        // Set simulation config if possible.
        ini.get("Config", "record_events_until", sim.config.record_events_until);
        ini.get("Config", "record_obs_until", sim.config.record_obs_until);
        ini.get("SSE", "$reoz_track_last_n_events", sim.config.sse.rreoz_track_last_n_events);
    }
    // qes::Program<> main_program = qes::from_file(qes_file);

    qes::Program<> main_program; //add to this
    std::set<int> allStabs;
    std::set<int> innerStabs;
    std::set<int> zStabs;
    std::set<int> xStabs;
    std::set<int> selective_stabs;


    //initial params
    int d = 3;
    int numEvents = 0;
    int numMeasurements = 0;
    std::map<int, int> latestMeasurement;
    std::map<int, int> previousMeasurement;
    std::vector<int> firstRoundMeasurements;
    int roundNum = 0;
    int numRounds = 10*d;

    std::map<int, std::vector<int>> p_to_d_map = sim.create_p_to_d_map(d);

    std::map<int, std::set<int>> p_to_d;

    //nice to play around with this cache idea
    std::map<int, int> sse_cache; //stab to int cache 

    std::set<int> hold_out; //hold these stabs when choosing for next round


    // std::queue<std::pair<int, int>> sse_cache; 
    //queue so we can add to front, remove from back when size >= selective_stabs amount


    //cache size will only be num of selective_stabs large, so that's all we can augment with
    //evict oldest lol

    //GENERALIZE THIS SHII - generalized

    for (int i = pow(d,2); i < 2*(pow(d,2)); i++) {
        std::set<int> curr;
        for (auto q : p_to_d_map[i]) {
            curr.insert(q);
        }
        p_to_d[i] = curr;
    }

    std::map<int, std::set<int>> surrounding_stabs;

    
    for (int i = pow(d,2); i < 2*(pow(d,2)); i++) {
        std::set<int> curr_data = p_to_d[i]; //check where else has this

        std::set<int> surrounding;

        for (int j = pow(d,2); j < 2*(pow(d,2)); j++) {
            if (i != j) {
                //use set intersection - see if any data qubits in common
                std::set<int> intersect; //common elements in here

                std::set_intersection(
                    curr_data.begin(), curr_data.end(), p_to_d[j].begin(), p_to_d[j].end(), std::inserter(intersect, intersect.begin())
                );
                
                //if any common data qubit, add j
                if (!intersect.empty()) {
                    surrounding.insert(j);
                }
            }

            if (j < (pow(d,2) + 2*(d-1))) {
            //outer
                if (j%2 == 1) {
                    surrounding.insert(j+1);
                } else {
                    surrounding.insert(j-1); //in pairz
                }
            }
        }
        if (i < (pow(d,2) + 2*(d-1))) {
            //outer
                if (i%2 == 1) {
                    surrounding.insert(i+1);
                } else {
                    surrounding.insert(i-1); //in pairz
                }
            }
        //now added for all data q of curr
        surrounding_stabs[i] = surrounding;
    }

    // for (const auto& [key, value_set] : surrounding_stabs) {
    //     std::cout << "Key: " << key << ", Values: {";
    //     for (const auto& val : value_set) {
    //         std::cout << val << " ";
    //     }
    //     std::cout << "}\n";
    // }



    // std::set<int> selective_stabs = {9,10,11,12,13,14,15,16};
    
    //init to 0 for prologue instruction number
    int instruction_number = 0;

    std::tie(allStabs, innerStabs, zStabs, xStabs, main_program, instruction_number) = sim.create_masked_prologue(d, main_program, instruction_number);

    // Setup error model.
    tables::ErrorAndTiming et;
    et *= p*1000;

    const uint64_t n = get_number_of_qubits(main_program);
    tables::populate(n, sim.config.errors, sim.config.timing, et);

    std::map<int, int> q_to_m;
    std::map<int, int> previous_q_to_m;
    qes::Program<> currProgram;

    bool error = false; //presence of error (if no error, just fix to 13,14,15,16 for now, let's see)

    std::random_device rand; //rand
    std::mt19937 gen(rand()); //seed
    // std::uniform_int_distribution<> dist(pow(d,2), 2*(pow(d,2)-1)); //d^2 to 2(d^2-1) inclusive

    std::vector<int> weights;
    for (int i = pow(d,2); i <= 2*(pow(d,2)-1); i++) {
        if (i >= (pow(d,2) + 2*(d-1))) {
            //inner stab
            weights.push_back(1); // regular as inner
        } else {
            weights.push_back(1); // Regular weight
        }
    }

    // std::discrete_distribution<> dist(weights.begin(), weights.end());

    std::uniform_int_distribution<int> dist(pow(d,2), 2*((pow(d,2))-1));

    std::uniform_real_distribution<double> dist2(0.20, 0.80);

    bool canskip = false;

    int roundSkipped = 0;

    int totalRoundsSkipped = 0;


    std::ofstream outFileMeta("round_measure_deets.txt");

    std::string selective_policy_used;

    int rSinceErr = 0;

    int dontusefixed = 0;

    std::set<int> save_selective_stabs; //to save for the "redundancy adding round"

    std::srand(std::time(0));
    
    int num_fixed_pols = 0;
    int num_rand_pols = 0;

    //file for skipped rounds info

    std::ofstream outFileSkippedRounds("skipped_rounds_data.txt");

    for (int r = 0; r < 10*d; r++) {

        double random_p = dist2(gen); //get a random_p per round now

        int num_stabs = int((pow(d,2)-1) * (random_p)); //make it int to get num of total stabs based on p

        outFileMeta << "round " << r << std::endl;
        outFileMeta << "we measured " << num_stabs << " stabilizers " << std::endl;


        //if we r holding out more than num_stabs, then just clear hold out and skip the round?

        if (r % 2  == 0) {
            canskip = false;
        } else {
            canskip = true;
        }

        // canskip = false; //do all rn
       
        if (canskip) {
            roundSkipped += 1;
            totalRoundsSkipped += 1;
            //PUT THE TIMING ERROR HERE
        
            if (roundSkipped == 1) {
                canskip = false;
                roundSkipped = 0; //reset, so now we allowed to skip more (double chance)
            }

            //cant just add timing error - so we gotta create a noiseless skipped round circuit
            //then roll in custom error model with from_qes func edits 
            
            //k so here - let's save the instruction number before skipped round to a file
            outFileSkippedRounds << instruction_number << std::endl;

            
            // canskip = false;
            std::cout<<"we skipping this " << r << std::endl;
            continue;
        }

        
        // if (hold_out.size() >= num_stabs-1) {
        //     hold_out.clear();
        //     std::cout<<"we skipping this " << r << std::endl;

        //     continue;
        // }

        hold_out.clear();

        if (!error) {
            rSinceErr++; //inc rounds since error
        }

        // for (auto stab : selective_stabs) {
        //     std::cout << stab << ", ";
        // }

        // std::cout << std::endl;

        // int random_bit = (random_p < 0.50) ? 0 : 1;


        if (!error || dontusefixed%2 == 0) {
            selective_policy_used = "random";
            // outFileMeta << "random in round " << r << std::endl;
            //get a random RANDOM P mask -- code from the mask prog in full_system_inl
            
            selective_stabs.clear();

            while (selective_stabs.size() < num_stabs) {
                //needa do less than not != cuz it might go over by like 1 but chill
                //unique, so will get w out repeats
                
                int random_stab = dist(gen); //got rid of + pow(d,2)

                // std::cout<<random_stab<<std::endl;

                while (hold_out.contains(random_stab)) {
                    // for (auto i : hold_out) {
                    //     std::cout << i << std::endl;
                    // }
                    // std::cout<< "hold out in play!!" << random_stab << std::endl;

                    std::uniform_int_distribution<int> dist(pow(d,2), 2*((pow(d,2))-1));
                    random_stab = dist(gen);
                }

                //do the outer thing tho (if even outer add the -1, if odd outer add the +1)
                if (random_stab < (pow(d,2) + 2*(d-1))) {
                    if (random_stab%2==0) {
                        selective_stabs.insert(random_stab-1);
                    } else {
                        selective_stabs.insert(random_stab+1);
                    }
                }

                selective_stabs.insert(random_stab);

            }
            
            // std::cout<<"no error in round " << r <<std::endl;
        } else {
            //if error before, now set error to false for this round
            // if (!error) {
            //     selective_stabs = save_selective_stabs;
            //     selective_policy_used = "fixed cuz error 2 rounds ago";
            // }
            error = false; //for next time, and lets use the regular selective stabs (surrounding) thing
        }

        if (selective_policy_used == "random") {
            num_rand_pols++;
        } else {
            num_fixed_pols++;
        }

        outFileMeta << "selective policy used: " << selective_policy_used << std::endl;

        if (r%(2*d) == 0) {
            // selective_stabs = {9,10,11,12,13,14,15,16}; // 3
            // selective_stabs = {25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48}; // 5
            // selective_stabs = {49,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48}; // 5

            //  {49,50,51,52,53,54,62,64,66,67,69,71,74,76,78,79,81,83,86,88,90,91,93,95,55,56,57,58,59,60,61,63,65,68,70,72,73,75,77,80,82,84,85,87,89,92,94,96,55,56,57,58,59,60,61,63,65,68,70,72,73,75,77,80,82,84,85,87,89,92,94,96,49,50,51,52,53,54,62,64,66,67,69,71,74,76,78,79,81,83,86,88,90,91,93,95}; //7

            for (int i = pow(d,2); i <= 2*(pow(d,2)-1); i++){
                selective_stabs.insert(i);
            }
            //for 9 - nah?!

        }
        
        outFileMeta << "we measuring these stabs: ";

        for (auto stab : selective_stabs) {
            outFileMeta << stab << ", ";
        }

        outFileMeta << std::endl;

        // selective_stabs = {25,26,29,30,31,45,38,39,47,35};

        // for (auto i : qontra::FrameSimulator::measures) {
        //     std:: cout << "measure is " << i << std::endl;
        // }

        selective_stabs = {9,10,11,12,13,14,15,16}; //set to all

        int disconnectRound = 0;

        pp.get("disconnect-round", disconnectRound, true); //ok
        
        std::tie(main_program, numEvents, numMeasurements, roundNum, previousMeasurement, latestMeasurement, firstRoundMeasurements, currProgram, instruction_number) = 
            sim.create_masked_se_round(d, main_program, numEvents, numMeasurements, roundNum, latestMeasurement, previousMeasurement, selective_stabs, allStabs, innerStabs, xStabs, zStabs, firstRoundMeasurements, instruction_number, disconnectRound);

        
        // int index2 = 0;
        // for (int i = 0; i < allStabs.size(); i++) {
        //     //for all stabs right?
        //     std::cout << qontra::FrameSimulator::measures[index2++] << std::endl;
        //     // q_to_m[i] = qontra::FrameSimulator::measures[index2++];
        // }

        

        // std::cout<<"AFTER"<<std::endl;
        histogram_t<uint64_t> shots_hist = sim.run_program<FrameSimulator>(main_program, 1); //1 shot - for currprogram ya?

        //--checking obs to get p(error)

        //running above so now check obs

        // std::cout << "round " << r << " obs " << std::endl;
        // std::cout << "num measures is " << numMeasurements << std::endl;

        // size_t num_rows = sim.syndrome_table.num_major_bits_padded();
        // size_t num_cols = sim.syndrome_table.num_minor_bits_padded();

        // for (size_t row = 0; row < num_rows; row++) {
        //     for (size_t col = 0; col < num_cols; col++) {
        //         if (sim.syndrome_table[row][col] == 1) {
        //             std::cout << "hi " << row << " " << col << std::endl;
        //         }
        //     }
        //     // std::cout << std::endl;
        // }

        // std::cout << "num error" << sim.logical_errors << std::endl;

        // std::cout << "round " << r << " obs copy " << std::endl;

        // num_rows = sim.observable_table_cpy.num_major_bits_padded();
        // num_cols = sim.observable_table_cpy.num_minor_bits_padded();

        // for (size_t row = 0; row < num_rows; row++) {
        //     for (size_t col = 0; col < num_cols; col++) {
        //         std::cout << (sim.observable_table_cpy[row][col] ? '1' : '0') << " ";
        //     }
        //     std::cout << std::endl;
        // }


        // std::cout << qontra::FrameSimulator::x_table << std::endl;

        // ----- cache logic -----
        int amount_stabs = selective_stabs.size();

        //compare with cache: if even # stabs, then take out of selective stabs

        int index = 0;

        // for (auto i : qontra::FrameSimulator::measures) {
        //     std:: cout << "measure is " << i << std::endl;
        // }

        
        for (auto stab : selective_stabs) {

            int countOfFlips = 0;
            if (qontra::FrameSimulator::measures[index++] == 1) {
                // std::cout << " got a flip " << std::endl;

                //then this stab is flip so 
                q_to_m[stab] = 1;

                countOfFlips++; //if curr flipped then +1 flips type
            } else {
                q_to_m[stab] = 0;
            }

            std::set<int> neighbors = surrounding_stabs[stab];

            for (const auto& pair : sse_cache) {
                if (neighbors.contains(pair.first) && pair.second == 1) {
                    //if this is a neighbor key and its a 1, then add 1 to the count
                    countOfFlips++;
                }
            }

            if (countOfFlips%2 == 0) {
                //skipping logic nice, cuz if no error on this AND neighbor then we'll say don't choose this in next random stab selection
                hold_out.insert(stab);
            }
        }


        index = 0;
        for (auto stab : selective_stabs) {
            sse_cache[stab] = qontra::FrameSimulator::measures[index++];
        }
        //cache has new rounds measures now

        std::set<int> flipped_stabs;

        // for (const auto& pair : q_to_m) {
        //     std::cout << pair.first << " -> " << pair.second << "\n";
        // }

        // for (const auto& pair : previous_q_to_m) {
        //     std::cout << pair.first << " -> " << pair.second << "\n";
        // }


        for (auto event : allStabs) {
            if (q_to_m[event] != previous_q_to_m[event]) {
                //add that event to the flipped stabs set
                // std::cout<<"HIII the flipped was " << event <<std::endl;
                flipped_stabs.insert(event);
                // std::cout<<event<<std::endl;
                // std::cout<<"pass"<<std::endl;
                dontusefixed++;
                error = true;
                rSinceErr = 0; //reset this
            }
        }

        previous_q_to_m = q_to_m; //set right?

        //save the selective_stabs:

        // save_selective_stabs = selective_stabs;

        //selective stabs shud be the flipped_stabs + their surrounding stabs

        if (flipped_stabs.size() >= (d-1)) {
            //min cuz of how the outer ones work
            selective_policy_used = "only flipped";
            selective_stabs = flipped_stabs;
        } else {
            selective_policy_used = "flipped + surrounding";
            selective_stabs.clear(); //recreate
            for (auto stab : flipped_stabs) {
                selective_stabs.insert(stab);
                for (auto neighbor : surrounding_stabs[stab]) {
                    // std::cout<<stab<< "WITH " << neighbor<<std::endl;
                    selective_stabs.insert(neighbor);
                }
            }
        }
        
        // for (auto i : selective_stabs) {
        //     std::cout<<i<<std::endl;
        // }
        // std::cout<<"pass"<<std::endl;

        // canskip = !error;
    }

    outFileMeta << "random policy used times: " << num_rand_pols << std::endl;
    outFileMeta << "fixed policy used times: " << num_fixed_pols << std::endl;


    //surrounding stabs logics fricked idk if its reading measures properly

    main_program = sim.create_masked_epilogue(d, main_program, numEvents, numMeasurements, roundNum, latestMeasurement, previousMeasurement, selective_stabs, allStabs, innerStabs, xStabs, zStabs, firstRoundMeasurements, instruction_number);

    //after epilogue obs shud become one??

    // std::cout << "num error" << sim.logical_errors << std::endl;

    // std::cout << "Statistic Data: ";
    // for (size_t i = 0; i < sim.logical_errors.size(); i++) {
    //     std::cout << sim.logical_errors.at(i) << " ";
    // }
    // std::cout << std::endl;

    outFileSkippedRounds.close();

    std::ofstream outFile("sse_cache.txt");

    for (auto instr : main_program) {
        outFile <<instr << std::endl;
    }

    outFile.close();

    std::cout << "we skipped total rounds: "<< totalRoundsSkipped << std::endl;

    outFileMeta.close();


    // histogram_t<double> norm_hist = histogram_normalize(shots_hist);

    // if (world_rank == 0) {
    //     std::cout << "Shots Histogram ------------------------\n"
    //             << shots_hist << "\n"
    //             << "Probability Histogram ------------------\n"
    //             << norm_hist << std::endl;
    // }





    MPI_Finalize();
    return 0;
}