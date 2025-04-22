//write slw decoder code here

#include "qontra/decoder/sliding_window.h"

#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <fstream>

namespace qontra {

    // SLWDecoder::SLWDecoder(DetailedStimCircuit circ, size_t window_size, size_t step_size):
    //                         Decoder(circ), w(window_size), s(step_size), b(w-s), inner_decoder(circ), total_correction(circ.count_observables())
    //                 {
    //                     //empty constructor body is chill
    //                 }

    

    //using template from mwpm decode_error
    Decoder::result_t
    SLWDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {

        total_correction.clear(); //clear it init

        const size_t num_dets = circuit.count_detectors();

        std::vector<Decoder::assign_t> error_assignments; //NEEDA ASSIGN THESE - its chill not needed

        timer.clk_start(); //timer for result

        size_t lo = 0;

        std::string before = total_correction.str();

        while (lo < num_dets) {
            //slide window n decode each individual

            // std::cout << "lo is " << lo << std::endl;

            // std::cout << "num_dets is " << num_dets << std::endl;

            const size_t hi = std::min(lo + w, num_dets); //for final one need to do min right, cuz dont wanna go past

            //decode window
            decode_window(lo, hi, syndrome);

            //prop syndrome
            propagate_syndrome(lo, hi, syndrome); //based on boundaries logic

            lo += s; //next window based on step size
            // std::cout << "total AFTER correciton vec THIS WINDOW is " << total_correction.str() << std::endl;
        }

        fp_t t = (fp_t)timer.clk_end(); //end clock



        //pack time, correction, error_assignments (dont need) into result object
        return {
            t, total_correction, error_assignments
        };
    }

    //ok finally

    //indiv window decoding func

    void SLWDecoder::decode_window(size_t lo, size_t hi, stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
        //u have whole syndrome but u j looking at window

        stim::simd_bits<SIMD_WIDTH> curr_window(hi - lo);
        curr_window.clear();
        //populate window with syndrome values
        for (int i = 0; i < (hi-lo); i++) {
            curr_window[i] = syndrome[lo + i]; //start at lo in syndrome
        }

        //decode just that piece with MWPM now

        Decoder::result_t window_res = inner_decoder.decode_error(curr_window); //pass only this syndrome in

        //merge this corr into total correction

        stim::simd_bits<SIMD_WIDTH> curr_correction = window_res.corr; //.corr for correction

        std::cout << "curr corr MWPM is " << curr_correction << std::endl;

        if (curr_correction.popcnt() != 0) {
            //any ones, then we xor the msb
            //mwpm based on this syndrome window said to flip obs
            total_correction[0] = !total_correction[0];
        }

        // for (int i = 0; i < (hi-lo); i++) {
        //     if (curr_correction == 1) {
        //         //then xor it INTO THE MSB
                
        //     }

        // }

    }

    void SLWDecoder::propagate_syndrome(size_t lo, size_t hi, stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {

        //ok so here's what we doing!!!

        //so lprop is lo+s and that is the overlapping boundary between the core and the buffer. 
        //anything from lo to lo+s is in the core so it gets committed and zeroed, and 
        //anything past is the buffer so it stays tentatively to be updated next time. 
        //we ensure consistency by hard setting the lo+s as the end of the core this time, 
        //so next time it is checked as the beginning of the buffer, that correction isn't recomputed

        if (lo + s >= hi) {
            return; //cuz last window u dont needa prop anything
        }

        const size_t prop = lo + s;

        syndrome[prop] ^= total_correction[prop];
    }







} //qontra type