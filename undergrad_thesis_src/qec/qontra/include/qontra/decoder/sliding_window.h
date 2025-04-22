//header for sliding window decoder


#ifndef SLW_DECODER_H
#define SLW_DECODER_H

#include <stim.h>
#include "qontra/decoder.h"
#include "qontra/graph/decoding_graph.h"
#include "/Users/aryan/Desktop/quantumResearch/qec/qontra/include/qontra/decoder/pymatching.h"
#include "stim/mem/simd_bits.h"
#include "stim/mem/simd_bits_range_ref.h"



namespace qontra {

class SLWDecoder : public Decoder {
public:
    //constructor for sliding window decoder
    SLWDecoder(DetailedStimCircuit circ, size_t window_size, size_t step_size)
        :Decoder(circ, graph::DecodingGraph::Mode::DO_NOT_BUILD), //that shud be mode right it said sum in readme
        w(window_size), s(step_size), b(w -s), inner_decoder(circ), total_correction(circ.count_observables()) //pass in all these params to init vars
    {}


    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) override;

private:
  //vars that we need

  const size_t w; //window
  const size_t s; //step
  const size_t b; //buffer


  //inner decoder, we want PyMatching so instantiate that

  PyMatching inner_decoder;

  stim::simd_bits<SIMD_WIDTH> total_correction;


  void decode_window(size_t lo, size_t hi, stim::simd_bits_range_ref<SIMD_WIDTH> syndrome); //decode each window

  void propagate_syndrome(size_t lo, size_t hi, stim::simd_bits_range_ref<SIMD_WIDTH> syndrome); 
  //prop syndrome so consistent corrections bw windows

};

}   // qontra

#endif
