#include "genome.h"

void repeatRegionInference(SEQ &sequence, const int &k,
                           const float &repPorCutoff, const int &stepLen,
                           const bool &hpc,
                           std::vector<std::pair<int, int>> &repRegions) {

  for (size_t i = 0; i < sequence.seq.length() - stepLen; i += stepLen) {
    std::string thisSeq = sequence.seq.substr(i, stepLen);
    KC kCounter(thisSeq, k); // kmer counter
    if (hpc) {
      SEQ thisSeq_S(thisSeq, "sub");
      thisSeq_S.compressSeq();          // HPC
      kCounter.seq = thisSeq_S.compSeq; // kmer counter
    }

    bool rep = kCounter.checkRepeative(repPorCutoff);

    if (rep) {
      std::pair<int, int> region(i, i + stepLen);
      if (!repRegions.empty()) {
        if (repRegions.back().second - i < 10) {
          repRegions.back().second = i + stepLen;
        } else
          repRegions.push_back(region);
      } else
        repRegions.push_back(region);
    }
  }
}