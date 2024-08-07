#ifndef GENOME_H
#define GENOME_H
#include <iostream>
#include <utility> 
#include "sequenceUtils.h"
#include "monomer.h"

void repeatRegionInference(SEQ &sequence, const int &k,
                           const float &repPorCutoff, const int &stepLen,
                           const bool &hpc,
                           std::vector<std::pair<int, int>> &repRegions);

#endif
