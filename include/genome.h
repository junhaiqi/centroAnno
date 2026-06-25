#ifndef GENOME_H
#define GENOME_H
#include <iostream>
#include <utility> 
#include "sequenceUtils.h"
#include "monomer.h"
#include "monoRefine.h"
#include "hor.h"
#include <stdexcept>

void repeatRegionInference(SEQ &sequence, const int &k,
                           const float &repPorCutoff, const int &stepLen,
                           const bool &hpc,
                           std::vector<std::pair<int, int>> &repRegions, const int &maxRegionLength);

void refineDemcomposeGenome(const SEQ &longSequence, vector<string> &newMonoTems,
                      vector<MonomerAlignment> &newDemcompBlockRes,
                      const int &threadNum, const bool &hpc, const int &k,
                      const float &fpsCutoff, const float &repPorCutoff);

void findGoodDemcomRegion(const vector<MonomerAlignment> &newDemcompBlockRes,
                         vector<region> &goodRegions, const int &minRegionLength,
                         const float &thresCutoff);

#endif
