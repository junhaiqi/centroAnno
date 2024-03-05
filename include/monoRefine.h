#ifndef MONOREFINE_H
#define MONOREFINE_H
#include "monomer.h"
#include "sampleDBSCAN.h"
#include "lightWightMSA.h"
#include "moutput.h"
#include <cmath>

typedef robin_hood::unordered_map<int, std::vector<std::string>> clusters;
void sampleClustering(const SEQ &sequence,
                      const vector<MonomerAlignment> &demcompBlockRes,
                      vector<int> &clusterLabels,
                      vector<string> &consensuses, // new monomer template.
                      size_t threadNum, const float epsilon = 0.90,
                      const int minSamples = 5);

void refineDemcompose(const SEQ &longSequence,
                      vector<string> &newMonoTems,
                      vector<MonomerAlignment> &newDemcompBlockRes,
                      const int &threadNum,
                      const bool &hpc, const int &k,
                      const float &fpsCutoff,
                      const float &repPorCutoff,
                      const std::string &demResfilename,
                      const std::string &inferedMonosfilename);

typedef struct REGION {
  int start;
  int end;
  int blockStartID;
  int blockEndID;
  REGION() = default;
  REGION(int start_, int end_, int blockStartID_, int blockEndID_)
      : start(start_), end(end_), blockStartID(blockStartID_),
        blockEndID(blockEndID_) {}
} region;

void findBadDemcomRegion(const vector<MonomerAlignment> &newDemcompBlockRes, 
                         vector<region> &badRegions,
                         const int &minRegionLength = 5000,
                         const float &thresCutoff = 0.85);

void findMinBadRegion(vector<region> &badRegions, region &minRegion, int &minLen);

inline void correctPosOfBadRegionDemRes(vector<MonomerAlignment> &thisRes, const int &offset){
    for (MonomerAlignment &m : thisRes){
        m.start_pos += offset;
        m.end_pos += offset;
    }
}

void updateIdentityByEditDist(const SEQ &longSequence,
                              vector<string> &newMonoTems,
                              vector<MonomerAlignment> &newDemcompBlockRes);

void checkBadBlockNums(const SEQ &sequence,
                       const std::vector<MonomerAlignment> &newDemcompBlockRes,
                       bool &badSignal, std::vector<std::string> &badBlocks,
                       std::vector<std::string> &monoTemsForBadBlocks,
                       const size_t &threadNum, const float epsilon = 0.9,
                       const int minSamples = 10);

inline void sortMap(const robin_hood::unordered_map<int, int> &horLenMap,
                    std::vector<std::pair<int, int>> &keyValuePairs) {

    // std::vector<std::pair<int, int>> keyValuePairs;
    keyValuePairs.reserve(horLenMap.size());
    for (const auto &kvp : horLenMap) {
        keyValuePairs.emplace_back(kvp.first, kvp.second);
    }

    std::sort(
        keyValuePairs.begin(), keyValuePairs.end(),
        [](const std::pair<int, int> &lhs, const std::pair<int, int> &rhs) {
          return lhs.second > rhs.second;
        });

    for (const auto &kvp : keyValuePairs) {
        std::cout << "Key: " << kvp.first << ", Value: " << kvp.second
                  << std::endl;
    }
}

inline float calculateMean(const std::vector<std::string> &sequences) {
    double sum = 0.0;
    for (const std::string &seq : sequences) {
        sum += seq.length();
    }
    return sum / sequences.size();
}

inline float calculateStandardDeviation(const std::vector<std::string> &sequences, const float &mean) {
    float sum = 0.0;
    for (const std::string &seq : sequences) {
        float diff = seq.length() - mean;
        sum += diff * diff;
    }
    return std::sqrt(sum / sequences.size());
}

inline void removeOutliers(std::vector<std::string> &sequences) {
    float mean = calculateMean(sequences);
    float stddev = calculateStandardDeviation(sequences, mean);
    float upperThreshold = mean + 3.0 * stddev; // 3-sigma
    float lowerThreshold = mean - 3.0 * stddev;

    auto iter = std::remove_if(sequences.begin(), sequences.end(),
                               [upperThreshold, lowerThreshold](const std::string &seq) {
                                 return (seq.length() > upperThreshold ||
                                         seq.length() < lowerThreshold);
                               });

    sequences.erase(iter, sequences.end());
}

#endif