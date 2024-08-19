#include "genome.h"


void repeatRegionInference(SEQ &sequence, const int &k,
                           const float &repPorCutoff, const int &stepLen,
                           const bool &hpc,
                           std::vector<std::pair<int, int>> &repRegions, const int &maxRegionLength) {
  
  // const int maxRegionLength = 500000;  // Define the maximum allowed length for a region
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
        if (repRegions.back().second - i < 10000) {
          repRegions.back().second = i + stepLen;
        } else
          repRegions.push_back(region);
      } else
        repRegions.push_back(region);
      
      while (repRegions.back().second - repRegions.back().first > maxRegionLength) {
          int start = repRegions.back().first;
          int end = start + maxRegionLength;
          repRegions.back().first = end;  // Adjust the start of the current region
          repRegions.insert(repRegions.end() - 1, std::make_pair(start, end));
      }
    }
  }
}

// decompose a sequence into blocks which have high indentity with our infered
// monomer templates, and find the failed blocks to redecompose.
void refineDemcomposeGenome(const SEQ &longSequence, vector<string> &newMonoTems,
                      vector<MonomerAlignment> &newDemcompBlockRes,
                      const int &threadNum, const bool &hpc, const int &k,
                      const float &fpsCutoff, const float &repPorCutoff) {

  if (newMonoTems.empty()) {
    std::cout << longSequence.seqName << " (length: " << longSequence.seq.size()
              << ") not a repeative seuence (because of low non-unique k-mers)."
              << std::endl;
    return;
  }

  std::vector<SEQ> monos;
  monos.resize(newMonoTems.size());
  for (size_t i = 0; i < monos.size(); i++) {
    SEQ mono(newMonoTems[i], to_string(i));
    monos[i] = mono;
  }
  demcomposer2(longSequence, monos, newDemcompBlockRes, 1);
  vector<region> badRegions;
  findBadDemcomRegion(newDemcompBlockRes, badRegions);
  printf("Number of bad regions: %zu\n", badRegions.size());

  if (!badRegions.empty()) {
    for (size_t t = 0; t < 1;
         t++) { // Iterate once, for speed, we can adjust it.
      robin_hood::unordered_map<int, region> hashRegions;
      robin_hood::unordered_map<int, vector<MonomerAlignment>>
          maln; // region.start map to its new decomposed result.
      for (const region &r : badRegions) {
        hashRegions[r.start] = r;
      }
      vector<SEQ> monosUsedForBadRegion;
      printf("Parse bad regions...\n");
      while (!badRegions.empty()) {
        region minRegion;
        int minLen;
        findMinBadRegion(badRegions, minRegion, minLen);
        const int lenForInferMonos = minLen > 300000 ? 300000 : minLen;
        SEQ badRegionSeq(longSequence.seq.substr(minRegion.start, minLen),
                         longSequence.seqName);
        SEQ segForMonosInfer(badRegionSeq.seq.substr(0, lenForInferMonos),
                             badRegionSeq.seqName);
        std::vector<std::string> tempMonoTemsForBadRegions;
        vector<SEQ> monosForBadRegion;
        monoTempInference(segForMonosInfer, k, fpsCutoff, 4,
                          tempMonoTemsForBadRegions, monosForBadRegion, hpc,
                          repPorCutoff);
        if (!monosForBadRegion.empty()) {
          std::vector<MonomerAlignment> tempDemcompBlockRes;
          demcomposer2(segForMonosInfer, monosForBadRegion, tempDemcompBlockRes,
                       1);
          std::vector<std::string> refinedMonoTemplates;
          vector<int> clusterLabels;
          printf("Clustering for refined monomer template...\n");
          sampleClustering(segForMonosInfer, tempDemcompBlockRes, clusterLabels,
                           refinedMonoTemplates, threadNum);
          printf("Refined templates have been generated, redecompose...\n");
          size_t currentMonosNum = newMonoTems.size();
          for (size_t n = 0; n < refinedMonoTemplates.size(); n++) {
            SEQ thisTem(refinedMonoTemplates[n],
                        to_string(currentMonosNum + n));
            monosUsedForBadRegion.push_back(thisTem);
          }
          newMonoTems.insert(newMonoTems.end(), refinedMonoTemplates.begin(),
                             refinedMonoTemplates.end());
          std::vector<MonomerAlignment> reDemcompBlockResForThisBadRegion;
          demcomposer2(badRegionSeq, monosUsedForBadRegion,
                       reDemcompBlockResForThisBadRegion, 1);
          correctPosOfBadRegionDemRes(reDemcompBlockResForThisBadRegion,
                                      minRegion.start);
          maln[minRegion.start] = reDemcompBlockResForThisBadRegion;
          if (!badRegions.empty()) {
            for (auto r = badRegions.rbegin(); r != badRegions.rend(); ++r) {
              SEQ otherBadRegionSeq(
                  longSequence.seq.substr((*r.base()).start,
                                          (*r.base()).end - (*r.base()).start),
                  longSequence.seqName);
              vector<MonomerAlignment> otherRegionDecomRes;
              demcomposer2(otherBadRegionSeq, monosUsedForBadRegion,
                           otherRegionDecomRes, 1);
              vector<region> thisSonBadRegions;
              findBadDemcomRegion(otherRegionDecomRes, thisSonBadRegions);
              if (thisSonBadRegions.empty()) {
                correctPosOfBadRegionDemRes(otherRegionDecomRes,
                                            (*r.base()).start);
                maln[(*r.base()).start] = otherRegionDecomRes;
                badRegions.erase(std::next(r).base());
              }
            }
          }
        }
      }
      // update the decompose result
      vector<MonomerAlignment> finalRefinedNewDemcompBlockRes;
      size_t currentBadRegionBound = 0;
      for (size_t i = 0; i < newDemcompBlockRes.size(); i++) {
        if (currentBadRegionBound > newDemcompBlockRes[i].start_pos) {
          continue;
        }
        if (maln.find(newDemcompBlockRes[i].start_pos) == maln.end()) {
          finalRefinedNewDemcompBlockRes.push_back(newDemcompBlockRes[i]);
        } else {
          vector<MonomerAlignment> currentDemcompBlockRes =
              maln[newDemcompBlockRes[i].start_pos];
          for (const MonomerAlignment &m : currentDemcompBlockRes) {
            finalRefinedNewDemcompBlockRes.push_back(m);
          }
          currentBadRegionBound =
              hashRegions[newDemcompBlockRes[i].start_pos].end;
        }
      }
      newDemcompBlockRes = finalRefinedNewDemcompBlockRes;
      badRegions = {};
      findBadDemcomRegion(newDemcompBlockRes, badRegions);
    }
  }

  bool badSignal = false;
  std::vector<std::string> badBlocks = {};
  std::vector<std::string> monoTemsForBadBlocks = {};
  checkBadBlockNums(longSequence, newDemcompBlockRes, badSignal, badBlocks,
                    monoTemsForBadBlocks, threadNum);

  if (badSignal) {
    printf("Too many bad blocks were found, reinfer monomer templates and "
           "redecompose...\n");
    
    newMonoTems.insert(newMonoTems.end(), monoTemsForBadBlocks.begin(),
                       monoTemsForBadBlocks.end());

    removeOutliers(newMonoTems);
    // std::cout << "removeOutliers pass" << "\n";
    monos.resize(newMonoTems.size());
    newDemcompBlockRes = {};
    for (size_t i = 0; i < monos.size(); i++) {
      SEQ mono(newMonoTems[i], to_string(i));
      monos[i] = mono;
    }
    demcomposer2(longSequence, monos, newDemcompBlockRes, 1);
    // std::cout << "demcomposer2 pass" << "\n";
    updateIdentityByEditDist(longSequence, newMonoTems, newDemcompBlockRes);
    // std::cout << "updateIdentityByEditDist pass" << "\n";
    // writeDecomposeRes(newDemcompBlockRes, demResfilename);
    // std::cout << "writeDecomposeRes pass" << "\n";
    // writeInferedMonos(newMonoTems, inferedMonosfilename);
    // std::cout << "writeInferedMonos pass" << "\n";
    return;
  }

  // use those refined templates to recompose try to improve some accuracy
  monos = {};
  newDemcompBlockRes = {};
  removeOutliers(newMonoTems);
  monos.resize(newMonoTems.size());
  for (size_t i = 0; i < monos.size(); i++) {
    SEQ mono(newMonoTems[i], to_string(i));
    monos[i] = mono;
  }
  demcomposer2(longSequence, monos, newDemcompBlockRes, 1);
  updateIdentityByEditDist(longSequence, newMonoTems, newDemcompBlockRes);
  // writeDecomposeRes(newDemcompBlockRes, demResfilename);
  // writeInferedMonos(newMonoTems, inferedMonosfilename);
}

void findGoodDemcomRegion(const vector<MonomerAlignment> &newDemcompBlockRes,
                         vector<region> &goodRegions, const int &minRegionLength,
                         const float &thresCutoff) {

  size_t firstPos = 0;
  size_t endPos = 0;
  size_t startID = 0;
  int tmp = 0;
  
  for (size_t i = 0; i < newDemcompBlockRes.size(); i++) {
    if (i < tmp)
      continue;
    
    float simScore = newDemcompBlockRes[i].identity;

    if (simScore > thresCutoff) {
      firstPos = newDemcompBlockRes[i].start_pos;
      startID = i;
      endPos = newDemcompBlockRes[i].end_pos; // Initialize endPos here

      for (size_t j = i + 1; j < newDemcompBlockRes.size(); j++) {
        simScore = newDemcompBlockRes[j].identity;
        if (simScore < thresCutoff || j == newDemcompBlockRes.size() - 1) {
          if (simScore >= thresCutoff) {
            endPos = newDemcompBlockRes[j].end_pos;
            tmp = j + 1;
          } else {
            tmp = j;
          }

          int goodLen = endPos - firstPos;
          if (goodLen > minRegionLength) {
            if (!goodRegions.empty() && firstPos - goodRegions.back().end < 2000) {
              goodRegions.back().end = endPos;
              goodRegions.back().blockEndID = tmp - 1;
            } else {
              region goodRegion(firstPos, endPos, startID, tmp - 1);
              goodRegions.push_back(goodRegion);
            }
          }
          break;
        } else {
          endPos = newDemcompBlockRes[j].end_pos; // Continuously update endPos
        }
      }
    }
  }
}

