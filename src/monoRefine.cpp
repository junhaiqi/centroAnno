#include "monoRefine.h"
#include <stdexcept>

template <typename T> static void print2DMatrix(vector<vector<T>> Matrix) {
  for (size_t i = 0; i < Matrix.size(); i++) {
    for (size_t j = 0; j < Matrix.size(); j++) {
      cout << Matrix[i][j] << " ";
    }
    cout << endl;
  }
}

// clustering by dbscan to infer possible monomer template.
void sampleClustering(const SEQ &sequence,
                      const vector<MonomerAlignment> &demcompBlockRes,
                      vector<int> &clusterLabels,
                      vector<string> &consensuses, // new templates
                      size_t threadNum, const float epsilon,
                      const int minSamples) {

  vector<string> blocks(demcompBlockRes.size(), ""); // extract all segments.
  size_t idx = 0;
  for (const MonomerAlignment &m : demcompBlockRes) {
    if (m.identity / static_cast<float>(m.end_pos - m.start_pos + 1) >
        0.5) { // 2024.1.18 Qi add this condition in order to avoid abnormal
               // blocks from participating in the clustering process.
      blocks[idx] = sequence.seq.substr(m.start_pos, m.end_pos - m.start_pos + 1);
      idx++;
    }
  }
  if (idx == 0) {
    printf("%s may not contain repeat unit! (can't inference monomer templates "
           "by clustering)\n",
           sequence.seqName);
    exit(-1);
  }
  blocks.resize(idx);
  vector<vector<float>> idenMatrix(demcompBlockRes.size(),
                                   vector<float>(demcompBlockRes.size(), 1));
#pragma omp parallel for num_threads(                                          \
    threadNum) // it can be replaced by our cudaEdit in future.
  for (size_t i = 0; i < demcompBlockRes.size() - 1; i++) {
    for (size_t j = i + 1; j < demcompBlockRes.size(); j++) {
      int edist = 0;
      size_t l = std::max(blocks[i].size(), blocks[j].size());
      calEditDist(blocks[i], blocks[j], edist);
      float identity = 1 - (static_cast<float>(edist) / l);
      idenMatrix[i][j] = identity;
      idenMatrix[j][i] = identity;
    }
  }

  // print2DMatrix(idenMatrix);
  DBSCAN dbscan(epsilon, minSamples, idenMatrix);
  dbscan.runDBSCAN();
  clusterLabels = dbscan.getLabels();
  clusters clusRes;
  idx = 0;
  for (const int &label : clusterLabels) {
    clusRes[label].push_back(blocks[idx]);
    idx++;
  }

  consensuses.resize(clusRes.size());
  idx = 0;
  for (const auto &key : clusRes) {
    std::string consensus;
    msaAndConsensus(key.second, consensus);
    if (!consensus.empty()) {
      consensuses[idx] = consensus;
      idx++;
    }
  }
  consensuses.resize(idx);

  removeOutliers(consensuses);

  std::cout << "Clustering is over, there are " << consensuses.size()
            << " templates."
            << " There are " << idx - consensuses.size()
            << " templates are filtered." << std::endl;
}

// decompose a sequence into blocks which have high indentity with our infered
// monomer templates, and find the failed blocks to redecompose.
void refineDemcompose(const SEQ &longSequence, vector<string> &newMonoTems,
                      vector<MonomerAlignment> &newDemcompBlockRes,
                      const int &threadNum, const bool &hpc, const int &k,
                      const float &fpsCutoff, const float &repPorCutoff,
                      const std::string &demResfilename,
                      const std::string &inferedMonosfilename) {

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
    writeDecomposeRes(newDemcompBlockRes, demResfilename);
    // std::cout << "writeDecomposeRes pass" << "\n";
    writeInferedMonos(newMonoTems, inferedMonosfilename);
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
  writeDecomposeRes(newDemcompBlockRes, demResfilename);
  writeInferedMonos(newMonoTems, inferedMonosfilename);
}

// find bad regions which have low identity with our monomer template.
void findBadDemcomRegion(const vector<MonomerAlignment> &newDemcompBlockRes,
                         vector<region> &badRegions, const int &minRegionLength,
                         const float &thresCutoff) {

  size_t firstPos = 0;
  size_t endPos = 0;
  int tmp = 0;
  for (size_t i = 0; i < newDemcompBlockRes.size(); i++) {
    if (i < tmp)
      continue;
    float simScore =
        newDemcompBlockRes[i].identity /
        (newDemcompBlockRes[i].end_pos - newDemcompBlockRes[i].start_pos);

    if (simScore < thresCutoff) {
      firstPos = newDemcompBlockRes[i].start_pos;
      for (size_t j = i + 1; j < newDemcompBlockRes.size(); j++) {
        simScore =
            newDemcompBlockRes[j].identity /
            (newDemcompBlockRes[j].end_pos - newDemcompBlockRes[j].start_pos);

        if (simScore > thresCutoff) {
          tmp = j;
          endPos = newDemcompBlockRes[j].end_pos;
          int badLen = endPos - firstPos;
          if (badLen > minRegionLength) {
            region badRegion(firstPos, endPos, i, j);
            badRegions.push_back(badRegion);
          }
          break;
        }
      }
    }
  }
}

// find the minimal bad region and remove it from vector
void findMinBadRegion(vector<region> &badRegions, region &minRegion,
                      int &minLen) {
  int miniRegionLen = badRegions[0].end - badRegions[0].start;
  size_t idx = 0;
  size_t miniIdx = 0;
  for (const region &r : badRegions) {
    const int thisLen = r.end - r.start;
    if (thisLen <= miniRegionLen) {
      minRegion = r;
      miniIdx = idx;
      miniRegionLen = thisLen;
    }
    idx++;
  }
  minLen = miniRegionLen;
  badRegions.erase(badRegions.begin() + miniIdx);
}

void updateIdentityByEditDist(const SEQ &longSequence,
                              vector<string> &newMonoTems,
                              vector<MonomerAlignment> &newDemcompBlockRes) {

  for (MonomerAlignment &m : newDemcompBlockRes) {
    std::string thisMonomer;
    if (m.monomer_name.find("'") != std::string::npos) {
      const int idx =
          std::stoi(m.monomer_name.substr(0, m.monomer_name.length() - 1));
      thisMonomer = reverse_complement(newMonoTems[idx]);
    } else {
      const int idx = std::stoi(m.monomer_name);
      thisMonomer = newMonoTems[idx];
    }
    if (m.start_pos < longSequence.seq.size() && m.start_pos + m.end_pos - m.start_pos + 1 < longSequence.seq.size()) {
      std::string seqBlock =
      longSequence.seq.substr(m.start_pos, m.end_pos - m.start_pos + 1);
      int edist = 0;
      size_t l = std::max(seqBlock.size(), thisMonomer.size());
      calEditDist(seqBlock, thisMonomer, edist);
      float identity = 1 - (static_cast<float>(edist) / l);
      m.identity = identity;
    }
    else  m.identity = 0;
  }
}

const float cutoffRetioOfBad = 0.05; // 0.05 is large enough.
void checkBadBlockNums(const SEQ &sequence,
                       const std::vector<MonomerAlignment> &newDemcompBlockRes,
                       bool &badSignal, std::vector<std::string> &badBlocks,
                       std::vector<std::string> &monoTemsForBadBlocks,
                       const size_t &threadNum, const float epsilon,
                       const int minSamples) {

  for (const MonomerAlignment &m : newDemcompBlockRes) {
    const float esIdentity =
        m.identity / static_cast<float>(m.end_pos - m.start_pos + 1);
    if (esIdentity < 0.90) { // 0.92
      if (m.start_pos < sequence.seq.size() && m.start_pos + m.end_pos - m.start_pos + 1 < sequence.seq.size()) {
        badBlocks.push_back(
        sequence.seq.substr(m.start_pos, m.end_pos - m.start_pos + 1));
      }
    } // it indicates a bad block might be caused by insufficient templates.
  }

  if (badBlocks.size() / static_cast<float>(newDemcompBlockRes.size()) >
      cutoffRetioOfBad) {
    badSignal = true;
    vector<vector<float>> idenMatrix(badBlocks.size(),
                                     vector<float>(badBlocks.size(), 1));
#pragma omp parallel for num_threads(                                          \
    threadNum) // it can be replaced by our cudaEdit in future.
        for (size_t i = 0; i < badBlocks.size() - 1; i++) {
          for (size_t j = i + 1; j < badBlocks.size(); j++) {
            int edist = 0;
            size_t l = std::max(badBlocks[i].size(), badBlocks[j].size());
            calEditDist(badBlocks[i], badBlocks[j], edist);
            float identity = 1 - (static_cast<float>(edist) / l);
            idenMatrix[i][j] = identity;
            idenMatrix[j][i] = identity;
          }
        }
        DBSCAN dbscan(epsilon, minSamples, idenMatrix);
        dbscan.runDBSCAN();
        const std::vector<int> clusterLabels = dbscan.getLabels();
        clusters clusRes;
        size_t idx = 0;
        for (const int &label : clusterLabels) {
          clusRes[label].push_back(badBlocks[idx]);
          idx++;
        }

        monoTemsForBadBlocks.resize(clusRes.size());
        idx = 0;
        for (const auto &key : clusRes) {
          std::string consensus;
          msaAndConsensus(key.second, consensus);
          if (!consensus.empty()) {
            monoTemsForBadBlocks[idx] = consensus;
            idx++;
          }
        }
        monoTemsForBadBlocks.resize(idx);
    } else {
        badSignal = false;
    }
}