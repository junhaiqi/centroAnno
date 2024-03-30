#include "hor.h"

void horTest(const std::string &filename, const std::string &monoFaFile,
             const std::string &seqFaFile, std::vector<std::string> &HORs) {

  std::vector<std::string> monoNames;
  std::vector<std::string> monoStrs = readFastaSequences(monoFaFile);
  std::vector<std::string> seqStrs = readFastaSequences(seqFaFile);
  // cout << seqStrs[0] << "\n";
  SEQ thisTestSEQ(seqStrs[0], "test");
  readDecomposeFileOfMonoNames(filename, monoNames);

  std::vector<std::string> testMonos(monoNames.begin(), monoNames.begin() + 20);
  // inferHORs(monoNames, thisTestSEQ, monoStrs, HORs);
  // cout << filename << "\n" << endl;
}

void inferHORs(
    const vector<MonomerAlignment> &monomerDecompBlockRes,
    const std::vector<std::string> &monos, const std::string &sequence,
    std::vector<HORINFO> &horDecompBlockRes,
    std::vector<std::string> &cleanedHORs,    // format is name, like 2_1_3...
    std::vector<std::string> &cleanedHORSeqs, // format is base sequence
    const int maxMonosOfHOR) {

  std::vector<std::string> monoNames(monomerDecompBlockRes.size());
  std::vector<std::string> HORs;
  for (size_t i = 0; i < monomerDecompBlockRes.size(); ++i) {
    monoNames[i] = monomerDecompBlockRes[i].monomer_name;
  }

  // std::cout << "1-loop pass" << std::endl;

  const int highBound = monoNames.size();
  const int ourMax =
      maxMonosOfHOR < highBound - 1 ? maxMonosOfHOR : highBound - 1;

  std::vector<std::vector<bool>> diagonals(ourMax);

  for (size_t i = 0; i < diagonals.size(); ++i) {
    std::vector<bool> thisdiagonal(highBound - i - 1, 0);
    for (size_t j = 0; j < thisdiagonal.size(); ++j) {
      if (monoNames[i + j + 1] == monoNames[j])
        thisdiagonal[j] = 1;
    }
    diagonals[i] = thisdiagonal;
  }

  // std::cout << "2-loop pass" << std::endl;

  std::vector<int> oneVec = countOneVec(diagonals);

  // std::cout << "countOneVec pass" << std::endl;

  int maxPosition = static_cast<int>(max_element(oneVec.begin(), oneVec.end()) -
                                     oneVec.begin());

  std::vector<int> factors;
  findFactor(maxPosition + 1, factors);

  // std::cout << "findFactor pass" << std::endl;

  int maxContinousOneNum1 = 0, maxEndPos1 = 0;
  countMaxContinousOne(diagonals[maxPosition], maxContinousOneNum1, maxEndPos1);
  int thisMaxNum = maxContinousOneNum1;
  int majorHORLen = maxPosition + 1;
  for (const int &factor : factors) {
    if (static_cast<float>(oneVec[factor - 1]) / oneVec[maxPosition] > 0.9) {
      int maxContinousOneNum2 = 0, maxEndPos2 = 0;
      countMaxContinousOne(diagonals[factor - 1], maxContinousOneNum2,
                           maxEndPos2);
      if (maxContinousOneNum2 > thisMaxNum) {
        majorHORLen = factor;
        thisMaxNum = maxContinousOneNum2;
      }
    }
  }

  std::cout << "Most frequent HOR length (monomer sequence) is: " << majorHORLen
            << "\n";

  std::vector<int> horRecorder(highBound, -1);
  findContiousAndDecomposeHORs(diagonals[majorHORLen - 1], horRecorder,
                               majorHORLen - 1, 0, monoNames, HORs);

  // std::cout << "findContiousAndDecomposeHORs pass" << std::endl;

  std::vector<std::string> horBaseSeqs;
  getHORBaseSeqs(HORs, monos, horBaseSeqs);

  // std::cout << "getHORBaseSeqs pass" << std::endl;

  std::vector<int> clusteringRes;
  horSmallClustering(horBaseSeqs, 0, clusteringRes);

  // std::cout << "horSmallClustering pass" << std::endl;

  if (HORs.empty()) {

    std::cout << "The HOR feature of sequence "
              << monomerDecompBlockRes[0].read_name << " was not detected.\n";
    return;
  }

  int maxIdx = *max_element(clusteringRes.begin(), clusteringRes.end());

  for (size_t i = 0; i < maxIdx + 1; ++i) {
    for (size_t j = 0; j < clusteringRes.size(); ++j) {
      if (clusteringRes[j] == i) {
        cleanedHORs.push_back(HORs[j]);
        cleanedHORSeqs.push_back(horBaseSeqs[j]);
        break;
      }
    }
  }

  // std::cout << "loop pass" << std::endl;

  for (size_t i = 0; i < horRecorder.size(); ++i) {
    if (horRecorder[i] != -1) {
      int idx =
          findFirstOccurrence(clusteringRes, clusteringRes[horRecorder[i]]);
      horRecorder[i] = clusteringRes[idx];
    }
  }

  // std::cout << "findFirstOccurrence pass" << std::endl;

  slidingWindowSloveIsolatedHORs(cleanedHORSeqs, monoNames, monos, 0,
                                 majorHORLen, horRecorder);

  // std::cout << "slidingWindowSloveIsolatedHORs pass" << std::endl;

  // check completed. 2024.1.28
  if (std::count(horRecorder.begin(), horRecorder.end(), -1) != 0) {
    for (size_t i = 1; i < diagonals.size(); ++i) {
      if (i != (majorHORLen - 1)) {
        const int currentHORNum = cleanedHORs.size();
        assert(currentHORNum > 0);
        std::vector<std::string> tempHORs;
        int t = i;
        findContiousAndDecomposeHORs(diagonals[i], horRecorder, t,
                                     currentHORNum, monoNames, tempHORs);

        if (tempHORs.empty())
          continue;

        std::vector<std::string> tempHORBaseSeqs;
        getHORBaseSeqs(tempHORs, monos, tempHORBaseSeqs);

        std::vector<int> tempClusteringRes;
        horSmallClustering(tempHORBaseSeqs, currentHORNum, tempClusteringRes);

        // std::cout << "clusters:" << tempClusteringRes.size() << std::endl;

        assert(*min_element(tempClusteringRes.begin(),
                            tempClusteringRes.end()) >= currentHORNum);

        std::vector<std::string>
            tempCleanedHORs; // format is name, like 2_1_3...
        std::vector<std::string> tempCleanedHORSeqs; // format is base sequence

        int maxIdx =
            *max_element(tempClusteringRes.begin(), tempClusteringRes.end());

        int minIdx =
            *min_element(tempClusteringRes.begin(), tempClusteringRes.end());

        for (size_t x = minIdx; x < maxIdx + 1; ++x) {
          for (size_t j = 0; j < tempClusteringRes.size(); ++j) {
            if (tempClusteringRes[j] == x) {
              tempCleanedHORs.push_back(tempHORs[j]);
              tempCleanedHORSeqs.push_back(tempHORBaseSeqs[j]);
              break;
            }
          }
        }

        for (size_t x = 0; x < horRecorder.size(); ++x) {
          if (horRecorder[x] >= currentHORNum) {
            int idx = findFirstOccurrence(
                tempClusteringRes,
                tempClusteringRes[horRecorder[x] - currentHORNum]);
            horRecorder[x] = tempClusteringRes[idx];
          }
        }

        cleanedHORs.insert(cleanedHORs.end(), tempCleanedHORs.begin(),
                           tempCleanedHORs.end());

        cleanedHORSeqs.insert(cleanedHORSeqs.end(), tempCleanedHORSeqs.begin(),
                              tempCleanedHORSeqs.end());

        slidingWindowSloveIsolatedHORs(cleanedHORSeqs, monoNames, monos, 0,
                                       t + 1, horRecorder);
      }
    }
  }

  // std::cout << "big loop pass" << std::endl;

  size_t end = 0;
  assert(horRecorder.size() == monomerDecompBlockRes.size());
  for (size_t i = 0; i < horRecorder.size(); ++i) {
    if (i < end)
      continue;

    if (horRecorder[i] == -1) {
      for (size_t j = i; j < horRecorder.size(); ++j) {
        if (horRecorder[j] != -1) {
          end = j;
          break;
        }
      }
      int tempEnd = end; // come to end;
      if (end > i)
        tempEnd--;
      else {
        end = monoNames.size();
        tempEnd = end - 1;
      }

      assert(end <= monoNames.size());
      std::string horName = getAHOR(monoNames, i, end);
      std::string thisHORBaseSeq = getAHORBaseSeq(horName, monos);
      std::string subSeq =
          sequence.substr(monomerDecompBlockRes[i].start_pos,
                          monomerDecompBlockRes[tempEnd].end_pos -
                              monomerDecompBlockRes[i].start_pos + 1);
      assert(subSeq.length() > 0);
      int edist;
      size_t l = std::max(subSeq.length(), thisHORBaseSeq.length());
      calEditDist(subSeq, thisHORBaseSeq, edist);
      float identity = 1 - (static_cast<float>(edist) / l);
      HORINFO m(horName, monomerDecompBlockRes[i].read_name,
                monomerDecompBlockRes[i].start_pos,
                monomerDecompBlockRes[tempEnd].end_pos, identity, "low");
      horDecompBlockRes.push_back(m);
    } else {
      std::vector<std::string> thisHORVec =
          parseMonoStr(cleanedHORs[horRecorder[i]]);

      int thisHORLen = thisHORVec.size();
      end =
          i + thisHORLen > monoNames.size() ? monoNames.size() : i + thisHORLen;
      std::string testThisHOR = getAHOR(monoNames, i, end);
      std::string horName = cleanedHORs[horRecorder[i]];
      std::string subSeq =
          sequence.substr(monomerDecompBlockRes[i].start_pos,
                          monomerDecompBlockRes[end - 1].end_pos -
                              monomerDecompBlockRes[i].start_pos + 1);

      int edist;
      size_t l =
          std::max(subSeq.length(), cleanedHORSeqs[horRecorder[i]].length());
      calEditDist(subSeq, cleanedHORSeqs[horRecorder[i]], edist);
      float identity = 1 - (static_cast<float>(edist) / l);
      HORINFO m(horName, monomerDecompBlockRes[i].read_name,
                monomerDecompBlockRes[i].start_pos,
                monomerDecompBlockRes[end - 1].end_pos, identity, "high");
      horDecompBlockRes.push_back(m);
    }
  }
}

void findContiousAndDecomposeHORs(const std::vector<bool> &diagonal,
                                  std::vector<int> &horRecorder,
                                  const int &regionLenBound,
                                  const int &idxOffset,
                                  const std::vector<std::string> &monoNames,
                                  std::vector<std::string> &HORs) {

  int thisContinousOneNum = 0;
  const int offset = regionLenBound + 1; // regionLenBound = idx[dia]
  int end = 0;
  for (size_t i = 0; i < diagonal.size(); ++i) {
    if (i < end || horRecorder[i] != -1) {
      thisContinousOneNum = 0;
      continue;
    }

    if (diagonal[i] == 1) { // horRecorder[i+offset] corresponding
                            // to a positon in monomer sequenece
      thisContinousOneNum++;
    } else {
      if (thisContinousOneNum > regionLenBound) { // a HOR region is found
        std::vector<std::string> monoRegion(
            i + monoNames.begin() - thisContinousOneNum, i + monoNames.begin());
        const std::string thisHOR = getAHOR(monoRegion, 0, regionLenBound + 1);
        HORs.push_back(thisHOR);
        const int horIdx = findFirstOccurrence(HORs, thisHOR);
        std::vector<int> maxAnnoteRegion(horRecorder.begin() + i -
                                             thisContinousOneNum,
                                         horRecorder.begin() + i + offset);
        int mightAnnotedLen = countElementOccurrences(maxAnnoteRegion, -1) -
                              (countElementOccurrences(maxAnnoteRegion, -1) %
                               (regionLenBound + 1));
        for (size_t j = i - thisContinousOneNum;
             j < i - thisContinousOneNum + mightAnnotedLen; ++j) {
          horRecorder[j] = horIdx + idxOffset;
        }
      }
      thisContinousOneNum = 0;
      end = i;
    }
  }

  if (thisContinousOneNum == diagonal.size()) { // prefect results.
    const std::string thisHOR = getAHOR(monoNames, 0, offset);
    HORs.push_back(thisHOR);
    for (int i = 0; i < horRecorder.size(); ++i) {
      horRecorder[i] = idxOffset;
    }
    // return;
  }

  else if (thisContinousOneNum > regionLenBound) {
    size_t i = diagonal.size() - 1;
    std::vector<std::string> monoRegion(
        i + monoNames.begin() - thisContinousOneNum, i + monoNames.begin());
    const std::string thisHOR = getAHOR(monoRegion, 0, regionLenBound + 1);
    HORs.push_back(thisHOR);
    const int horIdx = findFirstOccurrence(HORs, thisHOR);
    std::vector<int> maxAnnoteRegion(horRecorder.begin() + i -
                                         thisContinousOneNum,
                                     horRecorder.begin() + i + offset);
    int mightAnnotedLen =
        countElementOccurrences(maxAnnoteRegion, -1) -
        (countElementOccurrences(maxAnnoteRegion, -1) % (regionLenBound + 1));

    for (size_t j = i - thisContinousOneNum;
         j < i - thisContinousOneNum + mightAnnotedLen; ++j) {
      horRecorder[j] = horIdx + idxOffset;
    }
  }
}

// 2024.1.28, unused.
void refineDecomposeHORs(const std::vector<std::string> &monoNames,
                         const std::vector<std::string> &HORs,
                         const std::vector<std::string> &monos,
                         std::vector<MonomerAlignment> &horDemcompBlockRes,
                         const SEQ &sequence) {

  std::set<std::string> thisSet(HORs.begin(), HORs.end());
  std::vector<std::string> filteredHORs;
  filteredHORs.assign(thisSet.begin(), thisSet.end());
  std::vector<SEQ> HORSEQs;
  for (const std::string &hor : filteredHORs) {
    std::vector<std::string> monosV = parseMonoStr(hor);
    std::string horStr = "";
    for (const std::string &m : monosV) {
      std::string thisMonomer;
      if (m.find("'") != std::string::npos) {
        const int idx = std::stoi(m.substr(0, m.length() - 1));
        thisMonomer = reverse_complement(monos[idx]);
      } else {
        const int idx = std::stoi(m);
        thisMonomer = monos[idx];
      }
      horStr += thisMonomer;
    }
    SEQ thisHOR(horStr, hor);
    HORSEQs.push_back(thisHOR);
  }
  demcomposer2(sequence, HORSEQs, horDemcompBlockRes, 1);
}

// 2024.1.28, unused.
void smallTryDecomposeUseMonoNameSeq(
    const std::vector<std::string> &monoNames,
    const std::vector<std::string> &HORs,
    std::vector<MonomerAlignment> &horDemcompBlockRes,
    const std::string &seqName) {

  std::set<std::string> thisSet(HORs.begin(), HORs.end());
  std::vector<std::string> filteredHORs;
  filteredHORs.assign(thisSet.begin(), thisSet.end());

  std::string monoStr = getAHOR(monoNames, 0, monoNames.size());
  SEQ decomposedSeq(monoStr, seqName);

  std::vector<SEQ> horVec;
  for (const std::string &hor : filteredHORs) {
    SEQ horSEQ(hor, hor);
    horVec.push_back(horSEQ);
  }
  cout << horVec.size() << "," << decomposedSeq.seq << "\n";
  demcomposeForHOR(decomposedSeq, horVec, horDemcompBlockRes, 1);

  // update identity
  for (MonomerAlignment &block : horDemcompBlockRes) {
    const std::string horName = block.monomer_name;
    const std::string thisMonoSeg(monoStr.begin() + block.start_pos,
                                  monoStr.begin() + block.end_pos + 1);
    cout << horName << "," << thisMonoSeg << "\n";
  }
}

void horSmallClustering(const std::vector<std::string> &horBaseSeqs,
                        const int &idxOffset, std::vector<int> &clusteringRes) {

  if (horBaseSeqs.empty())
    return;
  else if (horBaseSeqs.size() == 1)
    clusteringRes = {idxOffset};

  if (clusteringRes.empty()) {
    std::vector<int> temp(horBaseSeqs.size(), -1);
    clusteringRes = temp;
  }

  size_t clusterIdx = idxOffset;
  for (size_t i = 0; i < horBaseSeqs.size() - 1; ++i) {
    if (clusteringRes[i] == -1) {
      clusteringRes[i] = clusterIdx;
      for (size_t j = i + 1; j < horBaseSeqs.size(); ++j) {
        int edist;
        size_t l = std::max(horBaseSeqs[i].size(), horBaseSeqs[j].size());
        // std::cout << "calEdit-f\n";
        calEditDist(horBaseSeqs[j], horBaseSeqs[i], edist);
        // std::cout << "calEdit-e\n";
        float identity = 1 - (static_cast<float>(edist) / l);
        if (identity > 0.95)
          clusteringRes[j] = clusterIdx;
      }
      ++clusterIdx;
    }
  }

  if (clusteringRes.back() == -1)
    clusteringRes.back() = clusterIdx;
}

void slidingWindowSloveIsolatedHORs(
    const std::vector<std::string> &cleanedHORSeqs,
    const std::vector<std::string> &monoNames,
    const std::vector<std::string> &monos, const int &idxOffset,
    const int &windowSize, std::vector<int> &horRecorder) {

  assert(horRecorder.size() == monoNames.size());
  for (size_t i = 0; i < horRecorder.size() - windowSize + 1; ++i) {
    std::vector<int> seg(horRecorder.begin() + i,
                         horRecorder.begin() + i + windowSize);
    if (countElementOccurrences(seg, -1) ==
        windowSize) { // a unannoted segment with size equal to $windowSize
      std::vector<std::string> monoSegs(monoNames.begin() + i,
                                        monoNames.begin() + i + windowSize);
      std::string hor = getAHOR(monoSegs, 0, monoSegs.size());
      std::string thisSeq = getAHORBaseSeq(hor, monos);
      for (size_t j = 0; j < cleanedHORSeqs.size(); ++j) {
        int edist;
        size_t l = std::max(thisSeq.size(), cleanedHORSeqs[j].size());
        calEditDist(thisSeq, cleanedHORSeqs[j], edist);
        float identity = 1 - (static_cast<float>(edist) / l);
        if (identity > 0.95) {
          for (size_t k = 0; k < windowSize; ++k) {
            horRecorder[i + k] = (idxOffset + j);
          }
          break;
        }
      }
    }
  }
}