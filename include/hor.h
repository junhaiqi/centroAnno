#ifndef HOR_H
#define HOR_H
#include "monoRefine.h"
#include "moutput.h"
#include <set>

void horTest(const std::string &filename, const std::string &monoFaFile,
             const std::string &seqFaFile, std::vector<std::string> &HORs);

struct miniPair {
  int start;
  int end;
  miniPair() = default;
  miniPair(int s, int e) : start(s), end(e){};
};

struct HORINFO {
  std::string monomer_name;
  std::string read_name;
  int start_pos;
  int end_pos;
  float identity = 0;
  std::string confidenceLevel;

  HORINFO() = default;

  HORINFO(string monomer_name_, string read_name_, int start_pos_, int end_pos_,
          float identity_, string high_)
      : monomer_name(monomer_name_), read_name(read_name_),
        start_pos(start_pos_), end_pos(end_pos_), identity(identity_),
        confidenceLevel(high_) {}
};

inline void writeHORDecomposeRes(const std::vector<HORINFO> &horDecomoseRes,
                                 const std::string &filename) {

  std::ofstream outFile(filename, std::ios::app);
  if (!outFile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  outFile << "Sequence name,Mononer name,Start position,End position,Estimated "
             "Identity,Length,Confidence"
          << "\n";
  for (const auto &a : horDecomoseRes) {
    const size_t len = a.end_pos - a.start_pos + 1;
    std::string s =
        a.read_name + "," + a.monomer_name + "," + std::to_string(a.start_pos) +
        "," + std::to_string(a.end_pos) + "," + std::to_string(a.identity) +
        "," + std::to_string(len) + "," + a.confidenceLevel;
    outFile << s << "\n";
  }
}

void inferHORs(
    const vector<MonomerAlignment> &monomerDecompBlockRes,
    const std::vector<std::string> &monos, const std::string &sequence,
    std::vector<HORINFO> &horDecompBlockRes,
    std::vector<std::string> &cleanedHORs,    // format is name, like 2_1_3...
    std::vector<std::string> &cleanedHORSeqs, // format is base sequence
    const int maxMonosOfHOR = 50);

inline void printfDiag(const std::vector<vector<bool>> &diagonals) {
  for (size_t i = 0; i < diagonals.size(); ++i) {
    for (size_t j = 0; j < diagonals[i].size(); ++j) {
      std::cout << diagonals[i][j] << " ";
    }
    std::cout << "\n";
  }
}

inline void countOne(const std::vector<bool> &diagonal, int &oneNum) {
  int sum = 0;
  for (const bool &d : diagonal) {
    if (d == 1) {
      sum++;
    }
  }
  oneNum = sum;
}

inline void findFactor(const int &num, std::vector<int> &factors) {
  const int factorBound = num / 2 + 1;
  for (size_t i = 2; i < factorBound; i++) {
    if (num % i == 0) {
      factors.push_back(i);
    }
  }
}

inline std::vector<int> countOneVec(std::vector<std::vector<bool>> diagonals) {
  std::vector<int> res(diagonals.size());
  size_t idx = 0;
  for (const std::vector<bool> &dia : diagonals) {
    int oneNum = 0;
    countOne(dia, oneNum);
    res[idx] = oneNum;
    idx++;
  }
  return res;
}

inline void countLimMaxContinousOne(const std::vector<bool> &diagonal,
                                    int &limitedContinousOneNum,
                                    const int &bound) {

  int thisContinousOneNum = 0;
  limitedContinousOneNum = 0;
  for (size_t i = 0; i < diagonal.size(); ++i) {
    if (diagonal[i] == 1) {
      thisContinousOneNum++;
    } else {
      if (thisContinousOneNum > bound) {
        limitedContinousOneNum += thisContinousOneNum;
      }
      thisContinousOneNum = 0;
    }
  }
  if (thisContinousOneNum > bound)
    limitedContinousOneNum += thisContinousOneNum;
}

inline void countMaxContinousOne(const std::vector<bool> &diagonal,
                                 int &maxContinousOneNum, int &maxEndPos) {

  int lastMaxContinousOneNum = 0;
  int lastItem = diagonal[0];
  int thisMaxContinousOneNum = lastItem > 0 ? 1 : 0;
  for (size_t i = 1; i < diagonal.size(); ++i) {
    if (diagonal[i] == 1 && lastItem == 1)
      thisMaxContinousOneNum++;
    else {
      if (thisMaxContinousOneNum > lastMaxContinousOneNum) {
        lastMaxContinousOneNum = thisMaxContinousOneNum;
        maxEndPos = i - 1;
      }
      lastItem = diagonal[i];
      thisMaxContinousOneNum = 0;
      if (lastItem == 1)
        thisMaxContinousOneNum++;
    }
  }

  if (lastMaxContinousOneNum > thisMaxContinousOneNum) {
    maxContinousOneNum = lastMaxContinousOneNum;
  } else {
    maxContinousOneNum = thisMaxContinousOneNum;
    maxEndPos = diagonal.size() - 1;
  }
}

/**
 * @brief A important function for decomposing sequence to HOR blocks
 * @param diagonal Diagonal of monomer adjacency matrix, its index indicates the
 * length of HOR
 * @param horRecorder Recorder of hor, record the hor corresponding to each
 * position in the monomer sequence
 * @param regionLenBound Consecutive 1s in the matrix represent a region, this
 * is a minimum bound
 * @param monoNames A name vector which indicates a monomer sequence
 * @param HORs Record all HORs, HOR in a region which length larger than
 * "regionLenBound"
 */
void findContiousAndDecomposeHORs(const std::vector<bool> &diagonal,
                                  std::vector<int> &horRecorder,
                                  const int &regionLenBound,
                                  const int &idxOffset,
                                  const std::vector<std::string> &monoNames,
                                  std::vector<std::string> &HORs);

void refineDecomposeHORs(const std::vector<std::string> &monoNames,
                         const std::vector<std::string> &HORs,
                         const std::vector<std::string> &monos,
                         std::vector<MonomerAlignment> &horDemcompBlockRes,
                         const SEQ &sequence);

inline std::string getAHOR(const std::vector<std::string> &monoRegion,
                           const int &start, const int &end) {

  const std::vector<std::string> HORChars(monoRegion.begin() + start,
                                          monoRegion.begin() + end);

  std::string thisHOR = "";
  for (size_t i = 0; i < HORChars.size(); ++i) {
    if (i != HORChars.size() - 1)
      thisHOR = thisHOR + HORChars[i] + "_";
    else
      thisHOR = thisHOR + HORChars[i];
  }
  return thisHOR;
}

template <typename T>
inline int findFirstOccurrence(const std::vector<T> &vec, const T &target) {
  auto it = std::find(vec.begin(), vec.end(), target);
  if (it != vec.end()) {
    return std::distance(vec.begin(), it);
  } else {
    return -1;
  }
}

inline std::vector<std::string> parseMonoStr(const std::string &inputString) {
  std::vector<std::string> strs;
  std::string tempString = "";

  for (char c : inputString) {
    if (c != '_') {
      tempString += c;
    } else {
      if (!tempString.empty()) {
        strs.push_back(tempString);
        tempString.clear();
      }
    }
  }
  if (!tempString.empty()) {
    strs.push_back(tempString);
  }

  return strs;
}

inline void findMaxOneIdx(const std::vector<std::vector<bool>> &diagonals,
                          int &maxIdx, int &maxContinousOneIdx,
                          int &maxLContinousOneIdx, int &maxEndPos,
                          int &maxCnum) {
  int idx = 0, maxNum = -1, mcmNum = -1, mclmNum = -1;
  for (const std::vector<bool> &dia : diagonals) {
    int thisOneNum = 0;
    int thisPos;
    int thisContinousOneNum = 0;
    int thisLimitedContinousOneNum = 0;
    countOne(dia, thisOneNum);
    countMaxContinousOne(dia, thisContinousOneNum, thisPos);
    countLimMaxContinousOne(dia, thisLimitedContinousOneNum, idx + 2);
    if (thisOneNum > maxNum) {
      maxIdx = idx;
      maxNum = thisOneNum;
    }

    if (thisContinousOneNum > mcmNum) {
      maxContinousOneIdx = idx;
      mcmNum = thisContinousOneNum;
      maxEndPos = thisPos;
    }

    if (thisLimitedContinousOneNum > mclmNum) {
      maxLContinousOneIdx = idx;
      mclmNum = thisLimitedContinousOneNum;
    }

    idx++;
  }
  maxCnum = mcmNum;
}

void smallTryDecomposeUseMonoNameSeq(
    const std::vector<std::string> &monoNames,
    const std::vector<std::string> &HORs,
    std::vector<MonomerAlignment> &horDemcompBlockRes,
    const std::string &seqName);

inline void demcomposeForHOR(const SEQ &sequence,
                             std::vector<SEQ> &bestMonoTemplateList,
                             std::vector<MonomerAlignment> &demcompBlockRes,
                             const int &threadNum) {
  int ins = -1, del = -1, mismatch = -1, match = 1;
  std::vector<Seq> reads = {Seq(sequence.seqName, sequence.seq)};
  vector<Seq> monomers;
  for (const SEQ &mono : bestMonoTemplateList) {
    monomers.push_back(Seq(mono.seqName, mono.seq));
  }
  MonomersAligner monomers_aligner(monomers, ins, del, mismatch, match);
  monomers_aligner.AlignReadsSet(reads, demcompBlockRes, threadNum, 5000, -1,
                                 500);
}

inline std::string getAHORBaseSeq(const std::string &hor,
                                  const std::vector<std::string> &monos) {
  std::vector<std::string> monosV = parseMonoStr(hor);
  std::string horBaseSeq = "";
  for (const std::string &m : monosV) {
    std::string thisMonomer;
    if (m.find("'") != std::string::npos) {
      const int idx = std::stoi(m.substr(0, m.length() - 1));
      thisMonomer = reverse_complement(monos[idx]);
    } else {
      const int idx = std::stoi(m);
      thisMonomer = monos[idx];
    }
    horBaseSeq += thisMonomer;
  }
  return horBaseSeq;
}

inline void getHORBaseSeqs(const std::vector<std::string> &hors,
                           const std::vector<std::string> &monos,
                           std::vector<std::string> &horBaseSeqs) {

  for (const std::string &hor : hors) {
    std::string thisHORBaseSeq = getAHORBaseSeq(hor, monos);
    horBaseSeqs.push_back(thisHORBaseSeq);
  }
}

template <typename T>
inline int countElementOccurrences(const std::vector<T> &vec,
                                   const T &targetElement) {
  int occurrences = std::count(vec.begin(), vec.end(), targetElement);
  return occurrences;
}

void horSmallClustering(const std::vector<std::string> &horBaseSeqs,
                        const int &idxOffset, std::vector<int> &clusteringRes);

void slidingWindowSloveIsolatedHORs(
    const std::vector<std::string> &cleanedHORSeqs,
    const std::vector<std::string> &monoNames,
    const std::vector<std::string> &monos, const int &idxOffset,
    const int &windowSize, std::vector<int> &horRecorder);

#endif