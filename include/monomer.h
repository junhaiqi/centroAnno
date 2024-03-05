#ifndef MONOMER_H
#define MONOMER_H
#include <assert.h>
#include <map>
#include <set>
#include <cstring>
#include <sstream>
#include <iterator>
#include <numeric>
#include "sequenceUtils.h"
using namespace std;
/**
 * @brief Inference the repeat template from the repetitive sequence
 * @param sequence Repetitive sequence
 * @param k K-mer size
 * @param kIdenCf The cutoff of identity to indicates different types of k-mers
 * @param kCountCf The cutoff of count to limit sampled k-mers
 * @param monoSeqs All temporary monomers with differenct length
 * @return bestMonoTemplate Most likely template
 */

void monoTempInference(SEQ &sequence, const int &k, const float &kIdenCf,
                       const int &kCountCf, std::vector<std::string> &monoSeqs,
                       std::vector<SEQ> &bestMonoTemplateList,
                       const bool &hpc, const float &repPorCutoff);

std::vector<int> sortReturnIdxList(const std::vector<int> &v);

template <typename T> T findMode(const std::vector<T> &nums, int &count) {
  robin_hood::unordered_map<T, int> counts;

  // Count occurrences of each element
  for (const T &num : nums) {
    counts[num]++;
  }

  // Find the element with the maximum count
  T mode = T();
  int maxCount = 0;
  for (const auto &entry : counts) {
    if (entry.second > maxCount) {
      mode = entry.first;
      maxCount = entry.second;
    }
  }
  count = maxCount;
  return mode;
}

void evalMonomer(std::vector<std::string> &monoList,
                 float &maxMeanScoreOfMonomer, int &maxAvgIndex);

template <typename T> float calculateMean(const std::vector<T> &values) {
  if (values.empty()) {
    std::cerr << "Error: Empty vector." << std::endl;
    return 0.0;
  }
  T sum = std::accumulate(values.begin(), values.end(), 0.0);
  return sum / values.size();
}

void calLossList(const std::vector<int> &v, const int &num,
                 std::vector<float> &lossList);

void inferRepeatLen(const std::vector<int> &monoLenList, int &inferedIdx);
void printStrLen(std::vector<std::string> seqList);

struct ReadId {
  std::string name;
  int id = -1;

  ReadId(string name_) : name(name_) {}
  ReadId(string name_, int id_) : name(name_), id(id_) {}
};

struct Seq {
  ReadId read_id;
  std::string seq;

  Seq(string name_, string seq_) : read_id(ReadId(name_)), seq(seq_) {
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
  }

  Seq(string name_, string seq_, int id_)
      : read_id(ReadId(name_, id_)), seq(seq_) {}

  size_t size() { return seq.size(); }
};

struct MonomerAlignment {
  std::string monomer_name;
  std::string read_name;
  int start_pos;
  int end_pos;
  float identity = 0;
  bool best;

  MonomerAlignment() {}

  MonomerAlignment(string monomer_name_, string read_name_, int start_pos_,
                   int end_pos_, float identity_, bool best_)
      : monomer_name(monomer_name_), read_name(read_name_),
        start_pos(start_pos_), end_pos(end_pos_), identity(identity_),
        best(best_) {}
};

inline bool sortby1(const pair<int, vector<MonomerAlignment>> &a,
             const pair<int, vector<MonomerAlignment>> &b);

class MonomersAligner {

public:
  MonomersAligner(vector<Seq> &monomers, int ins = -1, int del = -1,
                  int mismatch = -1, int match = 1)
      : monomers_(monomers), ins_(ins), del_(del), mismatch_(mismatch),
        match_(match) {}

  void AlignReadsSet(vector<Seq> &reads, vector<MonomerAlignment> &batch, int threads, int part_size, int ed_thr,
                     int overlap);

  vector<MonomerAlignment> AlignPartClassicDP(Seq &read,
                                              std::vector<Seq> &monomers);
  void SaveBatch(vector<MonomerAlignment> &batch);

  vector<MonomerAlignment> PostProcessing(vector<MonomerAlignment> &batch);
  vector<MonomerAlignment> ourPostProcessing(vector<MonomerAlignment> &batch);
  ~MonomersAligner() {}

private:
  double MonomerEditDistance(Seq &monomer, Seq &read);

  std::vector<Seq> FilterMonomersForRead(Seq &read, int ed_thr);

  vector<Seq> monomers_;
  const int SAVE_STEP = 1;
  int ins_;
  int del_;
  int mismatch_;
  int match_;
};

string reverse_complement(const string &s);
void add_reverse_complement(vector<Seq> &monomers);

typedef robin_hood::unordered_map<std::string, int> name2Idx;
void demcomposer(const SEQ &sequence, std::vector<SEQ> &bestMonoTemplateList, std::vector<MonomerAlignment> &demcompBlockRes);

void demcomposer2(const SEQ &sequence, std::vector<SEQ> &bestMonoTemplateList,
                  std::vector<MonomerAlignment> &demcompBlockRes, const int &threadNum);
#endif