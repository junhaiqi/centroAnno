
#ifndef SEQUtils_H
#define SEQUtils_H
#include <stdio.h>
#include <zlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "robin_hood.h"
#include "edlib.h"

typedef robin_hood::unordered_map<std::string, std::vector<int>> counter;
typedef robin_hood::unordered_map<std::string, int> counter_t;

class SEQ {
public:
  std::string seq;
  std::string seqName;
  bool comp;
  std::vector<int> compIdxList;
  std::string compSeq;
  // int seqID;
  // std::string seqName;
  SEQ() = default;  
  SEQ(const std::string &sequence, const std::string &name);

  void compressSeq();
  std::string recoverSeg(int startPos, int endPos);

  ~SEQ() {}
};

typedef class KMERCOUNTER {
public:
  std::string seq;
  counter h;
  int k = 11;
  float repeativeScore = 0;
  std::vector<std::string> sortedKmers;
  std::vector<std::string> selectedKmers;
  KMERCOUNTER() = default;
  KMERCOUNTER(const std::string &sequence, const int &k);
  void countKmer(); // counting kmer.
  void sortKmers(); // sorted by the count.
  void samplingKmers(const float &idenCutoff,
                     const int &kCountCutoff); // select kmers to get monomer by FPS.
  std::vector<int> calKmerDistList(const std::string &kmer); // get kmer distance.
  bool checkRepeative(const float &cutoff); // check for repeatability
  void printCounter();
  void printKmers();

  ~KMERCOUNTER() {}
} KC;

void calEditDist(const std::string &query, const std::string &ref, int &edist);
void seqAndSetIdentity(const std::string &query,
                       const std::vector<std::string> &refs, float &Identity);
template <typename T1>
void printOneDimVector(const std::vector<T1> &thisV){
  for (int i = 0; i < thisV.size(); i++)
    std::cout << thisV[i] << std::endl;
}
#endif