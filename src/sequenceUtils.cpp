#include "sequenceUtils.h"

static int cmpBy2(const std::pair<std::string, int> &x,
                  const std::pair<std::string, int> &y) {
  return x.second > y.second;
} // sort by the second element of pair.

void calEditDist(const std::string &query, const std::string &ref, int &edist) {

  EdlibAlignResult result =
      edlibAlign(query.c_str(), query.length(), ref.c_str(), ref.length(),
                 edlibDefaultAlignConfig());
  if (result.status == EDLIB_STATUS_OK) {
    edist = result.editDistance;
  }
  edlibFreeAlignResult(result);
  
}

void seqAndSetIdentity(const std::string &query,
                       const std::vector<std::string> &refs, float &Identity) {
  // Note: In our problem, the length of query and ref is same.
  int l = query.length();
  float maxScore = 0;
  for (size_t i = 0; i < refs.size(); i++) {
    int edist = 0;
    calEditDist(query, refs[i], edist);
    float thisScore = 1.0 - (static_cast<float>(edist) / l);
    if (maxScore < thisScore)
      maxScore = thisScore;
  }
  Identity = maxScore;
}

SEQ::SEQ(const std::string &sequence, const std::string &name) : seq(sequence), seqName(name), comp(false) {std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);}

void SEQ::compressSeq() {
  compSeq = seq[0];
  char lastChar = seq[0];
  int compNum = 0;
  compIdxList.resize(seq.size());
  int t = 0;
  for (char currentChar : seq) {
    if (currentChar == lastChar) {
      compNum += 1;
    } else {
      compIdxList[t] = compNum;
      lastChar = currentChar;
      compNum = 1;
      compSeq += currentChar;
      t += 1;
    }
  }
  compIdxList[t] = compNum;
  compIdxList.resize(compSeq.size());
  comp = true;
}

std::string SEQ::recoverSeg(int startPos, int endPos) {
  std::string rSeq;
  std::vector<int> idxList(compIdxList.begin() + startPos,
                           compIdxList.begin() + endPos);
  std::string compSeg(compSeq.begin() + startPos, compSeq.begin() + endPos);

  for (size_t i = 0; i < compSeg.size(); ++i) {
    rSeq.append(idxList[i], compSeg[i]);
  }

  return rSeq;
}

KC::KMERCOUNTER(const std::string &sequence, const int &k)
    : seq(sequence), k(k) {}

void KC::countKmer() {
  if (k > seq.length()){
    // printf("Too large k (too short sequence)!");
    return;
  }
  // std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
  for (size_t i = 0; i <= seq.length() - k; ++i) {
    auto kmer = seq.substr(i, k);
    if (kmer.find('N') < kmer.length())
      continue;
    h[kmer].push_back(i);
  }
}

void KC::sortKmers() {
  if (h.empty())
    countKmer();
  counter_t ht;
  for (auto &key : h) {
    ht[key.first] = key.second.size();
  }
  std::vector<std::pair<std::string, int>> tVector;
  for (auto &key : ht) {
    tVector.push_back(make_pair(key.first, key.second));
  }
  std::sort(tVector.begin(), tVector.end(), cmpBy2);
  sortedKmers.resize(tVector.size());
  for (size_t i = 0; i < tVector.size(); i++) {
    sortedKmers[i] = tVector[i].first;
  }
}

void KC::samplingKmers(const float &idenCutoff, const int &kCountCutoff) {
  if (sortedKmers.empty())
    sortKmers();
  std::vector<std::string> tempSeqList(sortedKmers);
  int counIdx = 0;
  for (size_t i = 0; i < tempSeqList.size(); i++) {
    if (h[tempSeqList[i]].size() == kCountCutoff) {
      counIdx = i;
      break;
    }
  }
  tempSeqList.resize(counIdx + 1);
  const std::string initSeq = tempSeqList[0];
  tempSeqList.erase(tempSeqList.begin());
  selectedKmers.push_back(initSeq);
  for (std::string &seq : tempSeqList) {
    float maxIden;
    seqAndSetIdentity(seq, selectedKmers, maxIden);
    if (maxIden < idenCutoff)
      selectedKmers.push_back(seq);
  }
}

std::vector<int> KC::calKmerDistList(const std::string &kmer) {
  if (h.empty())
    countKmer();
  std::vector<int> resList(h[kmer].size() - 1, 0);
  for (size_t i = 0; i < h[kmer].size() - 1; i++){
    resList[i] = h[kmer][i+1] - h[kmer][i];
  }
  // std::sort(resList.begin(), resList.end());
  return resList;
}

bool KC::checkRepeative(const float &cutoff){
  if (h.empty())
    countKmer();
  size_t totalNum = 0;
  size_t repNum = 0;
  for (auto &key : h){
    totalNum += key.second.size();
    if (key.second.size() > 1) {
      repNum += key.second.size();
    }
  }
  bool repeative = false;
  repeativeScore = static_cast<double>(repNum) / totalNum;
  if (repeativeScore > cutoff){
    repeative = true;
  };
  return repeative;
}

double KC::getRepeatRedio() {
  if (h.empty())
    countKmer();
  size_t totalNum = 0;
  size_t repNum = 0;
  for (auto &key : h){
    totalNum += key.second.size();
    if (key.second.size() > 1) {
      repNum += key.second.size();
    }
  }
  repeativeScore = static_cast<double>(repNum) / totalNum;
  return repeativeScore;
}

void KC::printCounter() {
  for (auto &key : h) {
    for (int pos : key.second)
      printf("%s, %d\n", key.first.c_str(), pos);
  }
}

void KC::printKmers() {
  for (size_t i = 0; i < sortedKmers.size(); i++) {
    std::cout << sortedKmers[i] << std::endl;
  }
}
