#include "monomer.h"
using namespace std;

size_t tolLen = 5;
void monoTempInference(SEQ &sequence, const int &k, const float &kIdenCf,
                       const int &kCountCf, std::vector<std::string> &monoSeqs,
                       std::vector<SEQ> &bestMonoTemplateList,
                       const bool &hpc, const float &repPorCutoff) {
  
  KC kCounter(sequence.seq, k); // kmer counter
  
  if (hpc){
    sequence.compressSeq();           // HPC
    kCounter.seq = sequence.compSeq; // kmer counter
  }

  bool rep = kCounter.checkRepeative(repPorCutoff);
  if (!rep){
    std::cout << sequence.seqName << " has a subsequence (length: " << sequence.seq.size() << ") not a repeative seuence (because of low non-unique k-mers)." << std::endl;
    return;
  }
  kCounter.samplingKmers(kIdenCf, kCountCf);
  if (kCounter.selectedKmers.empty()) {
    printf("This sequence may not contain repeat unit!\n");
    return;
  }
  for (std::string &kmer : kCounter.selectedKmers) {
    std::vector<int> posList = kCounter.h[kmer];
    std::vector<int> diffValueList = kCounter.calKmerDistList(kmer);
    std::vector<int> sortedIdxList = sortReturnIdxList(diffValueList);
    std::sort(diffValueList.begin(), diffValueList.end());
    int lastItem = 0;
    for (size_t i = 0; i < diffValueList.size(); i++) {
      if (diffValueList[i] > lastItem) {
        int repeatMightMaxLen = diffValueList[i] + tolLen;
        size_t suppRepeatLenCount = 0;
        for (size_t j = i; j < diffValueList.size(); j++) {
          if (diffValueList[j] < repeatMightMaxLen) {
            suppRepeatLenCount += 1;
          }
        }
        lastItem = diffValueList[i] + tolLen;
        size_t suppCf =
            10 < diffValueList.size() * 0.1 ? 10 : diffValueList.size() * 0.1;
        if (suppRepeatLenCount > suppCf) {
          std::vector<std::string> monoList;
          std::vector<int> subDistList(diffValueList.begin() + i,
                                       diffValueList.begin() + i +
                                           suppRepeatLenCount);
          int maxCount = 0;
          int mode = findMode(subDistList, maxCount);
          if (maxCount >= 2 &&
              mode > 10) { // We find the repeats which length larger 10.
            auto it =
                std::find(diffValueList.begin(), diffValueList.end(), mode);
            int firstModePos = std::distance(diffValueList.begin(), it);
            std::vector<int> posIdx(sortedIdxList.begin() + firstModePos,
                                    sortedIdxList.begin() + firstModePos +
                                        maxCount);
            std::vector<std::vector<int>> posPairList(maxCount,
                                                      std::vector<int>(2));
            int tmp = 0;
            for (std::vector<int> &posPair : posPairList) {
              assert(posList.size() > (posIdx[tmp] + 1));
              posPair[0] = posList[ posIdx[tmp] ];
              posPair[1] = posList[ posIdx[tmp] + 1 ];
              std::string momomer;
              if (hpc) {
                momomer = sequence.compSeq.substr(posPair[0], mode);
              }
              else{
                momomer = sequence.seq.substr(posPair[0], mode);
              }
              monoList.push_back(momomer);
              tmp += 1;
            }
            float maxMeanScoreOfMonomer;
            int maxAvgIndex;
            
            if (monoList.size() < 2)
              return;

            evalMonomer(monoList, maxMeanScoreOfMonomer, maxAvgIndex);

            if (maxMeanScoreOfMonomer > 0.9) {
              if (hpc){
                assert(posPairList.size() >  maxAvgIndex);
                std::string rSeq = sequence.recoverSeg(
                posPairList[maxAvgIndex][0], posPairList[maxAvgIndex][1]);
                monoSeqs.push_back(rSeq);
                break;
              }
              else {
                monoSeqs.push_back(sequence.seq.substr(posPairList[maxAvgIndex][0], posPairList[maxAvgIndex][1] - posPairList[maxAvgIndex][0]));
              }
            }
          }
        }
      }
    }
  }
  if (monoSeqs.empty()){
    return;
  }
  std::vector<int> monoLenList(monoSeqs.size(), 0);
  for (size_t i = 0; i < monoSeqs.size(); i++) {
    monoLenList[i] = monoSeqs[i].length();
  }
  int repeatIdx = 0;
  inferRepeatLen(monoLenList, repeatIdx);
  int bestLen = monoLenList[repeatIdx];
  size_t idx = 0;
  for (const std::string mono : monoSeqs) {
    if (mono.length() == bestLen) {
      SEQ monoClass(mono, std::to_string(idx));
      bestMonoTemplateList.push_back(monoClass);
      idx++;
    }
  }
}

void demcomposer2(const SEQ &sequence, std::vector<SEQ> &bestMonoTemplateList,
                  std::vector<MonomerAlignment> &demcompBlockRes, const int &threadNum){
  int ins = -1, del = -1, mismatch = -1, match = 1;
  // int ins = -1, del = 0, mismatch = 0, match = 1;
  std::vector<Seq> reads = {Seq(sequence.seqName, sequence.seq)};
  vector<Seq> monomers;
  for (const SEQ &mono : bestMonoTemplateList){
    monomers.push_back(Seq(mono.seqName, mono.seq));
  }
  add_reverse_complement(monomers);
  MonomersAligner monomers_aligner(monomers, ins, del, mismatch, match);
  monomers_aligner.AlignReadsSet(reads, demcompBlockRes, threadNum, 5000, -1, 500);
}

void demcomposer(const SEQ &sequence, std::vector<SEQ> &bestMonoTemplateList,
                 std::vector<MonomerAlignment> &demcompBlockRes) {
  int ins = -1, del = -1, mismatch = -1, match = 1;
  Seq read = Seq(sequence.seqName, sequence.seq);
  // name2Idx monoMap;
  vector<Seq> monomers;
  size_t idx = 0;
  for (const SEQ &mono : bestMonoTemplateList) {
    monomers.push_back(Seq(mono.seqName, mono.seq));
    // monoMap[mono.seqName] = idx;
    idx++;
  }
  add_reverse_complement(monomers);
  MonomersAligner monomers_aligner(monomers, ins, del, mismatch, match);
  demcompBlockRes = monomers_aligner.AlignPartClassicDP(read, monomers);
}

std::vector<int> sortReturnIdxList(const std::vector<int> &v) {
  std::vector<int> idxList(v.size());
  for (size_t i = 0; i < idxList.size(); ++i) {
    idxList[i] = i;
  }
  std::sort(idxList.begin(), idxList.end(),
            [&v](int i1, int i2) { return v[i1] < v[i2]; });
  return idxList;
}

void evalMonomer(std::vector<std::string> &monoList,
                 float &maxMeanScoreOfMonomer, int &maxAvgIndex) {
  std::vector<std::vector<float>> idenMatrix(
      monoList.size(), std::vector<float>(monoList.size(), 1));
  int l = monoList[0].length();
  // #pragma omp parallel for num_threads(2)
  for (size_t i = 0; i < monoList.size() - 1; i++) {
    for (size_t j = i + 1; j < monoList.size(); j++) {
      int edist = 0;
      calEditDist(monoList[i], monoList[j], edist);
      float identity = 1 - (static_cast<float>(edist) / l);
      idenMatrix[i][j] = identity;
      idenMatrix[j][i] = identity;
    }
  }

  std::vector<float> meanScoreList(monoList.size(), 0);
  // #pragma omp parallel for num_threads(2)
  for (size_t i = 0; i < monoList.size(); i++) {
    meanScoreList[i] = calculateMean(idenMatrix[i]);
  }

  maxMeanScoreOfMonomer =
      *max_element(meanScoreList.begin(), meanScoreList.end());
  maxAvgIndex = max_element(meanScoreList.begin(), meanScoreList.end()) -
                meanScoreList.begin();
}

void printStrLen(std::vector<std::string> seqList) {
  for (std::string s : seqList) {
    std::cout << s << ", length: " << s.length() << std::endl;
  }
}

void calLossList(const std::vector<int> &v, const int &num,
                 std::vector<float> &lossList) {

  for (int item : v) {
    if (item > num - num * 0.1) {
      float loss =
          static_cast<float>(item) / num - static_cast<float>(item / num);
      if (loss != 0) {
        float temp =
            static_cast<float>(item) / num - static_cast<float>(item / num);
        loss = temp < std::abs(temp - 1.0) ? temp : std::abs(temp - 1.0);
      }
      lossList.push_back(loss);
    } else {
      lossList.push_back(2.0);
    }
  }
}

void inferRepeatLen(const std::vector<int> &monoLenList, int &inferedIdx) {
  std::vector<int> matchList;
  // std::vector<float> aveLossList;

  for (int num : monoLenList) {
    std::vector<float> lossList;
    calLossList(monoLenList, num, lossList);
    int matchNum = std::count_if(lossList.begin(), lossList.end(),
                                 [](float item) { return item < 0.1; });
    // float aveLoss = std::accumulate(lossList.begin(), lossList.end(), 0.0) /
    // lossList.size(); aveLossList.push_back(aveLoss);
    matchList.push_back(matchNum);
  }

  auto maxMatchIdx = std::max_element(matchList.begin(), matchList.end());
  // auto minAveLossIdx = std::min_element(aveLossList.begin(),
  // aveLossList.end());

  int idx = std::distance(matchList.begin(), maxMatchIdx);
  // int idx2 = std::distance(aveLossList.begin(), minAveLossIdx);

  inferedIdx = idx;
}

void MonomersAligner::AlignReadsSet(vector<Seq> &reads, vector<MonomerAlignment> &batch, int threads,
                                    int part_size, int ed_thr,
                                    int overlap = 500) {
  vector<Seq> new_reads;
  vector<int> save_steps;
  for (const auto &r : reads) {
    int cnt = 0;
    // cout << r.seq.size() << endl;
    for (size_t i = 0; i < r.seq.size(); i += part_size) {
      if ((int)r.seq.size() - i >= overlap || r.seq.size() < overlap) {
        Seq seq = Seq(r.read_id.name,
                      r.seq.substr(i, min(part_size + overlap,
                                          static_cast<int>(r.seq.size() - i))),
                      i);
        new_reads.push_back(seq);
        ++cnt;
      }
    }
    save_steps.push_back(cnt);
  }
  // cerr << "Prepared reads\n";

  size_t start = 0, p = 0;
  int step = threads * 2;
  vector<pair<int, vector<MonomerAlignment>>> subbatches;
  for (size_t i = 0; i < new_reads.size(); i += step) {
#pragma omp parallel for num_threads(threads)
    for (size_t j = i; j < min(i + step, new_reads.size()); ++j) {
      std::vector<MonomerAlignment> aln;
      if (ed_thr > -1) {
        std::vector<Seq> filter_monomers =
            FilterMonomersForRead(new_reads[j], ed_thr);
        aln = AlignPartClassicDP(new_reads[j], filter_monomers);
      } else {
        aln = AlignPartClassicDP(new_reads[j], monomers_);
      }

#pragma omp critical(aligner)
      { subbatches.push_back(pair<int, vector<MonomerAlignment>>(j, aln)); }
    }
    sort(subbatches.begin() + i,
         subbatches.begin() + min(i + step, new_reads.size()), sortby1);
    while (p < save_steps.size() &&
           start + save_steps[p] <= subbatches.size()) {
      // vector<MonomerAlignment> batch;
      for (size_t j = start; j < start + save_steps[p]; ++j) {
        int read_index = subbatches[j].first;
        for (auto a : subbatches[j].second) {
          MonomerAlignment new_m_aln(
              a.monomer_name, a.read_name,
              new_reads[read_index].read_id.id + a.start_pos,
              new_reads[read_index].read_id.id + a.end_pos, a.identity, a.best);
          batch.push_back(new_m_aln);
        }
      }
      // cerr << (p + 1) * 100 / save_steps.size() << "%: Aligned "
      //      << batch[0].read_name << endl;
      batch = PostProcessing(batch);
      // SaveBatch(batch);
      start += save_steps[p];
      ++p;
    }
  }
  batch = ourPostProcessing(batch);
}

double MonomersAligner::MonomerEditDistance(Seq &monomer, Seq &read) {
  EdlibAlignResult result = edlibAlign(
      monomer.seq.c_str(), monomer.seq.size(), read.seq.c_str(),
      read.seq.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
  double res = result.editDistance;
  edlibFreeAlignResult(result);
  return res;
}

std::vector<Seq> MonomersAligner::FilterMonomersForRead(Seq &read, int ed_thr) {
  std::vector<Seq> monomers_for_read;
  std::vector<std::pair<double, int>> mn_edit;
  for (size_t i = 0; i < monomers_.size(); ++i) {
    mn_edit.push_back(
        std::make_pair(MonomerEditDistance(monomers_[i], read), i));
  }
  std::sort(mn_edit.begin(), mn_edit.end());
  monomers_for_read.push_back(monomers_[mn_edit[0].second]);
  for (size_t i = 1; i < mn_edit.size(); ++i) {
    if (mn_edit[i].first <= ed_thr) {
      monomers_for_read.push_back(monomers_[mn_edit[i].second]);
    }
  }
  return monomers_for_read;
}

vector<MonomerAlignment>
MonomersAligner::AlignPartClassicDP(Seq &read, std::vector<Seq> &monomers) {
  int ins = ins_;
  int del = del_;
  int match = match_;
  int mismatch = mismatch_;
  int INF = -1000000;
  int monomers_num = (int)monomers.size();
  vector<vector<vector<long long>>> dp(read.seq.size());
  // cout << dp.size() << endl;
  for (size_t i = 0; i < read.seq.size(); ++i) {
    for (const auto &m : monomers) {
      dp[i].push_back(vector<long long>(m.seq.size()));
      for (size_t k = 0; k < m.seq.size(); ++k) {
        dp[i][dp[i].size() - 1][k] = INF;
      }
    }
    dp[i].push_back(vector<long long>(1));
    dp[i][monomers_num][0] = INF;
  }

  for (size_t j = 0; j < monomers.size(); ++j) {
    Seq m = monomers[j];
    if (m.seq[0] == read.seq[0]) {
      dp[0][j][0] = match;
    } else {
      dp[0][j][0] = mismatch;
    }
    for (size_t k = 1; k < m.seq.size(); ++k) {
      long long mm_score = monomers[j].seq[k] == read.seq[0] ? match : mismatch;
      dp[0][j][k] =
          max(dp[0][j][k - 1] + del, (long long)(del * (k - 1) + mm_score));
    }
  }
  for (size_t i = 1; i < read.seq.size(); ++i) {
    for (size_t j = 0; j < monomers.size(); ++j) {
      dp[i][monomers_num][0] =
          max(dp[i][monomers_num][0], dp[i - 1][j][monomers[j].size() - 1]);
    }
    for (size_t j = 0; j < monomers.size(); ++j) {
      for (size_t k = 0; k < monomers[j].size(); ++k) {
        long long score = INF;
        int mm_score = monomers[j].seq[k] == read.seq[i] ? match : mismatch;
        if (dp[i][monomers_num][0] > INF) {
          score = max(score, dp[i][monomers_num][0] + mm_score +
                                 static_cast<long long>(k * del));
        }
        if (k > 0) {
          if (dp[i - 1][j][k - 1] > INF) {
            score = max(score, dp[i - 1][j][k - 1] + mm_score);
          }
          if (dp[i - 1][j][k] > INF) {
            score = max(score, dp[i - 1][j][k] + ins);
          }
          if (dp[i][j][k - 1] > INF) {
            score = max(score, dp[i][j][k - 1] + del);
          }
        }
        dp[i][j][k] = score;
      }
    }
  }
  int max_score = INF;
  int best_m = monomers_num;
  for (size_t j = 0; j < monomers.size(); ++j) {
    if (max_score < dp[read.seq.size() - 1][j][monomers[j].size() - 1]) {
      max_score = dp[read.seq.size() - 1][j][monomers[j].size() - 1];
      best_m = j;
    }
  }
  vector<MonomerAlignment> ans;
  long long i = read.seq.size() - 1;
  long long j = best_m;
  long long k = dp[i][j].size() - 1;
  bool monomer_changed = true;
  MonomerAlignment cur_aln;
  while (i >= 0) {
    if (k == static_cast<long long>(dp[i][j].size() - 1) && j != monomers_num &&
        monomer_changed) {
      cur_aln = MonomerAlignment(monomers[j].read_id.name, read.read_id.name, i,
                                 i, dp[i][j][k], true);
      monomer_changed = false;
    }
    if (j == monomers_num) {
      if (i != 0) {
        for (size_t p = 0; p < dp[i - 1].size(); ++p) {
          if (dp[i - 1][p][dp[i - 1][p].size() - 1] == dp[i][j][k]) {
            --i;
            j = p;
            k = dp[i][j].size() - 1;
            break;
          }
        }
      } else {
        --i;
      }
    } else {
      if (k != 0 && dp[i][j][k] == dp[i][j][k - 1] + del) {
        --k;
      } else {
        if (i != 0 && dp[i][j][k] == dp[i - 1][j][k] + ins) {
          --i;
        } else {
          int mm_score = monomers[j].seq[k] == read.seq[i] ? match : mismatch;
          if (i != 0 && k != 0 &&
              dp[i][j][k] == dp[i - 1][j][k - 1] + mm_score) {
            --i;
            --k;
          } else {
            monomer_changed = true;
            if (i != 0 &&
                dp[i][monomers_num][0] + k * del + mm_score == dp[i][j][k]) {
              cur_aln.start_pos = i;
              cur_aln.identity = cur_aln.identity - dp[i][monomers_num][0];
              ans.push_back(cur_aln);
              j = monomers_num;
              k = 0;
            } else {
              cur_aln.start_pos = i;
              ans.push_back(cur_aln);
              --i;
            }
          }
        }
      }
    }
  }
  reverse(ans.begin(), ans.end());
  return ans;
}

void MonomersAligner::SaveBatch(vector<MonomerAlignment> &batch) {
  // int prev_end = 0;
  for (auto a : batch) {
    string s = a.read_name + "\t" + a.monomer_name + "\t" +
               to_string(a.start_pos) + "\t" + to_string(a.end_pos) + "\t" +
               to_string(a.identity) + "\t" +
               //  to_string(a.start_pos - prev_end) + "\t" +
               to_string(a.end_pos - a.start_pos);
    // prev_end = a.end_pos;
    cout << s << "\n";
  }
}

vector<MonomerAlignment>
MonomersAligner::PostProcessing(vector<MonomerAlignment> &batch) {
  vector<MonomerAlignment> res;
  size_t i = 0;
  while (i < batch.size()) {
    for (size_t j = i + 1; j < min(i + 7, batch.size()); ++j) {
      if ((batch[i].end_pos - batch[j].start_pos) * 2 >
          (batch[j].end_pos - batch[j].start_pos)) {
        res.push_back(batch[i]);
        i = j + 1;
        break;
      }
    }
    if (i < batch.size())
      res.push_back(batch[i]);
    ++i;
  }
  return res;
}

vector<MonomerAlignment>
MonomersAligner::ourPostProcessing(vector<MonomerAlignment> &batch) {
  vector<MonomerAlignment> res;
  int lastPos = 0;
  for (size_t i = 0; i < batch.size(); i++){
    if (batch[i].start_pos > lastPos){
      res.push_back(batch[i]);
      lastPos = batch[i].end_pos;
    }
  } 
  return res;
}

string reverse_complement(const string &s) {
  string res = "";
  map<char, char> rc = {
      {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}, {'N', 'N'}};
  for (int i = (int)s.size() - 1; i >= 0; --i) {
    try {
      res += rc.at(s[i]);
    } catch (std::out_of_range &e) {
      cerr << e.what() << std::endl;
      exit(-1);
    }
  }
  return res;
}

void add_reverse_complement(vector<Seq> &monomers) {
  vector<Seq> rev_c_monomers;
  for (auto s : monomers) {
    rev_c_monomers.push_back(
        Seq(s.read_id.name + "'", reverse_complement(s.seq)));
  }
  monomers.insert(monomers.end(), rev_c_monomers.begin(), rev_c_monomers.end());
  return;
}

inline bool sortby1(const pair<int, vector<MonomerAlignment>> &a,
                    const pair<int, vector<MonomerAlignment>> &b) {
  return (a.first < b.first);
}
