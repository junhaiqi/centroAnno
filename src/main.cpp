#include "hor.h"
#include "ketopt.h"
#include "kseq.h"
#include "lightWightMSA.h"
#include "monoRefine.h"
#include "monomer.h"
#include "moutput.h"
#include "sequenceUtils.h"
#include "version_centroAnno.h"
#include "genome.h"

KSEQ_INIT(gzFile, gzread)

static void loadMonoTemps(const char *faPath, std::vector<SEQ> &bestMonos) {
  gzFile fp;
  kseq_t *ks;
  if ((fp = gzopen(faPath, "r")) == 0) {
    printf("%s not found!\n", faPath);
    exit(0);
  }

  ks = kseq_init(fp);
  while (kseq_read(ks) >= 0) {
    std::string seq = ks->seq.s;
    std::string name = ks->name.s;
    bestMonos.push_back(SEQ(seq, name));
  }
}

static void depSeq2Monoblock(const char *fn, int &windowSize,
                             const float &repCutoff, const bool &hpc,
                             const int &k, const int &threadNum,
                             const float &fpsCutoff, const float &epsilon,
                             const std::string &outDir,
                             const char *monoTemFa = "",
                             const int maxHORLen = 50,
                             const int lengthCutoff = 5000,
                             const bool onlyScan = false) {
  gzFile fp;
  kseq_t *ks;
  if ((fp = gzopen(fn, "r")) == 0)
    return;
  ks = kseq_init(fp);
  if (onlyScan) {
    while (kseq_read(ks) >= 0) {
      std::string seq = ks->seq.s;
      std::string name = ks->name.s;
      if (seq.length() < lengthCutoff) {
        continue;
      }
      windowSize = seq.length() > windowSize ? windowSize : seq.length();
      SEQ sequence(seq.substr(0, windowSize), name);
      KC kCounter(sequence.seq, k); // kmer counter

      if (hpc) {

        sequence.compressSeq();          // HPC
        kCounter.seq = sequence.compSeq; // kmer counter
      }

      double rep = kCounter.getRepeatRedio();
      if (rep >= repCutoff) {
        std::cout << ">" << name << "\n";
        std::cout << seq << "\n";
      }
    }
    return;
  }
  while (kseq_read(ks) >= 0) {
    std::string seq = ks->seq.s;
    std::string name = ks->name.s;
    std::cout << name << " Length: " << seq.length() << "\n";
    if (seq.length() < lengthCutoff) {
      std::cout << "The length of " << name << " is less than " << lengthCutoff << " and is not suitable for repeat annotation." << std::endl;
      continue;
    }
    windowSize = seq.length() > windowSize ? windowSize : seq.length();
    SEQ mySeq(seq.substr(0, windowSize),
              name); // Extract seuence to infer monomer template.

    std::vector<std::string> monoSeqs;
    std::vector<SEQ> bestMonos;

    if (monoTemFa == "") {
      monoTempInference(mySeq, k, fpsCutoff, 4, monoSeqs, bestMonos, hpc,
                        repCutoff);
      if (bestMonos.empty()) {
        std::cout << name << " may not be a reapeating sequence. (empty inferred templates)\n";
        continue;
      }
    }

    else {
      loadMonoTemps(monoTemFa, bestMonos);
    }

    std::vector<MonomerAlignment> demcompBlockRes;
    demcomposer2(mySeq, bestMonos, demcompBlockRes, 1);
    // NOTE: the above code needs to be modified when our method is applied to
    // more general sequences.

    vector<int> clusterLabels;
    vector<std::string> refinedBestMonomerTems; // Final monomer templates.
    sampleClustering(mySeq, demcompBlockRes, clusterLabels,
                     refinedBestMonomerTems, threadNum, epsilon);

    const SEQ longSequence(seq, name);
    vector<MonomerAlignment> newDemcompBlockRes;
    const std::string demResFileName =
        outDir + "/" + name + "_decomposedResult.csv";
    const std::string horDemResFileName =
        outDir + "/" + name + "_horDecomposedResult.csv";
    const std::string monoTemsFileName =
        outDir + "/" + name + "_monomerTemplates.fa";
    const std::string horFileName = outDir + "/" + name + "_HORs.fa";
    refineDemcompose(longSequence, refinedBestMonomerTems, newDemcompBlockRes,
                     threadNum, hpc, k, fpsCutoff, repCutoff, demResFileName,
                     monoTemsFileName);

    if (newDemcompBlockRes.size() > 10) {
      std::vector<HORINFO> horDecompBlockRes;
      std::vector<std::string> cleanedHORs;    // format is name, like 2_1_3...
      std::vector<std::string> cleanedHORSeqs; // format is base sequence
      inferHORs(newDemcompBlockRes, refinedBestMonomerTems, longSequence.seq,
                horDecompBlockRes, cleanedHORs, cleanedHORSeqs, maxHORLen);
      writeHORDecomposeRes(horDecompBlockRes, horDemResFileName);
      writeInferedHORs(cleanedHORs, cleanedHORSeqs, horFileName);
    }
  }
  kseq_destroy(ks);
  gzclose(fp);
}

void annoGenome(const char *fn, int &windowSize, const float &repCutoff,
                const bool &hpc, const int &k, const int &threadNum,
                const float &fpsCutoff, const float &epsilon,
                const std::string &outDir, const char *monoTemFa = "",
                const int maxHORLen = 50, const int lengthCutoff = 5000,
                const bool onlyScan = false, const int maxRegionLength = 500000, const int minRegionLength = 100, const float genomeCutoff = 0.8) {
  gzFile fp;
  kseq_t *ks;
  if ((fp = gzopen(fn, "r")) == 0)
    return;
  ks = kseq_init(fp);

  while (kseq_read(ks) >= 0) {
    std::string seq = ks->seq.s;
    std::string name = ks->name.s;
    std::cout << name << " Length: " << seq.length() << "\n";
    if (seq.length() < lengthCutoff) {
      std::cout << "The length of " << name << " is less than " << lengthCutoff << " and is not suitable for repeat annotation." << std::endl;
      continue;
    }

    std::vector<std::pair<int, int>> repRegions;
    SEQ sequence(seq, name);
    std::cout << "Detecting tandem repeat regions in the " << name << "...\n";
    repeatRegionInference(sequence, k, repCutoff, lengthCutoff, hpc,
                          repRegions, maxRegionLength);
    std::cout << "Detecting end.\n";
    std::cout << "There are " << repRegions.size() << " regions.\n";
    std::cout << "\n";

    const std::string demResFileName =
        outDir + "/" + name + "_decomposedResult.csv";
    const std::string horDemResFileName =
        outDir + "/" + name + "_horDecomposedResult.csv";
    const std::string monoTemsFileName =
        outDir + "/" + name + "_monomerTemplates.fa";
    const std::string horFileName = outDir + "/" + name + "_HORs.fa";

    const int tempWindowSize = windowSize;
    for (size_t i = 0; i < repRegions.size(); ++i) {
      std::cout << "***************************" << "Region " << i << ":Pos:" << repRegions[i].first << ":" << repRegions[i].second << "***************************\n";
      std::string thisSeq = seq.substr( repRegions[i].first, repRegions[i].second - repRegions[i].first );
      std::cout << "Repeat Region Length:" << thisSeq.length() << "\n";
      const SEQ longSequence(thisSeq, name);
      windowSize = thisSeq.length() > tempWindowSize ? tempWindowSize : thisSeq.length();
      // std::cout << "Window Length For :" << windowSize << "\n";
      SEQ mySeq(thisSeq.substr(0, windowSize),
                name); // Extract seuence to infer monomer template.

      std::vector<std::string> monoSeqs;
      std::vector<SEQ> bestMonos;

      if (monoTemFa == "") {
        monoTempInference(mySeq, k, fpsCutoff, 4, monoSeqs, bestMonos, hpc,
                          repCutoff);
        if (bestMonos.empty()) {
          std::cout << name << " may not be a reapeating sequence. (empty inferred templates)\n";
          continue;
        }
      }

      else {
        loadMonoTemps(monoTemFa, bestMonos);
      }

      std::vector<MonomerAlignment> demcompBlockRes;
      demcomposer2(mySeq, bestMonos, demcompBlockRes, 1);
      // NOTE: the above code needs to be modified when our method is applied to
      // more general sequences.

      vector<int> clusterLabels;
      vector<std::string> refinedBestMonomerTems; // Final monomer templates.
      sampleClustering(mySeq, demcompBlockRes, clusterLabels,
                       refinedBestMonomerTems, threadNum, epsilon);
      
      vector<MonomerAlignment> newDemcompBlockRes;
      refineDemcomposeGenome(longSequence, refinedBestMonomerTems, newDemcompBlockRes,
                       threadNum, hpc, k, fpsCutoff, repCutoff);

      std::vector<region> goodRegions;
      findGoodDemcomRegion(newDemcompBlockRes, goodRegions, minRegionLength, genomeCutoff);
      writeInferedMonosForGenome(refinedBestMonomerTems, monoTemsFileName, repRegions[i].first, repRegions[i].second, name);

      for (size_t j = 0; j < goodRegions.size(); ++j) {
        const int regionSize = goodRegions[j].end - goodRegions[j].start;
        std::vector<MonomerAlignment> goodBlockRegion( newDemcompBlockRes.begin() + goodRegions[j].blockStartID, 
                                      newDemcompBlockRes.begin() + goodRegions[j].blockEndID);
        

        writeDecomposeResForGenome(goodBlockRegion, demResFileName, repRegions[i].first, repRegions[i].second);

        if (regionSize > lengthCutoff && goodBlockRegion.size() > 10) {
          std::vector<HORINFO> horDecompBlockRes;
          std::vector<std::string> cleanedHORs;    // format is name, like 2_1_3...
          std::vector<std::string> cleanedHORSeqs; // format is base sequence
          inferHORs(goodBlockRegion, refinedBestMonomerTems, longSequence.seq,
                    horDecompBlockRes, cleanedHORs, cleanedHORSeqs, maxHORLen);

          writeHORDecomposeResForGenome(horDecompBlockRes, horDemResFileName, repRegions[i].first, repRegions[i].second);
          writeInferedHORsForGenome(cleanedHORs, cleanedHORSeqs, horFileName, repRegions[i].first, repRegions[i].second, name);
        }
      }

      std::cout << "***************************" << "Region " << i << ":Pos:" << repRegions[i].first << ":" << repRegions[i].second << "***************************\n";
      std::cout << "\n";
    }
  }

  kseq_destroy(ks);
  gzclose(fp);
}

int main(int argc, char *argv[]) {
  ketopt_t o = KETOPT_INIT;
  int32_t c;
  std::string outDir;
  int k = 11;
  float fpsCutoff = 0.6;
  float repCutoff = 0.2;
  int windowSize = 500000;
  bool hpc = 1;
  int threadNum = 8;
  float epsilon = 0.95;
  char *monoTemFaFile = "";
  int maxHORLen = 50;
  int lengthCutoff = 1000;
  int maxRegionLength = 1000000;
  int minRegionLength = 100;
  float genomeCutOff = 0.8;
  bool onlyScan = 0;
  // bool genome = 0;
  std::string annoMode = "anno-sat-asm";
  // std::vector<std::string> annoModeLst = {"anno-sat-asm", "anno-read", "anno-asm"};
  while ((c = ketopt(&o, argc, argv, 1, "o:k:f:r:w:c:e:t:m:x:M:L:S:G:A:N:F:", 0)) >= 0) {
    if (c == 'o') {
      if (o.arg == 0) {
        fprintf(stderr, "Error: -o requires an argument.\n");
        return 1;
      }
      outDir = o.arg;
    } else if (c == 'k')
      k = atoi(o.arg);
    else if (c == 'f') {
      fpsCutoff = atof(o.arg);
      if (fpsCutoff < 0) {
        fprintf(
            stderr,
            "Error: -f requires an argument larger than 0 and less than 1.\n");
        return 1;
      } else if (fpsCutoff > 1) {
        fprintf(stderr, "Warning: -f requires an argument larger than 0 and "
                        "less than 1.\n");
      }
    }

    else if (c == 'r') {
      repCutoff = atof(o.arg);
      if (repCutoff < 0 || repCutoff > 0.95) {
        fprintf(stderr, "Error: -r requires an argument larger than 0 and less "
                        "than 0.95.\n");
        return 1;
      }
    }

    else if (c == 'w') {
      windowSize = atoi(o.arg);
      if (windowSize < 10000) {
        fprintf(stderr, "Warning: the value of -w is recommended to be greater "
                        "than 10000.\n");
        return 1;
      }
    }

    else if (c == 'e') {
      epsilon = atof(o.arg);
      if (epsilon < 0.8) {
        fprintf(stderr, "Warning: the value of -e is recommended to be greater "
                        "than 0.8.\n");
        return 1;
      } else if (epsilon < 0) {
        fprintf(
            stderr,
            "Error: -e requires an argument larger than 0 and less than 1.\n");
        return 1;
      }
    }

    else if (c == 'c')
      hpc = atoi(o.arg);

    else if (c == 'x')
    {
      annoMode = o.arg;
      if(annoMode != "anno-sat-asm" && annoMode != "anno-asm" && annoMode != "anno-read")
      {
        std::cerr << "-x requires an argument is anno-sat-asm, anno-asm, or anno-read\n";
        return 1;
      }
    }

    else if (c == 't')
      threadNum = atoi(o.arg);

    else if (c == 'm')
      monoTemFaFile = o.arg;

    else if (c == 'M')
      maxHORLen = atoi(o.arg);
    
    else if (c == 'L') {
      lengthCutoff = atoi(o.arg);
    
      if (lengthCutoff < 100) {
        fprintf(stderr, "Warning: the value of -L is recommended to be greater "
                        "than 100.\n");
        return 1;
      }
    }

    else if (c == 'S')
      onlyScan = atoi(o.arg);
      
    // else if (c == 'G')
    //   genome = atoi(o.arg);

    else if (c == 'A')
      maxRegionLength = atoi(o.arg);

    else if (c == 'N')
      minRegionLength = atoi(o.arg);

    else if (c == 'F')
      genomeCutOff = atof(o.arg);
      
  }
  if (outDir.empty() || argc - o.ind < 1) {
    fprintf(stderr, "Version %d.%d.%d\n", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
    fprintf(stderr, "Usage: %s [Options:] <in.fa>\n", argv[0]);
    fprintf(stderr, "Options:\n");
    fprintf(stderr,
            "  -o STR     Specify the output folder [required parameters]\n");
    fprintf(stderr, "  -m STR     Specify the monomer template file with fasta "
                    "type [default = None]\n");
    fprintf(stderr, "  -k INT     Specify the k-mer size [default = 11]\n");
    fprintf(stderr, "  -f FLOAT   Specify the fps cutoff [default = 0.6]\n");
    fprintf(stderr,
            "  -r FLOAT   Specify the repeat redio cutoff [default = 0.2]\n");
    fprintf(stderr, "  -w INT     Specify the window size for infering "
                    "templates [default = 500000]\n");
    fprintf(stderr, "  -c BOOL    Specify the homopolymer compression (1: yes, 0: no) "
                    "[default = 1]\n");
    fprintf(stderr, "  -e FLOAT   Specify the indentity cutoff for DBSCAN "
                    "[default = 0.95]\n");
    fprintf(stderr, "  -x STR   Specify the annotated data type (anno-sat-asm: annotate centromeric alpha-satellite sequence (HiCAT/HORmon-like input), anno-asm: annotate chromosome/assembly, or anno-read: annotate sequencing reads)\n"
      "[default is anno-sat-asm]\n");
    fprintf(stderr, "  -t INT     Specify the number of threads for template "
                    "inference [default = 8]\n");
    fprintf(stderr, "  -M INT     Specify the maximum number of monomers that "
                    "a HOR can contain [default = 50]\n");
    fprintf(stderr, "  -L INT     Specify the length cutoff that "
                    "the annotated sequence needs to meet [default = 1000]\n");
    fprintf(stderr, "  -A INT     Specify the maxinum length cutoff that "
                    "the annotated region in the genome needs to meet for speed [default = 1000000]\n");
    fprintf(stderr, "  -N INT     Specify the mininum length cutoff that "
                    "the annotated region in the genome needs to meet for accuracy [default = 100]\n");
    fprintf(stderr, "  -F FLOAT   Specify the indentity cutoff for genome annotation "
                    "[default = 0.8]\n");
    fprintf(stderr, "  -S BOOL    Specify the repeated sequences are scanned "
                    "out without annotation (1: yes, 0: no) [default = 0]\n");
    // fprintf(stderr, "  -G BOOL    Specify the tendem repeat of chromosome/assembly are annotated (1: yes, 0: no) "
    //                 "[default = 0]\n");
    fprintf(stderr, "  example command: %s -o test test.fa\n", argv[0]);
    return 1;
  }

  if (annoMode == "anno-sat-asm")
  {
    createDirectory(outDir);
    depSeq2Monoblock(argv[o.ind], windowSize, repCutoff, hpc, k, threadNum,
                     fpsCutoff, epsilon, outDir, monoTemFaFile, maxHORLen,
                     lengthCutoff, onlyScan);
  }

  else if (annoMode == "anno-asm")
  {
    createDirectory(outDir);
    annoGenome(argv[o.ind], windowSize, repCutoff, hpc, k, threadNum, fpsCutoff,
               epsilon, outDir, monoTemFaFile, maxHORLen, lengthCutoff,
               onlyScan, maxRegionLength, minRegionLength, genomeCutOff);
  }

  else // anno-read
  {
    createDirectory(outDir);
    k = 9;
    hpc = 0;
    depSeq2Monoblock(argv[o.ind], windowSize, repCutoff, hpc, k, threadNum,
                     fpsCutoff, epsilon, outDir, monoTemFaFile, maxHORLen,
                     lengthCutoff, onlyScan);
  }

  return 0;
}
