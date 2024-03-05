#include "hor.h"
#include "ketopt.h"
#include "kseq.h"
#include "lightWightMSA.h"
#include "monoRefine.h"
#include "monomer.h"
#include "moutput.h"
#include "sequenceUtils.h"


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
                             const int maxHORLen = 50) {
  gzFile fp;
  kseq_t *ks;
  if ((fp = gzopen(fn, "r")) == 0)
    return;
  ks = kseq_init(fp);
  while (kseq_read(ks) >= 0) {
    std::string seq = ks->seq.s;
    std::string name = ks->name.s;
    windowSize = seq.length() > windowSize ? windowSize : seq.length();
    SEQ mySeq(seq.substr(0, windowSize),
              name); // Extract seuence to infer monomer template.

    std::vector<std::string> monoSeqs;
    std::vector<SEQ> bestMonos;

    if (monoTemFa == "") {
      monoTempInference(mySeq, k, fpsCutoff, 4, monoSeqs, bestMonos, hpc,
                        repCutoff);
      if (bestMonos.empty()) {
        std::cout << name << " may not be a reapeating sequence.\n";
        return;
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

    std::vector<HORINFO> horDecompBlockRes;
    std::vector<std::string> cleanedHORs;    // format is name, like 2_1_3...
    std::vector<std::string> cleanedHORSeqs; // format is base sequence
    // std::cout << "refine pass\n";
    inferHORs(newDemcompBlockRes, refinedBestMonomerTems, longSequence.seq,
              horDecompBlockRes, cleanedHORs, cleanedHORSeqs, maxHORLen);

    writeHORDecomposeRes(horDecompBlockRes, horDemResFileName);
    writeInferedHORs(cleanedHORs, cleanedHORSeqs, horFileName);
  }
  kseq_destroy(ks);
  gzclose(fp);
}

int main(int argc, char *argv[]) {
  ketopt_t o = KETOPT_INIT;
  int32_t c;
  std::string outDir;
  int k = 13;
  float fpsCutoff = 0.6;
  float repCutoff = 0.2;
  int windowSize = 500000;
  bool hpc = true;
  int threadNum = 8;
  float epsilon = 0.95;
  char *monoTemFaFile = "";
  int maxHORLen = 50;
  while ((c = ketopt(&o, argc, argv, 1, "o:k:f:r:w:c:e:t:m:M:", 0)) >= 0) {
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
      hpc = false;

    else if (c == 't')
      threadNum = atoi(o.arg);

    else if (c == 'm')
      monoTemFaFile = o.arg;

    else if (c == 'M')
      maxHORLen = atoi(o.arg);
  }
  if (outDir.empty() || argc - o.ind < 1) {
    fprintf(stderr, "Usage: %s [Options:] <in.fa>\n", argv[0]);
    fprintf(stderr, "Options:\n");
    fprintf(stderr,
            "  -o STR     Specify the output folder [required parameters]\n");
    fprintf(stderr, "  -m STR     Specify the monomer template file with fasta "
                    "type [default = None]\n");
    fprintf(stderr, "  -k INT     Specify the k-mer size [default = 13]\n");
    fprintf(stderr, "  -f FLOAT   Specify the fps cutoff [default = 0.6]\n");
    fprintf(stderr,
            "  -r FLOAT   Specify the repeat redio cutoff [default = 0.2]\n");
    fprintf(stderr, "  -w INT     Specify the window size for infering "
                    "templates [default = 500000]\n");
    fprintf(stderr, "  -c BOOL    Specify closing the homopolymer compression "
                    "[default = false]\n");
    fprintf(stderr, "  -e FLOAT   Specify the indentity cutoff for DBSCAN "
                    "[default = 0.95]\n");
    fprintf(stderr, "  -t INT     Specify the number of threads for template "
                    "inference [default = 8]\n");
    fprintf(stderr, "  -M INT     Specify the maximum number of monomers that "
                    "a HOR can contain [default = 50]\n");
    fprintf(stderr, "  example command: %s -o test test.fa\n", argv[0]);
    return 1;
  }

  // const std::string file =
  // "testData/ourBaselineFirstDecomPose/chm13_cen4/decomposedResult.csv";
  std::vector<std::string> files = {
      "testData/ourBaselineFirstDecomPose/chm13_cen1/"
      "chr1:121796048-126300487_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen2/"
      "chr2:92333543-94673023_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen3/"
      "chr3:91738002-96415026_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen4/"
      "chr4:49705154-55199795_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen5/"
      "chr5:47039134-49596625_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen6/"
      "chr6:58286706-61058390_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen7/"
      "chr7:60414372-63714499_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen8/"
      "chr8:44215832-46325080_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen9/"
      "chr9:44951775-47582595_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen10/"
      "chr10:39633793-41664589_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen11/"
      "chr11:51061948-54413484_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen12/"
      "chr12:34620838-37202490_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen13/"
      "chr13:15547593-17498291_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen14/"
      "chr14:10092112-12708411_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen15/"
      "chr15:16678794-17694466_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen16/"
      "chr16:35854528-37793352_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen17/"
      "chr17:23892419-27486939_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen18/"
      "chr18:15971633-20740248_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen19/"
      "chr19:25832447-29749519_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen20/"
      "chr20:26925852-29099655_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen21/"
      "chr21:10962853-11303831_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cen22/"
      "chr22:12788180-15711065_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cenX/"
      "chrX:57819763-60927195_decomposedResult.csv",
      "testData/ourBaselineFirstDecomPose/chm13_cenY/"
      "ChrY:10565750-10883085_decomposedResult.csv"};

  // std::vector<std::string> HORs;
  // horTest("testData/ourBaselineFirstDecomPose/chm13_cenX/chrX:57819763-60927195_decomposedResult.csv",
  //         "testData/ourBaselineFirstDecomPose/chm13_cenX/chrX:57819763-60927195_monomerTemplates.fa",
  //         "../human_CHM13_2.0_centro_regions/chm13_cenX.fasta",
  //         HORs);

  // HORs = {};
  // horTest("testData/ourBaselineFirstDecomPose/chm13_cen5/chr5:47039134-49596625_decomposedResult.csv",
  //         "testData/ourBaselineFirstDecomPose/chm13_cen5/chr5:47039134-49596625_monomerTemplates.fa",
  //         "../human_CHM13_2.0_centro_regions/chm13_cen5.fasta",
  //         HORs);

  // HORs = {};
  // horTest("testData/ourBaselineFirstDecomPose/chm13_cen1/chr1:121796048-126300487_decomposedResult.csv",
  //         "testData/ourBaselineFirstDecomPose/chm13_cen1/chr1:121796048-126300487_monomerTemplates.fa",
  //         "../human_CHM13_2.0_centro_regions/chm13_cen1.fasta",
  //         HORs);
  // horTest("testData/ourBaselineFirstDecomPose/chm13_cen18/chr18:15971633-20740248_decomposedResult.csv");

  // for (const std::string &file : files){
  //   std::vector<std::string> HORs;
  //   horTest(file, HORs);
  //   break;
  // }

  // std::cout << k << std::endl;
  // std::cout << maxHORLen << std::endl;

  createDirectory(outDir);
  depSeq2Monoblock(argv[o.ind], windowSize, repCutoff, hpc, k, threadNum,
                   fpsCutoff, epsilon, outDir, monoTemFaFile, maxHORLen);

  /*
  main code:
    createDirectory(outDir);
    depSeq2Monoblock(argv[o.ind], windowSize, repCutoff, hpc, k, threadNum,
                      fpsCutoff, outDir);
  */

  /*
  test code:
    const std::string demResFile = outDir + "/decomposedResult.csv";
    const std::string monoTemplatesFile = outDir + "/monomerTemplates.fa";
    testSEQ(argv[o.ind], demResFile, monoTemplatesFile);
  */
  return 0;
}
