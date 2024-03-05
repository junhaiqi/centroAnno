#ifndef SAMPLEDBSCAN_H
#define SAMPLEDBSCAN_H

#include <iostream>
#include <vector>

class DBSCAN {
public:
  DBSCAN(float eps, int minPts,
         const std::vector<std::vector<float>> &similarityMatrix)
      : eps(eps), minPts(minPts), similarityMatrix(similarityMatrix),
        dataSize(similarityMatrix.size()), currentCluster(0) {}

  void runDBSCAN();

  const std::vector<int> &getLabels() const { return labels; }

private:
  float eps;
  int minPts;
  const std::vector<std::vector<float>> &similarityMatrix;
  int dataSize;
  std::vector<bool> corePoints;
  std::vector<int> labels;
  int currentCluster;
  int countSimilarNeighbors(const int &point);
  void expandCluster(const int &point);
};

#endif // !SAMPLEDBSCAN
