#include "sampleDBSCAN.h"

void DBSCAN::runDBSCAN() {
  corePoints.resize(dataSize, false);
  labels.resize(dataSize, -1);

  for (int i = 0; i < dataSize; ++i) {
    if (countSimilarNeighbors(i) >= minPts) {
      corePoints[i] = true;
    }
  }

  for (int i = 0; i < dataSize; ++i) {
    if (labels[i] == -1 && corePoints[i]) {
      expandCluster(i);
      ++currentCluster;
    }
  }
}

int DBSCAN::countSimilarNeighbors(const int &point) {
  int count = 0;
  for (int i = 0; i < dataSize; ++i) {
    if (similarityMatrix[point][i] >= eps) {
      ++count;
    }
  }
  return count;
}

void DBSCAN::expandCluster(const int &point) {
  labels[point] = currentCluster;

  for (int i = 0; i < dataSize; ++i) {
    if (similarityMatrix[point][i] >= eps && labels[i] == -1) {
      labels[i] = currentCluster;
      if (corePoints[i]) {
        expandCluster(i);
      }
    }
  }
}