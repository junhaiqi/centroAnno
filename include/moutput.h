#ifndef MOUTPUT_H
#define MOUTPUT_H
#include <fstream>
#include <sys/stat.h>
#include "hor.h"

static bool createDirectory(const std::string &folderName) {
#ifdef _WIN32
  return _mkdir(folderName.c_str()) == 0;
#else
  return mkdir(folderName.c_str(), 0777) == 0;
#endif
}

inline void writeDecomposeRes(const std::vector<MonomerAlignment> &decomoseRes,
                              const std::string &filename) {

  std::ofstream outFile(filename, std::ios::app);
  if (!outFile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  outFile << "sequence name,mononer name,start position,end position,estimated "
             "identity,length"
          << "\n";
  for (const auto &a : decomoseRes) {
    const size_t len = a.end_pos - a.start_pos + 1;
    std::string s = a.read_name + "," + a.monomer_name + "," +
                    std::to_string(a.start_pos) + "," +
                    std::to_string(a.end_pos) + "," +
                    std::to_string(a.identity) + "," + std::to_string(len);
    outFile << s << "\n";
  }
//   outFile.close();
}

inline void writeDecomposeResForGenome(const std::vector<MonomerAlignment> &decomoseRes,
                              const std::string &filename, const int &startPos, const int &endPos) {

  std::ofstream outFile(filename, std::ios::app);
  if (!outFile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  // outFile << "sequence name\tmononer name\tstart position,end position\testimated "
  //            "identity\tlength"
  //         << "\n";
  for (const auto &a : decomoseRes) {
    const size_t len = a.end_pos - a.start_pos + 1;
    std::string s = a.read_name + "," + "Pos:" + std::to_string(startPos) + ":" + std::to_string(endPos) + "_" + a.monomer_name + "," +
                    std::to_string(a.start_pos + startPos) + "," +
                    std::to_string(a.end_pos + startPos) + "," +
                    std::to_string(a.identity) + "," + std::to_string(len);
    outFile << s << "\n";
  }
//   outFile.close();
}

inline void writeInferedMonos(const std::vector<std::string> &inferedMonos,
                              const std::string &filename) {

  std::ofstream outFile(filename);
  if (!outFile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  for (size_t i = 0; i < inferedMonos.size(); i++) {
    std::string monoID = '>' + std::to_string(i);
    outFile << monoID << "\n";
    outFile << inferedMonos[i] << "\n";
  }
  outFile.close();
}

inline void writeInferedMonosForGenome(const std::vector<std::string> &inferedMonos,
                              const std::string &filename, const int &startPos, const int &endPos, const std::string &seqName) {

  std::ofstream outFile(filename, std::ios::app);
  if (!outFile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  for (size_t i = 0; i < inferedMonos.size(); i++) {
    std::string monoID = '>' + seqName + ':' + std::to_string(startPos) + ':' + std::to_string(endPos) + '_' + std::to_string(i);
    outFile << monoID << "\n";
    outFile << inferedMonos[i] << "\n";
  }
  outFile.close();
}

inline void writeInferedHORs(const std::vector<std::string> &horNames,
                             const std::vector<std::string> &horSeqs,
                             const std::string &filename) {

  std::ofstream outFile(filename);
  if (!outFile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  for (size_t i = 0; i < horNames.size(); i++) {
    std::string monoID = '>' + horNames[i];
    outFile << monoID << "\n";
    outFile << horSeqs[i] << "\n";
  }
  outFile.close();
}

inline void writeInferedHORsForGenome(const std::vector<std::string> &horNames,
                             const std::vector<std::string> &horSeqs,
                             const std::string &filename, const int &startPos, const int &endPos, const std::string &seqName) {

  std::ofstream outFile(filename, std::ios::app);
  if (!outFile.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  for (size_t i = 0; i < horNames.size(); i++) {
    std::string monoID = '>' + seqName + ':' + std::to_string(startPos) + ':' + std::to_string(endPos) + ':' + horNames[i];
    outFile << monoID << "\n";
    outFile << horSeqs[i] << "\n";
  }
  outFile.close();
}

inline void readDecomposeFileOfMonoNames(const std::string &filename,
                                         std::vector<std::string> &monoNames) {

  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening the file." << std::endl;
    return;
  }

  std::string line;
  std::vector<std::string> tokens;

  // Read and ignore the first line
  std::getline(file, line);

  while (std::getline(file, line)) {
    std::istringstream ss(line);
    std::string token;

    // Use getline to split the line by commas
    std::getline(ss, token, ',');

    std::getline(ss, token, ',');
    monoNames.push_back(token);
  }

  file.close();
}

inline std::vector<std::string> readFastaSequences(const std::string &filename) {
  std::vector<std::string> sequences;
  std::ifstream file(filename);

  if (file.is_open()) {
    std::string line;
    std::string currentSequence = "";

    while (std::getline(file, line)) {
      if (line.empty()) {
        continue; // Skip empty lines
      }

      if (line[0] == '>' && !currentSequence.empty()) {
        // Skip header lines
        sequences.push_back(currentSequence);
        currentSequence = "";
        continue;
      } else {
        // Append to the current sequence
        if (line[0] != '>')
            currentSequence += line;
      }
    }

    // Save the last sequence in the file
    if (!currentSequence.empty()) {
      sequences.push_back(currentSequence);
    }

    file.close();
  } else {
    std::cerr << "Error opening file: " << filename << std::endl;
  }

  return sequences;
}

#endif