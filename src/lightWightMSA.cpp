#include "lightWightMSA.h"

void msaAndConsensus(const std::vector<std::string> &seqList,
                     std::string &consensus) {

  auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3); // linear gaps

  spoa::Graph graph{};

  for (const std::string &it : seqList) {
    auto alignment = alignment_engine->Align(it, graph);
    graph.AddAlignment(alignment, it);
  }

  consensus = graph.GenerateConsensus();

//   std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
//             << consensus << std::endl;

//   auto msa = graph.GenerateMultipleSequenceAlignment();

//   for (const auto &it : msa) {
//     std::cerr << it << std::endl;
//   }
}