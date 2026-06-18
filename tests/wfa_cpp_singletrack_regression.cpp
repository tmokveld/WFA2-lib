/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 */

#include <cstdio>

#include "bindings/cpp/WFAligner.hpp"

class InspectableGapAffine : public wfa::WFAlignerGapAffine {
public:
  InspectableGapAffine(
      const wfa::WFAligner::AlignmentScope alignment_scope,
      const wfa::WFAligner::MemoryModel memory_model) :
        wfa::WFAlignerGapAffine(4,6,2,alignment_scope,memory_model) {}

  wf_heuristic_strategy heuristic_strategy() const {
    return wfAligner->heuristic.strategy;
  }
};

int main() {
  InspectableGapAffine high(wfa::WFAligner::Alignment,wfa::WFAligner::MemoryHigh);
  InspectableGapAffine singletrack(wfa::WFAligner::Alignment,wfa::WFAligner::MemorySingletrack);

  if (high.heuristic_strategy() != wf_heuristic_wfadaptive) {
    std::fprintf(stderr,"high-memory default heuristic was not preserved\n");
    return 1;
  }
  if (singletrack.heuristic_strategy() != high.heuristic_strategy()) {
    std::fprintf(stderr,
        "Singletrack default heuristic mismatch: observed=%lu expected=%lu\n",
        (unsigned long)singletrack.heuristic_strategy(),
        (unsigned long)high.heuristic_strategy());
    return 1;
  }
  return 0;
}
