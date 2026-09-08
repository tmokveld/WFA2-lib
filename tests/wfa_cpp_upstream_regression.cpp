/*
 * Wavefront Alignment Algorithms
 * Regression coverage for upstream issues #126 and #129.
 */
#include <cstdio>
#include <cstring>
#include <string>

#include "bindings/cpp/WFAligner.hpp"

static int check_cigar_capacity() {
  wfa::WFAlignerEdit aligner(wfa::WFAligner::Alignment);
  if (aligner.alignEnd2End("A",1,"A",1) != wfa::WFAligner::StatusAlgCompleted ||
      aligner.getCIGAR(true) != "1=") {
    std::fprintf(stderr,"Single-operation SAM CIGAR failed\n");
    return 1;
  }
  return 0;
}

static int check_alignment_reuse() {
  wfa::WFAlignerEdit aligner(wfa::WFAligner::Alignment);
  if (aligner.alignEndsFree("T",1,0,0,"AAAAT",5,4,0) !=
      wfa::WFAligner::StatusAlgCompleted || aligner.getAlignmentScore() != 0) {
    std::fprintf(stderr,"Ends-free alignment failed\n");
    return 1;
  }
  if (aligner.alignEnd2End("AAAAA",5,"TTTTT",5) !=
      wfa::WFAligner::StatusAlgCompleted || aligner.getAlignmentScore() != 5 ||
      aligner.getAlignment() != "XXXXX") {
    std::fprintf(stderr,"End-to-end alignment inherited free ends: score=%d\n",
        aligner.getAlignmentScore());
    return 1;
  }
  return 0;
}

static int check_biwfa_reuse() {
  const std::string pattern(512,'A'), text(512,'T');
  for (auto scope : {wfa::WFAligner::Alignment, wfa::WFAligner::Score}) {
    for (bool free_pattern : {false,true}) {
      wfa::WFAlignerGapAffine aligner(4,6,2,scope,wfa::WFAligner::MemoryUltralow);
      const auto initial_status = free_pattern ?
          aligner.alignEndsFree("AAAATGGGG",9,4,4,"T",1,0,0) :
          aligner.alignEndsFree("T",1,0,0,"AAAATGGGG",9,4,4);
      if (initial_status != wfa::WFAligner::StatusAlgCompleted ||
          aligner.getAlignmentScore() != 0) return 1;
      if (aligner.alignEnd2End(pattern,text) != wfa::WFAligner::StatusAlgCompleted ||
          aligner.getAlignmentScore() != -2048) {
        std::fprintf(stderr,"BiWFA inherited free ends: scope=%d score=%d\n",
            scope,aligner.getAlignmentScore());
        return 1;
      }
    }
  }
  return 0;
}

int main(int argc, char** argv) {
  if (argc == 2 && std::strcmp(argv[1],"cigar") == 0) return check_cigar_capacity();
  if (argc == 2 && std::strcmp(argv[1],"reuse") == 0) return check_alignment_reuse();
  return check_cigar_capacity() | check_alignment_reuse() | check_biwfa_reuse();
}
