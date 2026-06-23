/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 */

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "wavefront/wavefront_align.h"

typedef void (*configure_heuristic_fn)(wavefront_aligner_t* const);

static const char* const gap_pattern =
    "ACGTACGTGATACCGT"
    "TTTTGGGGCCCCAA"
    "GATTACAGGCTAATCG"
    "CGTATGCAACGTAGCT";
static const char* const gap_text =
    "ACGTACGTGATACCGT"
    "GATTACAGGCTAATCG"
    "AAAACCCCGGGGTT"
    "CGTATGCAACGTAGCT";

static void configure_wfadaptive(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_heuristic_wfadaptive(wf_aligner,10,50,1);
}

static void configure_banded_static(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_heuristic_banded_static(wf_aligner,-16,16);
}

static int check_cigar_coherent(
    const char* const label,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const cigar_t* const cigar) {
  int pattern_pos = 0, text_pos = 0, i;
  for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
    switch (cigar->operations[i]) {
      case 'M':
        if (pattern_pos >= pattern_length || text_pos >= text_length ||
            pattern[pattern_pos] != text[text_pos]) {
          fprintf(stderr,"%s CIGAR match incoherent at (%d,%d)\n",
              label,pattern_pos,text_pos);
          return 1;
        }
        ++pattern_pos;
        ++text_pos;
        break;
      case 'X':
        if (pattern_pos >= pattern_length || text_pos >= text_length ||
            pattern[pattern_pos] == text[text_pos]) {
          fprintf(stderr,"%s CIGAR mismatch incoherent at (%d,%d)\n",
              label,pattern_pos,text_pos);
          return 1;
        }
        ++pattern_pos;
        ++text_pos;
        break;
      case 'I':
        if (text_pos >= text_length) {
          fprintf(stderr,"%s CIGAR insertion exceeds text length at %d\n",
              label,text_pos);
          return 1;
        }
        ++text_pos;
        break;
      case 'D':
        if (pattern_pos >= pattern_length) {
          fprintf(stderr,"%s CIGAR deletion exceeds pattern length at %d\n",
              label,pattern_pos);
          return 1;
        }
        ++pattern_pos;
        break;
      default:
        fprintf(stderr,"%s CIGAR has unknown operation '%c'\n",
            label,cigar->operations[i]);
        return 1;
    }
  }
  if (pattern_pos != pattern_length || text_pos != text_length) {
    fprintf(stderr,
        "%s CIGAR length mismatch: traced=(%d,%d) expected=(%d,%d)\n",
        label,pattern_pos,text_pos,pattern_length,text_length);
    return 1;
  }
  if (pattern_pos != cigar->end_v || text_pos != cigar->end_h) {
    fprintf(stderr,
        "%s CIGAR endpoint mismatch: traced=(%d,%d) stored=(%d,%d)\n",
        label,pattern_pos,text_pos,cigar->end_v,cigar->end_h);
    return 1;
  }
  return 0;
}

static int compute_cigar_score(
    const wavefront_aligner_t* const wf_aligner) {
  switch (wf_aligner->penalties.distance_metric) {
    case gap_affine:
      return cigar_score_gap_affine(
          wf_aligner->cigar,&wf_aligner->penalties.affine_penalties);
    case gap_affine_2p:
      return cigar_score_gap_affine2p(
          wf_aligner->cigar,&wf_aligner->penalties.affine2p_penalties);
    default:
      return wf_aligner->cigar->score;
  }
}

static int check_case(
    const char* const name,
    const distance_metric_t distance_metric,
    configure_heuristic_fn configure_heuristic) {
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = distance_metric;
  attributes.memory_mode = wavefront_memory_high;
  attributes.alignment_scope = compute_alignment;
  attributes.heuristic.strategy = wf_heuristic_none;
  attributes.system.check_alignment_correct = true;

  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  configure_heuristic(wf_aligner);
  const int pattern_length = (int)strlen(gap_pattern);
  const int text_length = (int)strlen(gap_text);
  const int status = wavefront_align(
      wf_aligner,gap_pattern,pattern_length,gap_text,text_length);

  int failed = 0;
  if (status != WF_STATUS_ALG_COMPLETED ||
      wf_aligner->align_status.status != WF_STATUS_ALG_COMPLETED) {
    fprintf(stderr,"%s status mismatch: return/status=%d/%d\n",
        name,status,wf_aligner->align_status.status);
    failed = 1;
  } else if (cigar_is_null(wf_aligner->cigar)) {
    fprintf(stderr,"%s produced a null CIGAR\n",name);
    failed = 1;
  } else {
    failed |= check_cigar_coherent(
        name,gap_pattern,pattern_length,gap_text,text_length,
        wf_aligner->cigar);
    const int cigar_score = compute_cigar_score(wf_aligner);
    if (cigar_score != wf_aligner->cigar->score) {
      fprintf(stderr,"%s stored/CIGAR score mismatch: %d/%d\n",
          name,wf_aligner->cigar->score,cigar_score);
      failed = 1;
    }
  }

  wavefront_aligner_delete(wf_aligner);
  return failed;
}

int main(void) {
  int failed = 0;
  failed |= check_case("affine-wfadaptive",gap_affine,configure_wfadaptive);
  failed |= check_case("affine-banded-static",gap_affine,configure_banded_static);
  failed |= check_case("affine2p-wfadaptive",gap_affine_2p,configure_wfadaptive);
  failed |= check_case(
      "affine2p-banded-static",gap_affine_2p,configure_banded_static);
  return failed;
}
