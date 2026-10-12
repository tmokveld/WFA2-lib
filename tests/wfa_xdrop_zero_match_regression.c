/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 */

#include <stdbool.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "wavefront/wavefront_align.h"

typedef struct {
  int status;
  int score;
  int end_v;
  int end_h;
  bool cigar_valid;
} xdrop_result_t;

static const char pattern[] =
    "ATGAGATTCGCCCCCTGACTTGGGACCAGAAGTCTTATACAAGTACGGAGATGGTAAGTTTGCT"
    "GTTCGACTGAACGCCTCCTCGATTAGATCCTCGTCCGGTCG";
static const char text[] =
    "ATGAGATCGCCCCCTGACTTGGGACCAGAAGTCTTATACAAGTACGGAGATGGTAAGTTTGCTG"
    "TTCGACTGAACGCTTCCTCGATTAGATCCTCGTCCTACAC";

static xdrop_result_t align_xdrop(
    const distance_metric_t distance_metric,
    const int match,
    const int xdrop) {
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = distance_metric;
  attributes.memory_mode = wavefront_memory_high;
  attributes.alignment_scope = compute_alignment;
  attributes.heuristic.strategy = wf_heuristic_none;
  attributes.system.check_alignment_correct = true;
  if (distance_metric == gap_affine) {
    attributes.affine_penalties.match = match;
    attributes.affine_penalties.mismatch = 4;
    attributes.affine_penalties.gap_opening = 6;
    attributes.affine_penalties.gap_extension = 2;
  } else {
    attributes.affine2p_penalties.match = match;
    attributes.affine2p_penalties.mismatch = 4;
    attributes.affine2p_penalties.gap_opening1 = 6;
    attributes.affine2p_penalties.gap_extension1 = 2;
    attributes.affine2p_penalties.gap_opening2 = 24;
    attributes.affine2p_penalties.gap_extension2 = 1;
  }
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  wavefront_aligner_set_alignment_extension(wf_aligner);
  if (xdrop >= 0) {
    wavefront_aligner_set_heuristic_xdrop(wf_aligner,xdrop,1);
  }
  xdrop_result_t result;
  result.status = wavefront_align(
      wf_aligner,pattern,(int)strlen(pattern),text,(int)strlen(text));
  result.end_v = wf_aligner->cigar->end_v;
  result.end_h = wf_aligner->cigar->end_h;
  result.cigar_valid = !cigar_is_null(wf_aligner->cigar);
  if (!result.cigar_valid) {
    result.score = INT_MIN;
  } else if (distance_metric == gap_affine) {
    result.score = cigar_score_gap_affine(
        wf_aligner->cigar,&wf_aligner->penalties.affine_penalties);
  } else {
    result.score = cigar_score_gap_affine2p(
        wf_aligner->cigar,&wf_aligner->penalties.affine2p_penalties);
  }
  result.cigar_valid = result.cigar_valid && cigar_check_alignment(
          stderr,pattern,result.end_v,text,result.end_h,
          wf_aligner->cigar,true);
  wavefront_aligner_delete(wf_aligner);
  return result;
}

static int check_xdrop_matches_baseline(
    const char* const label,
    const distance_metric_t distance_metric,
    const int match,
    const int xdrop_threshold) {
  const xdrop_result_t baseline = align_xdrop(distance_metric,match,-1);
  const xdrop_result_t xdrop =
      align_xdrop(distance_metric,match,xdrop_threshold);
  const bool baseline_ok =
      (baseline.status == WF_STATUS_ALG_COMPLETED ||
       baseline.status == WF_STATUS_ALG_PARTIAL) &&
      baseline.cigar_valid;
  const bool xdrop_ok =
      xdrop.status == baseline.status &&
      xdrop.score == baseline.score &&
      xdrop.end_v == baseline.end_v &&
      xdrop.end_h == baseline.end_h &&
      xdrop.cigar_valid;
  if (!baseline_ok || !xdrop_ok) {
    fprintf(stderr,
        "%s X-drop mismatch: baseline=(%d,%d,%d,%d,%d) "
        "xdrop=(%d,%d,%d,%d,%d)\n",
        label,
        baseline.status,baseline.score,baseline.end_v,baseline.end_h,
        baseline.cigar_valid,
        xdrop.status,xdrop.score,xdrop.end_v,xdrop.end_h,xdrop.cigar_valid);
  }
  return !baseline_ok || !xdrop_ok;
}

int main(void) {
  int failed = 0;
  failed |= check_xdrop_matches_baseline(
      "M=0 affine",gap_affine,0,10);
  failed |= check_xdrop_matches_baseline(
      "M=0 affine-2p",gap_affine_2p,0,10);
  failed |= check_xdrop_matches_baseline(
      "M=-1 affine",gap_affine,-1,100);
  failed |= check_xdrop_matches_baseline(
      "M=-1 affine-2p",gap_affine_2p,-1,100);
  return failed;
}
