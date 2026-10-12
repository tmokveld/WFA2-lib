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

static xdrop_result_t align_xdrop_mode(
    const distance_metric_t distance_metric,
    const wavefront_memory_t memory_mode,
    const int match,
    const int xdrop) {
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = distance_metric;
  attributes.memory_mode = memory_mode;
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

static xdrop_result_t align_xdrop(
    const distance_metric_t distance_metric,
    const int match,
    const int xdrop) {
  return align_xdrop_mode(distance_metric,wavefront_memory_high,match,xdrop);
}

static int check_xdrop_exhausted(
    const char* const label,
    const distance_metric_t distance_metric,
    const wavefront_memory_t memory_mode,
    const int match,
    const int xdrop_threshold) {
  const xdrop_result_t xdrop =
      align_xdrop_mode(distance_metric,memory_mode,match,xdrop_threshold);
  const bool xdrop_ok =
      xdrop.status == WF_STATUS_ALG_PARTIAL &&
      xdrop.cigar_valid &&
      xdrop.end_v > 0 && xdrop.end_h > 0;
  if (!xdrop_ok) {
    fprintf(stderr,
        "%s (memory mode %d) X-drop exhausted every diagonal without an alignment: "
        "(%d,%d,%d,%d,%d)\n",
        label,memory_mode,
        xdrop.status,xdrop.score,xdrop.end_v,xdrop.end_h,xdrop.cigar_valid);
  }
  return !xdrop_ok;
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
  const wavefront_memory_t memory_modes[] = {
      wavefront_memory_high,wavefront_memory_med,
      wavefront_memory_low,wavefront_memory_singletrack};
  int i;
  for (i=0;i<4;++i) {
    failed |= check_xdrop_exhausted(
        "M=0 affine",gap_affine,memory_modes[i],0,1);
    failed |= check_xdrop_exhausted(
        "M=-1 affine",gap_affine,memory_modes[i],-1,10);
    failed |= check_xdrop_exhausted(
        "M=-1 affine-2p",gap_affine_2p,memory_modes[i],-1,10);
  }
  return failed;
}
