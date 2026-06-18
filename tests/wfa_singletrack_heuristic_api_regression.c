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
typedef void (*configure_span_fn)(wavefront_aligner_t* const);

static const char* const default_pattern =
    "CACGATCAGGAAGCTGCATCCTTATCTAGGGTGTGTGATCGACAGCTTGGCTGCTCGTCCCACATATGAGTTGCAAGTAACGAATACGTTTCTACGAACC";
static const char* const default_text =
    "CACGATCAGGAAGCATGCATCCTTATCTAGGGTGTGTGATCGACAGCTTGGCTGCTCGTCCCACATATGAGTTGCAAGTAACGAATACGTTTCTACGAACC";
static const char* const endsfree_pattern =
    "AAAACCCCGGGG";
static const char* const endsfree_text =
    "TTTTAAAACCCCGGGGAAAA";
static const char* const extension_pattern =
    "AAAACCCCGGGGTTTT";
static const char* const extension_text =
    "AAAACCCCGGGGAAAA";
static const char* const zdrop_pattern =
    "CACGATCAGGAAGCTGCATCCTTATCTAGGGTGTGTGATCGACAGCTTGGCTGCTCGTCCCACATATGAGTTGCAAGTAACGAATACGTTTCTACGAACC";
static const char* const zdrop_text =
    "CTAGAATCAGAAGCATGCCTCCTATCCAGGGTGTGTGATCACAGCTCGCTGCTCGTCCCCATAATGAGTTCACAGTACGAATAAATTTCTACAAACC";

static wavefront_aligner_t* new_aligner(
    const wavefront_memory_t memory_mode,
    const distance_metric_t distance_metric) {
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = distance_metric;
  attributes.memory_mode = memory_mode;
  attributes.alignment_scope = compute_alignment;
  attributes.heuristic.strategy = wf_heuristic_none;
  attributes.system.check_alignment_correct = true;
  return wavefront_aligner_new(&attributes);
}

static int check_positive_penalty_invariants(
    const char* const label,
    wavefront_aligner_attr_t attributes) {
  attributes.memory_mode = wavefront_memory_singletrack;
  attributes.alignment_scope = compute_alignment;
  attributes.heuristic.strategy = wf_heuristic_none;

  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  int failed = 0;

  if (penalties->mismatch <= 0) {
    fprintf(stderr,"%s internal mismatch is not positive: %d\n",
        label,penalties->mismatch);
    failed = 1;
  }
  if (penalties->gap_extension1 <= 0) {
    fprintf(stderr,"%s internal gap-extension1 is not positive: %d\n",
        label,penalties->gap_extension1);
    failed = 1;
  }
  if (attributes.distance_metric == gap_affine_2p &&
      penalties->gap_extension2 <= 0) {
    fprintf(stderr,"%s internal gap-extension2 is not positive: %d\n",
        label,penalties->gap_extension2);
    failed = 1;
  }

  wavefront_aligner_delete(wf_aligner);
  return failed;
}

static int check_singletrack_penalty_invariants(void) {
  int failed = 0;
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;

  attributes.distance_metric = gap_affine;
  attributes.affine_penalties.match = 0;
  attributes.affine_penalties.mismatch = 1;
  attributes.affine_penalties.gap_opening = 0;
  attributes.affine_penalties.gap_extension = 1;
  failed |= check_positive_penalty_invariants("affine",attributes);

  attributes.affine_penalties.match = -5;
  attributes.affine_penalties.mismatch = 1;
  attributes.affine_penalties.gap_opening = 0;
  attributes.affine_penalties.gap_extension = 1;
  failed |= check_positive_penalty_invariants("affine-negative-match",attributes);

  attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = gap_affine_2p;
  attributes.affine2p_penalties.match = 0;
  attributes.affine2p_penalties.mismatch = 1;
  attributes.affine2p_penalties.gap_opening1 = 0;
  attributes.affine2p_penalties.gap_extension1 = 1;
  attributes.affine2p_penalties.gap_opening2 = 0;
  attributes.affine2p_penalties.gap_extension2 = 1;
  failed |= check_positive_penalty_invariants("affine2p",attributes);

  attributes.affine2p_penalties.match = -5;
  attributes.affine2p_penalties.mismatch = 1;
  attributes.affine2p_penalties.gap_opening1 = 0;
  attributes.affine2p_penalties.gap_extension1 = 1;
  attributes.affine2p_penalties.gap_opening2 = 0;
  attributes.affine2p_penalties.gap_extension2 = 1;
  failed |= check_positive_penalty_invariants(
      "affine2p-negative-match",attributes);

  return failed;
}

static void configure_end_to_end(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_alignment_end_to_end(wf_aligner);
}

static void configure_endsfree(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_alignment_free_ends(wf_aligner,4,4,4,4);
}

static void configure_extension(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_alignment_extension(wf_aligner);
}

static int check_heuristic(
    const char* const label,
    wavefront_aligner_t* const wf_aligner,
    const wf_heuristic_strategy expected_strategy,
    const int expected_min_wavefront_length,
    const int expected_max_distance_threshold,
    const int expected_steps_between_cutoffs) {
  if (wf_aligner->heuristic.strategy != expected_strategy) {
    fprintf(stderr,"%s strategy mismatch: observed=%lu expected=%lu\n",
        label,
        (unsigned long)wf_aligner->heuristic.strategy,
        (unsigned long)expected_strategy);
    return 1;
  }
  if (wf_aligner->heuristic.min_wavefront_length != expected_min_wavefront_length ||
      wf_aligner->heuristic.max_distance_threshold != expected_max_distance_threshold ||
      wf_aligner->heuristic.steps_between_cutoffs != expected_steps_between_cutoffs) {
    fprintf(stderr,
        "%s parameter mismatch: observed=(%d,%d,%d) expected=(%d,%d,%d)\n",
        label,
        wf_aligner->heuristic.min_wavefront_length,
        wf_aligner->heuristic.max_distance_threshold,
        wf_aligner->heuristic.steps_between_cutoffs,
        expected_min_wavefront_length,
        expected_max_distance_threshold,
        expected_steps_between_cutoffs);
    return 1;
  }
  return 0;
}

static void configure_wfadaptive(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_heuristic_wfadaptive(wf_aligner,10,50,1);
}

static int check_wfadaptive(
    const char* const label,
    wavefront_aligner_t* const wf_aligner) {
  return check_heuristic(label,wf_aligner,wf_heuristic_wfadaptive,10,50,1);
}

static void configure_wfmash(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_heuristic_wfmash(wf_aligner,10,50,1);
}

static int check_wfmash(
    const char* const label,
    wavefront_aligner_t* const wf_aligner) {
  return check_heuristic(label,wf_aligner,wf_heuristic_wfmash,10,50,1);
}

static void configure_combined_wfa_wfmash(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_heuristic_wfadaptive(wf_aligner,8,40,3);
  wavefront_aligner_set_heuristic_wfmash(wf_aligner,10,50,1);
}

static int check_combined_wfa_wfmash(
    const char* const label,
    wavefront_aligner_t* const wf_aligner) {
  return check_heuristic(
      label,wf_aligner,wf_heuristic_wfadaptive|wf_heuristic_wfmash,10,50,1);
}

static void configure_xdrop(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_heuristic_xdrop(wf_aligner,100,1);
}

static int check_xdrop(
    const char* const label,
    wavefront_aligner_t* const wf_aligner) {
  if (wf_aligner->heuristic.strategy != wf_heuristic_xdrop) {
    fprintf(stderr,"%s xdrop strategy mismatch: observed=%lu expected=%lu\n",
        label,
        (unsigned long)wf_aligner->heuristic.strategy,
        (unsigned long)wf_heuristic_xdrop);
    return 1;
  }
  if (wf_aligner->heuristic.xdrop != 100 ||
      wf_aligner->heuristic.steps_between_cutoffs != 1) {
    fprintf(stderr,"%s xdrop parameter mismatch: observed=(%d,%d) expected=(100,1)\n",
        label,
        wf_aligner->heuristic.xdrop,
        wf_aligner->heuristic.steps_between_cutoffs);
    return 1;
  }
  return 0;
}

static void configure_zdrop(
    wavefront_aligner_t* const wf_aligner) {
  wavefront_aligner_set_heuristic_zdrop(wf_aligner,10,1);
}

static int check_zdrop(
    const char* const label,
    wavefront_aligner_t* const wf_aligner) {
  if (wf_aligner->heuristic.strategy != wf_heuristic_zdrop) {
    fprintf(stderr,"%s zdrop strategy mismatch: observed=%lu expected=%lu\n",
        label,
        (unsigned long)wf_aligner->heuristic.strategy,
        (unsigned long)wf_heuristic_zdrop);
    return 1;
  }
  if (wf_aligner->heuristic.zdrop != 10 ||
      wf_aligner->heuristic.steps_between_cutoffs != 1) {
    fprintf(stderr,"%s zdrop parameter mismatch: observed=(%d,%d) expected=(10,1)\n",
        label,
        wf_aligner->heuristic.zdrop,
        wf_aligner->heuristic.steps_between_cutoffs);
    return 1;
  }
  return 0;
}

static int check_cigar_coherent(
    const char* const label,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const cigar_t* const cigar,
    const bool require_full_length,
    const bool check_stored_end_position) {
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
  if (check_stored_end_position &&
      (pattern_pos != cigar->end_v || text_pos != cigar->end_h)) {
    fprintf(stderr,
        "%s CIGAR endpoint mismatch: traced=(%d,%d) stored=(%d,%d)\n",
        label,pattern_pos,text_pos,cigar->end_v,cigar->end_h);
    return 1;
  }
  if (require_full_length &&
      (pattern_pos != pattern_length || text_pos != text_length)) {
    fprintf(stderr,
        "%s CIGAR length mismatch: traced=(%d,%d) expected=(%d,%d)\n",
        label,pattern_pos,text_pos,pattern_length,text_length);
    return 1;
  }
  return 0;
}

static int check_alignment_result(
    const char* const label,
    const wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const char* const text,
    const bool check_stored_end_position) {
  if (wf_aligner->align_status.status != WF_STATUS_ALG_COMPLETED &&
      wf_aligner->align_status.status != WF_STATUS_ALG_PARTIAL) {
    return 0;
  }
  if (cigar_is_null(wf_aligner->cigar)) {
    fprintf(stderr,"%s produced a null CIGAR\n",label);
    return 1;
  }
  const int failed = check_cigar_coherent(
      label,
      pattern,(int)strlen(pattern),text,(int)strlen(text),
      wf_aligner->cigar,
      wf_aligner->align_status.status == WF_STATUS_ALG_COMPLETED,
      check_stored_end_position);
  if (failed) {
    fprintf(stderr,"%s produced an incoherent CIGAR\n",label);
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

static int check_case_sequences(
    const char* const name,
    const distance_metric_t distance_metric,
    configure_heuristic_fn configure_heuristic,
    int (*check_config)(const char* const,wavefront_aligner_t* const),
    const char* const pattern,
    const char* const text,
    const bool expect_dropped,
    configure_span_fn configure_span,
    const bool check_end_position) {
  wavefront_aligner_t* const high = new_aligner(wavefront_memory_high,distance_metric);
  wavefront_aligner_t* const singletrack = new_aligner(wavefront_memory_singletrack,distance_metric);
  configure_span(high);
  configure_span(singletrack);
  configure_heuristic(high);
  configure_heuristic(singletrack);

  int failed = 0;
  failed |= check_config("high",high);
  failed |= check_config("singletrack",singletrack);

  const int high_status = wavefront_align(
      high,pattern,(int)strlen(pattern),text,(int)strlen(text));
  const int singletrack_status = wavefront_align(
      singletrack,pattern,(int)strlen(pattern),text,(int)strlen(text));

  if (high_status != singletrack_status ||
      high->align_status.status != singletrack->align_status.status) {
    fprintf(stderr,
        "%s status mismatch: high return/status=%d/%d singletrack return/status=%d/%d\n",
        name,high_status,high->align_status.status,
        singletrack_status,singletrack->align_status.status);
    failed = 1;
  }
  if (high->align_status.score != singletrack->align_status.score) {
    fprintf(stderr,"%s align-status score mismatch: high=%d singletrack=%d\n",
        name,high->align_status.score,singletrack->align_status.score);
    failed = 1;
  }
  if (high->align_status.dropped != singletrack->align_status.dropped) {
    fprintf(stderr,"%s dropped mismatch: high=%d singletrack=%d\n",
        name,high->align_status.dropped,singletrack->align_status.dropped);
    failed = 1;
  }
  if (expect_dropped &&
      (!high->align_status.dropped || !singletrack->align_status.dropped)) {
    fprintf(stderr,"%s did not exercise dropped alignment path: high=%d singletrack=%d\n",
        name,high->align_status.dropped,singletrack->align_status.dropped);
    failed = 1;
  }
  if (high->cigar->score != singletrack->cigar->score) {
    fprintf(stderr,"%s score mismatch: high=%d singletrack=%d\n",
        name,high->cigar->score,singletrack->cigar->score);
    failed = 1;
  }
  const int high_cigar_score = compute_cigar_score(high);
  const int singletrack_cigar_score = compute_cigar_score(singletrack);
  if (high_cigar_score != singletrack_cigar_score) {
    fprintf(stderr,"%s CIGAR score mismatch: high=%d singletrack=%d\n",
        name,high_cigar_score,singletrack_cigar_score);
    failed = 1;
  }
  if (check_end_position &&
      high->align_status.status == WF_STATUS_ALG_COMPLETED &&
      (high_cigar_score != high->cigar->score ||
       singletrack_cigar_score != singletrack->cigar->score)) {
    fprintf(stderr,
        "%s stored/CIGAR score mismatch: high=%d/%d singletrack=%d/%d\n",
        name,high->cigar->score,high_cigar_score,
        singletrack->cigar->score,singletrack_cigar_score);
    failed = 1;
  }
  if (check_end_position &&
      (high->cigar->end_v != singletrack->cigar->end_v ||
       high->cigar->end_h != singletrack->cigar->end_h)) {
    fprintf(stderr,"%s end-position mismatch: high=(%d,%d) singletrack=(%d,%d)\n",
        name,high->cigar->end_v,high->cigar->end_h,
        singletrack->cigar->end_v,singletrack->cigar->end_h);
    failed = 1;
  }
  failed |= check_alignment_result("high",high,pattern,text,check_end_position);
  failed |= check_alignment_result(
      "singletrack",singletrack,pattern,text,check_end_position);

  wavefront_aligner_delete(high);
  wavefront_aligner_delete(singletrack);
  return failed;
}

static int check_case(
    const char* const name,
    const distance_metric_t distance_metric,
    configure_heuristic_fn configure_heuristic,
    int (*check_config)(const char* const,wavefront_aligner_t* const)) {
  return check_case_sequences(
      name,distance_metric,configure_heuristic,check_config,
      default_pattern,default_text,false,configure_end_to_end,true);
}

static int check_span_case(
    const char* const name,
    const distance_metric_t distance_metric,
    configure_heuristic_fn configure_heuristic,
    int (*check_config)(const char* const,wavefront_aligner_t* const),
    configure_span_fn configure_span,
    const char* const pattern,
    const char* const text) {
  return check_case_sequences(
      name,distance_metric,configure_heuristic,check_config,
      pattern,text,false,configure_span,false);
}

static int check_span_suite(
    const char* const span_name,
    configure_span_fn configure_span,
    const char* const pattern,
    const char* const text) {
  char name[128];
  int failed = 0;

#define CHECK_SPAN_CASE(label,distance,configure,check) \
  snprintf(name,sizeof(name),"%s-%s",span_name,label); \
  failed |= check_span_case(name,distance,configure,check, \
      configure_span,pattern,text)

  CHECK_SPAN_CASE(
      "affine-wfadaptive",gap_affine,
      configure_wfadaptive,check_wfadaptive);
  CHECK_SPAN_CASE(
      "affine-wfmash",gap_affine,
      configure_wfmash,check_wfmash);
  CHECK_SPAN_CASE(
      "affine-xdrop",gap_affine,
      configure_xdrop,check_xdrop);
  CHECK_SPAN_CASE(
      "affine-zdrop",gap_affine,
      configure_zdrop,check_zdrop);
  CHECK_SPAN_CASE(
      "affine2p-wfadaptive",gap_affine_2p,
      configure_wfadaptive,check_wfadaptive);
  CHECK_SPAN_CASE(
      "affine2p-wfmash",gap_affine_2p,
      configure_wfmash,check_wfmash);
  CHECK_SPAN_CASE(
      "affine2p-xdrop",gap_affine_2p,
      configure_xdrop,check_xdrop);
  CHECK_SPAN_CASE(
      "affine2p-zdrop",gap_affine_2p,
      configure_zdrop,check_zdrop);

#undef CHECK_SPAN_CASE

  return failed;
}

int main(void) {
  int failed = 0;
  failed |= check_singletrack_penalty_invariants();
  failed |= check_case("wfmash",gap_affine,configure_wfmash,check_wfmash);
  failed |= check_case(
      "combined-wfa-wfmash",gap_affine,
      configure_combined_wfa_wfmash,check_combined_wfa_wfmash);
  failed |= check_case("xdrop",gap_affine,configure_xdrop,check_xdrop);
  failed |= check_case_sequences(
      "zdrop-partial",gap_affine,configure_zdrop,check_zdrop,
      zdrop_pattern,zdrop_text,true,configure_end_to_end,true);
  failed |= check_span_suite(
      "endsfree",configure_endsfree,endsfree_pattern,endsfree_text);
  failed |= check_span_suite(
      "extension",configure_extension,extension_pattern,extension_text);
  return failed;
}
