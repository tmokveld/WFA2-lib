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
#include "wavefront/wavefront_unialign.h"

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

static const char* const endsfree_pattern =
    "AGATAAACGCCCTTGCCGCTTTCACGCATCGTAATCCCTGAGCCCTTTCTTAATTACATGAATGTCACCCCTCGTGAATCCGTATTCTTAGTACCTGTA"
    "TCCGAATGGAATCTACACCGGGCATCGGTACCGATCAGTTTCTCGGCGTTCCT";
static const char* const endsfree_text =
    "CCTGACGCGATCTAGGTAGAACATGGGATTTGGGATGAATTGGCCCCGACAGCTAACCCGGGGCCAGTTCTTTATGGCGTCTAACGAATTGAGTTTGCAA"
    "TAAAATCAGATTAGACGGGCAAAAGCTCGTCTAAAAGAAAAAGGACCTTATGTTGCCTCGGTCCCAAGCATTCTCTTCTGTCGGCCTCGGGGAAAATCTC"
    "GCGCGATAATATGCTGCC";

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
    const cigar_t* const cigar,
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
  if (pattern_pos != pattern_length || text_pos != text_length) {
    fprintf(stderr,
        "%s CIGAR length mismatch: traced=(%d,%d) expected=(%d,%d)\n",
        label,pattern_pos,text_pos,pattern_length,text_length);
    return 1;
  }
  if (check_stored_end_position &&
      (pattern_pos != cigar->end_v || text_pos != cigar->end_h)) {
    fprintf(stderr,
        "%s CIGAR endpoint mismatch: traced=(%d,%d) stored=(%d,%d)\n",
        label,pattern_pos,text_pos,cigar->end_v,cigar->end_h);
    return 1;
  }
  return 0;
}

static int compute_cigar_score(
    const wavefront_aligner_t* const wf_aligner,
    const cigar_t* const cigar) {
  switch (wf_aligner->penalties.distance_metric) {
    case gap_affine:
      return cigar_score_gap_affine(
          cigar,&wf_aligner->penalties.affine_penalties);
    case gap_affine_2p:
      return cigar_score_gap_affine2p(
          cigar,&wf_aligner->penalties.affine2p_penalties);
    default:
      return cigar->score;
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
        wf_aligner->cigar,true);
    const int cigar_score = compute_cigar_score(wf_aligner,wf_aligner->cigar);
    if (cigar_score != wf_aligner->cigar->score) {
      fprintf(stderr,"%s stored/CIGAR score mismatch: %d/%d\n",
          name,wf_aligner->cigar->score,cigar_score);
      failed = 1;
    }
  }

  wavefront_aligner_delete(wf_aligner);
  return failed;
}

static int check_endsfree_result(
    const char* const label,
    const wavefront_aligner_t* const wf_aligner,
    const int status,
    const bool transpose) {
  if (status != WF_STATUS_ALG_COMPLETED ||
      wf_aligner->align_status.status != WF_STATUS_ALG_COMPLETED ||
      wf_aligner->align_status.dropped) {
    fprintf(stderr,"%s status mismatch: return/status/dropped=%d/%d/%d\n",
        label,status,wf_aligner->align_status.status,
        wf_aligner->align_status.dropped);
    return 1;
  }
  const cigar_t* const cigar = wf_aligner->cigar;
  if (cigar_is_null(cigar)) {
    fprintf(stderr,"%s produced a null CIGAR\n",label);
    return 1;
  }
  const char* const pattern = transpose ? endsfree_text : endsfree_pattern;
  const char* const text = transpose ? endsfree_pattern : endsfree_text;
  const int pattern_length = (int)strlen(pattern);
  const int text_length = (int)strlen(text);
  if (check_cigar_coherent(
          label,pattern,pattern_length,text,text_length,cigar,false)) {
    return 1;
  }
  if ((!transpose &&
       (cigar->end_v != 152 || cigar->end_h < 0 || cigar->end_h > 218)) ||
      (transpose &&
       (cigar->end_h != 152 || cigar->end_v < 0 || cigar->end_v > 218))) {
    fprintf(stderr,"%s invalid ends-free endpoint (%d,%d)\n",
        label,cigar->end_v,cigar->end_h);
    return 1;
  }
  int v = 0, h = 0, endpoint_index = -1, i;
  for (i=cigar->begin_offset;i<=cigar->end_offset;++i) {
    if (v == cigar->end_v && h == cigar->end_h) {
      endpoint_index = i;
      break;
    }
    if (i == cigar->end_offset) break;
    const char op = cigar->operations[i];
    if (op != 'I') ++v;
    if (op != 'D') ++h;
  }
  if (endpoint_index < 0) {
    fprintf(stderr,"%s endpoint (%d,%d) is not on the CIGAR path\n",
        label,cigar->end_v,cigar->end_h);
    return 1;
  }
  const char suffix_op = transpose ? 'D' : 'I';
  for (i=endpoint_index;i<cigar->end_offset;++i) {
    if (cigar->operations[i] != suffix_op) {
      fprintf(stderr,"%s free suffix has operation '%c', expected '%c'\n",
          label,cigar->operations[i],suffix_op);
      return 1;
    }
  }
  const int suffix_length = transpose ?
      pattern_length - cigar->end_v : text_length - cigar->end_h;
  if (cigar->end_offset - endpoint_index != suffix_length) {
    fprintf(stderr,"%s free suffix does not consume the remaining sequence\n",
        label);
    return 1;
  }
  cigar_t prefix = *cigar;
  prefix.end_offset = endpoint_index;
  const int prefix_score = compute_cigar_score(wf_aligner,&prefix);
  if (prefix_score != cigar->score) {
    fprintf(stderr,"%s stored/prefix score mismatch: %d/%d\n",
        label,cigar->score,prefix_score);
    return 1;
  }
  if (wf_aligner->penalties.distance_metric == gap_affine &&
      wf_aligner->heuristic.strategy == wf_heuristic_none &&
      cigar->score != -298) {
    fprintf(stderr,"%s exact score mismatch: observed=%d expected=-298\n",
        label,cigar->score);
    return 1;
  }
  return 0;
}

static int check_endsfree_case(
    const distance_metric_t distance_metric,
    const bool wfadaptive,
    const bool transpose) {
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = distance_metric;
  attributes.alignment_scope = compute_alignment;
  attributes.alignment_form.span = alignment_endsfree;
  attributes.alignment_form.pattern_begin_free = 0;
  attributes.alignment_form.text_begin_free = 0;
  attributes.alignment_form.pattern_end_free = transpose ? 218 : 0;
  attributes.alignment_form.text_end_free = transpose ? 0 : 218;
  if (!wfadaptive) attributes.heuristic.strategy = wf_heuristic_none;
  attributes.memory_mode = wavefront_memory_high;
  wavefront_aligner_t* const high = wavefront_aligner_new(&attributes);
  attributes.memory_mode = wavefront_memory_singletrack;
  wavefront_aligner_t* const singletrack = wavefront_aligner_new(&attributes);
  const char* const pattern = transpose ? endsfree_text : endsfree_pattern;
  const char* const text = transpose ? endsfree_pattern : endsfree_text;
  const int pattern_length = (int)strlen(pattern);
  const int text_length = (int)strlen(text);
  const int high_status =
      wavefront_align(high,pattern,pattern_length,text,text_length);
  const int singletrack_status =
      wavefront_align(singletrack,pattern,pattern_length,text,text_length);
  char high_label[80], singletrack_label[80];
  const char* const metric =
      distance_metric == gap_affine ? "affine" : "affine2p";
  const char* const heuristic = wfadaptive ? "wfadaptive" : "none";
  const char* const orientation = transpose ? "transposed" : "original";
  snprintf(high_label,sizeof(high_label),"endsfree-%s-%s-%s-high",
      metric,heuristic,orientation);
  snprintf(singletrack_label,sizeof(singletrack_label),
      "endsfree-%s-%s-%s-singletrack",metric,heuristic,orientation);
  int failed = check_endsfree_result(high_label,high,high_status,transpose);
  failed |= check_endsfree_result(
      singletrack_label,singletrack,singletrack_status,transpose);
  if (high->cigar->score != singletrack->cigar->score ||
      high->align_status.score != singletrack->align_status.score) {
    fprintf(stderr,
        "%s memory-mode score mismatch: stored=%d/%d wavefront=%d/%d\n",
        high_label,high->cigar->score,singletrack->cigar->score,
        high->align_status.score,singletrack->align_status.score);
    failed = 1;
  }
  wavefront_aligner_delete(high);
  wavefront_aligner_delete(singletrack);
  return failed;
}

static int check_singletrack_traceback_failure(
    const distance_metric_t distance_metric) {
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = distance_metric;
  attributes.memory_mode = wavefront_memory_singletrack;
  attributes.alignment_scope = compute_alignment;
  attributes.alignment_form.span = alignment_end2end;
  attributes.heuristic.strategy = wf_heuristic_none;
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  const char* const label = distance_metric == gap_affine ?
      "affine-singletrack-failure" : "affine2p-singletrack-failure";
  int failed = 0;
  int status = wavefront_align(wf_aligner,"A",1,"C",1);
  if (status != WF_STATUS_ALG_COMPLETED ||
      wf_aligner->align_status.status != WF_STATUS_ALG_COMPLETED ||
      cigar_is_null(wf_aligner->cigar) || wf_aligner->cigar->score != -4) {
    fprintf(stderr,"%s initial alignment failed: return/status/score=%d/%d/%d\n",
        label,status,wf_aligner->align_status.status,wf_aligner->cigar->score);
    failed = 1;
    goto cleanup;
  }
  failed |= check_cigar_coherent(label,"A",1,"C",1,wf_aligner->cigar,true);
  if (failed) goto cleanup;

  const int saved_score = wf_aligner->align_status.score;
  const int saved_k = wf_aligner->alignment_end_pos.k;
  const wf_offset_t saved_offset = wf_aligner->alignment_end_pos.offset;
  // Create an untraceable retained endpoint, not a natural one-base failure.
  wf_aligner->wf_components.mwavefronts[0]->offsets[0] = WAVEFRONT_OFFSET_NULL;
  status = wavefront_unialign(wf_aligner);
  const cigar_t* const cigar = wf_aligner->cigar;
  if (status != WF_STATUS_UNATTAINABLE ||
      wf_aligner->align_status.status != WF_STATUS_UNATTAINABLE ||
      !cigar_is_null(cigar) || cigar->score != INT32_MIN ||
      cigar->end_v != -1 || cigar->end_h != -1 || cigar->cigar_length != 0 ||
      wf_aligner->align_status.dropped ||
      wf_aligner->align_status.score != saved_score ||
      wf_aligner->alignment_end_pos.k != saved_k ||
      wf_aligner->alignment_end_pos.offset != saved_offset) {
    fprintf(stderr,
        "%s missing predecessor did not return a cleared unattainable result\n",
        label);
    failed = 1;
  }
  status = wavefront_align(wf_aligner,"A",1,"C",1);
  if (status != WF_STATUS_ALG_COMPLETED ||
      wf_aligner->align_status.status != WF_STATUS_ALG_COMPLETED ||
      cigar_is_null(cigar) || cigar->score != -4) {
    fprintf(stderr,"%s reuse failed: return/status/score=%d/%d/%d\n",
        label,status,wf_aligner->align_status.status,cigar->score);
    failed = 1;
  } else {
    failed |= check_cigar_coherent(label,"A",1,"C",1,cigar,true);
  }
cleanup:
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
  failed |= check_endsfree_case(gap_affine,true,false);
  failed |= check_endsfree_case(gap_affine,true,true);
  failed |= check_endsfree_case(gap_affine,false,false);
  failed |= check_endsfree_case(gap_affine,false,true);
  failed |= check_endsfree_case(gap_affine_2p,true,false);
  failed |= check_endsfree_case(gap_affine_2p,true,true);
  failed |= check_endsfree_case(gap_affine_2p,false,false);
  failed |= check_endsfree_case(gap_affine_2p,false,true);
  failed |= check_singletrack_traceback_failure(gap_affine);
  failed |= check_singletrack_traceback_failure(gap_affine_2p);
  return failed;
}
