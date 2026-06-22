/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 */

#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wavefront/wavefront_align.h"

typedef struct {
  int v;
  int h;
} cigar_coord_t;

typedef struct {
  const char* label;
  int pattern_begin_free;
  int pattern_end_free;
  int text_begin_free;
  int text_end_free;
} endsfree_case_t;

static void build_divergent_pair(
    char* const pattern,
    char* const text,
    const int length) {
  static const char bases[] = "ACGT";
  int i;
  for (i=0;i<length;++i) {
    pattern[i] = bases[i%4];
    text[i] = pattern[i];
    if (i%3 == 0) {
      text[i] = bases[(i+1)%4];
    }
  }
  pattern[length] = '\0';
  text[length] = '\0';
}

static wavefront_aligner_t* new_affine_biwfa_aligner(
    const alignment_scope_t scope,
    const wavefront_memory_t memory_mode,
    const endsfree_case_t* const config) {
  wavefront_aligner_attr_t attrs = wavefront_aligner_attr_default;
  attrs.memory_mode = memory_mode;
  attrs.alignment_scope = scope;
  attrs.distance_metric = gap_affine;
  attrs.affine_penalties.match = 0;
  attrs.affine_penalties.mismatch = 4;
  attrs.affine_penalties.gap_opening = 6;
  attrs.affine_penalties.gap_extension = 2;
  attrs.score_only_recover_endpoints = true;
  if (config != NULL) {
    attrs.alignment_form.span = alignment_endsfree;
    attrs.alignment_form.extension = false;
    attrs.alignment_form.pattern_begin_free = config->pattern_begin_free;
    attrs.alignment_form.pattern_end_free = config->pattern_end_free;
    attrs.alignment_form.text_begin_free = config->text_begin_free;
    attrs.alignment_form.text_end_free = config->text_end_free;
  }
  return wavefront_aligner_new(&attrs);
}

static wavefront_aligner_t* new_affine_reward_aligner(
    const alignment_scope_t scope,
    const endsfree_case_t* const config) {
  wavefront_aligner_attr_t attrs = wavefront_aligner_attr_default;
  attrs.memory_mode = wavefront_memory_high;
  attrs.alignment_scope = scope;
  attrs.distance_metric = gap_affine;
  attrs.affine_penalties.match = -1;
  attrs.affine_penalties.mismatch = 4;
  attrs.affine_penalties.gap_opening = 6;
  attrs.affine_penalties.gap_extension = 2;
  attrs.alignment_form.span = alignment_endsfree;
  attrs.alignment_form.extension = false;
  attrs.alignment_form.pattern_begin_free = config->pattern_begin_free;
  attrs.alignment_form.pattern_end_free = config->pattern_end_free;
  attrs.alignment_form.text_begin_free = config->text_begin_free;
  attrs.alignment_form.text_end_free = config->text_end_free;
  return wavefront_aligner_new(&attrs);
}

static int align_or_fail(
    const char* const label,
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const char* const text) {
  const int status = wavefront_align(
      wf_aligner,pattern,(int)strlen(pattern),text,(int)strlen(text));
  if (status != WF_STATUS_ALG_COMPLETED ||
      wf_aligner->align_status.status != WF_STATUS_ALG_COMPLETED) {
    fprintf(stderr,"%s: expected completed alignment, return/status=%d/%d\n",
        label,status,wf_aligner->align_status.status);
    return 1;
  }
  return 0;
}

static bool legal_start(
    const cigar_coord_t coord,
    const endsfree_case_t* const config) {
  return (coord.h == 0 && coord.v <= config->pattern_begin_free) ||
         (coord.v == 0 && coord.h <= config->text_begin_free);
}

static bool legal_end(
    const cigar_coord_t coord,
    const int pattern_length,
    const int text_length,
    const endsfree_case_t* const config) {
  const int pattern_left = pattern_length - coord.v;
  const int text_left = text_length - coord.h;
  return (coord.v == pattern_length &&
          text_left >= 0 &&
          text_left <= config->text_end_free) ||
         (coord.h == text_length &&
          pattern_left >= 0 &&
          pattern_left <= config->pattern_end_free);
}

static int trace_cigar(
    const char* const label,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const cigar_t* const cigar,
    cigar_coord_t** const coords_out,
    int* const num_ops_out) {
  const int num_ops = cigar->end_offset - cigar->begin_offset;
  cigar_coord_t* const coords = malloc((num_ops+1)*sizeof(cigar_coord_t));
  if (coords == NULL) {
    fprintf(stderr,"%s: unable to allocate CIGAR coordinates\n",label);
    return 1;
  }
  int v = 0, h = 0, i;
  coords[0].v = 0;
  coords[0].h = 0;
  for (i=0;i<num_ops;++i) {
    const char op = cigar->operations[cigar->begin_offset+i];
    switch (op) {
      case 'M':
        if (v >= pattern_length || h >= text_length || pattern[v] != text[h]) {
          fprintf(stderr,"%s: incoherent match at (%d,%d)\n",label,v,h);
          free(coords);
          return 1;
        }
        ++v; ++h;
        break;
      case 'X':
        if (v >= pattern_length || h >= text_length || pattern[v] == text[h]) {
          fprintf(stderr,"%s: incoherent mismatch at (%d,%d)\n",label,v,h);
          free(coords);
          return 1;
        }
        ++v; ++h;
        break;
      case 'I':
        if (h >= text_length) {
          fprintf(stderr,"%s: insertion overruns text at %d\n",label,h);
          free(coords);
          return 1;
        }
        ++h;
        break;
      case 'D':
        if (v >= pattern_length) {
          fprintf(stderr,"%s: deletion overruns pattern at %d\n",label,v);
          free(coords);
          return 1;
        }
        ++v;
        break;
      default:
        fprintf(stderr,"%s: unknown CIGAR operation '%c'\n",label,op);
        free(coords);
        return 1;
    }
    coords[i+1].v = v;
    coords[i+1].h = h;
  }
  if (v != pattern_length || h != text_length) {
    fprintf(stderr,"%s: CIGAR traces to (%d,%d), expected (%d,%d)\n",
        label,v,h,pattern_length,text_length);
    free(coords);
    return 1;
  }
  *coords_out = coords;
  *num_ops_out = num_ops;
  return 0;
}

static int find_endpoint_index(
    const char* const label,
    const cigar_t* const cigar,
    const cigar_coord_t* const coords,
    const int num_ops,
    int* const endpoint_index) {
  int i;
  for (i=0;i<=num_ops;++i) {
    if (coords[i].v == cigar->end_v && coords[i].h == cigar->end_h) {
      *endpoint_index = i;
      return 0;
    }
  }
  fprintf(stderr,"%s: stored endpoint (%d,%d) is not on the CIGAR path\n",
      label,cigar->end_v,cigar->end_h);
  return 1;
}

static int check_terminal_suffix(
    const char* const label,
    const cigar_t* const cigar,
    const cigar_coord_t* const coords,
    const int endpoint_index,
    const int num_ops,
    const int pattern_length,
    const int text_length) {
  if (coords[endpoint_index].v == pattern_length &&
      coords[endpoint_index].h == text_length) {
    if (endpoint_index != num_ops) {
      fprintf(stderr,"%s: full endpoint has trailing CIGAR operations\n",label);
      return 1;
    }
    return 0;
  }
  const char expected_op =
      (coords[endpoint_index].v == pattern_length) ? 'I' :
      (coords[endpoint_index].h == text_length) ? 'D' : '\0';
  if (expected_op == '\0') {
    fprintf(stderr,"%s: endpoint (%d,%d) is not on a semiglobal boundary\n",
        label,coords[endpoint_index].v,coords[endpoint_index].h);
    return 1;
  }
  int i;
  for (i=endpoint_index;i<num_ops;++i) {
    const char op = cigar->operations[cigar->begin_offset+i];
    if (op != expected_op) {
      fprintf(stderr,"%s: trailing free suffix has operation '%c', expected '%c'\n",
          label,op,expected_op);
      return 1;
    }
  }
  return 0;
}

static int score_affine_segment(
    const cigar_t* const cigar,
    const int begin_index,
    const int end_index) {
  char last_op = '\0';
  int score = 0, i;
  for (i=begin_index;i<end_index;++i) {
    const char op = cigar->operations[cigar->begin_offset+i];
    switch (op) {
      case 'M':
        break;
      case 'X':
        score -= 4;
        break;
      case 'I':
      case 'D':
        score -= 2;
        if (op != last_op) score -= 6;
        break;
      default:
        fprintf(stderr,"score: unknown CIGAR operation '%c'\n",op);
        exit(1);
    }
    last_op = op;
  }
  return score;
}

static int check_endpoint_score(
    const char* const label,
    const cigar_t* const cigar,
    const cigar_coord_t* const coords,
    const int endpoint_index,
    const endsfree_case_t* const config,
    const int expected_score) {
  int best_score = INT_MIN, i;
  for (i=0;i<=endpoint_index;++i) {
    if (!legal_start(coords[i],config)) continue;
    const int score = score_affine_segment(cigar,i,endpoint_index);
    if (score > best_score) best_score = score;
  }
  if (best_score != expected_score) {
    fprintf(stderr,"%s: endpoint/CIGAR score mismatch observed=%d expected=%d "
        "endpoint=(%d,%d)\n",
        label,best_score,expected_score,cigar->end_v,cigar->end_h);
    return 1;
  }
  return 0;
}

static int check_biwfa_endpoint_case(
    const endsfree_case_t* const config,
    const char* const pattern,
    const char* const text) {
  int failed = 0;
  wavefront_aligner_t* const oracle = new_affine_biwfa_aligner(
      compute_score,wavefront_memory_high,config);
  wavefront_aligner_t* const biwfa = new_affine_biwfa_aligner(
      compute_alignment,wavefront_memory_ultralow,config);
  failed |= align_or_fail(config->label,oracle,pattern,text);
  failed |= align_or_fail(config->label,biwfa,pattern,text);
  if (failed) {
    wavefront_aligner_delete(oracle);
    wavefront_aligner_delete(biwfa);
    return 1;
  }

  const int pattern_length = (int)strlen(pattern);
  const int text_length = (int)strlen(text);
  cigar_t* const cigar = biwfa->cigar;
  cigar_coord_t* coords = NULL;
  int num_ops = 0, endpoint_index = -1;
  failed |= trace_cigar(
      config->label,pattern,pattern_length,text,text_length,
      cigar,&coords,&num_ops);
  if (!failed) {
    failed |= find_endpoint_index(
        config->label,cigar,coords,num_ops,&endpoint_index);
  }
  if (!failed && !legal_end(
      coords[endpoint_index],pattern_length,text_length,config)) {
    fprintf(stderr,"%s: endpoint (%d,%d) is not legal for free ends (%d,%d,%d,%d)\n",
        config->label,cigar->end_v,cigar->end_h,
        config->pattern_begin_free,config->pattern_end_free,
        config->text_begin_free,config->text_end_free);
    failed = 1;
  }
  if (!failed) {
    failed |= check_terminal_suffix(
        config->label,cigar,coords,endpoint_index,num_ops,
        pattern_length,text_length);
  }
  if (!failed && cigar->score != oracle->cigar->score) {
    fprintf(stderr,"%s: stored score mismatch observed=%d expected=%d\n",
        config->label,cigar->score,oracle->cigar->score);
    failed = 1;
  }
  if (!failed) {
    failed |= check_endpoint_score(
        config->label,cigar,coords,endpoint_index,config,
        oracle->cigar->score);
  }

  free(coords);
  wavefront_aligner_delete(oracle);
  wavefront_aligner_delete(biwfa);
  return failed;
}

static int check_score_only_endpoint(
    const char* const label,
    const wavefront_memory_t memory_mode,
    const endsfree_case_t* const config,
    const char* const pattern,
    const char* const text,
    const int expected_score,
    const int expected_end_v,
    const int expected_end_h) {
  wavefront_aligner_t* const wf_aligner = new_affine_biwfa_aligner(
      compute_score,memory_mode,config);
  int failed = align_or_fail(label,wf_aligner,pattern,text);
  if (!failed && wf_aligner->cigar->score != expected_score) {
    fprintf(stderr,"%s: score-only score mismatch observed=%d expected=%d\n",
        label,wf_aligner->cigar->score,expected_score);
    failed = 1;
  }
  if (!failed &&
      (wf_aligner->cigar->end_v != expected_end_v ||
       wf_aligner->cigar->end_h != expected_end_h)) {
    fprintf(stderr,"%s: score-only endpoint mismatch observed=(%d,%d) "
        "expected=(%d,%d)\n",
        label,wf_aligner->cigar->end_v,wf_aligner->cigar->end_h,
        expected_end_v,expected_end_h);
    failed = 1;
  }
  wavefront_aligner_delete(wf_aligner);
  return failed;
}

static int check_score_only_endpoint_case(
    const endsfree_case_t* const config,
    const char* const pattern,
    const char* const text,
    const int expected_end_v,
    const int expected_end_h) {
  int failed = 0;
  wavefront_aligner_t* const oracle = new_affine_biwfa_aligner(
      compute_score,wavefront_memory_high,config);
  failed |= align_or_fail(config->label,oracle,pattern,text);
  const int expected_score = oracle->cigar->score;
  wavefront_aligner_delete(oracle);
  if (failed) return failed;

  char label[128];
  snprintf(label,sizeof(label),"%s-score-high",config->label);
  failed |= check_score_only_endpoint(
      label,wavefront_memory_high,config,pattern,text,
      expected_score,expected_end_v,expected_end_h);
  snprintf(label,sizeof(label),"%s-score-ultralow",config->label);
  failed |= check_score_only_endpoint(
      label,wavefront_memory_ultralow,config,pattern,text,
      expected_score,expected_end_v,expected_end_h);
  return failed;
}

static int check_unidirectional_reward_score_case(
    const endsfree_case_t* const config,
    const char* const pattern,
    const char* const text) {
  int failed = 0;
  wavefront_aligner_t* const score_only = new_affine_reward_aligner(
      compute_score,config);
  wavefront_aligner_t* const alignment = new_affine_reward_aligner(
      compute_alignment,config);
  failed |= align_or_fail(config->label,score_only,pattern,text);
  failed |= align_or_fail(config->label,alignment,pattern,text);
  if (!failed && score_only->cigar->score != alignment->cigar->score) {
    fprintf(stderr,"%s-reward: score-only score mismatch observed=%d "
        "expected=%d\n",
        config->label,score_only->cigar->score,alignment->cigar->score);
    failed = 1;
  }
  if (!failed &&
      (score_only->cigar->end_v != alignment->cigar->end_v ||
       score_only->cigar->end_h != alignment->cigar->end_h)) {
    fprintf(stderr,"%s-reward: score-only endpoint mismatch observed=(%d,%d) "
        "expected=(%d,%d)\n",
        config->label,
        score_only->cigar->end_v,score_only->cigar->end_h,
        alignment->cigar->end_v,alignment->cigar->end_h);
    failed = 1;
  }
  wavefront_aligner_delete(score_only);
  wavefront_aligner_delete(alignment);
  return failed;
}

static int check_global_endpoint(void) {
  char pattern[256], text[256];
  build_divergent_pair(pattern,text,180);
  wavefront_aligner_t* const biwfa = new_affine_biwfa_aligner(
      compute_alignment,wavefront_memory_ultralow,NULL);
  int failed = align_or_fail("global",biwfa,pattern,text);
  if (!failed &&
      (biwfa->cigar->end_v != (int)strlen(pattern) ||
       biwfa->cigar->end_h != (int)strlen(text))) {
    fprintf(stderr,"global: expected full endpoint (%lu,%lu), observed (%d,%d)\n",
        (unsigned long)strlen(pattern),(unsigned long)strlen(text),
        biwfa->cigar->end_v,biwfa->cigar->end_h);
    failed = 1;
  }
  wavefront_aligner_delete(biwfa);
  return failed;
}

int main(void) {
  char pattern[1024], text[1024], core_pattern[1024], core_text[1024];
  build_divergent_pair(core_pattern,core_text,900);

  int failed = 0;
  snprintf(pattern,sizeof(pattern),"%s",core_pattern);
  snprintf(text,sizeof(text),"%sTTTTTTTT",core_text);
  const endsfree_case_t text_end = {
      "text-end",0,0,0,8
  };
  failed |= check_biwfa_endpoint_case(&text_end,pattern,text);
  failed |= check_score_only_endpoint_case(&text_end,pattern,text,900,900);
  failed |= check_unidirectional_reward_score_case(&text_end,pattern,text);

  snprintf(pattern,sizeof(pattern),"%sGGGGGGGG",core_pattern);
  snprintf(text,sizeof(text),"%s",core_text);
  const endsfree_case_t pattern_end = {
      "pattern-end",0,8,0,0
  };
  failed |= check_biwfa_endpoint_case(&pattern_end,pattern,text);
  failed |= check_score_only_endpoint_case(&pattern_end,pattern,text,900,900);
  failed |= check_unidirectional_reward_score_case(&pattern_end,pattern,text);

  snprintf(pattern,sizeof(pattern),"CCCC%s",core_pattern);
  snprintf(text,sizeof(text),"%sAAAA",core_text);
  const endsfree_case_t combined = {
      "pattern-begin-text-end",4,0,0,4
  };
  failed |= check_biwfa_endpoint_case(&combined,pattern,text);
  failed |= check_score_only_endpoint_case(&combined,pattern,text,904,900);
  failed |= check_unidirectional_reward_score_case(&combined,pattern,text);

  const endsfree_case_t empty_text_pattern_end = {
      "empty-text-pattern-end",0,1,0,0
  };
  failed |= check_biwfa_endpoint_case(&empty_text_pattern_end,"A","");
  failed |= check_score_only_endpoint_case(&empty_text_pattern_end,"A","",0,0);

  const endsfree_case_t empty_pattern_text_end = {
      "empty-pattern-text-end",0,0,0,1
  };
  failed |= check_biwfa_endpoint_case(&empty_pattern_text_end,"","A");
  failed |= check_score_only_endpoint_case(&empty_pattern_text_end,"","A",0,0);

  const endsfree_case_t both_empty = {
      "both-empty",0,0,0,0
  };
  failed |= check_biwfa_endpoint_case(&both_empty,"","");
  failed |= check_score_only_endpoint_case(&both_empty,"","",0,0);

  failed |= check_global_endpoint();
  return failed;
}
