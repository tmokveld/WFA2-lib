/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "benchmark/benchmark_gap_affine.h"
#include "gap_affine/affine_matrix.h"
#include "gap_affine/swg.h"
#include "system/mm_allocator.h"
#include "wavefront/wavefront_align.h"

static int captured_score = INT32_MIN;
static int captured_calls = 0;

void benchmark_print_output(
    align_input_t* const align_input,
    const distance_metric_t distance_metric,
    const bool score_only,
    cigar_t* const cigar) {
  (void)align_input;
  (void)distance_metric;
  (void)score_only;
  captured_score = cigar->score;
  ++captured_calls;
}

void benchmark_record_wfa_memory(
    align_input_t* const align_input) {
  (void)align_input;
}

void benchmark_check_alignment(
    align_input_t* const align_input,
    cigar_t* const cigar_computed) {
  (void)align_input;
  (void)cigar_computed;
}

int wavefront_align(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  (void)wf_aligner;
  (void)pattern;
  (void)pattern_length;
  (void)text;
  (void)text_length;
  return 0;
}

int wavefront_align_lambda(
    wavefront_aligner_t* const wf_aligner,
    alignment_match_funct_t const match_funct,
    void* match_funct_arguments,
    const int pattern_length,
    const int text_length) {
  (void)wf_aligner;
  (void)match_funct;
  (void)match_funct_arguments;
  (void)pattern_length;
  (void)text_length;
  return 0;
}

static int direct_swg_endsfree_score(
    mm_allocator_t* const mm_allocator,
    affine_penalties_t* const penalties,
    const char* const pattern,
    const char* const text,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free) {
  const int pattern_length = (int)strlen(pattern);
  const int text_length = (int)strlen(text);
  affine_matrix_t affine_matrix;
  affine_matrix_allocate(
      &affine_matrix,pattern_length+1,text_length+1,mm_allocator);
  cigar_t* const cigar = cigar_new(pattern_length+text_length+8);

  swg_align_endsfree(
      &affine_matrix,penalties,
      pattern,pattern_length,
      text,text_length,
      pattern_begin_free,pattern_end_free,
      text_begin_free,text_end_free,cigar);
  const int score = cigar->score;

  cigar_free(cigar);
  affine_matrix_free(&affine_matrix,mm_allocator);
  return score;
}

static int check_wrapper_score(
    const char* const label,
    const char* const pattern,
    const char* const text,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free) {
  affine_penalties_t penalties = {
    .match = 0,
    .mismatch = 4,
    .gap_opening = 6,
    .gap_extension = 2,
  };
  mm_allocator_t* const mm_allocator = mm_allocator_new(1<<20);
  const int expected_score = direct_swg_endsfree_score(
      mm_allocator,&penalties,pattern,text,
      pattern_begin_free,pattern_end_free,
      text_begin_free,text_end_free);

  align_input_t align_input;
  memset(&align_input,0,sizeof(align_input));
  align_input.pattern = (char*)pattern;
  align_input.pattern_length = (int)strlen(pattern);
  align_input.text = (char*)text;
  align_input.text_length = (int)strlen(text);
  align_input.pattern_begin_free = pattern_begin_free;
  align_input.pattern_end_free = pattern_end_free;
  align_input.text_begin_free = text_begin_free;
  align_input.text_end_free = text_end_free;
  align_input.output_file = stdout;
  align_input.mm_allocator = mm_allocator;

  captured_score = INT32_MIN;
  captured_calls = 0;
  benchmark_gap_affine_swg_endsfree(&align_input,&penalties);

  const int failed =
      captured_calls != 1 || captured_score != expected_score;
  if (failed) {
    fprintf(stderr,
        "%s: expected score %d, captured %d, calls %d\n",
        label,expected_score,captured_score,captured_calls);
  }

  mm_allocator_delete(mm_allocator);
  return failed;
}

int main(void) {
  int failed = 0;

  failed |= check_wrapper_score(
      "pattern-end-only","A","",0,1,0,0);
  failed |= check_wrapper_score(
      "symmetric-pattern-free","A","",1,1,0,0);
  failed |= check_wrapper_score(
      "global-like","A","",0,0,0,0);

  return failed;
}
