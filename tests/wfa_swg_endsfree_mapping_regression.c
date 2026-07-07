/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 */

#include <stdio.h>
#include <string.h>

#include "system/mm_allocator.h"
#include "gap_affine/affine_matrix.h"
#include "gap_affine/swg.h"

static int swg_endsfree_score(
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

static int check_swg_score(
    const char* const label,
    const char* const pattern,
    const char* const text,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free,
    const int expected_score) {
  affine_penalties_t penalties = {
    .match = 0,
    .mismatch = 4,
    .gap_opening = 6,
    .gap_extension = 2,
  };
  mm_allocator_t* const mm_allocator = mm_allocator_new(1<<20);
  const int observed_score = swg_endsfree_score(
      mm_allocator,&penalties,pattern,text,
      pattern_begin_free,pattern_end_free,
      text_begin_free,text_end_free);

  const int failed = observed_score != expected_score;
  if (failed) {
    fprintf(stderr,
        "%s: expected score %d, observed %d\n",
        label,expected_score,observed_score);
  }

  mm_allocator_delete(mm_allocator);
  return failed;
}

int main(void) {
  int failed = 0;

  failed |= check_swg_score(
      "pattern-end-only","A","",0,1,0,0,0);
  failed |= check_swg_score(
      "symmetric-pattern-free","A","",1,1,0,0,0);
  failed |= check_swg_score(
      "global-like","A","",0,0,0,0,8);

  return failed;
}
