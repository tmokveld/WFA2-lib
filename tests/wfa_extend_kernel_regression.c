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
#include "wavefront/wavefront.h"
#include "wavefront/wavefront_extend_kernels.h"
#include "wavefront/wavefront_sequences.h"

static wf_offset_t extend_offset(
    const char* const pattern,
    const char* const text) {
  wavefront_aligner_t wf_aligner;
  wavefront_t wavefront;
  memset(&wf_aligner,0,sizeof(wf_aligner));

  mm_allocator_t* const mm_allocator = mm_allocator_new(4096);
  wavefront_sequences_allocate(&wf_aligner.sequences);
  wavefront_sequences_init_ascii(&wf_aligner.sequences,
      pattern,(int)strlen(pattern),
      text,(int)strlen(text),
      false);

  wavefront_allocate(&wavefront,1,false,mm_allocator);
  wavefront_init(&wavefront,0,0);
  wavefront_set_limits(&wavefront,0,0);
  wavefront.offsets[0] = 0;

  wavefront_extend_matches_packed_end2end(&wf_aligner,&wavefront,0,0);
  const wf_offset_t offset = wavefront.offsets[0];

  wavefront_free(&wavefront,mm_allocator);
  wavefront_sequences_free(&wf_aligner.sequences);
  mm_allocator_delete(mm_allocator);
  return offset;
}

static int check_mismatch_window(
    const int window_begin,
    const int window_end) {
  int failed = 0;
  char pattern[129];
  char text[129];

  for (int mismatch_pos=window_begin;mismatch_pos<=window_end;++mismatch_pos) {
    memset(pattern,'A',128);
    memset(text,'A',128);
    pattern[128] = '\0';
    text[128] = '\0';
    text[mismatch_pos] = 'C';

    const wf_offset_t offset = extend_offset(pattern,text);
    if (offset != mismatch_pos) {
      fprintf(stderr,
          "extend kernel stopped at offset %d with mismatch at %d\n",
          offset,mismatch_pos);
      failed = 1;
    }
  }
  return failed;
}

static int check_sentinel_termination(void) {
  const char pattern[] =
      "ACGTACGTACGTACGTACGTACGTACGTACGT"
      "ACGTACGTACGTACGTACGTACGTACGTACGT"
      "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const char text[] =
      "ACGTACGTACGTACGTACGTACGTACGTACGT"
      "ACGTACGTACGTACGTACGTACGTACGTACGT"
      "ACGTACGTACGTACGTACGTACGTACGTACGT";
  const wf_offset_t offset = extend_offset(pattern,text);
  if (offset != 96) {
    fprintf(stderr,
        "extend kernel did not stop at sentinels after exact match: offset=%d\n",
        offset);
    return 1;
  }
  return 0;
}

int main(void) {
  int failed = 0;
  failed |= check_mismatch_window(0,16);
  failed |= check_mismatch_window(64,80);
  failed |= check_sentinel_termination();
  return failed;
}
