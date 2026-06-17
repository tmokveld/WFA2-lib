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
#include <stdlib.h>
#include <string.h>

#include "wavefront/wavefront_align.h"

static uint8_t encode_base(const char base) {
  switch (base) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:
      fprintf(stderr,"Unexpected base '%c'\n",base);
      exit(1);
  }
}

static uint8_t* pack2bits(const char* const sequence) {
  const size_t sequence_length = strlen(sequence);
  const size_t packed_length = (sequence_length + 3) / 4;
  uint8_t* const packed = calloc(packed_length,sizeof(uint8_t));
  if (packed == NULL) {
    fprintf(stderr,"Could not allocate packed sequence\n");
    exit(1);
  }
  size_t i;
  for (i=0;i<sequence_length;++i) {
    packed[i/4] |= encode_base(sequence[i]) << (2*(i%4));
  }
  return packed;
}

static int align_ascii(
    const char* const pattern,
    const char* const text,
    const wavefront_memory_t memory_mode) {
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = edit;
  attributes.heuristic.strategy = wf_heuristic_none;
  attributes.memory_mode = memory_mode;

  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  wavefront_align(wf_aligner,
      pattern,(int)strlen(pattern),
      text,(int)strlen(text));
  const int score = wf_aligner->cigar->score;
  wavefront_aligner_delete(wf_aligner);
  return score;
}

static int align_packed2bits(
    const char* const pattern,
    const char* const text,
    const wavefront_memory_t memory_mode) {
  uint8_t* const packed_pattern = pack2bits(pattern);
  uint8_t* const packed_text = pack2bits(text);

  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.distance_metric = edit;
  attributes.heuristic.strategy = wf_heuristic_none;
  attributes.memory_mode = memory_mode;

  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  wavefront_align_packed2bits(wf_aligner,
      packed_pattern,(int)strlen(pattern),
      packed_text,(int)strlen(text));
  const int score = wf_aligner->cigar->score;

  wavefront_aligner_delete(wf_aligner);
  free(packed_pattern);
  free(packed_text);
  return score;
}

static int check_case(
    const char* const pattern,
    const char* const text,
    const wavefront_memory_t memory_mode) {
  const int ascii_score = align_ascii(pattern,text,memory_mode);
  const int packed2bits_score = align_packed2bits(pattern,text,memory_mode);
  if (ascii_score != packed2bits_score) {
    fprintf(stderr,
        "Packed-2bit score mismatch for pattern='%s' text='%s' mode=%d: ascii=%d packed2bits=%d\n",
        pattern,text,memory_mode,ascii_score,packed2bits_score);
    return 1;
  }
  return 0;
}

int main(void) {
  int failed = 0;
  failed |= check_case("AAAATTTT","AAAACCCC",wavefront_memory_high);
  failed |= check_case("AAAATTTT","AAAACCCC",wavefront_memory_ultralow);
  failed |= check_case("ACGTA","ACGTC",wavefront_memory_high);
  failed |= check_case("ACGTA","ACGTC",wavefront_memory_ultralow);
  failed |= check_case("ACGTACGTA","ACGTACGTC",wavefront_memory_high);
  failed |= check_case("ACGTACGTA","ACGTACGTC",wavefront_memory_ultralow);
  return failed;
}
