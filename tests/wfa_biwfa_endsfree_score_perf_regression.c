/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "wavefront/wavefront_align.h"

#define SEQUENCE_LENGTH 10000
#define NUM_REPEATS 24
#define NUM_SAMPLES 3
#define MAX_ENDSFREE_GLOBAL_RATIO 3.0

static double now_seconds(void) {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
}

static void build_pair(
    char* const pattern,
    char* const text,
    const int length) {
  static const char bases[] = "ACGT";
  int i;
  for (i=0;i<length;++i) {
    pattern[i] = bases[i%4];
    text[i] = pattern[i];
    if (i%47 == 0) text[i] = bases[(i+1)%4];
    if (i%113 == 0) text[i] = bases[(i+2)%4];
  }
  pattern[length] = '\0';
  text[length] = '\0';
}

static wavefront_aligner_t* new_aligner(
    const wavefront_memory_t memory_mode,
    const bool ends_free) {
  wavefront_aligner_attr_t attrs = wavefront_aligner_attr_default;
  attrs.memory_mode = memory_mode;
  attrs.alignment_scope = compute_score;
  attrs.distance_metric = gap_affine;
  attrs.affine_penalties.match = 0;
  attrs.affine_penalties.mismatch = 4;
  attrs.affine_penalties.gap_opening = 6;
  attrs.affine_penalties.gap_extension = 2;
  attrs.heuristic.strategy = wf_heuristic_none;
  if (ends_free) {
    attrs.alignment_form.span = alignment_endsfree;
    attrs.alignment_form.pattern_begin_free = 0;
    attrs.alignment_form.pattern_end_free = 500;
    attrs.alignment_form.text_begin_free = 0;
    attrs.alignment_form.text_end_free = 500;
  }
  return wavefront_aligner_new(&attrs);
}

static int align_once(
    wavefront_aligner_t* const aligner,
    const char* const pattern,
    const char* const text,
    int* const score_out) {
  const int status = wavefront_align(
      aligner,pattern,SEQUENCE_LENGTH,text,SEQUENCE_LENGTH);
  if (status != WF_STATUS_ALG_COMPLETED ||
      aligner->align_status.status != WF_STATUS_ALG_COMPLETED) {
    fprintf(stderr,"alignment failed: return/status=%d/%d\n",
        status,aligner->align_status.status);
    return 1;
  }
  if (score_out != NULL) *score_out = aligner->cigar->score;
  return 0;
}

static double time_batch(
    wavefront_aligner_t* const aligner,
    const char* const pattern,
    const char* const text) {
  double best = 0.0;
  int sample;
  for (sample=0;sample<NUM_SAMPLES;++sample) {
    const double begin = now_seconds();
    int i;
    for (i=0;i<NUM_REPEATS;++i) {
      if (align_once(aligner,pattern,text,NULL)) return -1.0;
    }
    const double elapsed = now_seconds() - begin;
    if (best == 0.0 || elapsed < best) best = elapsed;
  }
  return best;
}

int main(void) {
  char* const pattern = malloc(SEQUENCE_LENGTH+1);
  char* const text = malloc(SEQUENCE_LENGTH+1);
  if (pattern == NULL || text == NULL) {
    fprintf(stderr,"unable to allocate sequences\n");
    free(pattern);
    free(text);
    return 1;
  }
  build_pair(pattern,text,SEQUENCE_LENGTH);

  wavefront_aligner_t* const oracle =
      new_aligner(wavefront_memory_high,true);
  wavefront_aligner_t* const global =
      new_aligner(wavefront_memory_ultralow,false);
  wavefront_aligner_t* const endsfree =
      new_aligner(wavefront_memory_ultralow,true);

  int oracle_score = 0, endsfree_score = 0;
  int failed = align_once(oracle,pattern,text,&oracle_score);
  failed |= align_once(endsfree,pattern,text,&endsfree_score);
  if (!failed && oracle_score != endsfree_score) {
    fprintf(stderr,"score mismatch: oracle=%d endsfree-ultralow=%d\n",
        oracle_score,endsfree_score);
    failed = 1;
  }

  if (!failed) {
    const double global_time = time_batch(global,pattern,text);
    const double endsfree_time = time_batch(endsfree,pattern,text);
    const double ratio = (global_time > 0.0) ?
        endsfree_time / global_time : MAX_ENDSFREE_GLOBAL_RATIO + 1.0;
    if (global_time <= 0.0 || endsfree_time <= 0.0 ||
        ratio > MAX_ENDSFREE_GLOBAL_RATIO) {
      fprintf(stderr,"ends-free score-only too slow: global=%.6fs "
          "endsfree=%.6fs ratio=%.2f max=%.2f\n",
          global_time,endsfree_time,ratio,MAX_ENDSFREE_GLOBAL_RATIO);
      failed = 1;
    }
  }

  wavefront_aligner_delete(oracle);
  wavefront_aligner_delete(global);
  wavefront_aligner_delete(endsfree);
  free(pattern);
  free(text);
  return failed;
}
