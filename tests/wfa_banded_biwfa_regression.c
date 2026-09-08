/*
 * Wavefront Alignment Algorithms
 * Regression coverage for static BiWFA band coordinates and adaptive movement.
 */
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "wavefront/wavefront_align.h"

static int check_static_band(
    const char* const pattern,
    const char* const text,
    const int min_k,
    const int max_k,
    const int text_end_free) {
  const int pattern_length = (int)strlen(pattern);
  const int text_length = (int)strlen(text);
  int expected_score = 0;
  const wavefront_memory_t modes[] = {
      wavefront_memory_high, wavefront_memory_ultralow,
      wavefront_memory_singletrack};
  for (int mode=0;mode<3;++mode) {
    const int scopes = modes[mode] == wavefront_memory_singletrack ? 1 : 2;
    for (int score_only=0;score_only<scopes;++score_only) {
      wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
      attr.distance_metric = gap_affine_2p;
      attr.affine2p_penalties = (affine2p_penalties_t){0,5,8,2,24,1};
      attr.memory_mode = modes[mode];
      attr.alignment_scope = score_only ? compute_score : compute_alignment;
      attr.score_only_recover_endpoints = true;
      attr.heuristic.strategy = wf_heuristic_none;
      wavefront_aligner_t* const aligner = wavefront_aligner_new(&attr);
      wavefront_aligner_set_heuristic_banded_static(aligner,min_k,max_k);
      if (text_end_free > 0) {
        wavefront_aligner_set_alignment_free_ends(aligner,0,0,0,text_end_free);
      }
      const int status = wavefront_align(aligner,pattern,pattern_length,text,text_length);
      if (mode == 0 && !score_only) expected_score = aligner->cigar->score;
      int failed = status != WF_STATUS_ALG_COMPLETED ||
          aligner->cigar->score != expected_score;
      if (!score_only && status == WF_STATUS_ALG_COMPLETED) {
        failed |= !cigar_check_alignment(stderr,pattern,pattern_length,
            text,text_length,aligner->cigar,true);
      }
      if (score_only && text_end_free > 0) {
        failed |= aligner->cigar->end_v != pattern_length ||
            aligner->cigar->end_h < text_length-text_end_free ||
            aligner->cigar->end_h > text_length;
      }
      if (failed) {
        fprintf(stderr,"Band [%d,%d] memory=%d score_only=%d end_free=%d: "
            "status=%d score=%d expected=%d endpoint=(%d,%d)\n",
            min_k,max_k,modes[mode],score_only,text_end_free,status,
            aligner->cigar->score,expected_score,
            aligner->cigar->end_v,aligner->cigar->end_h);
      }
      wavefront_aligner_delete(aligner);
      if (failed) return 1;
    }
  }
  return 0;
}

static int check_unreachable_band(void) {
  char pattern[128], text[256];
  memset(pattern,'A',sizeof(pattern));
  memset(text,'C',sizeof(text));
  const distance_metric_t distances[] = {indel,edit};
  const wavefront_memory_t modes[] = {
      wavefront_memory_high, wavefront_memory_ultralow};
  for (int distance=0;distance<2;++distance) {
    for (int mode=0;mode<2;++mode) {
      const int expected_status = modes[mode] == wavefront_memory_ultralow ?
          WF_STATUS_UNATTAINABLE : WF_STATUS_ALG_PARTIAL;
      for (int score_only=0;score_only<2;++score_only) {
        wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
        attr.distance_metric = distances[distance];
        attr.memory_mode = modes[mode];
        attr.alignment_scope = score_only ? compute_score : compute_alignment;
        attr.heuristic.strategy = wf_heuristic_none;
        attr.system.max_alignment_steps = 1024;
        wavefront_aligner_t* const aligner = wavefront_aligner_new(&attr);
        wavefront_aligner_set_heuristic_banded_static(aligner,-4,4);
        for (int reverse=0;reverse<2;++reverse) {
          const int status = reverse ?
              wavefront_align(aligner,text,sizeof(text),pattern,sizeof(pattern)) :
              wavefront_align(aligner,pattern,sizeof(pattern),text,sizeof(text));
          if (status != expected_status) {
            fprintf(stderr,"Unreachable band: distance=%d memory=%d "
                "score_only=%d reverse=%d status=%d\n",
                distances[distance],modes[mode],score_only,reverse,status);
            wavefront_aligner_delete(aligner);
            return 1;
          }
        }
        wavefront_aligner_delete(aligner);
      }
    }
  }
  return 0;
}

static int check_wide_band(
    const char* const pattern,
    const char* const text,
    const int expected_score) {
  const int pattern_length = (int)strlen(pattern);
  const int text_length = (int)strlen(text);
  for (int score_only=0;score_only<2;++score_only) {
    wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
    attr.distance_metric = gap_affine;
    attr.affine_penalties = (affine_penalties_t){0,4,6,2};
    attr.memory_mode = wavefront_memory_ultralow;
    attr.alignment_scope = score_only ? compute_score : compute_alignment;
    attr.heuristic.strategy = wf_heuristic_none;
    wavefront_aligner_t* const aligner = wavefront_aligner_new(&attr);
    wavefront_aligner_set_heuristic_banded_static(aligner,INT_MIN,INT_MAX);
    const int status = wavefront_align(aligner,pattern,pattern_length,text,text_length);
    int failed = status != WF_STATUS_ALG_COMPLETED ||
        aligner->cigar->score != expected_score;
    if (!score_only && status == WF_STATUS_ALG_COMPLETED) {
      failed |= !cigar_check_alignment(stderr,pattern,pattern_length,
          text,text_length,aligner->cigar,true);
    }
    if (failed) {
      fprintf(stderr,"Wide band: score_only=%d status=%d score=%d expected=%d\n",
          score_only,status,aligner->cigar->score,expected_score);
    }
    wavefront_aligner_delete(aligner);
    if (failed) return 1;
  }
  return 0;
}

static int check_adaptive_movement(void) {
  char pattern[129], text[149];
  for (int i=0;i<128;++i) pattern[i] = "ACGT"[i%4];
  pattern[128] = '\0';
  memset(text,'T',20);
  memcpy(text+20,pattern,sizeof(pattern));
  const wavefront_memory_t modes[] = {
      wavefront_memory_high, wavefront_memory_singletrack};
  for (int mode=0;mode<2;++mode) {
    wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
    attr.memory_mode = modes[mode];
    attr.heuristic.strategy = wf_heuristic_none;
    wavefront_aligner_t* const aligner = wavefront_aligner_new(&attr);
    wavefront_aligner_set_heuristic_banded_adaptive(aligner,-4,4,1);
    const int status = wavefront_align(aligner,pattern,128,text,148);
    const int failed = status != WF_STATUS_ALG_COMPLETED ||
        !cigar_check_alignment(stderr,pattern,128,text,148,aligner->cigar,true);
    if (failed) fprintf(stderr,"Adaptive band could not move to diagonal 20 (memory=%d)\n",modes[mode]);
    wavefront_aligner_delete(aligner);
    if (failed) return 1;
  }
  return 0;
}

int main(void) {
  char pattern[257], text[273];
  unsigned int state = 42;
  for (int i=0;i<256;++i) {
    state = 1664525u*state + 1013904223u;
    pattern[i] = "ACGT"[state>>30];
  }
  pattern[256] = '\0';
  memset(text,'T',16);
  memcpy(text+16,pattern,240);
  text[256] = '\0';
  for (int i=24;i<240;i+=7) text[i] = text[i]=='A' ? 'C' : 'A';
  int failed = check_static_band(pattern,text,-4,20,0);
  failed |= check_static_band(text,pattern,-20,4,0);
  failed |= check_static_band(pattern,text,-4,4,0);
  // End-free score-only recovery uses a separate recursive/base path.
  memset(text+256,'G',16);
  text[272] = '\0';
  failed |= check_static_band(pattern,text,-4,20,16);
  failed |= check_adaptive_movement();
  failed |= check_unreachable_band();
  char wide_pattern[1025], wide_text[1025];
  memset(wide_pattern,'A',512);
  memset(wide_text,'C',512);
  wide_pattern[512] = wide_text[512] = '\0';
  failed |= check_wide_band(wide_pattern,wide_text,-2048);
  // Force recursive/base sub-alignments with nonzero diagonal shifts.
  memset(wide_pattern+512,'C',512);
  memset(wide_text,'G',128);
  memcpy(wide_text+128,wide_pattern,896);
  wide_pattern[1024] = wide_text[1024] = '\0';
  failed |= check_wide_band(wide_pattern,wide_text,-524);
  failed |= check_wide_band(wide_text,wide_pattern,-524);
  return failed;
}
