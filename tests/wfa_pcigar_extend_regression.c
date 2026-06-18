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

int pcigar_unpack_extend(
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    int v,
    int h,
    char* cigar_buffer);

static int check_no_remaining_bases(void) {
  const char pattern[8] = { 0 };
  const char text[8] = { 0 };
  char cigar[16];
  memset(cigar,'?',sizeof(cigar));

  const int matches = pcigar_unpack_extend(pattern,0,text,0,0,0,cigar);
  if (matches != 0) {
    fprintf(stderr,
        "pcigar_unpack_extend returned %d matches with no bases remaining\n",
        matches);
    return 1;
  }
  return 0;
}

static int check_exact_block_match(void) {
  const char pattern[] = "ACGTACGT";
  const char text[] = "ACGTACGT";
  char cigar[16];
  memset(cigar,'?',sizeof(cigar));

  const int matches = pcigar_unpack_extend(pattern,8,text,8,0,0,cigar);
  if (matches != 8 || memcmp(cigar,"MMMMMMMM",8) != 0) {
    fprintf(stderr,
        "pcigar_unpack_extend failed an exact 8-base block match: matches=%d\n",
        matches);
    return 1;
  }
  return 0;
}

static int check_mismatch_inside_block(void) {
  const char pattern[] = "ACGTACGT";
  const char text[] = "ACGTTCGT";
  char cigar[16];
  memset(cigar,'?',sizeof(cigar));

  const int matches = pcigar_unpack_extend(pattern,8,text,8,0,0,cigar);
  if (matches != 4 || memcmp(cigar,"MMMM",4) != 0 || cigar[4] != '?') {
    fprintf(stderr,
        "pcigar_unpack_extend failed a mismatch inside a block: "
        "matches=%d cigar[4]=%c\n",
        matches,cigar[4]);
    return 1;
  }
  return 0;
}

int main(void) {
  int failed = 0;
  failed |= check_no_remaining_bases();
  failed |= check_exact_block_match();
  failed |= check_mismatch_inside_block();
  return failed;
}
