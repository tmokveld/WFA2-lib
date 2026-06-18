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

#include "alignment/cigar.h"

static int check_sam_cigar(
    const char* const operations,
    const bool show_mismatches,
    const char* const expected) {
  cigar_t* const cigar = cigar_new((int)strlen(operations));
  strcpy(cigar->operations,operations);
  cigar->begin_offset = 0;
  cigar->end_offset = (int)strlen(operations);

  char buffer[64];
  cigar_sprint_SAM_CIGAR(buffer,cigar,show_mismatches);

  const int failed = strcmp(buffer,expected) != 0;
  if (failed) {
    fprintf(stderr,
        "Unexpected SAM CIGAR for operations='%s' show_mismatches=%d: "
        "expected '%s', got '%s'\n",
        operations,show_mismatches,expected,buffer);
  }

  cigar_free(cigar);
  return failed;
}

int main(void) {
  int failed = 0;

  failed |= check_sam_cigar("X",false,"1M");
  failed |= check_sam_cigar("XM",false,"2M");
  failed |= check_sam_cigar("XXM",false,"3M");
  failed |= check_sam_cigar("XIX",false,"1M1I1M");

  failed |= check_sam_cigar("XM",true,"1X1=");

  return failed;
}
