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

#include "utils/vector.h"

static void trigger_out_of_range(void) {
  vector_t* const vector = vector_new(0,uint32_t);
  vector_insert(vector,1,uint32_t);
  (void) vector_get_elm(vector,1,uint32_t);
  vector_delete(vector);
}

static void trigger_type_mismatch(void) {
  vector_t* const vector = vector_new(0,uint32_t);
  vector_insert(vector,1,uint32_t);
  (void) vector_get_elm(vector,0,uint64_t);
  vector_delete(vector);
}

static void trigger_empty_last(void) {
  vector_t* const vector = vector_new(0,uint32_t);
  (void) vector_get_last_elm(vector,uint32_t);
  vector_delete(vector);
}

int main(
    int argc,
    char** argv) {
  if (argc != 2) {
    fprintf(stderr,"Usage: %s out_of_range|type_mismatch|empty_last\n",argv[0]);
    return 2;
  }
  if (strcmp(argv[1],"out_of_range") == 0) {
    trigger_out_of_range();
  } else if (strcmp(argv[1],"type_mismatch") == 0) {
    trigger_type_mismatch();
  } else if (strcmp(argv[1],"empty_last") == 0) {
    trigger_empty_last();
  } else {
    fprintf(stderr,"Unknown vector debug regression mode '%s'\n",argv[1]);
    return 2;
  }
  return 0;
}
