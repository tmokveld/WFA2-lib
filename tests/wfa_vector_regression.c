/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

#include "utils/vector.h"

static int check_zero_capacity_append(void) {
  vector_t* const vector = vector_new(0,uint64_t);
  if (vector->memory != NULL || vector->used != 0 ||
      vector->elements_allocated != 0 ||
      vector->element_size != sizeof(uint64_t)) {
    fprintf(stderr,"Unexpected zero-capacity vector state\n");
    vector_delete(vector);
    return 1;
  }

  vector_insert(vector,42,uint64_t);
  if (vector->memory == NULL || vector_get_used(vector) != 1 ||
      vector->elements_allocated < 1 ||
      *vector_get_elm(vector,0,uint64_t) != 42) {
    fprintf(stderr,"Append after zero-capacity construction failed\n");
    vector_delete(vector);
    return 1;
  }

  vector_delete(vector);
  return 0;
}

static int check_reserve_from_null(void) {
  vector_t* const vector = vector_new(0,uint32_t);
  vector_reserve(vector,4,true);

  if (vector->memory == NULL || vector->elements_allocated < 4) {
    fprintf(stderr,"Reserve from NULL memory failed\n");
    vector_delete(vector);
    return 1;
  }

  uint32_t* const memory = vector_get_mem(vector,uint32_t);
  for (uint64_t i=0;i<4;++i) {
    if (memory[i] != 0) {
      fprintf(stderr,"Reserve zero_mem did not clear unused element %" PRIu64 "\n",i);
      vector_delete(vector);
      return 1;
    }
  }

  vector_delete(vector);
  return 0;
}

static int check_self_copy(void) {
  vector_t* const vector = vector_new(0,uint64_t);
  vector_insert(vector,11,uint64_t);
  vector_insert(vector,22,uint64_t);

  vector_copy(vector,vector);
  if (vector_get_used(vector) != 2 ||
      *vector_get_elm(vector,0,uint64_t) != 11 ||
      *vector_get_elm(vector,1,uint64_t) != 22) {
    fprintf(stderr,"Self-copy changed vector contents\n");
    vector_delete(vector);
    return 1;
  }

  vector_delete(vector);
  return 0;
}

static int check_if_else_macros(void) {
  vector_t* const vector = vector_new(0,uint32_t);
  int condition = 1;

  if (condition)
    vector_insert(vector,7,uint32_t);
  else
    vector_insert(vector,8,uint32_t);

  if (vector_get_used(vector) != 1 ||
      *vector_get_elm(vector,0,uint32_t) != 7) {
    fprintf(stderr,"vector_insert failed in if/else context\n");
    vector_delete(vector);
    return 1;
  }

  if (condition)
    vector_prepare(vector,2,uint16_t);
  else
    vector_clear(vector);

  if (vector_get_used(vector) != 0 ||
      vector->element_size != sizeof(uint16_t)) {
    fprintf(stderr,"vector_prepare failed in if/else context\n");
    vector_delete(vector);
    return 1;
  }

  uint16_t* element = NULL;
  if (condition)
    vector_alloc_new(vector,uint16_t,element);
  else
    vector_clear(vector);
  *element = 13;

  vector_reserve(vector,4,false);
  vector_push_unsafe(vector,21,uint16_t);
  vector_alloc_new_unsafe(vector,uint16_t,element);
  *element = 34;

  if (vector_get_used(vector) != 3 ||
      *vector_get_elm(vector,0,uint16_t) != 13 ||
      *vector_get_elm(vector,1,uint16_t) != 21 ||
      *vector_get_elm(vector,2,uint16_t) != 34) {
    fprintf(stderr,"Append helper contents are incorrect\n");
    vector_delete(vector);
    return 1;
  }

  vector_delete(vector);
  return 0;
}

int main(void) {
  int failed = 0;
  failed |= check_zero_capacity_append();
  failed |= check_reserve_from_null();
  failed |= check_self_copy();
  failed |= check_if_else_macros();
  return failed;
}
