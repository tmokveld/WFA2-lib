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
#include <stdlib.h>
#include <time.h>

#include "utils/vector.h"

#ifndef vector_push_unsafe
#define vector_push_unsafe(vector,element,type) do { \
  vector_get_mem((vector),type)[(vector)->used++] = (element); \
} while (0)
#endif

static double elapsed_seconds(
    const clock_t begin,
    const clock_t end) {
  return (double)(end-begin)/(double)CLOCKS_PER_SEC;
}

static uint64_t checksum_vector(
    vector_t* const vector) {
  uint64_t checksum = 0;
  uint64_t* const memory = vector_get_mem(vector,uint64_t);
  for (uint64_t i=0;i<vector_get_used(vector);++i) {
    checksum += memory[i];
  }
  return checksum;
}

static uint64_t bench_vector_insert(
    const uint64_t num_elements,
    double* const seconds) {
  vector_t* const vector = vector_new(0,uint64_t);
  const clock_t begin = clock();
  for (uint64_t i=0;i<num_elements;++i) {
    vector_insert(vector,i,uint64_t);
  }
  const clock_t end = clock();
  *seconds = elapsed_seconds(begin,end);
  const uint64_t checksum = checksum_vector(vector);
  vector_delete(vector);
  return checksum;
}

static uint64_t bench_vector_alloc_new(
    const uint64_t num_elements,
    double* const seconds) {
  vector_t* const vector = vector_new(0,uint64_t);
  const clock_t begin = clock();
  for (uint64_t i=0;i<num_elements;++i) {
    uint64_t* element = NULL;
    vector_alloc_new(vector,uint64_t,element);
    *element = i;
  }
  const clock_t end = clock();
  *seconds = elapsed_seconds(begin,end);
  const uint64_t checksum = checksum_vector(vector);
  vector_delete(vector);
  return checksum;
}

static uint64_t bench_push_unsafe(
    const uint64_t num_elements,
    double* const seconds) {
  vector_t* const vector = vector_new(0,uint64_t);
  vector_reserve(vector,num_elements,false);
  const clock_t begin = clock();
  for (uint64_t i=0;i<num_elements;++i) {
    vector_push_unsafe(vector,i,uint64_t);
  }
  const clock_t end = clock();
  *seconds = elapsed_seconds(begin,end);
  const uint64_t checksum = checksum_vector(vector);
  vector_delete(vector);
  return checksum;
}

static uint64_t bench_raw_pointer(
    const uint64_t num_elements,
    double* const seconds) {
  vector_t* const vector = vector_new(0,uint64_t);
  vector_reserve(vector,num_elements,false);
  uint64_t* pointer = vector_get_mem(vector,uint64_t)+vector->used;
  const clock_t begin = clock();
  for (uint64_t i=0;i<num_elements;++i) {
    *pointer++ = i;
  }
  vector->used = (uint64_t)(pointer-vector_get_mem(vector,uint64_t));
  const clock_t end = clock();
  *seconds = elapsed_seconds(begin,end);
  const uint64_t checksum = checksum_vector(vector);
  vector_delete(vector);
  return checksum;
}

int main(
    int argc,
    char** argv) {
  const uint64_t num_elements =
      argc > 1 ? (uint64_t) strtoull(argv[1],NULL,10) : 1000000;

  double insert_seconds = 0.0;
  double alloc_new_seconds = 0.0;
  double unsafe_seconds = 0.0;
  double raw_seconds = 0.0;
  const uint64_t insert_checksum =
      bench_vector_insert(num_elements,&insert_seconds);
  const uint64_t alloc_new_checksum =
      bench_vector_alloc_new(num_elements,&alloc_new_seconds);
  const uint64_t unsafe_checksum =
      bench_push_unsafe(num_elements,&unsafe_seconds);
  const uint64_t raw_checksum =
      bench_raw_pointer(num_elements,&raw_seconds);

  printf("elements,method,seconds,checksum\n");
  printf("%" PRIu64 ",vector_insert,%.6f,%" PRIu64 "\n",
      num_elements,insert_seconds,insert_checksum);
  printf("%" PRIu64 ",vector_alloc_new,%.6f,%" PRIu64 "\n",
      num_elements,alloc_new_seconds,alloc_new_checksum);
  printf("%" PRIu64 ",vector_push_unsafe,%.6f,%" PRIu64 "\n",
      num_elements,unsafe_seconds,unsafe_checksum);
  printf("%" PRIu64 ",raw_pointer,%.6f,%" PRIu64 "\n",
      num_elements,raw_seconds,raw_checksum);
  return insert_checksum != alloc_new_checksum ||
      insert_checksum != unsafe_checksum ||
      insert_checksum != raw_checksum;
}
