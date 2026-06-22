/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include "wavefront/wavefront_align.h"

static const char* const negative_match_message =
    "[WFA] BiWFA ends-free with negative match rewards (match < 0) "
    "requires aligned-length-aware breakpoint scoring and is not implemented yet\n";
static const char* const extension_message =
    "[WFA] BiWFA extension is not implemented yet (let me know and I'll add it)\n";

typedef void (*child_test_fn)(void);

static wavefront_aligner_attr_t make_biwfa_attrs(
    const alignment_scope_t scope,
    const int match) {
  wavefront_aligner_attr_t attrs = wavefront_aligner_attr_default;
  attrs.memory_mode = wavefront_memory_ultralow;
  attrs.alignment_scope = scope;
  attrs.distance_metric = gap_affine;
  attrs.affine_penalties.match = match;
  attrs.affine_penalties.mismatch = 1;
  attrs.affine_penalties.gap_opening = 2;
  attrs.affine_penalties.gap_extension = 1;
  return attrs;
}

static void set_attr_endsfree(
    wavefront_aligner_attr_t* const attrs,
    const int pattern_begin_free,
    const int pattern_end_free,
    const int text_begin_free,
    const int text_end_free) {
  attrs->alignment_form.span = alignment_endsfree;
  attrs->alignment_form.extension = false;
  attrs->alignment_form.pattern_begin_free = pattern_begin_free;
  attrs->alignment_form.pattern_end_free = pattern_end_free;
  attrs->alignment_form.text_begin_free = text_begin_free;
  attrs->alignment_form.text_end_free = text_end_free;
}

static int read_stderr(
    const int fd,
    char* const buffer,
    const size_t capacity) {
  size_t used = 0;
  while (used+1 < capacity) {
    const ssize_t count = read(fd,buffer+used,capacity-used-1);
    if (count > 0) {
      used += (size_t)count;
    } else if (count == 0) {
      break;
    } else if (errno != EINTR) {
      perror("read");
      return 1;
    }
  }
  buffer[used] = '\0';
  return 0;
}

static int run_child_test(
    const char* const label,
    child_test_fn const test_fn,
    const int expected_exit_status,
    const char* const expected_stderr) {
  int pipe_fds[2];
  if (pipe(pipe_fds) != 0) {
    perror("pipe");
    return 1;
  }
  const pid_t pid = fork();
  if (pid < 0) {
    perror("fork");
    close(pipe_fds[0]);
    close(pipe_fds[1]);
    return 1;
  }
  if (pid == 0) {
    close(pipe_fds[0]);
    if (dup2(pipe_fds[1],STDERR_FILENO) < 0) {
      perror("dup2");
      _exit(127);
    }
    close(pipe_fds[1]);
    test_fn();
    fflush(stderr);
    _exit(0);
  }

  close(pipe_fds[1]);
  char stderr_buffer[4096];
  const int read_failed = read_stderr(
      pipe_fds[0],stderr_buffer,sizeof(stderr_buffer));
  close(pipe_fds[0]);

  int status = 0;
  if (waitpid(pid,&status,0) < 0) {
    perror("waitpid");
    return 1;
  }
  if (read_failed) return 1;

  if (!WIFEXITED(status) ||
      WEXITSTATUS(status) != expected_exit_status) {
    fprintf(stderr,"%s: expected exit %d, observed status 0x%x\n",
        label,expected_exit_status,status);
    fprintf(stderr,"%s: stderr: %s",label,stderr_buffer);
    return 1;
  }
  if (expected_stderr != NULL &&
      strcmp(stderr_buffer,expected_stderr) != 0) {
    fprintf(stderr,"%s: unexpected stderr\nExpected: %sObserved: %s",
        label,expected_stderr,stderr_buffer);
    return 1;
  }
  return 0;
}

static void construct_negative_endsfree_score_only(void) {
  wavefront_aligner_attr_t attrs = make_biwfa_attrs(compute_score,-5);
  set_attr_endsfree(&attrs,1,0,0,0);
  wavefront_aligner_new(&attrs);
}

static void construct_negative_endsfree_full_alignment(void) {
  wavefront_aligner_attr_t attrs = make_biwfa_attrs(compute_alignment,-5);
  set_attr_endsfree(&attrs,0,1,0,0);
  wavefront_aligner_new(&attrs);
}

static void align_negative_endsfree_score_only(void) {
  wavefront_aligner_attr_t attrs = make_biwfa_attrs(compute_score,-5);
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attrs);
  wavefront_aligner_set_alignment_free_ends(wf_aligner,0,0,1,0);
  wavefront_align(wf_aligner,"A",1,"A",1);
}

static void align_negative_endsfree_full_alignment(void) {
  wavefront_aligner_attr_t attrs = make_biwfa_attrs(compute_alignment,-5);
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attrs);
  wavefront_aligner_set_alignment_free_ends(wf_aligner,0,0,0,1);
  wavefront_align(wf_aligner,"A",1,"A",1);
}

static void construct_extension(void) {
  wavefront_aligner_attr_t attrs = make_biwfa_attrs(compute_alignment,0);
  attrs.alignment_form.span = alignment_endsfree;
  attrs.alignment_form.extension = true;
  wavefront_aligner_new(&attrs);
}

static void align_extension(void) {
  wavefront_aligner_attr_t attrs = make_biwfa_attrs(compute_alignment,0);
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attrs);
  wavefront_aligner_set_alignment_extension(wf_aligner);
  wavefront_align(wf_aligner,"A",1,"A",1);
}

static void align_negative_zero_free_ends(void) {
  wavefront_aligner_attr_t attrs = make_biwfa_attrs(compute_score,-5);
  set_attr_endsfree(&attrs,0,0,0,0);
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attrs);
  const int status = wavefront_align(wf_aligner,"ACGT",4,"ACCT",4);
  if (status != WF_STATUS_ALG_COMPLETED) {
    fprintf(stderr,"global-like ends-free alignment returned status %d\n",status);
    exit(2);
  }
  wavefront_aligner_delete(wf_aligner);
}

int main(void) {
  int failed = 0;
  failed |= run_child_test(
      "construct-negative-score-only",
      construct_negative_endsfree_score_only,1,negative_match_message);
  failed |= run_child_test(
      "construct-negative-full-alignment",
      construct_negative_endsfree_full_alignment,1,negative_match_message);
  failed |= run_child_test(
      "align-negative-score-only",
      align_negative_endsfree_score_only,1,negative_match_message);
  failed |= run_child_test(
      "align-negative-full-alignment",
      align_negative_endsfree_full_alignment,1,negative_match_message);
  failed |= run_child_test(
      "construct-extension",
      construct_extension,1,extension_message);
  failed |= run_child_test(
      "align-extension",
      align_extension,1,extension_message);
  failed |= run_child_test(
      "align-negative-zero-free-ends",
      align_negative_zero_free_ends,0,NULL);
  return failed;
}
