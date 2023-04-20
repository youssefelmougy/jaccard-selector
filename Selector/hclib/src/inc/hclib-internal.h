/* Copyright (c) 2015, Rice University

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1.  Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.
3.  Neither the name of Rice University
     nor the names of its contributors may be used to endorse or
     promote products derived from this software without specific
     prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

/*
 * hclib-internal.h
 *  
 *      Author: Vivek Kumar (vivekk@rice.edu)
 *      Acknowledgments: https://wiki.rice.edu/confluence/display/HABANERO/People
 */

#ifndef HCLIB_INTERNAL_H_
#define HCLIB_INTERNAL_H_

#include <stdarg.h>
#include <stdint.h>
#include "hclib-tree.h"
#include "hclib-deque.h"
#include "hclib.h"
#include "litectx.h"
#include "hclib-locality-graph.h"

#define LOG_LEVEL_FATAL         1
#define LOG_LEVEL_WARN          2
#define LOG_LEVEL_INFO          3
#define LOG_LEVEL_DEBUG         4
#define LOG_LEVEL_TRACE         5

/* set the current log level */
#define LOG_LEVEL LOG_LEVEL_FATAL

#define WHEREARG __FILE__,__LINE__

#define LOG_(level, ...) if (level<=LOG_LEVEL) log_(WHEREARG, CURRENT_WS_INTERNAL, __VA_ARGS__);

/* We more or less mimic log4c without the ERROR level */
#define LOG_FATAL(...)  LOG_(LOG_LEVEL_FATAL, __VA_ARGS__)
#define LOG_WARN(...)   LOG_(LOG_LEVEL_WARN,  __VA_ARGS__)
#define LOG_INFO(...)   LOG_(LOG_LEVEL_INFO,  __VA_ARGS__)
#define LOG_DEBUG(...)  LOG_(LOG_LEVEL_DEBUG, __VA_ARGS__)
#define LOG_TRACE(...)  LOG_(LOG_LEVEL_TRACE, __VA_ARGS__)

/* log the msg using the fatal logger and abort the program */
#define log_die(... ) { LOG_FATAL(__VA_ARGS__); abort(); }
#define check_log_die(cond, ... ) if(cond) { log_die(__VA_ARGS__) }

#define CACHE_LINE_L1 8

// Default value of a promise datum
#define UNINITIALIZED_PROMISE_DATA_PTR NULL

// For waiting frontier (last element of the list)
#define SATISFIED_FUTURE_WAITLIST_PTR NULL
#define SENTINEL_FUTURE_WAITLIST_PTR ((void*) -1)

typedef struct {
    volatile int flag;
    int padding;
    void * pad[CACHE_LINE_L1 - 1];
} worker_done_t;

/*
 * Global context information for the HC runtime, shared by all worker threads.
 */
typedef struct hclib_context {
    struct _hclib_worker_state** workers; /* all workers */
    hclib_locality_graph *graph;
    hclib_worker_paths *worker_paths;
    int nworkers; /* # of worker threads created */
    int ncores; /* physical number of cores detected */
    /* a simple implementation of wait/wakeup condition */
    volatile int workers_wait_cond;
    worker_done_t *done_flags;
#ifdef HC_CUDA
    hclib_memory_tree_node *pinned_host_allocs;
    cudaStream_t stream;
#endif
} hclib_context;

#include "hclib-finish.h"

typedef struct _hclib_deque_t {
    /* The actual deque, WARNING: do not move declaration !
     * Other parts of the runtime rely on it being the first one. */
    hclib_internal_deque_t deque;
    struct _hclib_worker_state * ws;
    struct _hclib_deque_t *nnext;
    struct _hclib_deque_t *prev; /* the deque list of the worker */
    hclib_locale_t *locale;
} hclib_deque_t;

void log_(const char * file, int line, hclib_worker_state * ws, const char * format,
        ...);

// promise
int register_on_all_promise_dependencies(hclib_task_t *wrapper_task);
void try_schedule_async(hclib_task_t * async_task, hclib_worker_state *ws);

int static inline _hclib_promise_is_satisfied(hclib_promise_t *p) {
    return p->wait_list_head == SATISFIED_FUTURE_WAITLIST_PTR;
}

#endif /* HCLIB_INTERNAL_H_ */
