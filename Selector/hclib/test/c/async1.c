/* Copyright (c) 2013, Rice University

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

/**
 * DESC: Fork a bunch of asyncs in a top-level loop
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "hclib.h"

#define NB_ASYNC 127

int ran[NB_ASYNC];

void async_fct(void * arg) {
    int idx = *((int *) arg);
    assert(ran[idx] == -1);
    ran[idx] = idx;
}

void entrypoint(void *arg) {
    int i = 0;
    // This is ok to have these on stack because this
    // code is alive until the end of the program.
    int indices [NB_ASYNC];
    for (i = 0; i < NB_ASYNC; i++) {
        ran[i] = -1;
    }

    hclib_start_finish();

    for (i = 0; i < NB_ASYNC; i++) {
        indices[i] = i;
        //Note: Forcefully pass the address we want to write to as a void **
        hclib_async(async_fct, (void*) (indices+i), NO_FUTURE, 0,
                ANY_PLACE);
    }

    hclib_end_finish();

    printf("Call Finalize\n");
}

int main (int argc, char ** argv) {
    printf("Call Init\n");
    char const *deps[] = { "system" };
    hclib_launch(entrypoint, NULL, deps, 1);
    printf("Check results: ");
    int i = 0;
    while(i < NB_ASYNC) {
        assert(ran[i] == i);
        i++;
    }
    printf("OK\n");
    return 0;
}
