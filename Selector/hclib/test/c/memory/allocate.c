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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "hclib.h"

void entrypoint(void *arg) {
    int i;
    int nlocales = hclib_get_num_locales();
    hclib_locale_t *locales = hclib_get_all_locales();

    for (i = 0; i < nlocales; i++) {
        hclib_future_t *future = hclib_allocate_at(30, locales + i);
        void *allocation = hclib_future_wait(future);
        fprintf(stderr, "%d: allocation = %p\n", i, allocation);
        future = hclib_reallocate_at(allocation, 50, locales + i);
        void *reallocation = hclib_future_wait(future);
        fprintf(stderr, "%d: reallocation = %p\n", i, reallocation);
        future = hclib_memset_at(reallocation, 0, 50, locales + i);
        hclib_future_wait(future);

        hclib_free_at(reallocation, locales + i);
    }

    printf("Passed\n");
}

int main (int argc, char ** argv) {
    char const *deps[] = { "system" };
    hclib_launch(entrypoint, NULL, deps, 1);
    fprintf(stderr, "Finished launch\n");
    return 0;
}
