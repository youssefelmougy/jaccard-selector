/*
 *  RICE University
 *  Habanero Team
 *  
 *  This file is part of HC Test.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "hclib.h"

void async_fct(void * arg) {
    void ** argv = (void **) arg;
    int index = *((int *) argv[0]);
    hclib_future_t *dependent_on = (hclib_future_t *)argv[1];
    hclib_promise_t *put_on = (hclib_promise_t *)argv[2];

    int *prev = (int *)hclib_future_get(dependent_on);
    assert(prev == NO_DATUM);

    /* Check value set by predecessor */
    printf("Async %d putting\n", index);
    hclib_promise_put(put_on, NO_DATUM);
    free(argv);
}

void entrypoint(void *arg) {
    hclib_start_finish();
    int n = 5;
    int index = 0;
    // Create asyncs
    // Building 'n' NULL-terminated lists of a single promise each
    // n + 1 to prevent an out of bound error on the final put
    hclib_promise_t **promise_list = (hclib_promise_t **)malloc(sizeof(hclib_promise_t *) * (n + 1));

    for (index = 0 ; index <= n; index++) {
        promise_list[index] = hclib_promise_create();
    }

    for(index = n - 1; index >= 1; index--) {
        // Build async's arguments
        // Pass down the whole promise_list, and async uses index*2 to resolve promises it needs
        void ** argv = (void **)malloc(sizeof(void *) * 3);
        hclib_future_t *fut = hclib_get_future_for_promise(promise_list[index - 1]);

        argv[0] = malloc(sizeof(int));
        *((int *)argv[0]) = index;
        argv[1] = (void *)fut;
        argv[2] = (void *)promise_list[index];

        hclib_async(async_fct, argv, &fut, 1, ANY_PLACE);
    }

    printf("Putting in promise 0\n");
    hclib_promise_put(promise_list[0], NO_DATUM);
    hclib_end_finish();
    // freeing everything up
    for (index = 0 ; index <= n; index++) {
        hclib_promise_free(promise_list[index]);
    }
    free(promise_list);
}

/*
 * Create async await and enable them (by a put) in the 
 * reverse order they've been created.
 */
int main(int argc, char ** argv) {
    setbuf(stdout,NULL);
    char const *deps[] = { "system" };
    hclib_launch(entrypoint, NULL, deps, 1);
    printf("Exiting...\n");
    return 0;
}
