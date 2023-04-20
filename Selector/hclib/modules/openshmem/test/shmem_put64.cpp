#include "hclib_cpp.h"
#include "hclib_openshmem.h"
#include "hclib_system.h"

#include <iostream>

int main(int argc, char **argv) {
    const char *deps[] = { "system", "openshmem" };
    hclib::launch(deps, 2, [] {
        int pe = hclib::shmem_my_pe();
        std::cout << "Hello world from rank " << pe << std::endl;

        uint64_t *shared = (uint64_t *)hclib::shmem_malloc(sizeof(uint64_t));
        *shared = pe;

        hclib::shmem_barrier_all();

        uint64_t source = pe;
        int target = pe - 1;
        if (target == -1) target = hclib::shmem_n_pes() - 1;

        hclib::shmem_put64(shared, &source, 1, target);

        hclib::shmem_barrier_all();

        assert(*shared == ((source + 1) % hclib::shmem_n_pes()));
        fprintf(stderr, "PE %d done, got %lu\n", pe, *shared);

    });
    return 0;
}
