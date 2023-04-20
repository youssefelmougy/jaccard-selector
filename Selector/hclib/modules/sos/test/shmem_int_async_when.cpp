#include "hclib_cpp.h"
#include "hclib_sos.h"
#include "hclib_system.h"

#include <unistd.h>
#include <iostream>

int shared_var = 0;

int main(int argc, char **argv) {
    const char *deps[] = { "system", "sos" };
    hclib::launch(deps, 2, [] {
        const int pe = hclib::shmem_my_pe();
        std::cout << "Hello world from rank " << pe << std::endl;

        hclib::shmem_barrier_all();

        hclib::finish([pe] {

            if (pe < hclib::shmem_n_pes() - 1) {
                const unsigned long long start_wait = hclib_current_time_ns();
                hclib::shmem_int_async_when(&shared_var, SHMEM_CMP_EQ, pe + 1, [start_wait, pe] {
                        assert(shared_var == pe + 1);
                        const unsigned long long end_wait = hclib_current_time_ns();

                        std::cerr << "PE " << pe <<
                            " waited for " << ((end_wait - start_wait) / 1000000) <<
                            " ms" << std::endl;

                        sleep(1);

                        if (pe > 0) {
                            int pe_no = pe;
                            hclib::shmem_int_put(&shared_var, &pe_no, 1, pe - 1);
                        }

                    });
            } else {
                int pe_no = pe;
                hclib::shmem_int_put(&shared_var, &pe_no, 1, pe - 1);
            }
        });

        hclib::shmem_barrier_all();
    });
    return 0;
}
