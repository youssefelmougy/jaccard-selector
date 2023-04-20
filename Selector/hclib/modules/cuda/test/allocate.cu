#include "hclib_cpp.h"
#include "hclib_cuda.h"
#include "hclib_system.h"

#include <iostream>

int main(int argc, char **argv) {
    const char *deps[] = { "system", "cuda" };
    hclib::launch(deps, 2, [] {
        const int N = 1024;
        int gpu_id = hclib::get_gpu_locale_id();
        int num_gpus = hclib::get_num_gpu_locales();
        assert(num_gpus > 0);

        hclib::locale_t *gpu_locale = hclib::get_closest_gpu_locale();
        std::cout << "Hello world, gpu locale id = " << gpu_id <<
            ", closest gpu = " << gpu_locale << " (" <<
            hclib::get_gpu_name(gpu_locale) << ")" << std::endl;
        assert(gpu_locale);

        int ngpus;
        hclib::locale_t **gpu_locales = hclib::get_gpu_locales(&ngpus);
        std::cout << "Got " << ngpus << " GPU locale(s)" << std::endl;

        hclib::future_t<void*> *fut = hclib::allocate_at(N * sizeof(int), gpu_locale);
        void *alloc = fut->wait();
        assert(alloc);

        hclib::locale_t *cpu_locale = hclib::get_closest_cpu_locale();
        fut = hclib::allocate_at(N * sizeof(int), cpu_locale);
        int *host_alloc = (int *)fut->wait();
        assert(host_alloc);

        for (int i = 0; i < N; i++) host_alloc[i] = i;

        fut = hclib::async_copy(gpu_locale, alloc, cpu_locale,
                host_alloc, N * sizeof(int));
        fut->wait();

        fut = hclib::allocate_at(N * sizeof(int), cpu_locale);
        int *out_host_alloc = (int *)fut->wait();
        assert(out_host_alloc);

        fut = hclib::async_copy(cpu_locale, out_host_alloc, gpu_locale, alloc,
                N * sizeof(int));
        fut->wait();

        for (int i = 0; i < N; i++) {
            assert(out_host_alloc[i] == i);
        }

        fut = hclib::memset_at(alloc, 0x00, N * sizeof(int), gpu_locale);
        fut->wait();
        fut = hclib::async_copy(cpu_locale, out_host_alloc, gpu_locale, alloc,
                N * sizeof(int));
        fut->wait();

        for (int i = 0; i < N; i++) {
            assert(out_host_alloc[i] == 0);
        }
    });
    return 0;
}
