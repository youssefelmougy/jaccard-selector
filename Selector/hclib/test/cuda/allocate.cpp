#include <iostream>

#include "hclib_cpp.h"

int main(int argc, char **argv) {
    hclib::launch(&argc, argv, []() {
        hclib::place_t *root_pl = hclib::get_root_place();

        int num_toplevel;
        hclib::place_t **toplevel = hclib::get_children_of_place(root_pl,
                &num_toplevel);
        for (int i = 0; i < num_toplevel; i++) {
            if (toplevel[i]->type == NVGPU_PLACE) {
                std::cout << "Found GPU place" << std::endl;
                void *d_ptr = hclib::allocate_at<int>(toplevel[i], 10, 0);
                HASSERT(d_ptr);
                hclib::free_at(toplevel[i], d_ptr);
            } else {
                std::cout << "Found CPU place" << std::endl;
                void *h_ptr = hclib::allocate_at<int>(toplevel[i], 10, PHYSICAL);
                HASSERT(h_ptr);
                hclib::free_at(toplevel[i], h_ptr);
            }
        }
    });
    return 0;
}
