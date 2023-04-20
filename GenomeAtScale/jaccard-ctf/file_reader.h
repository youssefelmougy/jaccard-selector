#ifndef __GRAPH_AUX_H__
#define __GRAPH_AUX_H__

#include <ctf.hpp>
#include <float.h>
#define __STDC_FORMAT_MACROS 1
#include <inttypes.h>

uint64_t read_file(int myid, int ntask, const char *fpath, uint64_t **edge);
uint64_t read_file_mpiio(int myid, int ntask, const char *fpath, uint64_t **edge, char ***led);
void process_data(char **led, uint64_t ned, int myid, uint64_t **edge);
#endif

