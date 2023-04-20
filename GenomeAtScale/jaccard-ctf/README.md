## Build instructions

See the README file at https://github.com/youssefelmougy/jaccard-selector for setup instructions.

## Run instructions

See the README file at https://github.com/youssefelmougy/jaccard-selector for run instructions.

To test with default parameters, just run `./jaccard` or `mpirun -np 4 ./jaccard`.

To run, use e.g. `./jaccard -m 4000 -n 100 -p .01 -nbatch 10`, which would generate a 4000-by-100 k-mer bit matrix with 1% nonzeros, then compute a 100-by-100 similarity matrix by accumulating batches of 400 rows at a time.

