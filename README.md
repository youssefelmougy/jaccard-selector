# Asynchronous Distributed Actor-based Approach to Jaccard Similarity for Genome Comparisons

## Abstract

The computation of genome similarity is important in computational biology applications, and is assessed by calculating the Jaccard similarity of DNA sequencing sets. However, it's challenging to find solutions that can compute Jaccard similarity with the efficiency and scalability needed to fully utilize capabilities of modern HPC hardware. We introduce a novel approach for computing Jaccard similarity for genome comparisons, founded on an actor-based programming model. Our algorithm takes advantage of fine-grained asynchronous computations, distributed/shared memory, and the Fine-grained Asynchronous Bulk-Synchronous Parallelism execution model. Our performance results on the NERSC Perlmutter supercomputer demonstrate that this approach scales to 16,384 cores, showing an average of 4.94x and 5.5x improvement in execution time at the largest scale and relevant hardware performance monitors at medium scale compared to a state-of-the-art baseline. Our approach is also able to process much larger scale genomic datasets than this baseline.

The code for our algorithm can be found [here](https://github.com/youssefelmougy/jaccard-selector/blob/main/Selector/hclib/modules/bale_actor/jaccard-selector/jaccard_kmer_selector.cpp) (in `/hclib/modules/bale_actor/jaccard-selector/`).

## Installation Instructions

The following installation instructions are for the Perlmutter supercomputer at the National Energy Research Scientific Computing Center (NERSC). 

### Load the appropriate modules to prepare for setup

This loads the modules for both Selector and GenomeAtScale to prepare for setup.

```
source scripts/modules.sh
```

### First time setup and installation

This sets up and installs both the Selector and GenomeAtScale applications and their backend runtimes.

```
source scripts/setup.sh
```

## Running Instructions

The following running instructions are for the Perlmutter supercomputer at the National Energy Research Scientific Computing Center (NERSC). 

The run script (`/scripts/run.sh`) has 4 options:

```
   source /scripts/run.sh [selector | ctf | both] [1...inf] [1...inf] [0...5]
   
   [selector | ctf | both]             Selects which application (or both) to run
   [1...inf]                           Selects the number of cores for the run ([32...16384] used in the paper)
   [1...inf]                           Selects the number of nodes for the run ([1...512] used in the paper)
   [0...5]                             Selects the set of HWPC to collect (0:none, 1:L1DA/L1DM/L1IA/L1IM, 2:L2DR/L2DM/L2IR/L2IM, 3:TLBDM/TLBIM, 4:BRINS/BRMSP, 5:INS/CYC)
```

Note: when selecting the number of nodes for the run, please remember that GenomeAtScale uses 32 cores/node and Selector uses either 32 or 64 cores/node.

For example, `source /scripts/run.sh selector 1024 16 2` will run an experiment for the Selector application using 1024 cores on 16 nodes, collecting L2 cache statistics.

This will submit an sbatch file to the run queue at Perlmutter. At job completion, a `jaccard_selector.out` or `jaccard_ctf.out` or both will be created, showing the CMD output of the run. Moreover, if HWPC were collected, a directory with the structure `jaccard*+pat+*` will be created in `/Selector/hclib/modules/bale_actor/jaccard-selector/` or `/GenomeAtScale/jaccard-ctf/` or both. Please see the Output Interpretation section for instructions on how to understand these results.

## Output Interpretation

The following instructions are for understanding the results and relating them to the results found in the paper.

At the completion of each run, there are two outputs that are created:

```
jaccard_selector.out OR jaccard_ctf.out OR both          Output file from submitted job
jaccard_selector+pat+* OR jaccard+pat+* OR both          Output folder (in respective directory) from a CrayPat run if HWPC were collected
```

The `*.out` files contain the execution times of the run for the specific version. This result directly relates to Figure 2 (q) in the paper. An example output is shown below, where `0.06150 seconds` would be reported as the resulting value for the run.

```
...
Running jaccard on 128 threads
K-mer Matrix is 15000x5000 and has 15248 nonzeros.

Jaccard Similarity Matrix is 5000x5000 and has 12497374 values.

Running Jaccard Similarity K-mers (selector): 
  0.06150 seconds
...
```

The `jaccard*+pat+*` folders contain information dumped by the CrayPat profiler (for more information see https://docs.nersc.gov/tools/performance/craypat/). To generate human-readable content, we run `pat_report` on the respective directory. This will display information of interest for the specified HWPC in the run, and will directly relate to Figures 2 (a-p). An example output is shown below, where we can see the L1 cache statistics which would be reported as the resulting values for the run.

```
user@perlmutter: ~> pat_report $PWD/Selector/hclib/modules/bale_actor/jaccard-selector/jaccard_kmer_selector+pat+190420-8718377t
  CrayPat/X:  Version 23.02.0 Revision a53634a72  01/11/23 17:17:09

  Number of PEs (MPI ranks):   128

  Numbers of PEs per Node:      64  PEs on each of  2  Nodes

  Numbers of Threads per PE:     2

  Number of Cores per Socket:   64

  Execution start time:  Sun Mar 19 10:25:36 2023

  System name and speed:  nid004836  2.552 GHz (nominal)

  AMD   Milan                CPU  Family: 25  Model:  1  Stepping:  1

  Core Performance Boost:  256 PEs have CPB capability


  Current path to data file:
  /Selector/hclib/modules/bale_actor/jaccard-selector/jaccard_kmer_selector+pat+190420-8718377t   (RTS, 2 data files)

  ...
  ...

  Processing step 7 of 10
  Notes for table 5:
  ...
  ...
  ==============================================================================
  USER / #1.selector_jaccard
  ------------------------------------------------------------------------------
  Time%                                                  2.8% 
  Time                                               0.060836 secs
  Imb. Time                                          0.000013 secs
  Imb. Time%                                             0.0% 
  Calls                          16.438 /sec              1.0 calls
  PAPI_L1_DCM                     0.057G/sec    2,369,390.898 misses
  PAPI_L1_DCA                     2.252G/sec  110,478,052.633 refs
  Average Time per Call                              0.060836 secs
  CrayPat Overhead : Time          0.0%                       
  perf::PERF_COUNT_HW_CACHE_L1I:ACCESS              1,214,778 
  perf::PERF_COUNT_HW_CACHE_L1I:MISS                    5,868
  ==============================================================================

  ...
  ...

  Hardware performance counter events:
  PAPI_L1_DCM                           Level 1 data cache misses
  PAPI_L1_DCA                           Level 1 data cache accesses
  perf::PERF_COUNT_HW_CACHE_L1I:ACCESS  Undocumented counter
  perf::PERF_COUNT_HW_CACHE_L1I:MISS    Undocumented counter

  Estimated minimum instrumentation overhead per call of a traced function,
  which was subtracted from the data shown in this report
  (for raw data, use the option:  -s overhead=include):
      Time  0.114  microsecs

  Number of traced functions that were called:  7

  (To see the list, specify:  -s traced_functions=show)
user@perlmutter: ~> 
```

## Top-Level Directory Organization

The folder structure of this repository is as follows:

    .
    ├── Selector                                             # Contains files for the Actor-based runtime and the Jaccard k-mer Selector application
    │   ├── hclib                                            # Contains the HClib library and the  Actor-based runtime
    │   │   ├── ...                                         
    │   └── ─── modules                         
    │   │   │   ├── ...                             
    │   └── ─── ─── bale_actor                       
    │   │   │   │   ├── ...                                
    │   └── ─── ─── ─── jaccard-selector                     # Contains the Jaccard k-mer Selector application files
    │   │   │   │   ├── jaccard_kmer_selector.cpp            # Application code for Selector version of Jaccard similarity for genome comparisons
    │   │   │   │   ├── jaccard_kmer_locality_selector.cpp   # Application code for locality-aware Selector version of Jaccard similarity for genome comparisons
    │   │   │   │   ├── kmer_matrix.mtx                      # K-mer matrix file for evaluation
    │   └── ─── ─── ─── ...                             
    ├── GenomeAtScale                                        # Contains files for the CTF library and the GenomeAtScale application
    │   ├── ctf                                              # Contains the CTF library
    │   │   ├── ...                               
    │   ├── jaccard-ctf                                      # Contains the GenomeAtScale (jaccard-ctf) files
    │   │   ├── jaccard.cxx                                  # Application code for GenomeAtScale
    │   │   ├── kmer_matrix.zip                              # K-mer matrix files for evaluation
    │   └── ─── ...                                    
    ├── scripts                                              # Contains installation, running, and modules scripts and sample Perlmutter sbatch files
    │   ├── setup.sh                                         # Installation and build script for the system backends and application code for both the Selector application and the GenomeAtScale application
    │   ├── run.sh                                           # Run script for both the selector application and GenomeAtScale application
    │   ├── modules.sh                                       # Modules script to prepare for running experiments (only used following first time setup using setup.sh, has to be re-run everytime you login to a cluster/supercomputer)
    │   └── ...                                          
    └── README.md

## Citation

If you use our application in your work, please cite [our paper]().

> Youssef Elmougy, Akhiro Hayashi, and Vivek Sarkar. 2024. Asynchronous Distributed Actor-based Approach to Jaccard Similarity for Genome Comparisons.

Corresponding author: Youssef Elmougy ([yelmougy3@gatech.edu](mailto:yelmougy3@gatech.edu))

## Acknowledgement

This research is based upon work supported by the Office of the Director of National Intelligence (ODNI), Intelligence Advanced Research Projects Activity (IARPA), through the Advanced Graphical Intelligence Logical Computing Environment (AGILE) research program, under Army Research Office (ARO) contract number W911NF22C0083. The views and conclusions contained herein are those of the authors and should not be interpreted as necessarily representing the official policies or endorsements, either expressed or implied, of the ODNI, IARPA, or the U.S. Government.

