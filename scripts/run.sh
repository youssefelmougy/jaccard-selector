#!/bin/bash

# FORMAT: ./run.sh [selector | ctf | both] [1...inf] [1...inf] [0...5]
#
#           where
#                [selector | ctf | both]             Selects which application (or both) to run
#                [1...inf]                           Selects the number of cores for the run ([32...16384] used in the paper)
#                [1...inf]                           Selects the number of nodes for the run ([1...512] used in the paper)
#                [0...5]                             Selects the set of HWPC to collect (0:none, 1:L1DA/L1DM/L1IA/L1IM, 2:L2DR/L2DM/L2IR/L2IM, 3:TLBDM/TLBIM, 4:BRINS/BRMSP, 5:INS/CYC)


# Capture flags
#while getopts v:n:h:N: flag
#do
#    case "${flag}" in
#        v) application=${OPTARG};;
#        n) cores=${OPTARG};;
#        h) hw=${OPTARG};;
#        N) nodes=${OPTARG};;
#    esac
#done

#export CORES=$cores
#export NODES=$nodes
#export HWPC=$hw
application=$1
CORES=$2
NODES=$3
HWPC=$4
echo $HWPC > hwpc.txt

# Run sbatch files corresponding to application version from  user
if [ $application = "selector" ] ; then
  echo "Running evaluations of the Selector application using $CORES cores and $NODES nodes"
  sbatch -N $NODES -n $CORES $PWD/scripts/jaccard_selector.sbatch
elif [ $application = "ctf" ] ; then
  echo "Running evaluations of the GenomeAtScale application using $CORES cores and $NODES nodes"
  sbatch -N $NODES -n $CORES $PWD/scripts/jaccard_ctf.sbatch
else
  echo "Running evaluations of both the Selector application and the GenomeAtScale application using $CORES cores and $NODES nodes"
  sbatch -N $NODES -n $CORES $PWD/scripts/jaccard_selector.sbatch
  sbatch -N $NODES -n $CORES $PWD/scripts/jaccard_ctf.sbatch
fi

echo "---------------------------------------------------------------------------------------------------------------"
echo "Job(s) submitted! Please check job status for completion and see README for output interpretation instructions."
echo "---------------------------------------------------------------------------------------------------------------"


