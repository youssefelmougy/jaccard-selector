#!/bin/bash
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH -t 0:10:00
#SBATCH --output=jaccard_selector.out

# Initialize modules
# Always run modules.sh first to make sure everything is loaded
source $PWD/scripts/modules.sh

HWPC=$(< $PWD/hwpc.txt)

cd $PWD/Selector/hclib/modules/bale_actor/jaccard-selector/

# Export corresponding HWPC set specified by user
counters=" "
if [ $HWPC -eq 0 ] ; then
  unset PAT_RT_PERFCTR
elif [ $HWPC -eq 1 ] ; then
  export PAT_RT_PERFCTR="PAPI_L1_DCA,PAPI_L1_DCM,perf::PERF_COUNT_HW_CACHE_L1I:ACCESS,perf::PERF_COUNT_HW_CACHE_L1I:MISS"
  counters="L1 cache"
elif [ $HWPC -eq 2 ] ; then
  export PAT_RT_PERFCTR="PAPI_L2_DCR,PAPI_L2_DCM,PAPI_L2_ICR,PAPI_L2_ICM"
  counters="L2 cache"
elif [ $HWPC -eq 3 ] ; then
  export PAT_RT_PERFCTR="PAPI_TLB_DM,PAPI_TLB_IM"
  counters="TLB"
elif [ $HWPC -eq 4 ] ; then
  export PAT_RT_PERFCTR="PAPI_BR_INS,PAPI_BR_MSP"
  counters="Branches"
elif [ $HWPC -eq 5 ] ; then
  export PAT_RT_PERFCTR="PAPI_TOT_INS,PAPI_TOT_CYC"
  counters="Instructions and Cycles"
fi

# Compile application
rm -rf jaccard_kmer_selector jaccard_kmer_selector+pat
make jaccard_kmer_selector
pat_build jaccard_kmer_selector

# Run  application
if [ $HWPC -eq 0 ] ; then
  echo "Running Selector Jaccard k-mer application using $CORES cores on $NODES nodes (no HWPC collected)"
  srun -N $NODES -n $CORES -c 4 --cpu-bind=cores ./jaccard_kmer_selector
else
  echo "Running Selector Jaccard k-mer application using $CORES cores on $NODES nodes (HWPC for $counters collected)"
  srun -N $NODES -n $CORES -c 4 --cpu-bind=cores ./jaccard_kmer_selector+pat
fi

echo "-------------------------------------------------------------------------"
echo "Completed..."
echo "-------------------------------------------------------------------------"

