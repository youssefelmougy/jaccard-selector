include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.post.mak

SRUN ?= oshrun

TARGETS= jaccard_kmer_selector jaccard_kmer_locality_selector

all: $(TARGETS)

%: %.cpp
	$(CXX) -g -O3 -std=c++11 -DUSE_SHMEM=1 $(HCLIB_CFLAGS) $(HCLIB_LDFLAGS) -o $@ $^ $(HCLIB_LDLIBS) -lspmat -lconvey -lexstack -llibgetput -lhclib_bale_actor -lm

test:
	for target in $(TARGETS); do $(SRUN) ./$$target; done

clean:
	rm -f $(TARGETS)
