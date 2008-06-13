lib_common := lib/libseq.a
lib_dir := common/
psimscan_dir := psimscan/
psimscan := $(psimscan_dir)psimscan
nsimscan_dir := nsimscan/
nsimscan := $(nsimscan_dir)nsimscan
apps := $(psimscan) $(nsimscan)

.PHONY: all $(lib_common) $(apps)


all: $(apps)


$(psimscan) :
	$(MAKE) --directory=$(psimscan_dir)

$(nsimscan) :
	$(MAKE) --directory=$(nsimscan_dir)

$(lib_common) :
	$(MAKE) --directory=$(lib_dir)
	

$(apps): $(lib_common)

clean:
	for d in $(psimscan_dir) $(nsimscan_dir) $(lib_dir); \
	do \
	  $(MAKE) --directory=$$d clean ; \
	done \
