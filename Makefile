##############################################################################
## This software module is developed by SciDM (Scientific Data Management) in 1998-2015
## 
## This program is free software; you can redistribute, reuse,
## or modify it with no restriction, under the terms of the MIT License.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## 
## For any questions please contact Denis Kaznadzey at dkaznadzey@yahoo.com
##############################################################################

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
