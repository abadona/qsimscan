
##////////////////////////////////////////////////////////////////////////////
## This software module is developed by SCIDM team in 1999-2015.
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
##
## For any questions please contact SciDM team by email at team@scidm.org
##////////////////////////////////////////////////////////////////////////////

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
