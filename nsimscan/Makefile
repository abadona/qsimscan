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

.PHONY : clean all

UNAME:=$(shell uname)
ifneq (,$(findstring CYGWIN,$(UNAME)))
  EXESUFFIX:=.exe
else
  EXESUFFIX:=
endif

CC=$(CCPATH)g++
PIC := $(if $(filter $(shell uname -i),x86_64),-fPIC,)
OPTLEVEL := $(if $(OPTLEVEL),$(OPTLEVEL),-O2)
DEBUGLEVEL := $(if $(DEBUGLEVEL),$(DEBUGLEVEL),-g0)
CPPFLAGS=$(PIC) $(DEBUGLEVEL) $(OPTLEVEL) -D__x86__

TARGETNAME=nsimscan$(EXESUFFIX)
TARGETDIR=.
TARGET=$(TARGETDIR)/$(TARGETNAME)

SYSINCLUDE_DIRS=
SYSINCLUDE_DIRS_SPEC=$(addprefix -I,$(SYSINCLUDE_DIRS))

INCLUDE_DIRS=../common
INCLUDE_DIRS_SPEC=$(addprefix -I,$(INCLUDE_DIRS))

SYSLIBS=pthread
SYSLIBS_SPEC=$(addprefix -l,$(SYSLIBS))

SYSLIB_DIRS=
SYSLIB_DIRS_SPEC=$(addprefix -L,$(SYSLIB_DIRS))

LIBS=seq
LIBS_SPEC=$(addprefix -l,$(LIBS))

LIB_DIRS=../lib
LIB_DIRS_SPEC=$(addprefix -L,$(LIB_DIRS))

OBJS=nsimscan_params.o\
    nsimscan.o

DEPENDENCIES=$(patsubst %.o,%.dep,$(OBJS))

vpath %.h $(INCLUDE_DIRS)
vpath lib% $(LIB_DIRS)

%.o : %.cpp
	$(CC) $(CPPFLAGS) $(INCLUDE_DIRS_SPEC) $(SYSINCLUDE_DIRS_SPEC) -c -o $@ $<
	$(CC) $(CPPFLAGS) $(INCLUDE_DIRS_SPEC) $(SYSINCLUDE_DIRS_SPEC) -MM -MG $< >$*.dep

all : $(TARGETDIR) $(TARGET)

$(TARGET) : $(OBJS) $(LIBS_SPEC)
	$(CC) -o $@ $(OBJS) $(LIB_DIRS_SPEC) $(LIBS_SPEC) $(SYSLIB_DIRS_SPEC) $(SYSLIBS_SPEC)

-include $(DEPENDENCIES)

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)
	rm -f $(DEPENDENCIES)
