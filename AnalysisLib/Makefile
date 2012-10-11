#OS_DEFS SOEXE,SOFLAGS,DEFS
PACKAGE=$(shell basename `pwd`)
UNAME= $(shell uname)
ifeq ($(UNAME), Linux)
SOFLAGS     = -shared
DEFS        = -DOS_LINUX
LIB_EXT    := .so
endif
ifeq ($(UNAME), Darwin)
SOFLAGS     = -dynamiclib
DEFS        = -DOS_DARWIN
LIB_EXT    := .dylib
endif
#AR           = ar
#LD_SHARED    = g++

CXX	     = g++
CC           = gcc
ROOTCINT     = $(ROOTSYS)/bin/rootcint
DEFS		+= -DCERNLIB_TYPE -DCERNLIB_LINUX -DCERNLIB_UNIX -DHEP_SHORT_NAMES

CXX_SOURCES  = $(wildcard src/*.cc)
CXX_DICT     = $(wildcard include/*LinkDef.h)
CXX_DICTSRC  = $(CXX_DICT:include/%LinkDef.h=src/%Dict.cxx)
CXX_OBJECTS  = $(filter-out %Dict.o,$(CXX_SOURCES:src/%.cc=tmp/%.o))
CXX_DICTOBJ  = $(CXX_DICTSRC:src/%.cxx=tmp/*.o)

CXXFLAGS     = -fPIC -Wall
CXXFLAGS    += $(shell $(ROOTSYS)/bin/root-config --cflags)
CXXFLAGS    += -O2 -g -fmessage-length=0


LIBS        += -L$(BOOSTDIR)/lib
LIBS        += -L$(CLHEP_BASE_DIR)/lib -lCLHEP
LIBS        += -L$(G4INSTALL)/lib 
LIBS        += -L$(E14_TOP_DIR)/lib/so
LIBS        += $(shell $(ROOTSYS)/bin/root-config --glibs) -lEG
LIBS        += -lGsimData -lklong -lpi0 -lgamma -lrec2g -lgnana -lE14Fsim -lcsimap -lcluster -lGsimTools -lGsimPersistency
#LIBS        += -L$(ANALYSISLIB)/lib -lAnalysisLib

INCLUDEDIR  += .  ./include 
INCLUDEDIR  += $(ROOTSYS)/include 
INCLUDEDIR  += $(BOOSTDIR)/include 
INCLUDEDIR  += $(CLHEP_BASE_DIR)/include
INCLUDEDIR  += $(G4INSTALL)/include
INCLUDEDIR  += $(E14_TOP_DIR)/include 
#INCLUDEDIR  += $(ANALYSISLIB)/include

INCLUDES     = $(INCLUDEDIR:%=-I%)

COMPILE_CC  := $(CC)  -c $(DEFS) $(INCLUDES) $(CFLAGS) $(ARG)
COMPILE_CXX := $(CXX) -c $(DEFS) $(INCLUDES) $(CXXFLAGS) $(ARG)
DEPEND_CC   := $(CC) -MM $(DEFS) $(INCLUDES) $(filter-out -fPIC, $(CFLAGS))
DEPEND_CXX  := $(CC) -MM $(DEFS) $(INCLUDES) $(filter-out -fPIC, $(CXXFLAGS))
LINK_CC     := $(CC)
LINK_CXX    := $(CXX)

TARGET_SRC   = $(wildcard *.cc)
TARGET_BIN   = $(TARGET_SRC:%.cc=bin/%)
#TARGET_OBJS  = $(TARGET_SRC:%.cc=tmp/%.o)


all: rootDict makelib target

target: $(TARGET_BIN)
	@ basename `pwd`
makelib: $(CXX_OBJECTS)
	@ echo $(ROOTSYS)
	@ echo $(CXX_OBJECTS)	
	@[ -d ./lib ] || mkdir -p ./lib
	$(CXX) $(SOFLAGS) -o lib/lib$(PACKAGE)$(LIB_EXT) $(CXX_OBJECTS) $(CXX_DICTSRC) $(CXXFLAGS) $(INCLUDES) $(LIBS) $(ARG)
ifeq "$(LIB_EXT)" ".dylib"
	rm -f $(PWD)/lib/*.so	
	ln -s $(PWD)/lib/libConvLIB.dylib $(PWD)/lib/libConvLIB.so
endif

rootDict:
ifneq "$(shell ls ./include/*LinkDef.h 2>/dev/null)" ""
	@ echo making rootdicts
	@ for i in `ls ./include/*LinkDef.h | sed -n "s/\(.*\)LinkDef.h/\1.h/p"`;do \
	j=`basename $$i .h`; \
	$(ROOTCINT) -f src/$${j}Dict.cxx -c $(INCLUDES) ./include/$${j}.h ./include/$${j}LinkDef.h; \
	done;
	@ echo End making rootdicts
endif

#bin/%: tmp/%.o $(CXX_OBJECTS)
bin/%: tmp/%.o 
	@ echo making $(@F)
	@ [ -d ./bin ] || mkdir -p ./bin
	@ echo "##########################################################"	
	@ echo COMPILE $(@F)		
	@ echo "##########################################################"
	$(CXX) -o $@ $^ $(filter-out -fPIC,$(CXXFLAGS)) $(CXX_DICTSRC) $(INCLUDES) $(LIBS) -lAnalysisLib -L./lib
	@ echo "##########################################################"

tmp/%.o: src/%.c
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making $(@F) with $(<F) and $(^F)
	@ $(COMPILE_CC) -c $< $(CXXFLAGS) -o $@

tmp/%.d: src/%.c
	@ [-d ./tmp ] || mkdir -p ./tmp
	@ echo making dependencies of $(<F)
	@ $(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*\.o[ :]*/tmp\/$(@F) tmp\/&/g" > $@'

tmp/%.o: src/%.cc
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making $(@F) with $(<F) and $(^F)
	@ $(COMPILE_CXX) $< -o $@

tmp/%.d: src/%.cc
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making dependencies of $(<F)
	@ $(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*\.o[ :]*/tmp\/$(@F) tmp\/&/g" > $@'


tmp/%.o: %.cc
	@ echo $(INCLUDEDIR)
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making $(@F)
	$(COMPILE_CXX) $< -o $@

tmp/%.d: %.cc
	@ [ -d ./tmp ] || mkdir -p ./tmp
	@ echo making dependencies of $(@F)
	@ $(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*\.o[ :]*/tmp\/$(@F) tmp\/&/g" > $@'




.PHONY: clean
clean:
	rm -f $(TARGET_BIN)
	rm -f $(TARGET_OBJS)
	rm -f ./lib/*
	rm -f ./tmp/*.d
	rm -f ./tmp/*.o
	rm -f ./src/*Dict.*
