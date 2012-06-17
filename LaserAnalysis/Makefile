# -*- Makefile -*-

MAINS =LaserAnalysis

PACKAGE = $(shell basename `pwd`)
CXXFLAGS  += -I. -I$(ANALYSISLIB)/include 
EXTRALIBS += -L$(ANALYSISLIB)/lib -lAnalysisLib


INSTALL_LIBS = no	
INSTALL_BINS = no

include $(E14_CONFIG_DIR)/Makefile.bin
