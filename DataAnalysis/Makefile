# -*- Makefile -*-

MAINS =DataAnalysis
PACKAGE = $(shell basename `pwd`)
CXXFLAGS  += -I. -I$(ANALYSISLIB)/include -g -Wall
EXTRALIBS += -L$(ANALYSISLIB)/lib -lAnalysisLib

INSTALL_LIBS = no	
INSTALL_BINS = no

include $(E14_CONFIG_DIR)/Makefile.bin
