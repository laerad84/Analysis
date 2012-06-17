# -*- Makefile -*-

MAINS = E14_CALIBRATION_CLUSTER_BUILDER_V1
PACKAGE = $(shell basename `pwd`)
CXXFLAGS  += -I. -I$(ANALYSISLIB)/include 
EXTRALIBS += -L$(ANALYSISLIB)/lib -lAnalysisLib

INSTALL_LIBS = no	
INSTALL_BINS = no

include $(E14_CONFIG_DIR)/Makefile.bin
