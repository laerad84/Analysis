#
# Makefile for Prodction
#

CC = g++

CFLAGS = -O2 -Wall -fPIC

#PATHLIB = -L$(E14_TOP_DIR)/lib/so
#PATHINC = -I$(E14_TOP_DIR)/include

ROOTLIBS =$(shell root-config --glibs)
ROOTFLAGS=$(shell root-config --cflags)
#ROOTINC=-I$(ROOTSYS)/include

TARGET = libE14Data.so
SRCS = E14RawData.C E14DataReadVME.C E14Fill.C E14Mapper.C E14Calibrator.C SemiOfflineHist.C E14DataDict.C
OBJS = E14RawData.o E14DataReadVME.o E14Fill.o E14Mapper.o E14Calibrator.o SemiOfflineHist.o E14DataDict.o

all: $(TARGET)

.C.o:
	$(CC) $(CFLAGS) $(ROOTFLAGS) $(PATHINC) -c  $<

$(TARGET): $(OBJS) $(SRCS)
	$(CC) -shared -o $(TARGET) $(OBJS) $(PATHLIB) $(ROOTLIBS)
#	$(CC) -shared -o libE14Mapper.so E14Mapper.o $(PATHLIB) $(ROOTLIBS)
#	@rm -f $(OBJS)

E14DataDict.C: E14RawData.h E14DataReadVME.h E14Fill.h E14Mapper.h E14Calibrator.h SemiOfflineHist.h E14DataLinkdef.h
	rootcint -f $@ -c $(PATHINC) $(ROOTFLAGS) $^

clean: 
	rm -f *.so *.o *Dict* *~
