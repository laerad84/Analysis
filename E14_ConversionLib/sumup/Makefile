#
# Makefile for Prodction
#

CC = g++

CFLAGS = -O2 -Wall -fPIC

PATHLIB = -L$(E14_TOP_DIR)/lib/so
PATHINC = -I$(E14_TOP_DIR)/include

ROOTLIBS =$(shell root-config --glibs)
ROOTFLAGS=$(shell root-config --cflags)
#ROOTINC=-I$(ROOTSYS)/include

TARGET = libE14VMESum.so
SRCS = E14VMESum.C E14Mapper.C E14Calibrator.C EventTree.C SemiOnlinePlot.C CsIImage.C IDHandler.C E14DataDict.C
OBJS = E14VMESum.o E14Mapper.o E14Calibrator.o EventTree.o SemiOnlinePlot.o CsIImage.o IDHandler.o E14DataDict.o 

all: $(TARGET)

.C.o:
	$(CC) $(CFLAGS) $(ROOTFLAGS) $(PATHINC) -c  $<

$(TARGET): $(OBJS) $(SRCS)
	$(CC) -shared -o $(TARGET) $(OBJS) $(PATHLIB) $(ROOTLIBS) 
#	$(CC) -shared -o libE14Mapper.so E14Mapper.o $(PATHLIB) $(ROOTLIBS)
#	@rm -f $(OBJS)

E14DataDict.C: E14VMESum.h E14Mapper.h E14Calibrator.h EventTree.h SemiOnlinePlot.h CsIImage.h IDHandler.h E14DataLinkdef.h
	rootcint -f $@ -c $(PATHINC) $(ROOTFLAGS) $^

clean: 
	rm -f *.so *.o *Dict* *~
