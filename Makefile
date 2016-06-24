# Makefile for Jpsi2invi project
# Author: SHI Xin <shixin@ihep.ac.cn>
# Created: [2016-06-08 Wed 08:24]

BIN=./bin
SRC=./src

CC = g++
GCCFLAGS  = -Wall -g

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))

FLAGS=$(GCCFLAGS) $(ROOTCFLAGS) $(ROOTLIBS) -lHistPainter

PROG=runPlotter

LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)
	@echo "Build successful."


$(LIST): | $(BIN)


$(BIN): 
	mkdir -p $(BIN)


$(BIN)/runPlotter: $(SRC)/runPlotter.cc
	$(CC) $< $(FLAGS) -o $@

clean:
	rm -f $(BIN)/runPlotter


