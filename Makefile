#
# Makefile for PARAMETERIZED GRAPHLET DECOMPOSITION (PGD)
# Website: http://nesreenahmed.com/graphlets
#
# AUTHORS:
# Nesreen K. Ahmed (http://nesreenahmed.com)
# Ryan A. Rossi (http://ryanrossi.com)
#
# Copyright, 2012-2015
#

.KEEP_STATE:

ROOTDIR  = $(shell pwd)
OPSYS    = $(shell uname -s )

OPENMP_LIB ?= 1

all: pgd

CFLAGS = -O3 

ifeq ($(OPENMP_LIB),1)
        CFLAGS += -fopenmp
else
        CFLAGS += -Wno-unknown-pragmas
endif


CFLAGS += -ffast-math -funroll-loops -fno-strict-aliasing \
	 -fomit-frame-pointer -fexpensive-optimizations -funroll-loops \
	 -fmove-loop-invariants -fprefetch-loop-arrays -ftree-loop-optimize \
	 -ftree-vect-loop-version -ftree-vectorize
CXX          = g++


H_FILES     = graphlet.h
.cpp.o:
	$(CXX) $(CFLAGS) -c $<

IO_SRC 				= graphlet_utils.cpp graphlet_core.cpp graphlet_rand.cpp
GRAPHLET_MAIN		= graphlet_driver.cpp

OBJ_GRAPHLET	= $(GRAPHLET_MAIN:%.cpp=%.o) $(IO_SRC)
$(OBJ_GRAPHLET): $(H_FILES) Makefile
pgd: $(OBJ_GRAPHLET) $(H_FILES)
	$(CXX) $(CFLAGS) -o pgd $(OBJ_GRAPHLET)	

	
test:
	./pgd -f sample_graph.csv
	
test-4-clique:
	./pgd -f data/4-clique.txt --micro data/4-clique.edges

test-4-chordal-cycle:
	./pgd -f data/4-chordal-cycle.txt --micro data/4-chordal-cycle.txt
	
help:  
	@echo " "
	@echo "Operating System Detected: $(OPSYS) "
	@echo " "
	@echo "USAGE: "
	@echo "make help       To get this listing"
	@echo "make            To compile pgd library"
	@echo "make clean      Remove *.o and executable files"
	@echo "make list       List the compilers in current environment"
	@echo "make tar        Compress the graphlet package to graphlet.tgz"
	@echo "make docs	   Create documentation in 'docs' directory"
	@echo "make cleandocs  Clean 'docs' directory"
	@echo "make test	   Test example './pgd -f sample_graph.csv'"
	@echo " "

list:
	@echo
	@echo "OPSYS:       $(OPSYS)"
	@echo "ROOTDIR:     $(ROOTDIR)"
	@echo " "

docs:
	@echo "Creating documentation in 'docs' directory"
	@echo "NOTE: doxygen must be installed to generate docs"
	@echo "      "
	brew install doxygen # try to install doxygen if not installed (mac osx only)
	brew install graphviz
	doxygen Doxyfile
	
cleandoc:
	@echo "Removing 'docs' directory"
	rm -r docs


DL_DIR = download
DL_SRC_DIR = src
DL_EXEC_DIR = exec
DL_DATA_DIR = data
DL_FULL_DIR = full

TIME = `date +%Y-%m-%d_%H-%M`

tar:
	make clean; cd ../;  tar cvzf $(DL_DIR)/$(DL_SRC_DIR)/pgd_$(TIME).tgz pgd/
	make clean; cd ../;  zip -r $(DL_DIR)/$(DL_SRC_DIR)/pgd_$(TIME).zip pgd/
	make clean; cd ../;  7z a $(DL_DIR)/$(DL_SRC_DIR)/pgd_$(TIME).7z pgd/
	
tar-full: # src code AND data
	make clean; cd ../;  tar cvzf $(DL_DIR)/$(DL_FULL_DIR)/pgd_full_`date +%Y-%m-%d_%H-%M`.tgz pgd/ data/
	make clean; cd ../;  zip -r $(DL_DIR)/$(DL_FULL_DIR)/pgd_full_`date +%Y-%m-%d_%H-%M`.zip pgd/ data/
	make clean; cd ../;  7z a $(DL_DIR)/$(DL_FULL_DIR)/pgd_full_`date +%Y-%m-%d_%H-%M`.7z pgd/ data/
	
tar-exec: 
	make clean; cd ../;  tar cvzf $(DL_DIR)/$(DL_EXEC_DIR)/pgd_exec_`date +%Y-%m-%d_%H-%M`.tgz pgd
	make clean; cd ../;  zip -r $(DL_DIR)/$(DL_EXEC_DIR)/pgd_exec_`date +%Y-%m-%d_%H-%M`.zip pgd
	make clean; cd ../;  7z a $(DL_DIR)/$(DL_EXEC_DIR)/pgd_exec_`date +%Y-%m-%d_%H-%M`.7z pgd
	
clean:
	@echo "Removing *.o and executable files"
	rm -rf *.o pgd