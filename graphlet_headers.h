/**
 ============================================================================
 Name        : Parallel Parameterized Graphlet Decomposition (PGD) Library

 Author      : Nesreen K. Ahmed, (nesreen.k.ahmed@intel.com),
 		   	   Ryan A. Rossi   (rrossi@parc.com)

 Description : A general high-performance parallel framework for computing
               the graphlet decomposition. The library is designed to be fast
               for both large sparse graphs as well as dense graphs.

 Copyright (C) 2012-2015,
 Nesreen K. Ahmed (http://nesreenahmed.com), All rights reserved.

 Please cite the following paper:
	Nesreen K. Ahmed, Jennifer Neville, Ryan A. Rossi, Nick Duffield,
    Efficient Graphlet Counting for Large Networks, IEEE International
    Conference on Data Mining (ICDM), pages 10, 2015.
	Download PDF: http://www.nesreenahmed.com/publications/ahmed-et-al-icdm2015.pdf

	@inproceedings{ahmed2015icdm,
		title={Efficient Graphlet Counting for Large Networks},
		author={Nesreen K. Ahmed and Jennifer Neville and Ryan A. Rossi and Nick Duffield},
		booktitle={ICDM},
		pages={1--10},
		year={2015}
	}

 See http://nesreenahmed.com/graphlets for more information.
 ============================================================================
 */

#ifndef GRAPHLET_HEADERS_H_
#define GRAPHLET_HEADERS_H_

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <map>
#include <fstream>
#include <sstream>
#include <stdint.h>

#define schedule_type dynamic

#ifdef _OPENMP
#  include <omp.h>
#else
int omp_get_max_threads()       { return 1; }
void omp_set_num_threads(int)   {}
int omp_get_thread_num()        { return 0; }
#endif

using namespace std;

#endif
