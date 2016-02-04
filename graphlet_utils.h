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

#ifndef GRAPHLET_UTILS_H_
#define GRAPHLET_UTILS_H_

#ifdef WIN32

#else
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#endif

#include <cstddef>
#include <iostream>
#include "assert.h"
#include <errno.h>
#include <string>
#include <set>
#include <vector>
#include "graphlet_headers.h"
#include "graphlet_rand.h"

using namespace std;

bool fexists(const char *filename);
void usage(char *argv0);
double get_time();
double tic();
void toc(double & start);
string memory_usage();
void validate(bool condition, const string& msg);
void indent(int level);
void indent(int level, string str);
void print_line(int n = 80, string sym = "-");
string extract_filename(string fn, bool remove_ext = true);
string remove_file_extension(string fn);
int getdir (string dir, vector<string> &files);

inline
void write_buffer(ofstream & myfile, ostringstream & str_stream) {
	myfile << str_stream.str();
	str_stream.clear();
	str_stream.str("");
}

template<typename T>
void write_results_batch(std::vector<T> &data, string &filename, bool output_id = false,
		unsigned int buffer_size = 128000000) {
    ofstream myfile;
    char *fn = (char*)filename.c_str();
    myfile.open(fn);
    ostringstream str_stream;
    if (output_id) {
        for (long long e = 0; e < data.size(); e++) {
            str_stream << e << "," << data[e] << "\n";
            if (str_stream.str().size() > buffer_size) { write_buffer(myfile, str_stream); };
        }
    }
    else {
        for (long long e = 0; e < data.size(); e++) {
            str_stream << data[e] << "\n";
            if (str_stream.str().size() > buffer_size) { write_buffer(myfile, str_stream); };
        }
    }
    write_buffer(myfile, str_stream);
    myfile.close();
}

template<typename T>
void write_results(std::vector<T> &data, string filename, bool output_id = false) {
    ofstream myfile;
    char *fn = (char*)filename.c_str();
    myfile.open(fn);
    if (output_id) {
        for (long long e = 0; e < data.size(); e++) {
            myfile << e << "," << data[e] << "\n";
        }
    }
    else {
        for (long long e = 0; e < data.size(); e++) {
            myfile << data[e] << "\n";
        }
    }
    myfile.close();
}

template<typename T>
void write_results_line(std::vector<T> &data, string filename, string delim = "\t") {
    ofstream myfile;
    char *fn = (char*)filename.c_str();
    myfile.open(fn);
    for (int e=0; e<data.size(); e++) {
    	myfile << data[e] << "\t";
    }
    myfile << "\n";
    myfile.close();
}

void write_string(string &data, string filename);

template<typename T>
void bin_values(std::vector<T> &data, vector<int> & bin) {
    for (int v = 0; v < data.size(); ++v) {
        int val = data[v];
        bin[val]++;
    }
}
template<typename T>
void write_vector(std::vector<T> &data, string suffix);

int sample_rand_dist(int num_vals, vector<double> & dist_tmp);
#endif
