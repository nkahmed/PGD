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

#include "graphlet_utils.h"

using namespace std;

bool fexists(const char *filename) {
    ifstream ifile(filename);
#ifdef WIN32
    return ifile!=0;
#else
    return (bool)ifile;
#endif
}

#if defined(_WIN32) || defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
double get_time() {
    LARGE_INTEGER t, freq;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&freq);
    return double(t.QuadPart)/double(freq.QuadPart);
}
#else
#include <sys/types.h>
#include <sys/timeb.h>
#include <sys/time.h>
double get_time() {
    struct timeval t;
    gettimeofday(&t, 0);
    return (t.tv_sec*1.0 + t.tv_usec/1000000.0);
}
#endif

double tic() { return get_time(); }
void toc(double & start) {
    start = get_time() - start;
}

string memory_usage() {
    ostringstream mem;
    ifstream proc("/proc/self/status");
    string s;
    while(getline(proc, s), !proc.fail()) {
        if(s.substr(0, 6) == "VmSize") {
            mem << s;
            return mem.str();
        }
    }
    return mem.str();
}
void indent(int level = 0, string str = "") {
    for (int i = 0; i < level; i++) cout << "   ";
    cout << "(" << level << ") ";
}
void validate(bool condition, const string& msg) {
    if (!condition) {
        cerr << msg << endl;
        assert(condition);
    }
}
int getdir (string dir, vector<string> &files) {
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }
    while ((dirp = readdir(dp)) != NULL) {
        if (dirp->d_name != ".")
            files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}
void print_line(int n, string sym) {
    for (int i = 0; i < n; ++i) cout << sym;
    cout <<endl;
}

template<typename T>
void write_vector(std::vector<T> &data, string suffix, string fn) {
    string::size_type result;
    string::size_type beg;
    string name = "stats/";

    string rawname, fn_str;
    int last_index = fn.find_last_of(".");
    if (last_index != string::npos) {
        rawname = fn.substr(0, last_index);
        int beg_index = rawname.find_last_of("/");
        string fn_new = rawname.substr(beg_index,rawname.size()-1);
        fn_str = string("stats") + fn_new + suffix;
    }

    ofstream myfile;
    char *fileName = (char*)fn_str.c_str();
    myfile.open (fileName);

    int sum = 0;
    for (long long e = 0; e < data.size(); e++) {
        myfile << data[e] << "\n";
        sum += data[e];
    }
    myfile.close();
}
template void write_vector(std::vector<long long> &data, string suffix, string fn);
template void write_vector(std::vector<int> &data, string suffix, string fn);
template void write_vector(std::vector<double> &data, string suffix, string fn);
string extract_filename(string fn, bool remove_ext) {
    if (remove_ext) { fn = remove_file_extension(fn); }
    int idx = 0;
    idx = fn.find_last_of("/")+1;
    return fn.substr(idx, idx - fn.size());
}
string remove_file_extension(string fn) {
    return fn.substr(0,fn.size() - (fn.size()-fn.find_last_of(".")));
}

int sample_rand_dist(int num_vals, vector<double> & dist_tmp) {
	double rand_num = get_cust_rand();
	int i = 0;
	for (i = 0; i < num_vals; i++) {
		if (dist_tmp[i] > rand_num) { break; }
	}
	return i;
};
