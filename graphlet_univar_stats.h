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

#ifndef UNIVAR_STATS_H_
#define UNIVAR_STATS_H_

#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;
namespace graphlet {

    class univar_stats {
        private:
			double sum_sq_diff;
			unsigned long long sum;
			double variance;
			double std;
			double mean;
			double median;
			unsigned long long min;
			unsigned long long max;
			unsigned long long range;
			unsigned long long num_unique_values;
			double q1;
			double q3;
			double iqr; /* The interquartile range is often used to find outliers in data. */
			double ub;
			double lb;

		    bool static ascending_func(unsigned long long v,  unsigned long long u) { return (v < u); }; // smallest to largest (incr_bound)
		    bool static is_even(unsigned long long val) { return val % 2 == 0 ? true : false; }

        public:
			univar_stats(vector<unsigned long long> & data) { initialize(); compute_univariate_stats(data); };

			univar_stats() { initialize(); };
			void reset() { initialize(); }
			void initialize() {
				sum_sq_diff = sum = variance = std = mean = median = range = 0;
				num_unique_values = min = max = 0;
				lb = ub = q1 = q3 = iqr = 0;
			}

			void compute_univariate_stats(vector<unsigned long long> &data) {
				initialize();
				for (unsigned long long i = 0; i < data.size(); i++) {
					data[i] = +data[i];
					sum_sq_diff += data[i] * data[i];
					sum += data[i];
					if (data[i] > max) max = data[i];
					if (data[i] < min) min = data[i];
					num_unique_values++;
				}
				mean = sum / num_unique_values;
				variance = (sum_sq_diff / num_unique_values) - (mean * mean);
				std = sqrt(variance);
				range = max - min;

				vector<unsigned long long> arr = data; // deep copy
				sort(arr.begin(), arr.end(), ascending_func); // sort arr from smallest to largest

				unsigned long long med_idx = floor(num_unique_values / 2);
				unsigned long long quartiles_idx = floor(num_unique_values / 4); // plus/minus quartiles_idx to get Q1 and Q3
				unsigned long long Q1_idx = med_idx - quartiles_idx;
				unsigned long long Q3_idx = med_idx + quartiles_idx;

				median = arr[med_idx];
				if (is_even(arr.size())) { median = (arr[med_idx - 1] + arr[med_idx]) / 2; }

				q1 = arr[Q1_idx];
				q3 = arr[Q3_idx];
				iqr = q3 - q1;

				/// The interquartile range is often used to find outliers, that is, observations that fall below Q1 - 1.5(IQR) or above Q3 + 1.5(IQR)
				lb = q1 - (1.5 * iqr);
				ub = q3 + (1.5 * iqr);
			}

			double get_mean() { return mean; }
			unsigned long long get_min() { return min; }
			unsigned long long get_max() { return max; }
			double get_variance() { return variance; }
			double get_std() { return std; }

			string tostring(string delim="\n") {
				ostringstream os;
				os << "mean\t= "<< mean <<delim;
				os << "median\t= "<< median <<delim;
				os << "max\t= "<< max <<delim;
				os << "min\t= "<< min <<delim;
				os << "range\t= "<< range <<delim;
				os << "std\t= "<< std <<delim;
				os << "var\t= "<< variance <<delim;
				os << "iqr\t= "<< iqr <<delim;
				os << "q1\t= "<< q1 <<delim;
				os << "q3\t= "<< q3 <<delim;
				return os.str();
			}
    };
};
#endif /* UNIVAR_STATS_H_ */
