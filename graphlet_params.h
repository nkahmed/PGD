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

#ifndef GRAPHLET_PARAMS_H_
#define GRAPHLET_PARAMS_H_

#include "graphlet_headers.h"
#include "graphlet_utils.h"

using namespace std;

class params {
public:
	string algorithm;
	int adj_limit;
	double density_cutoff;
	bool graph_stats;
	bool verbose;
	bool help;
	string format;

	/**
	 * @brief Graph file given as input
	 */
	string graph;

	/**
	 * @brief Filename for storing the global MACRO graphlet counts,
	 * that is, the count for each graphlet motif.
	 */
	string macro_stats_filename;

	/**
	 * @brief Filename for storing the micro-level graphlet counts,
	 * that is, the motif counts of the various graphlets for each edge
	 */
	string micro_stats_filename; /// MICRO GRAPHLET FEATURES (Motif count for each edge)

	/** @brief Strategy that determines the order in which graphlet counts are computed for each edge in an edge-centric graphlet counting algorithm (node). */
	string ordering;
	bool is_small_to_large; /// direction of ordering (for edges/nodes): default is false
	/** @brief Strategy that determines how the neighbors of each vertex in the csc/csr sparse graph representation are ordered. */
	string ordering_csc_neighbor;
	bool is_small_to_large_csc_neighbor;
	/**
	 * @brief Graph data representation to use: adj = adjacency matrix structure, csc = compressed sparse column, hybrid = both
	 */
	string graph_representation;

	/**
	 * @brief Number of workers (cores, ...)
	 */
	int workers;

	/**
	 * @brief Batch of jobs to assign each worker
	 */
	int block_size;

	/**
	 * @brief Check if micro statistics should be computed
	 *
	 * @return Returns true if micro statistics should be computed otherwise returns false
	 */
	bool is_micro_stats() {
		if (micro_stats_filename == "") return false;
		return true;
	}

	bool is_macro_stats() {
		if (macro_stats_filename == "") return false;
		return true;
	}

	/**
	 * @brief Check if representation is manually set or if it should be determined automatically.
	 * Optimal representation should be used or if  is set manually or
	 *
	 * @return Returns true if set manually, otherwise returns false (determined automatically)
	 */
	bool is_representation_set_manually() {
		if (graph_representation == "") return false;
		return true;
	}

	/**
	 * @brief Initialize default input parameters
	 */
	void initialize() {
		algorithm = "exact";
		verbose = false;
		graph_stats = false;
		help = false;
		graph = "sample_graph.csv";
		macro_stats_filename = "";
		micro_stats_filename = "";
		workers = omp_get_max_threads();
		block_size = 64;
		ordering = "deg";
		is_small_to_large = false;
		ordering_csc_neighbor = "deg";
		is_small_to_large_csc_neighbor = false;
		graph_representation = "";
		adj_limit = 10000;
		density_cutoff = 0.80;
	}

	/**
	 * @brief Display usage information to help the user
	 *
	 * @param argv0
	 */
	void print_help(char *argv0) {
		const char *params =
				"Usage: %s -f path\n"
				"\t \n"
				"\t pgd options:\n"
				"\t=================================================================================\n"
				"\tParallel Parameterized Graphlet Decomposition (PGD) Library\n"
				"\t=================================================================================\n"
				"\t-f, --file,--graph              : Input GRAPH file for computing the graphlets (e.g., matrix market format, simple edge list). \n"
				"\t-a, --algorithm                 : Algorithm for the GRAPHLET DECOMPOSITION. Default: exact, approximate, etc.  \n"
				"\t---------------------------------------------------------------------------------\n"
				"\t-w, --workers                   : Number of PROCESSING UNITS (workers) for the algorithm to use (default = max). \n"
				"\t-b, --block_size                : Size of batch (number of jobs) dynamically assigned to the processing unit, that is, 1, 64, 512, etc.  Default: -b 64\n"
				"\t---------------------------------------------------------------------------------\n"
				"\t-o, --ordering                  : Strategy used to determine the order in which the edge/node graphlets are computed.\n"
				"\t                                  Default: -o degree (kcore, rand, natural/off, etc.).\n"
				"\t    --s2l                       : Direction of the ordering (default: largest to smallest).\n"
				"\t                                  Note to order from smallest to largest, set '--s2l'  \n"
				"\t-n, --neigh_ordering            : Strategy used to order the neighbors of each node. Default: degree (kcore, rand, natural, etc.)\n"
				"\t                                  Note only applicable if CSC/CSR is used (-r csc).\n"
				"\t    --s2l_neigh                 : Order direction for neighbor/csc ordering strategy\n"
				"\t                                  (e.g., --neigh_ordering degree --s2l_neigh, default: largest to smallest)\n"
				"\t---------------------------------------------------------------------------------\n"
				"\t-c, --counts,--macro            : Compute MACRO (GLOBAL) GRAPHLET FEATURES and save counts to file (e.g., --counts name.graphlets)\n"
				"\t-m, --micro                     : Compute MICRO (LOCAL) GRAPHLET FEATURES and save EDGE-by-MOTIF feature matrix (-m name.micro_graphlets)\n"
				"\t                                  Default: OFF. NOTE: Turn ON edge graphlet counting by specifying an output file, e.g., '-m name.micro_graphlets' \n"
				"\t---------------------------------------------------------------------------------\n"
				"\t-v, --verbose                   : Output additional details to the screen. \n"
				"\t-?, -h, --help                  : Print out this help menu. \n"
				"\n"
				"\n"
				"\tREPRESENTATION: Example: ./pgd -r csc (adj, etc.)\n"
				"\t=================================================================================\n"
				"\t-r,   --rep                     : Graph representation [adj, csc, hybrid, auto, etc].\n"
				"\t                                  Default: Auto select optimal. \n"
				"\t\t'adj'    - adjacency matrix   : dense n by n matrix, O(|V|^2) storage cost\n"
				"\t\t'csc'    - comp. sparse col   : large sparse graphs, O(|V|+|E|) storage cost\n"
				"\t\t'hybrid' -  csc + adj         : small sparse and dense graphs, O(|V|^2 + |V| + |E|) storage cost\n"
				"\t-l, --adj_limit                 : Threshold/limit for creating adj representation. Default: '-l 10000' (that is <10000 nodes).\n"
				"\n"
				"\n"
				"\tORDERING TECHNIQUES: Example: ./pgd -o degree (kcore, rand, etc.)\n"
				"\t=================================================================================\n"
				"\t'degree',   'deg'                    : O(|V|)\n"
				"\t'kcore',                             : O(|E|)\n"
				"\t'rand', 'random'                     : O(|V|)\n"
				"\t'off',  'natural'                    \n"
				"\n"
				"\t Other methods for ordering include: \n"
				"\t'kcore_degree', 'kcore_deg'          : O(|V|)\n"
				"\t'degree_vol',   'deg_vol'            : O(|E|)\n"
				"\t'kcore_vol',                         : O(|E|)\n"
				"\t'deg_kcore_vol'                      : O(|E|)\n"
				"\t------------------------------------------------------------------\n"
				"\tNOTE: Default ordering is kcore (degeneracy order). For natural order, use '-o off' or '-o natural'\n"
				"\n"
				"\tCopyright Nesreen K. Ahmed (http://nesreenahmed.com) and Ryan A. Rossi (http://ryanrossi.com).\n"
				"\tWebsite: http://nesreenahmed.com/graphlets for news and updates.\n"
				"\t";
		fprintf(stderr, params, argv0);
		exit(-1);
	}

	/**
	 * @brief Set default parameter values
	 */
	params() { initialize(); }

	/**
	 * @brief Handles the command line interface (CLI) and sets the appropriate input params
	 *
	 * @param argc
	 * @param argv
	 */
	params(int argc, char **argv) {
		if (argc <= 1) print_help(argv[0]);
		initialize();
		map<string, int*> int_flags;
		int_flags["-w"] = int_flags["--workers"] = int_flags["--threads"] = int_flags["--cores"] = &workers;
		int_flags["-b"] = int_flags["--block_size"] = &block_size;
		int_flags["-l"] = int_flags["--adj_limit"] = int_flags["--adj_size"] = &adj_limit;
		map<string, bool*> bool_flags;
		bool_flags["-?"] = bool_flags["-h"] = bool_flags["--help"] = &help;
		bool_flags["-v"] = bool_flags["--verbose"] = &verbose;

		/// direction of ordering (small to largest, or largest to smallest)
		bool_flags["--s2l"] = bool_flags["--small_to_large"] = &is_small_to_large;
		bool_flags["--s2l_neigh"] = bool_flags["--small_to_large_neigh"] = &is_small_to_large_csc_neighbor;
		map<string, double*> double_flags;
		map<string, string*> string_flags;
		string_flags["-f"] = string_flags["--graph"] = string_flags["--file"] = string_flags["-g"] = &graph;
		string_flags["-a"] = string_flags["--algorithm"] = string_flags["--alg"] = &algorithm;

		/** @brief MACRO/MICRO GRAPHLET COUNTS */
		string_flags["-c"] = string_flags["--counts"] = string_flags["--macro"] = &macro_stats_filename; // file to output local counts
		string_flags["-m"] = string_flags["--micro"] = &micro_stats_filename; // file to output local counts

		/**
		 * @brief Representation to use for storing the graph and performing computations
		 */
		string_flags["-r"] = string_flags["--rep"] = &graph_representation;

		/**
		 * @brief Edge order for graphlet algorithms, each edge corresponds to a parallel job, ... (Search/computation order)
		 * See is_small_to_large (for determining the direction of the ordering)
		 */
		string_flags["-o"] = string_flags["--ordering"] = string_flags["--order"] = &ordering;

		/**
		 * @brief Order technique used for ordering the neighbors of each vertex (in CSC/CSR)
		 *
		 */
		string_flags["-n"] = string_flags["--neigh_ordering"] = string_flags["--neigh_order"] = string_flags["--csc_ordering"] = string_flags["--neighbor_ordering"] = &ordering_csc_neighbor;

		validate(argc >= 2, "Error: graph must be supplied.");
		for (int i = 1; i < argc; ++i) {
			string flag(argv[i]);
			if (bool_flags.find(flag) != bool_flags.end()) { *bool_flags[flag] = true; continue; }
			if(i + 1 < argc) {
				if (int_flags.find(flag) != int_flags.end()) *int_flags[flag] = atoi(argv[++i]);
				else if (double_flags.find(flag) != double_flags.end())
					*double_flags[flag] = atof(argv[++i]);
				else if (string_flags.find(flag) != string_flags.end())
					*string_flags[flag] = string(argv[++i]);
				else { printf("Invalid argument %s, ignoring.\n",flag.c_str()); }
			}
		}

		if (*bool_flags["-h"] == true) {
			print_help(argv[0]);
			exit(-1);
		}
		if (workers <= 0) workers = 1;
		if (!fexists(graph.c_str())) {
			print_help(argv[0]);
			exit(-1);
		}
		FILE* fin = fopen(graph.c_str(), "r+t");
		if (fin == NULL) {
			print_help(argv[0]);
			exit(-1);
		}
		fclose(fin);
		omp_set_num_threads(workers);
		if (verbose) { print_params(); }
	}

	/**
	 * @brief Print the list of parameters.
	 * In short, displays the 'param' object including each parameter and its value
	 */
	void print_params() {
		print_line(80,"=");
		cout << "Parallel PARAMETERIZED GRAPHLET DECOMPOSITION (PGD) Library" <<endl;
		print_line(80,"=");
		cout << "file: " << graph.c_str() << endl;
		if (!fexists(graph.c_str()) ) { cout << "File not found!" << endl; return; }
		cout << "algorithm: " << algorithm <<endl;
		if (graph_representation=="" || graph_representation=="auto") cout << "graph representation: automatically select optimal." <<endl;
		else cout << "graph representation: " << graph_representation <<endl;
		cout << "adj limit: " << adj_limit <<endl;
		cout << "-----------ORDERING STRATEGIES------------------\n";
		if (is_small_to_large) cout << "ordering = " << ordering << ", smallest to largest  " <<endl;
		else cout << "ordering: " << ordering << ", largest to smallest" <<endl;
		if (is_small_to_large_csc_neighbor) cout << "csc/neighbor ordering: " << ordering_csc_neighbor << ", smallest to largest" <<endl;
		else cout << "csc/neighbor ordering: " << ordering_csc_neighbor << ", largest to smallest " <<endl;
		cout << "-----------PARALLEL SETTINGS--------------------\n";
		cout << "workers: " << workers <<endl;
		cout << "block size: " << block_size <<endl;
		cout << "-----------MICRO/MACRO SETTINGS-----------------\n";
		cout << "macro (per graph counts) output: " << macro_stats_filename <<endl;
		cout << "micro (per edge counts)  output: " << micro_stats_filename <<endl;
		print_line(80);
		cout <<endl;
	}
};
#endif
