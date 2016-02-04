/*!\mainpage Parallel Parameterized Graphlet Decomposition (PGD) Library
 *
 * \section INTRODUCTION
 *
 * A general high-performance parameterized parallel framework for computing the graphlet decomposition.
 * The library is designed to be fast for both large sparse graphs as well as dense graphs.
 * See the readme for a complete list of features.
 *
 * \section AUTHORS
 *
 * Nesreen K. Ahmed, (nesreen.k.ahmed@intel.com),<BR>
 * Ryan A. Rossi   (rrossi@parc.com)<BR>
 *
 *
 * \section DOWNLOAD
 *
 * <a href="http://nesreenahmed.com/graphlets/pgd.zip">PGD Library</a><BR>
 * <a href="http://networkrepository.com">Network Data Repository</a><BR>
 * See <a href="http://nesreenahmed.com/graphlets">http://nesreenahmed.com/graphlets</a> for more information.<BR>
 *
 *
 * \section PUBLICATIONS
 *
 * If used, cite the following paper:
 *
 * 		Nesreen K. Ahmed, Jennifer Neville, Ryan A. Rossi, Nick Duffield,
 *		Efficient Graphlet Counting for Large Networks, IEEE International
 * 		Conference on Data Mining (ICDM), pages 10, 2015.
 *
 * 		Download PDF: <a href="http://www.nesreenahmed.com/publications/ahmed-et-al-icdm2015.pdf">http://www.nesreenahmed.com/publications/ahmed-et-al-icdm2015.pdf</a><BR>
 *
 * 		@inproceedings{ahmed2015icdm,
 * 			title={Efficient Graphlet Counting for Large Networks},
 * 			author={Nesreen K. Ahmed and Jennifer Neville and Ryan A. Rossi and Nick Duffield},
 * 			booktitle={ICDM},
 * 			pages={1--10},
 * 			year={2015}
 * 		}
 *
 * \section COPYRIGHT & CONTACT
 *
 * Copyright (C) 2012-2015,<BR>
 * Nesreen K. Ahmed (http://nesreenahmed.com), All rights reserved.<BR><BR>
 *
 * Email Nesreen Ahmed at nesreen.k.ahmed [at] intel [dot] com.<BR>
 *
 */
#include "graphlet.h"

using namespace std;
using namespace graphlet;

int main(int argc, char *argv[]) {
	/** @brief parse command args */
	params p(argc, argv);
	if (p.verbose) cout << "graph file: " << p.graph <<endl;

	/** @brief read graph, optimize alg/data structs, etc. */
	graphlet_core G(p);

	G.compute_assortativity();
	cout << "r = " << G.get_assortativity() << endl; // Assortativity

	/** @brief Compute k-core decomposition (k-core numbers, and degen. ordering) */
	G.compute_cores();
	cout << "K = " << G.get_max_core() <<endl; // Max K-core

	/// Creates E_ordered where edges are ordered via a strategy (kcore, degree, etc.)
	double s = tic();
	G.sort_edges(p.ordering, p.is_small_to_large);

	if (p.verbose) cout << "edge/job ordering: " << p.ordering << ", time: " << get_time()-s <<endl;
	if (p.is_micro_stats()) { G.graphlet_decomposition_micro(p.workers); }
	else { G.graphlet_decomposition(p.workers); }
	toc(s);
	G.print_graphlet_counts();
	if (p.is_micro_stats()) G.print_micro_stats();
	cout << "graphlet decomposition time: " << s << " sec" <<endl;

	/// Save total counts of each motif (global macro stats) to the file specified by the user.
	if (p.is_macro_stats()) G.write_macro_stats(p.macro_stats_filename);
	if (p.is_micro_stats()) G.write_micro_stats(p.micro_stats_filename);
	G.print_GFD();
	G.print_connected_GFD();
	G.print_disconnected_GFD();
	if (p.is_micro_stats()) {
		univar_stats s;
		s.compute_univariate_stats(G.local_4_clique);
		cout << s.tostring() <<endl;
	}
	return 0;
}
