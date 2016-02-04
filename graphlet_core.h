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

#ifndef GRAPHLET_CORE_H_
#define GRAPHLET_CORE_H_

#ifdef WIN32
#else
#include <sys/time.h>
#include <unistd.h>
#endif
#include <stdlib.h>
#include <float.h>
#include <cstddef>
#include <iostream>
#include <limits>
#include "math.h"
#include "graphlet_rand.h"
#include "graphlet_headers.h"
#include "graphlet_utils.h"
#include "graphlet_vertex.h"
#include "graphlet_params.h"
#include "graphlet_univar_stats.h"

namespace graphlet {
    class graphlet_core {
        private:
	     /**
		* @brief Read mtx file and error-correcting
		*
		* Read MTX graph file and attempt to correct any problems that arise
		*
		* @param filename of the graph to read
		*/
    	int read_mtx(const string& filename);

    	/**
    	 * @brief Read generic edge list
    	 *
    	 * Read generic edge list with arbitrary delimiters (comma, space, tab, etc.)
    	 * File may also have edge weights and attributes, which are ignored if included.
    	 * May remove self-loops and other potential issues if necessary, as well as
    	 * additional preprocessing routines that enable fast and efficient computations
    	 *
    	 * @param filename of the graph to read
    	 */
    	void read_edge_list(const string& filename);

        public:
            /// number of jobs each worker is assigned (at a time). Must be a positive integer.
            int block_size;
            double density_threshold;

            /// stores the (row/col) indices of each nonzero (row/col ind),
            vector<int> edges;
            /// stores the position of the neighbors in edges
            vector<long long> vertices;
            /// input edge weights (if given)
            vector<double> wt;

            /**
             * @brief Map to get original ids
             *
             * NOTE: Vertex IDs are only mapped when needed.
             * For instance, if graph file uses string names for vertices instead of numeric IDS.
             * Also, graph files with certain problems (gaps in IDs, etc.) are remapped to enable more efficient computations
             */
            vector<int> vertex_lookup;

            /// Store vertex degrees
            vector<int> degree;

            /// Vertex degree stats
            int min_degree, max_degree, avg_degree;

            /// Properties of graph
            bool is_adj, is_weighted, is_gstats, verbose;
            string fn;

            /**
             * @brief Reset object variables
             */
            void initialize();

            /**
             * @brief Optimize graph operations (initial)
             *
             * @param adj_limit is the "Adjacency threshold" to use for determining
             * when to perform certain optimizations.
             * This includes optimizations with respect to the
             * 	+ data structure
             * 	+ algorithmic
             *  	+ and many other optimizations
             *
             *  Set rep to "auto" or "" to automatically select/tune appropriately
             */
            void optimize_graph_ops(int adj_limit, string rep = "auto");

            /**
             * @brief Check if adj structure exists.
             * Adj is used for graphs that are small or dense.
             *
             * @return Returns true if adj structure A exists, otherwise returns false
             */
            inline
			bool is_adj_exists() { return is_adj; }

            /** @brief adj. graph rep. using short int to allow for edge weights to be used as indicators. */
            vector< vector<short int> > A;
            void create_adj_mat();

            /**
             * @brief Allows the user to set the representation manually (for specific or special cases).
             * NOTE that if adj structure already exists, then we return immediately to avoid recomputing the data structure.
             *
             * @param rep is a string indicating the type of representation to use
             */
            void set_rep_manually(string &rep);

            ~graphlet_core();
            graphlet_core() {}
            graphlet_core(params & p);
            graphlet_core(const string& filename, string ordering_technique="degree", bool is_small_to_large=false, int adj_limit=10000, string rep="auto", int block_sz=64);
            graphlet_core(const string& filename);
            graphlet_core(bool graph_stats, const string& filename);
            graphlet_core(const string& filename, bool make_adj);
            graphlet_core(vector<long long> vs, vector<int> es) {
                edges = es;
                vertices = vs;
                vertex_degrees();
            }
            graphlet_core(int nverts, int nedges, std::pair<int,int>* edges);
            graphlet_core(int nverts, int* heads, int* tails);

            void preprocess(string ordering_csc_neighbor="degree", bool is_small_to_large="false", int adj_limit=10000, string rep="auto");

            /**
             * @brief   Functions for reading "almost" ANY graph given to us.
             *          Meant to be flexible, and assumes the user inputs a
             *          reasonable graph.
             *
             * @param filename
             */
            void read_graph(const string& filename);

            /**
             * @brief Detect delimiter automatically from the input graph file
             *
             * @param line read by the graph reader and used for determining the delimiter used
             */
            void detect_delim(string & line, string & delim);

            /**
             * @brief Detect weighted graph from input graph file
             *
             * @param line read by the graph reader to parse and ultimately determine if one or more edge weights or attributes exists
             * @param delim to use for parsing line (see @detect_delim)
             * @return true if graph is weighted or attributed otherwise returns false
             */
            bool detect_weighted_graph(string & line, string & delim);

            void get_token(int & v, string & line, string & delim, size_t & pos, size_t & prev);
            void get_token(double & weight, string & line, string & delim, size_t & pos, size_t & prev, bool & is_weighted_graph);

            /**
             * @brief Discard graphs with multiple identical edges.
             *
             * Multiple identical edges are removed which improves efficiency of PGD and
             * reduces the chance of errors arising from these identical edges.
             * This routine takes O(|E|) time.
             * Used by the graph readers that read arbitrary edge lists.
             * Further, it is only used if the graph file contains such identical edges.
		 */
            void remove_multiple_edges();

            /// Array of edges ordered via input strategy
            vector<Vertex> E_ordered;
            void sort_edges(string edge_ordering = "degree", bool is_small_to_large = false);
            /** @brief Create an ordered edge list from a set of sampled edges E_s such that |E_s| < |E| */
            void sort_edges(string ordering_technique, bool is_small_to_large, vector<unsigned long> & E_s);

            void get_ordering_value(int &ordering, long long &v, long long &u, long long &val);

            /**
             * @brief Get the neighbor (vertex id)
             * @param pos is the position of the neighbor in edges for a specific vertex
             */
            inline
			long long get_neighbor(int & pos) { return edges[pos]; }

            int num_vertices() { return vertices.size() - 1; }
            int num_edges() { return edges.size()/2; }
            vector <long long>* get_vertices(){ return &vertices; }
            vector<int>* get_edges(){ return &edges; }
            vector<int>* get_degree(){ return &degree; }
            vector< vector<short int> > get_adj() { return A; }
            vector<int> get_edges_array() { return edges; }
            vector<long long> get_vertices_array() { return vertices; };
            vector<long long> e_v, e_u;

            /**
             * @brief Compute degree of each vertex and store it in an array (degree).
             * Max and average degree of the input graph are also stored,
             * see get_max_degree() and get_avg_degree().
             *
             * Compute vertex degrees and store them in "degree" array.
             * Also, computes maximum degree and average degree.
             * See get_max_degree() and get_avg_degree().
             */
            void vertex_degrees();

            /**
             * @brief Density of the graph G
             *
             * @return Graph density
             */
            double density() { return (double)num_edges() / (num_vertices() * (num_vertices() - 1.0) / 2.0); }

            /**
             * @brief Computes density assuming n vertices and m edges (given as input)
             *
             * @param n is the number of vertices
             * @param m is the number of edges
             * @return Graph density resulting from n vertices and m edges given as input
             */
            double density(int n, long long m) { return (double) m / (n * (n - 1.0) / 2.0); }

            /**
             * @brief Get the maximum vertex degree
             *
             * @return The maximum vertex degree
             */
            int get_max_degree() { return max_degree; }

            /**
             * @brief Get the minimum vertex degree
             *
             * @return The minimum vertex degree
             */
            int get_min_degree() { return min_degree; }

            /**
             * @brief Get the average vertex degree
             *
             * @return The average vertex degree
             */
            double get_avg_degree() { return avg_degree; }

            /**
             * @brief Degree of a vertex v (Number of neighbors adjacency to v).
             *
             * @param v is the id of the vertex
             * @return degree or number of neighbors adjacent to v denoted as $|N(v)|$ or $d_{v}$
             */
            inline
            int get_degree(long long & v) { return vertices[v+1] - vertices[v]; }

            void get_progress(long long &e);
            void load_progress_bar(long long x, long long n, long long w = 50);

            /** @brief Create e_u and e_v arrays. Useful after csc neighbors are ordered, as e_u and e_v are ordered as well. */
            void create_edge_list_arrays();

            /** @brief triangle counts AND wedges */
            unsigned long long total_wedges;
            unsigned long long total_t, max_t_edge;
            unsigned long long num_triangles() { return total_t; }
            unsigned long long get_max_edge_triangles() { return max_t_edge; }
            double get_avg_edge_triangles() { return double(max_t_edge) / total_t; }
            unsigned long long get_triangle_count(long long e) { return tri[e]; }
            double get_clustering_coeff(long long e) {
            	return tri[e] / (get_degree(e_v[e]) + get_degree(e_u[e]));
            }

            /** @brief clustering coefficients */
            double mean_cc, global_cc;
            double get_cc_avg() { return mean_cc; }
            double get_cc_global() { return global_cc; }

            /// graphlet measures
            unsigned long long total_2_1edge;
            unsigned long long total_2_indep;

            unsigned long long total_3_tris;
            unsigned long long total_2_star;
            unsigned long long total_3_1edge;
            unsigned long long total_3_indep;

            unsigned long long total_4_clique;
            unsigned long long total_4_chordcycle;
            unsigned long long total_4_tailed_tris;
            unsigned long long total_4_cycle;
            unsigned long long total_3_star;
            unsigned long long total_4_path;

            unsigned long long total_4_2edge;
            unsigned long long total_4_1edge;
            unsigned long long total_4_2star;
            unsigned long long total_4_tri;
            unsigned long long total_4_indep;

            /** @brief micro-level graphlet counts */

            // connected (k=3) motifs
            /** @brief triangle counts (per edge). */
            vector<unsigned long long> tri;

            // connected (k=4) motifs
            vector<unsigned long long> local_4_clique;
            vector<unsigned long long> local_4_chordal_cycle;
            vector<unsigned long long> local_4_tailed_tris;
            vector<unsigned long long> local_4_cycle;
            vector<unsigned long long> local_3_star;
            vector<unsigned long long> local_4_path;

            /**
             * @brief Get number of 2-star motifs for $e \in E$ in o(1) time.
             * Prereq is the triangle count for that edge
             *
             * @param e is an edge id
             * @return Number of 2-star motifs for edge $e \in E$
             */
            inline
            unsigned long long get_2_star(long long & e, long long & deg_v, long long & deg_u) {
            	return (deg_v - tri[e]) + (deg_u - tri[e]) - 2;
            }

            /**
             * @brief Get number of 2-star motifs for $e \in E$ in o(1) time.
             * Prereq is the triangle count for that edge
             *
             * @param e is an edge id
             * @return Number of 2-star motifs for edge $e \in E$
             */
            inline
            unsigned long long get_2_star(long long & e) {
            	return (get_degree(e_v[e]) - tri[e]) + (get_degree(e_u[e]) - tri[e]) - 2;
            }


            /**
             * @brief This function gets the names of the k=4 graphlets.
             * Note that the graphlet names are in the same order as the
             * vector of graphlet stats provided by 'get_graphlet_size4_stats()'
             *
             * @return vector of the k=4 graphlet names/labels
             */
            vector<string> get_graphlet_size4_names() {
            	vector<string> g4_names;
            	g4_names.reserve(11);
            	g4_names.push_back("4-clique");
            	g4_names.push_back("4-chordal-cycle");
            	g4_names.push_back("4-tailed-tri");
            	g4_names.push_back("4-cycle");
            	g4_names.push_back("3-star");
            	g4_names.push_back("4-path");
            	// DISCONNECTED
            	g4_names.push_back("4-node-1-tri");
            	g4_names.push_back("4-node-2-star");
            	g4_names.push_back("4-node-2-edge");
            	g4_names.push_back("4-node-1-edge");
            	g4_names.push_back("4-node-indep");
            	return g4_names;
            }

            /**
             * @brief This function gets the stats of the k=4 graphlets.
             *
             * @return vector of the k=4 graphlet stats
             */
            vector<unsigned long long> get_graphlet_size4_stats() {
            	vector<unsigned long long> g4_stats;
            	g4_stats.reserve(11);
            	g4_stats.push_back(total_4_clique);
            	g4_stats.push_back(total_4_chordcycle);
            	g4_stats.push_back(total_4_tailed_tris);
            	g4_stats.push_back(total_4_cycle);
            	g4_stats.push_back(total_3_star);
            	g4_stats.push_back(total_4_path);
            	// DISCONNECTED
            	g4_stats.push_back(total_4_tri);
            	g4_stats.push_back(total_4_2star);
            	g4_stats.push_back(total_4_2edge);
            	g4_stats.push_back(total_4_1edge);
            	g4_stats.push_back(total_4_indep);
            	return g4_stats;
            }

            /**
             * @brief Arrays for storing distributions for GFD, Connected GFD, Disconnected GFD.
             */
            vector<double> GFD, connected_GFD, disconnected_GFD;

            /**
             * @brief Get methods for the various graphlet distirbutions.
             * Note if the distributions do not exist or are stale, then they are computed directly and the result is returned.
             *
             * @return Distribution(s) for the various graphlet patterns
             */
            vector<double> get_GFD() 				{ if (GFD.size()==0) compute_GFD(); return GFD; }
            vector<double> get_connected_GFD() 		{ if (connected_GFD.size()==0) compute_connected_GFD(); return connected_GFD; }
            vector<double> get_disconnected_GFD() 	{ if (disconnected_GFD.size()==0) compute_disconnected_GFD(); return disconnected_GFD; }

            /**
             * @brief compute graphlet frequency distribution (GFD) using both connected and disconnected motif counts
             *
             * @return string of the graphlet freq. distribution
             */
            string compute_GFD();

            /**
             * @brief compute GFD using ONLY _connected motifs_
             *
             * @return string of the graphlet freq. distribution
             */
            string compute_connected_GFD();

            /**
             * @brief compute GFD using ONLY _disconnected motifs_
             *
             * @return string of the graphlet freq. distribution
             */
            string compute_disconnected_GFD();

            void print_GFD();
            void print_connected_GFD();
            void print_disconnected_GFD();

            void write_GFD(string fn);
            void write_connected_GFD(string fn);
            void write_disconnected_GFD(string fn);

            /**
             * @brief Get the micro-level graphlet counts for each edge (or vertex) in the graph
             *
             * Example:
             * 	src,dst,triangle,2-star,3-node-1-edge,3-node-independent,4-clique,4-chordal-cycle,4-tailed-triangle, ...
             * 	1, 2, ...
             *
             * @param output_fn is the output filename containing the "micro-level" motif counts (per edge/vertex)
             */
            void write_micro_stats(string &filename, bool output_id = true, unsigned int buffer_size = 0);

             /**
              * @brief Get string consisting of the motif counts for each edge (node) in the graph
              * Example: string micro_stat_str = get_micro_stats();
              *
              * @param filename
              * @param output_id
              * @param buffer_size
              * @return String containing the motif counts for each edge
              */
             string get_micro_stats(bool output_id=true, string delim=",");

             /**
              * @brief Print the micro-level graphlet statistics, that is, motif counts for each edge in the graph
              *
              * @param output_id is a boolean flag indicating whether to include the edge id. Note that if true, then first two columns is the src and dst node ids
              * @param delim is the delimiter to use (default is comma)
              */
             void print_micro_stats(bool output_id=true, string delim=",");

             void write_macro_stats(string filename);

            /**
             * @brief Output counts of the various motif patterns
             */
            void debug_print_graphlet_counts(params & p, string prefix = "");

            string get_graphlet_counts_line(string delim=",");
            string get_graphlet_names_line(string delim=",");

            string basic_names_line(string delim=",", string prefix="", string suffix="");
            string basic_stats_line(string delim=",", string prefix="", string suffix="");

            /**
             * @brief Output counts of the various motif patterns
             */
            void print_graphlet_counts();

            /**
             * @brief Get a string consisting of the name and frequency of each graphlet motif.
             *
             * @return Returns a formatted string containing the name of each graphlet and its frequency.
             */
            string get_graphlet_counts();

            /**
             * @brief Get a string containing the graphlet frequency distribution (GFD)
             *
             * @return Returns a string of the graphlet frequency distribution (GFD)
             */
            string get_graphlet_frequency_distribution();

            /**
             * @brief Sort the neighbors of each vertex (in CSC/CSR)
             *
             * NOTE:  Function may be called by order_vertex_neighbors
             *
             */
            void bucket_sort_neighbors_parallel(vector<int> &bound, bool is_small_to_large = false);

            /**
             * @brief Neighbors of each vertex are ordered (in CSC/CSR)
             *
             * Function calls "bucket_sort_neighbors_parallel"
             *
             * @param p parameters given as input
             */
            void order_vertex_neighbors(string & ordering_technique, bool is_small_to_large = false);

            void compute_graph_measure(string &ordering_technique);

            void solve_graphlet_equations(vector< vector<unsigned long long> > & thread_n_count, vector<unsigned long long> & n_count,
            		vector<unsigned long long> & thread_tmp_triangles, vector<unsigned long long> & thread_tmp_3_star,
            		vector<unsigned long long> & thread_tmp_4_2edge, vector<unsigned long long> & thread_tmp_3_1edge,
            		long long & deg_v, long long & deg_u, unsigned long long & tri_count, unsigned long long & w_local_count, long long & m, long long & n);

            inline unsigned long long compute_4indep(unsigned long long n) {
//            	double double div_val = n/(4.0*3.0*2.0);
//            	total_4_indep = (unsigned long long)(div_val * (n-1) * (n-2) * (n-3));

//            	long double div_val = n/(4.0*3.0*2.0);
            	total_4_indep = (n*(n-1)*(n-2)*(n-3)) / (4*3*2);
            	total_4_indep = total_4_indep - (total_4_clique+total_4_chordcycle+total_4_cycle+total_4_path+total_4_tailed_tris+total_3_star);
            	total_4_indep = total_4_indep - (total_4_2edge+total_4_1edge+total_4_2star+total_4_tri);
            	return total_4_indep;
            }

            void graphlet_decomposition(int max_num_workers = omp_get_max_threads());
            void graphlet_decomposition_micro(int max_num_workers = omp_get_max_threads());

            /**
             * @brief ADJ GRAPH REP. is COMPUTED (Small/Dense graphs): Allows us to check edge existence in O(1) time.
             * Avoids the cost of creating (and resetting) a perfect hash on the neighbors of v (or u), that is, O(2|N(v)|)
             */
            void triangles_and_wedges_adj(long long & v, long long & u, vector<long long> & T_vu, unsigned long long & tri_count,
            		vector<long long> & W_u, unsigned long long & w_local_count, vector< vector<short int> > & adj_mat);
            void cycle_adj(unsigned long long & w_local_count, vector<long long> & W_u, unsigned long long & cycle4_count, long long & v, vector< vector<short int> > & adj_mat);
            void clique_adj(unsigned long long & tri_count, vector<long long> & T_vu, unsigned long long & clique4_count, long long & v, vector< vector<short int> > & adj_mat);

            /** @brief CSC SPARSE GRAPH REP. ONLY */
            void triangles_and_wedges(long long & v, long long & u, vector<long long> & T_vu, unsigned long long & tri_count,
            		vector<long long> & W_u, unsigned long long & w_local_count, vector<long long> & ind);
            void clique(unsigned long long & tri_count, vector<long long> & T_vu, unsigned long long & clique4_count, long long & v, vector<long long> & ind);
            void cycle(unsigned long long & w_local_count, vector<long long> & W_u, unsigned long long & cycle4_count, long long & v, vector<long long> & ind);


            void cycle_adj_micro(unsigned long long & w_local_count, vector<long long> & W_u,
            		unsigned long long & cycle4_count, long long & v, vector< vector<short int> > & adj_mat,
            		unsigned long long & tailed_tri_4_count);

            void cycle_micro(unsigned long long & w_local_count, vector<long long> & W_u,
            		unsigned long long & cycle4_count, long long & v, vector<long long> & ind,
            		unsigned long long & tailed_tri_4_count);

            /**
             * @brief Create perfect hash table for checking edge existence in O(1) time.
             * Takes only O(|N(v)|) time to create perfect hash table since we only need to mark the neighbors of v, which can be done in O(1) by directly indexing...
             */
            void mark_neighbors(long long & v, long long & u, vector<long long> & ind);
            /**
             * @brief Reset perfect hash table in O(|N(v)|) time
             */
            void reset_perfect_hash(long long & v, vector<long long> & ind);
            void reset_graphlet_counts();

            /**
             * @brief Check correctness of graphlet counts
             *
             * @param thread_n_count is an array containing the counts for each graphlet motif
             */
            bool test_graphlet_counts(vector< vector<unsigned long long> > & thread_n_count);

            /** @brief maximum k-core number of G*/
            int max_core;

            /** @brief k-core number of each vertex $v \in V$ */
            vector<int> kcore;

            /** @brief degeneracy order */
            vector<int> kcore_order;

            /**
             * @brief K-core number of each vertex
             *
             * @return
             */
            vector<int>* get_kcores() { return &kcore; }

            /**
             * @brief Get the k-core number of vertex v
             *
             * @param v is the vertex to use for finding the k-core number
             * @return Returns the k-core number of the vertex v (given as input)
             */
            int get_kcore_number(int v) { return kcore[v]; }

            /**
             * @brief Get the degereracy order
             *
             * @return Returns the array containing the degeneracy ordering
             */
            vector<int>* get_kcore_ordering() { return &kcore_order; }

            /**
             * @brief Get the maximum k-core number of the graph G denoted as K(G).
             *
             * @return Returns K(G) -- the maximum k-core number of the graph G
             */
            int get_max_core() { return max_core; }

            /**
             * @brief Compute k-core decomposition fast for large graphs.
             * Takes only $O(|E|)$ time and leverages bucket sort.
             * Many components are parallelized for best performance.
             *
             * This function computes the
             * 	+ k-core number for each vertex (stored in the kcore array)
             * 	+ degeneracy order (also known as k-core order)
             * 	+ maximum k-core number of G, that is, $K(G)$.
             *
             * NOTE: This function is a prereq for the following functions:
             * 	- get_kcore_number(),
             * 	- get_max_core(),
             * 	- get_kcores(),
             * 	- get_kcore_ordering()
             *
             */
            void compute_cores();

            /** @brief assortativity */
            double r;
            /**
             * @brief Get the assortativity coefficient r
             *
             * @return r is the assortativity coefficient
             */
            double get_assortativity() { return r; }

            /**
             * @brief Compute assortativity using a fast parallel edge-centric algorithm.
             *
             * This assortativity approach has the following advantages:
             * 	+ Computes assortivity in an edge-centric fashion
             * 	+ Parallel algorithm for computing assortivity fast for very large graphs
             *
             * NOTE: It also may be used for computing the assortativity of each edge or vertex in parallel.
             */
            void compute_assortativity();

            /**
             * @brief print basic stats,
             * basic_stats("", "\n"), or basic_stats("", ", ");
             */
            void basic_stats(string prefix = "", string suffix = "\n");

            /**
             * @brief Get the file extension from the path
             *
             * @param filename is the string containing the path/filename
             * @return extension of the path or filename
             */
            string get_file_extension(const string& filename);

            /**
             * @brief Get the filename (name of file and its extension) from the path
             *
             * @param s is the path "/path/to/filename.txt"
             * @return Name of file and extension, that is, any path information discarded
             */
            string get_filename_from_path(const string& s);
    };
}
#endif
