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

#include "graphlet_core.h"
#include <algorithm>

using namespace graphlet;
using namespace std;

void graphlet_core::initialize() {
    max_degree = 0;
    min_degree = 0;
    avg_degree = 0;
    max_core = 0;
    max_t_edge = 0;
    total_t = 0;
    fn = "";
    is_gstats = false;
    is_adj = false;
    verbose = false;
    block_size = 64;
    density_threshold = 0.8;
}
graphlet_core::~graphlet_core() { }

/**
 * @brief Constructor to initialize params, including the parallel params, and reads the graph
 *
 */
graphlet_core::graphlet_core(params & p) {
    initialize();
    fn = p.graph;
    is_gstats = p.graph_stats;
    block_size = p.block_size;
    read_graph(p.graph);

    preprocess(p.ordering_csc_neighbor, p.is_small_to_large_csc_neighbor, p.adj_limit, p.graph_representation);
}

/**
 * @brief Order neighbors of each vertex in csc/csr rep., create e_v and e_u arrays for fast lookup and optimize alg/data struct, etc.
 *
 * @param ordering_csc_neighbor
 * @param is_small_to_large
 * @param adj_limit
 * @param rep
 */
void graphlet_core::preprocess(string ordering_csc_neighbor, bool is_small_to_large, int adj_limit, string rep) {
    /// ordering neighbors of each vertex
    order_vertex_neighbors(ordering_csc_neighbor, is_small_to_large);

    /**
     * @brief e_v and e_u are used for PGD, |e_v| = |E|/2.
     * edge-csc (for fast computations on the edges)
     *
     * NOTE: Create ordered e_v and e_u arrays by utilizing the fact that csc neighbors are ordered
     */
	create_edge_list_arrays();

	/**
	 * @brief Automatically choose the optimal graph representation
	 */
    optimize_graph_ops(adj_limit, rep);
}

/**
 * @brief Initialize and read graph with various custom settings
 *
 * @param filename of graph to read
 */
graphlet_core::graphlet_core(const string& filename, string ordering_technique, bool is_small_to_large, int adj_limit, string rep, int block_sz) {
    initialize();
    fn = filename;
    block_size = block_sz;
    read_graph(filename);
    preprocess(ordering_technique, is_small_to_large, adj_limit, rep);
}


/**
 * @brief Initialize and read graph!
 *
 * @param filename of graph to read
 */
graphlet_core::graphlet_core(const string& filename) {
    initialize();
    fn = filename;
    read_graph(filename);
    preprocess();
}
graphlet_core::graphlet_core(bool graph_stats, const string& filename) {
    initialize();
    fn = filename;
    is_gstats = graph_stats; // construct e_u and e_v
    read_graph(filename);
    preprocess();
}
graphlet_core::graphlet_core(int nverts, int* heads, int* tails) {
    initialize();
    int num_vs = nverts, num_es = heads[nverts];
    vector<long long> V(nverts,0);
    vector<int> E;
    E.reserve(num_es);
    int start = 0;
    for (int i = 0; i < nverts; i++) {
        start = E.size();
        for (long long j = heads[i]; j < heads[i + 1]; j++ ) { E.push_back(tails[j]); }
        V[i] = start;
        V[i + 1] = E.size();
    }
    vertices = V;
    edges = E;
    vertex_degrees();
}

/**
 * pgd library interface
 *
 * Basic assumptions:
 *  - for undirected graphs, we assume one unique pair for each edge
    - vertex ids are assumed to start at 0
 */
graphlet_core::graphlet_core(int nverts, int nedges, std::pair<int,int>* edge_pairs) {
    initialize();
    int num_vs = nverts, num_es = nedges;
    map< int, vector<int> > vert_list;
    int v = 0, u = 0, self_edges = 0;

    for (int i = 0; i < num_es; ++i) {
        v = edge_pairs[i].first;
        u = edge_pairs[i].second;
        if (v == u)  self_edges++;
        else {
            vert_list[v].push_back(u);
            vert_list[u].push_back(v);
            if (is_gstats) {
                e_v.push_back(v);
                e_u.push_back(u);
            }
        }
    }

    vertices.push_back(edges.size());
    for (int i=0; i < vert_list.size(); i++) {
        edges.insert(edges.end(),vert_list[i].begin(),vert_list[i].end());
        vertices.push_back(edges.size());
    }

    vert_list.clear();
    vertex_degrees();
    if (verbose) cout << "self-loops: " << self_edges <<endl;
    vert_list.clear();
}

bool graphlet_core::detect_weighted_graph(string & line, string & delim) {
    int num_tokens = 0;
    string buf; // Have a buffer string
    stringstream ss(line); // Insert the string into a stream
    vector<string> tokens; // Create vector to hold our words

    while (ss >> buf) { tokens.push_back(buf); }
    if (verbose) printf("number of tokens in line = %lu \n", tokens.size());
    if (tokens.size() == 3) return true; // weighted graph (3rd column)
    return false;   // unweighted, only edge list with two columns
}

void graphlet_core::detect_delim(string & line, string & delim) {
    if (delim == "") {
        std::size_t prev = 0, pos;
        std::string tab_spaces("   ");
        if ((pos = line.find_first_of(',', prev)) != std::string::npos) {
            delim = ',';
        }
        else if ((pos = line.find_first_of('\t', prev)) != std::string::npos) {
            delim = '\t';
        }
        else if ((pos = line.find(tab_spaces)) != std::string::npos) {
            printf("found tab-spaces delimiter \n");
            delim = "   ";
        }
        else if ((pos = line.find_first_of(' ', prev)) != std::string::npos) {
            delim = ' ';
        }
    }

    if (delim == "") {
        if (get_file_extension(fn) == "csv")
            delim = ',';
        else if (get_file_extension(fn) == "tab")
            delim = '\t';
        else if (get_file_extension(fn) == "mtx")
            delim = ' ';
        else
            delim = ' ';
        if (verbose) cout << "[pgd]  no delimiter recognized, using \"" << delim.c_str()[0] << "\" as delimiter!" <<endl;
    }
    else
        if (verbose) cout << "[pgd]  detected \"" << delim << "\" as delimiter! " <<endl;
}

inline
void graphlet_core::get_token(int & v, string & line, string & delim, size_t & pos, size_t & prev) {
    if ((pos = line.find(delim, prev)) != std::string::npos) {
        if (pos > prev) {
            v = atoi(line.substr(prev, pos-prev).c_str());
        }
        prev = pos+1;
    }
    else if (prev < line.length())
        v = atoi(line.substr(prev, std::string::npos).c_str());
}

inline
void graphlet_core::get_token(double & weight, string & line, string & delim, size_t & pos, size_t & prev, bool & is_weighted_graph) {
    if ((pos = line.find(delim, prev)) != std::string::npos) {
        if (pos > prev) {
            weight = atof(line.substr(prev, pos-prev).c_str());
        }
        prev = pos+1;
    }
    else if (prev < line.length())
        weight = atof(line.substr(prev, std::string::npos).c_str());
}

/**
 * @brief Selects appropriate reader for graph
 *
 * @param filename or path of the graph to read
 */
void graphlet_core::read_graph(const string& filename) {
    is_gstats = true;
    fn = filename;
    double sec = get_time();
    string ext = get_file_extension(filename);

    if (verbose) cout << "[pgd: graph reader]  All graphs are assumed to be undirected" <<endl;
    if (verbose) cout << "[pgd: graph reader]  Self-loops and weights (if any) are discarded" <<endl;

    if (ext == "edges" || ext == "eg2" || ext == "txt" || ext == "csv") {
        if (verbose) cout << "[pgd: general graph reader]  reading the edge list" <<endl;
        read_edge_list(filename);
    }
    else if (ext == "mtx") {
        if (verbose) cout << "[pgd: mtx graph reader]  Assuming matrix is undirected, and upper-triangular " <<endl;
        read_mtx(filename);
    }
    else {
        if (verbose) cout << "[pgd: general graph reader] Reading the graph" <<endl;
        read_edge_list(filename);
    }
    if (verbose) cout << "Reading time " << get_time() - sec << endl;
    vertex_degrees();
    basic_stats();
}

/**
 * /brief Reads a general edge list, makes limited assumptions about the graph
 *
 * WEIGHTS: All weights are discarded, unless the graph is temporal
 * LABELS:  Vertices are relabeled, and the old ids are discarded, unless specified.
 */
void graphlet_core::read_edge_list(const string& filename) {
	cout << "[reading generic edge list: read_edge_list func]  filename: " << filename <<endl;
    map< int, vector<int> > vert_list;
    int v = -1, u = -1, num_es = 0, self_edges = 0;
    double weight;
    string delimiters = " ,\t", delim="", line="", token="";
    string graph_exif = "";

    ifstream file (filename.c_str());
    if (!file) { if (verbose) cout << filename << "File not found!" <<endl; return; }

    // check if vertex ids start at 0 or 1
    is_weighted = false;
    bool fix_start_idx = true, ignore_first_line = false;
    stringstream iss;

    // save graph info/comments at top of file
    while(std::getline(file, line) && (line[0] == '#' || line[0] == '%')) {
        graph_exif += line;
        if (line.find("MatrixMarket matrix coordinate pattern symmetric") != std::string::npos) {
            delim = ' ';
            ignore_first_line = true;
        }
    }

    int num_verts = 0, num_edges = 0;
    if (get_file_extension(filename) == ".mtx") {
        if (verbose) cout << "[pgd: graph reader]  mtx file detected!" <<endl;
        iss << line;
        int cols = 0;
        iss >> num_verts >> cols >> num_edges;
        if(num_verts!=cols) { cout<<"[pgd]  error: this is not a square matrix, attempting to proceed."<<endl; }
    }

    // detect the delim for reading the graph
    detect_delim(line,delim);

    // detect if line has three columns, third is assumed to be for weights
    is_weighted = detect_weighted_graph(line,delim);
    if (is_weighted)    printf("weighted graph detected \n");


    // handle the first line (find starting vertex id)
    if (line != "") {
        iss.clear();
        iss.str(line);
        iss >> v >> u; //>> weight;
        if (v == 0 || u == 0) { fix_start_idx = false; }
    }

    int max = 0; // largest vertex id (assumed to be ints)
    if (verbose) cout << "[pgd: graph reader]  reading a general edge list (limited assumptions)" <<endl;
    // find starting vertex id, compute the number of vertices to expect (since gaps in vertex ids are possible)
    while(std::getline(file, line)) {
        if (line != "") { // ensure line actually contains data
            iss << line;

            // ignore comments
            if (line[0] == '%' || line[0] == '#') continue;

            std::size_t prev = 0, pos;
            get_token(v,line,delim,pos,prev);
            get_token(u,line,delim,pos,prev);

            if (v == 0 || u == 0) { fix_start_idx = false; }
            if (v > max) max = v;
            if (u > max) max = u;
        }
    }
    if (verbose) cout << "[pgd: graph reader]  largest vertex id is " << max <<endl;

    file.close();
    if (verbose) {
	    if (fix_start_idx) cout << "[pgd: graph reader]  vertex ids from the file begin at 1" <<endl;
	    else cout << "[pgd: graph reader]  vertex ids begin at 0" <<endl;
    }

    ifstream fin (filename.c_str());
    if (!fin) { cout << filename << "Error: file not found!" <<endl; return; }

    int vertex_id = 0;
    vector<int> vert_lookup(max+1,-1);
    if (is_weighted) {
        while(std::getline(fin, line)) {
            if (line != "") { // ensure line actually contains data
                iss << line;

                // ignore comments
                if (line[0] == '%' || line[0] == '#') continue;

                std::size_t prev = 0, pos; // prev is last location in the line
                get_token(v,line,delim,pos,prev);
                get_token(u,line,delim,pos,prev);

                // get the weight (3rd column)
                get_token(weight,line,delim,pos,prev,is_weighted);
                if (fix_start_idx) {
                    v--;
                    u--;
                }
                if (v == u)  self_edges++;
                else {
                    if (vert_lookup[v] == -1) { // new vertex
                        vert_lookup[v] = vertex_id; // store the new id
                        vertex_id++;
                    }
                    v = vert_lookup[v]; // get the consistent vertex id

                    if (vert_lookup[u] == -1) { // new vertex
                        vert_lookup[u] = vertex_id; // store the new id
                        vertex_id++;
                    }
                    u = vert_lookup[u]; // get the consistent vertex id
                    vert_list[v].push_back(u);
                    vert_list[u].push_back(v);
                }
            }
        }
    }
    else { // unweighted graph (two columns)
        while(std::getline(fin, line)) {
            if (line != "") { // ensure line actually contains data
                iss << line;
                if (line[0] == '%' || line[0] == '#') continue;

                std::size_t prev = 0, pos; // prev is last location in the line
                get_token(v,line,delim,pos,prev);
                get_token(u,line,delim,pos,prev);
                if (fix_start_idx) {
                    v--;
                    u--;
                }
                if (v == u)  self_edges++;
                else {
                    if (vert_lookup[v] == -1) { // new vertex
                        vert_lookup[v] = vertex_id; // store the new id
                        vertex_id++;
                    }
                    v = vert_lookup[v]; // get the consistent vertex id

                    if (vert_lookup[u] == -1) { // new vertex
                        vert_lookup[u] = vertex_id; // store the new id
                        vertex_id++;
                    }
                    u = vert_lookup[u]; // get the consistent vertex id

                    vert_list[v].push_back(u);
                    vert_list[u].push_back(v);
                }
            }
        }
    }
    fin.close();
    vert_lookup.clear();

    if (verbose) cout << "vert_list size: " << vert_list.size() <<endl;

    vertices.push_back(edges.size());
    for (int i=0; i < vert_list.size(); i++) {
    	edges.insert(edges.end(),vert_list[i].begin(),vert_list[i].end());
    	vertices.push_back(edges.size());
    }
    vert_list.clear();
    bool delete_multiple_edges = true; /// todo: set via command line
    if (delete_multiple_edges) { remove_multiple_edges(); }
    if (verbose) cout << "self-loops: " << self_edges <<endl;
}

int graphlet_core::read_mtx(const string& filename) {
    string line = "";
    map< int, vector<int> > vert_list;
    int v = -1, u = -1, num_es = 0, self_edges = 0;
    string delim = " ", graph_exif = "";

    ifstream file (filename.c_str());
    if (!file) { cout << filename << "Error: file not found!" <<endl; return 0; }

    // check if vertex ids start at 0 or 1
    bool fix_start_idx = true, ignore_first_line = false;
    string token;
    stringstream iss;

    // save graph info/comments at top of file
    while(std::getline(file, line) && (line[0] == '#' || line[0] == '%')) { graph_exif += line; }

    int num_verts = 0, num_edges = 0, cols = 0;
    iss << line;
    iss >> num_verts >> cols >> num_edges;
    if(num_verts!=cols) cout<<"[pgd]  error: this is not a square matrix, attempting to proceed."<<endl;

    // handle the first line (find starting vertex id)
    if (line != "") {
        iss.clear();
        iss.str(line);
        iss >> v >> u;
        if (v == 0 || u == 0) { fix_start_idx = false; }
    }

    double value;
    if (verbose) cout << "[pgd: graph reader]  reading mtx file" <<endl;
    /*
     * @brief Assume there are no gaps in vertex ids (mtx standard), break if a zero is encountered
     * Find starting vertex id (break asap)
     */
    while(std::getline(file, line) && fix_start_idx) {
        if (line != "") { // ensure line actually contains data

            // ignore comments
            if (line[0] == '%' || line[0] == '#') continue;

            iss.clear();
            iss.str(line);
            iss >> v >> u >> value;

            v--;
            u--;
            if (v == u)  self_edges++;
            else {
                vert_list[v].push_back(u);
                vert_list[u].push_back(v);
            }
        }
    }
    vertices.push_back(edges.size());
    for (int i=0; i < vert_list.size(); i++) {
        edges.insert(edges.end(),vert_list[i].begin(),vert_list[i].end());
        vertices.push_back(edges.size());
    }
    vert_list.clear();
    if (verbose) cout << "self-loops: " << self_edges <<endl;
    return 1;
}

/**
 * @brief Specifically designed to be _fast_ for very large graphs
 * Impossible to store full adj of large sparse graphs, instead
 * we create a lookup table for each vertex, and build it on the fly,
 * using this info to mark and essentially remove the multiple edges
 */
void graphlet_core::remove_multiple_edges() {
    vector<int> ind(vertices.size(),0);
    vector<long long> vs(vertices.size(),0);
    vector<int> es;
    es.reserve(edges.size());

    int start = 0;
    for (int i = 0; i < vertices.size()-1; i++) {
        start = es.size();
        for (long long j = vertices[i]; j < vertices[i + 1]; j++) {
            int u = edges[j];
            if (ind[u] == 0) {
                es.push_back(edges[j]);
                ind[u] = 1;
            }
        }
        vs[i] = start;
        vs[i + 1] = es.size();
        for (long long j = vertices[i]; j < vertices[i + 1]; j++) { ind[edges[j]] = 0; }
    }
    if (verbose) cout << "[pgd: graph reader]  removed " << (edges.size() - es.size())/2 << " duplicate edges (multigraph)" <<endl;
    if (verbose) cout << "[remove multiple edges] " << edges.size() <<endl;
    vertices = vs;
    edges = es;
    vs.clear();
    es.clear();
}

/**
 * @brief Output basic graph stats such as |V|, |E|, density, max degree, avg degree
 * Examples: basic_stats("", "\n"), or basic_stats("", ", ");
 */
void graphlet_core::basic_stats(string prefix, string suffix) {
    cout << prefix << "|V|: " << num_vertices() <<suffix;
    cout << prefix << "|E|: " << num_edges() <<suffix;
    cout << prefix << "p: " << density() <<suffix;
    cout << prefix << "d_max: " << get_max_degree() <<suffix;
    cout << prefix << "d_avg: " << get_avg_degree() <<suffix;
}

string graphlet_core::basic_names_line(string delim, string prefix, string suffix) {
	ostringstream str_stream;
    str_stream << prefix << "|V|" << delim <<suffix;
    str_stream << prefix << "|E|" << delim <<suffix;
    str_stream << prefix << "p" << delim <<suffix;
    str_stream << prefix << "d_max" << delim <<suffix;
    str_stream << prefix << "d_avg" << delim <<suffix;
    if (max_core>0) str_stream << prefix << "K" << delim <<suffix;
    return str_stream.str();
}

string graphlet_core::basic_stats_line(string delim, string prefix, string suffix) {
	ostringstream str_stream;
    str_stream << prefix << num_vertices() << delim <<suffix;
    str_stream << prefix << num_edges() << delim <<suffix;
    str_stream << prefix << density() << delim <<suffix;
    str_stream << prefix << get_max_degree() << delim <<suffix;
    str_stream << prefix << get_avg_degree() << delim <<suffix;
    if (max_core>0) str_stream << prefix << max_core << delim <<suffix;
    return str_stream.str();
}

/**
 * @brief       adapt graph representation based on some simple statistics
 *              and memory requirements. This is meant for performance and
 *              flexibility.
 *
 * @todo        remove this from the user, simply decide when graph is read
 *              or created via other graph operations
 *
 * @param adj_limit Threshold for deciding if an adjacency matrix should be constructed for checking edge existence in O(1) time.
 */
void graphlet_core::optimize_graph_ops(int adj_limit, string rep) {
    if (num_vertices() < adj_limit && density() > density_threshold) {
        create_adj_mat();
        if (verbose) cout << "DENSE and/or SMALL graph detected.  Optimizing internal data structures" <<endl;
    }
    else {
        if (verbose) cout << "SPARSE graph detected.  Optimizing internal data structures" <<endl;
    }
    set_rep_manually(rep);
}

/**
 * @brief Creates an adjacency matrix for checking edge existence in O(1)
 *
 * Used for graphlet decomposition of DENSE and/or small networks (offers speedup for both types)
 *
 * \return An adjacency matrix "A" where A[i][j] = 1 if $(i,j) \in E$
 */
void graphlet_core::create_adj_mat() {
	if (is_adj==false) {
		is_adj = true;
		double sec = get_time();
		int size = num_vertices();
		A.resize(size);
		for (int i = 0; i < size; i++) { A[i].resize(size,0); }
		for (int i = 0; i < size; i++) {
			for (long long j = vertices[i]; j < vertices[i + 1]; j++ )
				A[i][edges[j]] = 1;
		}
		if (verbose) cout << "Created adjacency matrix in " << get_time() - sec << " seconds" <<endl;
	}
}

enum  { // representation
	AUTO, 	// 0
	CSC,        // 1
	ADJ,        // 2
	HYBRID,     // 3
};

/**
 * @brief Get the appropriate representation id/enum.
 * Note that auto indicates to automatically determine optimal (default).
 *
 * @param rep
 * @return
 */
static int get_representation_enum(string &rep) {
	if 	  	(rep == "auto" || rep == "")			return 0;
	else if (rep == "csc" || rep == "sparse")		return 1;
	else if (rep == "adj" || rep == "adjcency")		return 2;
	else if (rep == "hybrid" || rep == "csc_adj")	return 3;
	return 0;
}

/**
 * @brief Allows the user to set the representation manually (for specific or special cases).
 * NOTE that if adj structure already exists, then we return immediately to avoid recomputing the data structure.
 *
 * @param rep is a string indicating the type of representation to use
 */
void graphlet_core::set_rep_manually(string &rep) {
	switch(get_representation_enum(rep)) {
		case AUTO: {
			break;
		}
		case CSC: {
			is_adj = false;
			break;
		}
		case ADJ: {
			create_adj_mat();
			break;
		}
		case HYBRID: {
			create_adj_mat();
			break;
		}
		default: {
			break;
		}
	}
}

enum  { // order and pruning methods
	NATURAL,            // 0
	RAND,               // 1
	DEGREE,             // 2
	KCORE,              // 3
	KCORE_DEG,          // 4
	DEGREE_VOL,         // 5
	KCORE_VOL,          // 6
	DEGREE_KCORE_VOL,   // 7
	VAR,                // 8
	TRIANGLES,          // 9
	WEDGES,             // 10
	TRIANGLE_CORES,     // 11
	COLORING,           // 12
	TRIANGLE_VOL,       // 13
	TRIANGLE_CORE_VOL,  // 14
	TRIANGLES_ONLY,     // 15 (for global pruning)
	DYNAMIC_LARGEST_FIRST,          // 16
	TRIANGLE_CORE_MAX,              // 17
	DEGREE_TRIANGLES,               // 18
	KCORE_TRIANGLES,                // 19
	KCORE_DEG_TRI,                  // 20
	KCORE_TRIANGLE_VOL,             // 21
	DEGREE_KCORE_TRIANGLE_VOL,      // 22
	DIST_TWO_LARGEST_FIRST,         // 23
	DIST_TWO_DYNAMIC_LARGEST_FIRST, // 24
	DIST_TWO_SMALLEST_LAST,         // 25
	INCIDENCE_DEGREE,               // 26
	DIST_TWO_INCIDENCE_DEGREE,      // 27
	TRIANGLE_CORE_DIST_TWO,         // 28
	TRIANGLE_CORE_MIN,              // 29
	COMMON_NEIGHBORS,			  	// 30
};

static int get_ordering_enum(string &ordering) {
	if (ordering == "natural")                                                     	return 0;
	else if (ordering == "rand" || ordering == "random")                                return 1;
	else if (ordering == "deg" || ordering == "degree")                            	return 2;
	else if (ordering == "kcore" || ordering == "kcores")                               return 3;
	else if (ordering == "kcore_deg" || ordering == "kcore_degree")                     return 4;
	else if (ordering == "deg_vol" || ordering == "degree_vol")                         return 5;
	else if (ordering == "kcore_vol" || ordering == "kcore_volume")                     return 6;
	else if (ordering == "kcore_deg_vol" || ordering == "kcore_degree_vol")             return 7;
	else if (ordering == "var")                                                         return 8;

	else if (ordering == "triangles" || ordering == "tri" || ordering == "triangle")    return 9;
	else if (ordering == "wedges")                                                      return 10;
	else if (ordering == "tcores" || ordering == "triangle_cores" ||
			ordering == "triangle_core" || ordering == "tcore")                     return 11;
	else if (ordering == "coloring" || ordering == "greedy_coloring")                   return 12;
	else if (ordering == "tri_vol" || ordering == "triangle_vol")                       return 13;
	else if (ordering == "tcore_vol" || ordering == "tcore_volume" ||
			ordering == "triangle_core_vol")                                        return 14;
	else if (ordering == "triangles_only" || ordering == "triangle_only")               return 15;
	else if (ordering == "dynamic_largest_first" || ordering == "lfo" ||
			ordering =="LFO")  									return 16; // graph coloring
	else if (ordering == "triangle_core_max" || ordering == "tcore_max" ||
			ordering == "tcores_max")                                               return 17;
	else if (ordering == "degree_triangles" || ordering == "degree_tri" ||
			ordering == "deg_tri" || ordering == "deg_triangles")                   return 18;
	else if (ordering == "kcore_tri" || ordering == "kcore_triangles" ||
			ordering =="kcore_triangle")                                            return 19; // graph coloring
	else if (ordering == "kcore_deg_tri" || ordering == "kcore_degree_triangle" ||
			ordering == "kcore_degree_tri")                                         return 20;
	else if (ordering == "kcore_tri_vol" || ordering == "kcore_triangle_vol")           return 21;
	else if (ordering == "kcore_deg_tri_vol" || ordering == "kcore_degree_tri_vol")     return 22;
	else if (ordering == "dist_two_lfo" || ordering == "dist_two_largest_first")        return 23;
	else if (ordering == "dist_two_dlfo" || ordering == "dist_two_dynamic_lfo")         return 24;
	else if (ordering == "dist_two_slo" || ordering == "dist_two_smallest_last")        return 25;
	else if (ordering == "ido" || ordering == "incidence" ||
			ordering == "incidence_deg")                                            return 26;
	else if (ordering == "dist_two_ido" || ordering == "dist_two_incidence")            return 27;
	else if (ordering == "dist_two_tcore" || ordering == "dist_two_triangle_core")      return 28;
	else if (ordering == "tcore_min" || ordering == "triangle_core_min" ||
			ordering =="tri_core_min")                                              return 29; // graph coloring
	else if (ordering == "common" || ordering=="common_neighbors")				return 30;
	else return 0; // default natural
}

/**
 * @brief
 *
 * @param ordering enum from get_ordering_enum
 * @param v is a vertex id for edge e (src)
 * @param u is a vertex id for edge e (dst)
 * @param val that is modified and returned via reference
 */
inline
void graphlet_core::get_ordering_value(int &ordering, long long &v, long long &u, long long &val) {
	switch (ordering) {
		case NATURAL: {
			break;
		}
		case RAND: {
			val = get_cust_rand_int();
			break;
		}
		case DEGREE: {
			val = degree[u]+degree[v];
			break;
		}
		case KCORE: {
			val = kcore[u]+kcore[v];
			break;
		}
		case KCORE_DEG: {
			val = (degree[u]+degree[v]) * (kcore[u]+kcore[v]);
			break;
		}
		case DEGREE_VOL: {
			val = 0;
			for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
				val = val + vertices[edges[j] + 1] - vertices[edges[j]];
			}
			for (long long j = vertices[v]; j < vertices[v + 1]; j++) {
				val = val + vertices[edges[j] + 1] - vertices[edges[j]];
			}
			break;
		}
		case KCORE_VOL: {
			val = 0;
			for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
				val = val + kcore[edges[j]];
			}
			for (long long j = vertices[v]; j < vertices[v + 1]; j++) {
				val = val + kcore[edges[j]];
			}
			break;
		}
		case DEGREE_KCORE_VOL: {
			val = 0;
			for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
				val = val + ((vertices[edges[j] + 1] - vertices[edges[j]]) * kcore[edges[j]]);
			}
			for (long long j = vertices[v]; j < vertices[v + 1]; j++) {
				val = val + ((vertices[edges[j] + 1] - vertices[edges[j]]) * kcore[edges[j]]);
			}
			break;
		}
		case VAR: {
			val = (kcore[u] * ((int)degree[u]/kcore[u])) + (kcore[v] * ((int)degree[v]/kcore[v]));
			break;
		}
		default: {
			val = vertices[u + 1] - vertices[u];
			break;
		}
	}
}

/** @brief Create an ordered edge list from a set of sampled edges E_s such that |E_s| < |E| */
void graphlet_core::sort_edges(string ordering_technique, bool is_small_to_large, vector<unsigned long> & E_s) {
	int ordering = get_ordering_enum(ordering_technique);
	E_ordered.reserve(E_s.size());
	set_custom_seed(get_time());

	for (int e=0; e<E_s.size(); e++) {
		int edge_id = E_s[e];
		long long v = e_v[edge_id], u = e_u[edge_id];
		long long val = e;
		get_ordering_value(ordering, v, u, val);
		E_ordered.push_back(Vertex(edge_id,val));
	}
	if (is_small_to_large) { 	std::sort(E_ordered.begin(), E_ordered.end(), incr_bound); }
	else { std::sort(E_ordered.begin(), E_ordered.end(), decr_bound); }
}

/** @brief Create an ordered edge list from FULL set of edges E */
void graphlet_core::sort_edges(string ordering_technique, bool is_small_to_large) {
	int ordering = get_ordering_enum(ordering_technique);
	int m = e_v.size();
	E_ordered.reserve(m);
	set_custom_seed(get_time());
	for (int e=0; e<m; e++) {
		long long v = e_v[e], u = e_u[e];
		long long val = e;
		get_ordering_value(ordering, v, u, val);
		E_ordered.push_back(Vertex(e,val));
	}
	if (is_small_to_large) { 	std::sort(E_ordered.begin(), E_ordered.end(), incr_bound); }
	else { std::sort(E_ordered.begin(), E_ordered.end(), decr_bound); }
}

/**
 * @brief Compute the vertex degrees, and the max vertex degree.
 *
 * Returns the vertex degrees in the "degree" array as well as:
 * - the maximum vertex degree stored in "max_degree"
 * - average vertex degree "avg_degree"
 *
 */
void graphlet_core::vertex_degrees() {
    int n = vertices.size() - 1;
    degree.resize(n);
    int max_degree_tmp = vertices[1] - vertices[0];

    #pragma omp parallel for schedule(schedule_type,block_size) \
		reduction(max:max_degree_tmp)
    for (long long v=0; v<n; v++) {
        degree[v] = vertices[v+1] - vertices[v];
        if (max_degree_tmp < degree[v])  max_degree_tmp = degree[v];
    }
    max_degree = max_degree_tmp;
    avg_degree = (double)edges.size()/n;
    return;
}

/**
 * @brief Get the extension from a full filename given as input
 *
 * @param filename a filename
 * @return file extension, e.g., "filename.txt", then ".txt" is returned.
 */
string graphlet_core::get_file_extension(const string& filename) {
    string::size_type result;
    string fileExtension = "";
    result = filename.rfind('.', filename.size() - 1);
    if(result != string::npos)
        fileExtension = filename.substr(result+1);
    return fileExtension;
}

/**
 * @brief Gets the filename from an arbitrary path
 * @param s string containing the path of the file
 * @return filename consisting of the name and extension, that is, "file.txt"
 */
string graphlet_core::get_filename_from_path(const string& s) {
   char sep = '/';
#ifdef _WIN32
   sep = '\\';
#endif
   size_t i = s.rfind(sep, s.length( ));
   if (i != string::npos) {
      return(s.substr(i+1, s.length( ) - i));
   }
   return("");
}

/**
 * @brief Create e_u and e_v arrays
 *
 */
void graphlet_core::create_edge_list_arrays() {
	int n = vertices.size() - 1;
	long long m = edges.size();
	long long unique_edges = m/2;
	e_v.reserve(unique_edges+1);
	e_u.reserve(unique_edges+1);
	for (long long v=0; v<n; v++) {
		for (long long j=vertices[v]; j<vertices[v+1]; j++) {
			long long u = edges[j];
			if (v<u) {
				if ((vertices[v+1]-vertices[v]) < (vertices[u+1]-vertices[u])) {
					e_v.push_back(v);
					e_u.push_back(u);	
				} else {
					e_v.push_back(u);
					e_u.push_back(v);
				}
			}
		}
	}
}

/**
 * \brief Compute K-Core Decomposition
 *
 * Computes degeneracy ordering, as well as
 * 	- k-core number for each vertex
 * 	- maximum k-core number
 */
void graphlet_core::compute_cores() {
    long long j;
    int n = vertices.size(), d, i, start, num, md;
    int v, u, w, du, pu, pw, md_end;
    if (kcore.size()>0) return;
    vector <int> pos(n);
    if (kcore_order.size() > 0) {
        vector<int> tmp(n,0);
        kcore = tmp;
        kcore_order = tmp;
    }
    else {
        kcore_order.resize(n);
        kcore.resize(n);
    }

    md = 0;
	#pragma omp parallel for schedule(schedule_type,block_size) reduction(max:md)
    for (int v=1; v<n; v++) {
        kcore[v] = vertices[v] - vertices[v-1];
        if (kcore[v] > md)  md = kcore[v];
    }
    md_end = md+1;
    vector < int > bin(md_end,0);
    for (v=1; v < n; v++)  bin[kcore[v]]++;
    start = 1;
    for (d=0; d < md_end; d++) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }
    for (v=1; v<n; v++) {
        pos[v] = bin[kcore[v]];
        kcore_order[pos[v]] = v;
        bin[kcore[v]]++;
    }
    for (d=md; d > 1; d--)  bin[d] = bin[d-1];
    bin[0] = 1;
    for (i=1; i<n; i++) {
        v=kcore_order[i];
        for (j=vertices[v-1]; j<vertices[v]; j++) {
            u = edges[j] + 1;
            if (kcore[u] > kcore[v]) {
                du = kcore[u];   pu = pos[u];
                pw = bin[du];    w = kcore_order[pw];
                if (u != w) {
                    pos[u] = pw;   kcore_order[pu] = w;
                    pos[w] = pu;   kcore_order[pw] = u;
                }
                bin[du]++;   kcore[u]--;
            }
        }
    }
    max_core = 0;
    for (v = 0; v < n-1; v++) {
        kcore_order[v] = kcore_order[v+1]-1;
        if (kcore[kcore_order[v]] > max_core) max_core = kcore[kcore_order[v]];
    }
    bin.clear();
    pos.clear();
}

/**
 * @brief Compute Assortativity
 *
 * Parallel assortativity algorithm for computing the assortativity coefficient "r".
 * The algorithm is edge-centric making it extremely fast and can also be used for computing the assortativity of each edge.
 * \f$r(\alpha,\beta)=\frac{\sum_i (j^\alpha_i-\bar{j^\alpha})(k^\beta_i-\bar{k^\beta})}{ \sqrt{\sum_i (j^\alpha_i-\bar{j^\alpha})^2} \sqrt{\sum_i (k^\beta_i-\bar{k^\beta})^2} }.\f$
 *
 * NOTE: block_size is the number of "jobs" to assign to a worker (at a time).
 */
void graphlet_core::compute_assortativity() {
	long long n=vertices.size()-1, v, e, deg_u, deg_v;
	double mu = 0, jd = 0, psi = 0, tau = 0; r = 0;

	if (degree.size() == n) {
		#pragma omp parallel for schedule(schedule_type,block_size) reduction(+:mu,jd,psi)
		for (e=0; e<e_v.size(); e++) {
			mu += double(degree[e_v[e]] * degree[e_u[e]]);
			jd += (0.5 * double((degree[e_v[e]] + degree[e_u[e]])));
			psi += (0.5 * double((degree[e_v[e]]*degree[e_v[e]]) + (degree[e_u[e]]*degree[e_u[e]])));
		}
	}
	else {
		#pragma omp parallel for schedule(schedule_type,block_size) reduction(+:mu,jd,psi)
		for (e=0; e<e_v.size(); e++) {
			deg_v = vertices[e_v[e]+1] - vertices[e_v[e]];
			deg_u = vertices[e_u[e]+1] - vertices[e_u[e]];
			mu += double(deg_v * deg_u);
			jd += (0.5 * double((deg_v + deg_u)));
			psi += (0.5 * double((deg_v*deg_v) + (deg_u*deg_u)));
		}
	}
	double m = e_v.size() * 1.0;
	tau = (jd / m) * (jd / m);
	r = ((mu / m) - tau) / ((psi / m) - tau);
	if (isnan(r)) { r = 1; }
	if (r < -1 || r > 1.0) { r = 1; }
	if (m == 0) { r = 0.0; }
}

/**
 * @brief Create perfect hash table for checking edge existence in O(1) time.
 *
 * Takes only O(|N(v)|) time to create perfect hash table since we only need to mark the neighbors of v, which can be done in O(1) by directly indexing...
 *
 *\param v vertex to mark neighbors, taking $O(|N(v)|)$
 *\param u vertex id to ignore
 *\param ind a perfect hash table used to mark the neighbors of v
 */
inline
void graphlet_core::mark_neighbors(long long & v, long long & u, vector<long long> & ind) {
	for (long long i = vertices[v]; i < vertices[v+1]; ++i) {
		if (edges[i]==u) { continue; }
		ind[edges[i]] = 1;
	}
}

/**
 * @brief Reset perfect hash table in O(|N(v)|) time
 *
 * Resets the perfect hash table
 *
 * \param v vertex id used to get the neighbors for resetting the lookup table "ind". Required to get the neighbors for resetting ind.
 * \param ind array of size |V| used to check the existence of a neighbor in O(1) time.
 *
 */
inline
void graphlet_core::reset_perfect_hash(long long & v, vector<long long> & ind) {
	for (long long i = vertices[v]; i < vertices[v+1]; ++i) { ind[edges[i]] = 0; }
}

/**
 * @brief Compute number of triangles and 2-stars for vertex u.
 * Note 2-stars for v is computed directly (see solve_graphlet_equations).
 * Designed for CSC/CSR and Adjacency data structure.
 *
 * @param v
 * @param u
 * @param T_vu
 * @param tri_count
 * @param W_u
 * @param w_local_count
 * @param adj_mat
 */
inline
void graphlet_core::triangles_and_wedges_adj(long long & v, long long & u, vector<long long> & T_vu, unsigned long long & tri_count,
		vector<long long> & W_u, unsigned long long & w_local_count, vector< vector<short int> > & adj_mat) {
	for (long long j = vertices[u]; j < vertices[u+1]; ++j) {
		long long w = edges[j];
		if (w==v) { continue; }
		if (adj_mat[v][w] == 1) {
			adj_mat[v][w] = 3;
			T_vu[tri_count] = w;
			tri_count++;
		}
		else {
			W_u[w_local_count] = w;
			w_local_count++;
			adj_mat[v][w] = 2;
		}
	}
}

/**
 * @brief Compute number of triangles and 2-stars for vertex u.
 * Note 2-stars for v is computed directly (see solve_graphlet_equations).
 * Designed for CSC/CSR only.
 *
 * @param v
 * @param u
 * @param T_vu
 * @param tri_count
 * @param W_u
 * @param w_local_count
 * @param adj_mat
 */
inline
void graphlet_core::triangles_and_wedges(long long & v, long long & u, vector<long long> & T_vu, unsigned long long & tri_count,
		vector<long long> & W_u, unsigned long long & w_local_count, vector<long long> & ind) {
	for (long long j = vertices[u]; j < vertices[u+1]; ++j) {
		long long w = edges[j];
		if (w==v) { continue; }
		if (ind[w] == 1){
			ind[w] = 3;
			T_vu[tri_count] = w;
			tri_count++;
		}
		else {
			W_u[w_local_count] = w;
			w_local_count++;
			ind[w] = 2;
		}
	}
}

inline
void graphlet_core::cycle_adj(unsigned long long & w_local_count, vector<long long> & W_u,
		unsigned long long & cycle4_count, long long & v, vector< vector<short int> > & adj_mat) {
	for (long long j = 0; j < w_local_count; j++) {
		long long w = W_u[j];
		for (long long i = vertices[w]; i < vertices[w+1]; ++i) {
			if (adj_mat[v][edges[i]] == 1) cycle4_count++;
		}
		W_u[j] = 0;
	}
}

inline
void graphlet_core::cycle(unsigned long long & w_local_count, vector<long long> & W_u,
		unsigned long long & cycle4_count, long long & v, vector<long long> & ind) {
	for (long long j = 0; j < w_local_count; j++) {
		long long w = W_u[j];
		for (long long i = vertices[w]; i < vertices[w+1]; ++i) {
			if (ind[edges[i]]==1) cycle4_count++;
		}
		W_u[j] = 0;
	}
}

inline
void graphlet_core::cycle_adj_micro(unsigned long long & w_local_count, vector<long long> & W_u,
		unsigned long long & cycle4_count, long long & v, vector< vector<short int> > & adj_mat,
		unsigned long long & tailed_tri_4_count) {
	for (long long j = 0; j < w_local_count; j++) {
		long long w = W_u[j];
		for (long long i = vertices[w]; i < vertices[w+1]; ++i) {
			if (adj_mat[v][edges[i]] == 1) cycle4_count++;
			else if (adj_mat[v][edges[i]]==2) tailed_tri_4_count++;
		}
		W_u[j] = 0;
	}
}

inline
void graphlet_core::cycle_micro(unsigned long long & w_local_count, vector<long long> & W_u,
		unsigned long long & cycle4_count, long long & v, vector<long long> & ind,
		unsigned long long & tailed_tri_4_count) {
	for (long long j = 0; j < w_local_count; j++) {
		long long w = W_u[j];
		for (long long i = vertices[w]; i < vertices[w+1]; ++i) {
			if (ind[edges[i]]==1) cycle4_count++;
			else if (ind[edges[i]]==2) tailed_tri_4_count++;
		}
		W_u[j] = 0;
	}
}

inline
void graphlet_core::clique_adj(unsigned long long & tri_count, vector<long long> & T_vu,
		unsigned long long & clique4_count, long long & v, vector< vector<short int> > & adj_mat) {
	for (long long tr_i=0 ; tr_i < tri_count ; tr_i++) {
		long long w = T_vu[tr_i];
		for (long long i = vertices[w]; i < vertices[w+1]; ++i) {
			if (adj_mat[v][edges[i]] == 3) { clique4_count++; }
		}
		adj_mat[v][w] = 1;
		T_vu[tr_i] = 0;
	}
}

/**
 * @brief Count the total cliques (of size k) for a given edge (or vertex).
 * Specialized for CSR/CSC (Compressed sparse row) format.
 *
 * @param tri_count is the total triangles centered at (v,u)
 * @param T_vu is an array containing the
 * @param clique4_count is the total cliques
 * @param v is the vertex id to use for counting
 * @param ind is the perfect hash table for checking in O(1) time if edge (triangle, etc) exists
 */
inline
void graphlet_core::clique(unsigned long long & tri_count, vector<long long> & T_vu,
		unsigned long long & clique4_count, long long & v, vector<long long> & ind) {
	for (long long tr_i=0 ; tr_i < tri_count ; tr_i++) {
		long long w = T_vu[tr_i];
		for (long long i = vertices[w]; i < vertices[w+1]; ++i) {
			if (ind[edges[i]] == 3) clique4_count++;
		}
		ind[w] = 0;
		T_vu[tr_i] = 0;
	}
}

/**
 * @brief Reset the counts of all graphlet motifs to zero
 */
inline void graphlet_core::reset_graphlet_counts() {
	total_2_1edge = 0.0;
	total_2_indep = 0.0;
	// connected k=3 motifs
	total_3_tris = 0.0;
	total_2_star = 0.0;
	// disconnected k=3 motifs
	total_3_1edge = 0.0;
	total_3_indep = 0.0;
	// connected k=4 motifs
	total_4_clique = 0.0;
	total_4_chordcycle = 0.0;
	total_4_tailed_tris = 0.0;
	total_4_cycle = 0.0;
	total_3_star = 0.0;
	total_4_path = 0.0;
	// disconnected k=4 motifs
	total_4_2edge = 0.0;
	total_4_1edge = 0.0;
	total_4_2star = 0.0;
	total_4_tri = 0.0;
	total_4_indep = 0.0;
}

/**
 * @brief Check correctness of graphlet counts
 *
 * @param thread_n_count is an array containing the counts for each graphlet motif
 */
inline
bool graphlet_core::test_graphlet_counts(vector< vector<unsigned long long> > & thread_n_count) {
	unsigned long long ver_n1 = 3*total_3_star + 3*total_4_tailed_tris + 4*total_4_cycle + total_4_path + 5*total_4_chordcycle + 6*total_4_clique;
	unsigned long long ver_n2 = 3*total_3_star + total_4_tailed_tris + 4*total_4_cycle + total_4_path;
	unsigned long long ver_n3 = total_4_chordcycle + 6*total_4_clique;
	unsigned long long ver_n4 = 2*total_4_tailed_tris + 4*total_4_chordcycle;
	unsigned long long ver_n5 = total_4_path + 4*total_4_cycle;
	unsigned long long ver_n6 = 3*total_3_star + total_4_tailed_tris;
	unsigned long long ver_n7 = 2*total_4_tailed_tris + 4*total_4_cycle + 2*total_4_path + 4*total_4_chordcycle + 6*total_4_clique + 2*total_4_2edge;
	unsigned long long ver_n8 = 2*total_4_2edge + total_4_1edge;
	unsigned long long ver_n9 = 2*total_4_path + 2*total_4_2star;
	unsigned long long ver_n10 = total_4_tailed_tris + 3*total_4_tri;
	int total_verify = 0;
	total_verify = ((ver_n1-thread_n_count[0][1]) > 0) + ((ver_n2-thread_n_count[0][2]) > 0) + ((ver_n3-thread_n_count[0][3]) > 0);
	total_verify += ((ver_n4-thread_n_count[0][4]) > 0) + ((ver_n5-thread_n_count[0][5]) > 0) + ((ver_n6-thread_n_count[0][6]) > 0);
	total_verify += ((ver_n7-thread_n_count[0][7]) > 0) + ((ver_n8-thread_n_count[0][8]) > 0) + ((ver_n9-thread_n_count[0][9]) > 0);
	total_verify += ((ver_n10-thread_n_count[0][10]) > 0);
	if (total_verify != 0) {
		cout << "\tver_n1 = " << ver_n1-thread_n_count[0][1] <<endl;
		cout << "\tver_n2 = " << ver_n2-thread_n_count[0][2] <<endl;
		cout << "\tver_n3 = " << ver_n3-thread_n_count[0][3] <<endl;
		cout << "\tver_n4 = " << ver_n4-thread_n_count[0][4] <<endl;
		cout << "\tver_n5 = " << ver_n5-thread_n_count[0][5] <<endl;
		cout << "\tver_n6 = " << ver_n6-thread_n_count[0][6] <<endl;
		cout << "\tver_n7 = " << ver_n7-thread_n_count[0][7] <<endl;
		cout << "\tver_n8 = " << ver_n8-thread_n_count[0][8] <<endl;
		cout << "\tver_n9 = " << ver_n9-thread_n_count[0][9] <<endl;
		cout << "\tver_n10 = " << ver_n10-thread_n_count[0][10] <<endl;
		return false;
	}
	return true;
}

/**
 * @brief Direct computation of numerous graphlet counts by leveraging relationships
 * and solutions to a number of graphlet equations.
 * Resulting counts are obtained in O(1) time.
 *
 * @param thread_n_count
 * @param n_count
 * @param thread_tmp_triangles
 * @param thread_tmp_3_star
 * @param thread_tmp_4_2edge
 * @param thread_tmp_3_1edge
 * @param deg_v
 * @param deg_u
 * @param tri_count
 * @param w_local_count
 * @param m
 * @param n
 */
inline
void graphlet_core::solve_graphlet_equations(vector< vector<unsigned long long> > & thread_n_count, vector<unsigned long long> & n_count,
		vector<unsigned long long> & thread_tmp_triangles, vector<unsigned long long> & thread_tmp_3_star,
		vector<unsigned long long> & thread_tmp_4_2edge, vector<unsigned long long> & thread_tmp_3_1edge,
		long long & deg_v, long long & deg_u, unsigned long long & tri_count, unsigned long long & w_local_count, long long & m, long long & n) {

	unsigned long long star3_count = 0, disconn_nodes = 0;
	thread_tmp_triangles[omp_get_thread_num()] += tri_count;

	star3_count = deg_v - tri_count - 1;
	star3_count = star3_count + deg_u - tri_count - 1;
	thread_tmp_3_star[omp_get_thread_num()] += star3_count;

	thread_tmp_4_2edge[omp_get_thread_num()] += (m - deg_v - deg_u + 1);
	disconn_nodes = n - (star3_count + tri_count) - 2;
	thread_tmp_3_1edge[omp_get_thread_num()] += disconn_nodes;

	n_count[1] 	= 	(tri_count+star3_count)*(tri_count+star3_count-1) / 2.0;
	n_count[2] 	= 	(star3_count*(star3_count - 1 ) / 2.0);
	n_count[3] 	= 	(tri_count*(tri_count - 1 ) / 2.0);
	n_count[4] 	=	(tri_count*star3_count);
	n_count[5] 	= 	(deg_v - tri_count -1)*(deg_u - tri_count -1);
	n_count[6] 	= 	(deg_v - tri_count -1)*(deg_v - tri_count -2) /2;
	n_count[6] 	= 	n_count[6] + ((deg_u - tri_count -1)*(deg_u - tri_count -2) /2);
	n_count[7] 	= 	(m - deg_v - deg_u + 1);
	n_count[8] 	= 	(disconn_nodes*(disconn_nodes-1) /2.0);
	n_count[9] 	= 	star3_count*disconn_nodes;
	n_count[10] = 	tri_count*disconn_nodes;

	thread_n_count[omp_get_thread_num()][1] += n_count[1];
	thread_n_count[omp_get_thread_num()][2] += n_count[2];
	thread_n_count[omp_get_thread_num()][3] += n_count[3];
	thread_n_count[omp_get_thread_num()][4] += n_count[4];
	thread_n_count[omp_get_thread_num()][5] += n_count[5];
	thread_n_count[omp_get_thread_num()][6] += n_count[6];
	thread_n_count[omp_get_thread_num()][7] += n_count[7];
	thread_n_count[omp_get_thread_num()][8] += n_count[8];
	thread_n_count[omp_get_thread_num()][9] += n_count[9];
	thread_n_count[omp_get_thread_num()][10] += n_count[10];
}

void graphlet_core::graphlet_decomposition(int max_num_workers) {
	omp_set_num_threads(max_num_workers);
	long long e, v, u, w, r, i, j, tr_i, n = num_vertices(), m = num_edges();
	vector<long long> T_vu(max_degree+1,0), W_u(max_degree+1,0);
	vector<unsigned long long> thread_tmp_triangles(max_num_workers,0), thread_tmp_3_star(max_num_workers,0),
			thread_tmp_3_1edge(max_num_workers,0), thread_tmp_cliques(max_num_workers,0),
			thread_tmp_cycles(max_num_workers,0), thread_tmp_4_2edge(max_num_workers,0);
	vector< vector<unsigned long long> > thread_n_count(max_num_workers,vector<unsigned long long>(11,0));
	vector<unsigned long long> n_count(11,0);
	if (verbose) cout << "number of workers = " << max_num_workers-1 << ", |E|/m = " << m << ", |E_ordered| = " << E_ordered.size() <<endl;
	if (E_ordered.size()==0) {
		sort_edges();
	}
	double sec = tic();
	long long nnz = E_ordered.size();
	if (is_adj_exists()) {
		if (verbose) cout << "[PGD: DENSE and/or SMALL]  Using adj_matrix" <<endl;
		vector< vector<short int> > adj_mat = get_adj();
		sec = tic();
		#pragma omp parallel for schedule(schedule_type,block_size) \
			firstprivate(T_vu,W_u, n_count,adj_mat) private(v,u,w,r,i,j,tr_i,e)
		for (e=0; e<E_ordered.size(); e++) {
			int edge_id = E_ordered[e].get_id();
			long long v = e_v[edge_id], u = e_u[edge_id];
			long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];
			unsigned long long w_local_count = 0, tri_count = 0, clique4_count = 0, cycle4_count = 0;
			adj_mat[v][u] = 0;
			triangles_and_wedges_adj(v,u,T_vu,tri_count,W_u,w_local_count,adj_mat);
			solve_graphlet_equations(thread_n_count,n_count,thread_tmp_triangles,thread_tmp_3_star,thread_tmp_4_2edge,thread_tmp_3_1edge,
					deg_v,deg_u,tri_count,w_local_count,m,n);
			cycle_adj(w_local_count, W_u, cycle4_count, v, adj_mat);
			clique_adj(tri_count, T_vu, clique4_count, v, adj_mat);
			thread_tmp_cliques[omp_get_thread_num()] += clique4_count;
			thread_tmp_cycles[omp_get_thread_num()] += cycle4_count;
			adj_mat[v][u] = 1;
		}
		toc(sec);
	}
	else {
		if (verbose) cout << "[PGD: CSC Graph Rep.]  Using optimized sparse data structs" <<endl;
		vector<long long> ind(n, 0);
		sec = tic();
		#pragma omp parallel for schedule(schedule_type,block_size) \
			firstprivate(ind,T_vu,W_u,n_count) private(v,u,w,r,i,j,tr_i,e)
		for (e=0; e<E_ordered.size(); e++) {
			int edge_id = E_ordered[e].get_id();
			long long v = e_v[edge_id], u = e_u[edge_id];
			long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];
			unsigned long long w_local_count = 0, tri_count = 0, star3_count = 0, disconn_nodes = 0, clique4_count = 0, cycle4_count = 0;
			mark_neighbors(v,u,ind);
			triangles_and_wedges(v,u,T_vu,tri_count,W_u,w_local_count,ind);
			solve_graphlet_equations(thread_n_count,n_count,thread_tmp_triangles,thread_tmp_3_star,thread_tmp_4_2edge,thread_tmp_3_1edge,
					deg_v,deg_u,tri_count,w_local_count,m,n);
			cycle(w_local_count, W_u, cycle4_count, v, ind);
			clique(tri_count, T_vu, clique4_count, v, ind);
			thread_tmp_cliques[omp_get_thread_num()] += clique4_count;
			thread_tmp_cycles[omp_get_thread_num()] += cycle4_count;
			reset_perfect_hash(v,ind);
		}
		toc(sec);
		ind.clear();
	}

	reset_graphlet_counts();

	for (int tid = 1; tid < max_num_workers; ++tid) {
		thread_tmp_triangles[0] += thread_tmp_triangles[tid];
		thread_tmp_3_star[0] += thread_tmp_3_star[tid];
		thread_tmp_3_1edge[0] += thread_tmp_3_1edge[tid];
		thread_tmp_cliques[0] += thread_tmp_cliques[tid];
		thread_tmp_cycles[0] += thread_tmp_cycles[tid];

		thread_n_count[0][1] += thread_n_count[tid][1];
		thread_n_count[0][2] += thread_n_count[tid][2];
		thread_n_count[0][3] += thread_n_count[tid][3];
		thread_n_count[0][5] += thread_n_count[tid][5];
		thread_n_count[0][4] += thread_n_count[tid][4];
		thread_n_count[0][6] += thread_n_count[tid][6];
		thread_n_count[0][7] += thread_n_count[tid][7];
		thread_n_count[0][8] += thread_n_count[tid][8];
		thread_n_count[0][9] += thread_n_count[tid][9];
		thread_n_count[0][10] += thread_n_count[tid][10];
	}
	total_2_1edge = m;
	total_2_indep = (n*(n-1) / 2.0) - m;
	total_3_tris = thread_tmp_triangles[0] / 3.0;
	total_2_star = thread_tmp_3_star[0] / 2.0;
	total_3_1edge = thread_tmp_3_1edge[0];
	total_3_indep = (n*(n-1)*(n-2)) / (3*2);
	total_3_indep = total_3_indep - (total_3_tris+total_2_star+total_3_1edge);
	total_wedges = (3*total_3_tris) + total_2_star;
	global_cc = (3*(long double)total_3_tris)/(long double)total_wedges;
	total_4_clique = thread_tmp_cliques[0] / 6.0;
	total_4_chordcycle = thread_n_count[0][3] - (6*total_4_clique);
	total_4_cycle = thread_tmp_cycles[0] / 4.0;
	total_4_path = thread_n_count[0][5] - (4*total_4_cycle);
	total_4_tailed_tris = (thread_n_count[0][4] - (4*total_4_chordcycle)) / 2.0;
	total_3_star = (thread_n_count[0][6] - total_4_tailed_tris) / 3.0;
	total_4_2edge = thread_n_count[0][7];
	total_4_2edge = total_4_2edge - (6*total_4_clique) - (4*total_4_cycle) - (4*total_4_chordcycle) - (2*total_4_path) - (2*total_4_tailed_tris);
	total_4_2edge = total_4_2edge / 2.0;
	total_4_1edge = thread_n_count[0][8] - (2*total_4_2edge);
	total_4_2star = thread_n_count[0][9] - (2*total_4_path);
	total_4_2star = total_4_2star / 2.0;
	total_4_tri = thread_n_count[0][10] - total_4_tailed_tris;
	total_4_tri = total_4_tri / 3.0;
	total_4_indep = compute_4indep(n);
	if (verbose) { test_graphlet_counts(thread_n_count);}
}

/**
 * @brief Parameterized graphlet decomposition designed to be extremely fast and
 * space-efficient for computing micro-level graphlet statistics.
 * In short, the method computes the number of different graphlet motifs for each edge (or vertex)
 *
 * @param max_num_workers
 */
void graphlet_core::graphlet_decomposition_micro(int max_num_workers) {
	omp_set_num_threads(max_num_workers);
	long long e, v, u, w, r, i, j, tr_i, n = num_vertices(), m = num_edges();
	vector<long long> T_vu(max_degree+1,0), W_u(max_degree+1,0);
	vector<unsigned long long> thread_tmp_triangles(max_num_workers,0), thread_tmp_3_star(max_num_workers,0),
			thread_tmp_3_1edge(max_num_workers,0), thread_tmp_cliques(max_num_workers,0),
			thread_tmp_cycles(max_num_workers,0), thread_tmp_4_2edge(max_num_workers,0);
	vector< vector<unsigned long long> > thread_n_count(max_num_workers,vector<unsigned long long>(11,0));
	vector<unsigned long long> n_count(11,0);
	if (verbose) cout << "number of workers = " << max_num_workers-1 << ", |E|/m = " << m << ", |E_ordered| = " << E_ordered.size() <<endl;
	if (verbose) cout << "micro stats" <<endl;

	if (E_ordered.size()==0) {
		sort_edges();
	}
	tri.resize(m,0);
	local_4_clique.resize(m,0);
	local_4_chordal_cycle.resize(m,0);
	local_4_tailed_tris.resize(m,0);
	local_4_cycle.resize(m,0);
	local_3_star.resize(m,0);
	local_4_path.resize(m,0);

	if (is_adj_exists()) {
		if (verbose) cout << "[PGD: DENSE and/or SMALL]  Using adj_matrix" <<endl;
		vector< vector<short int> > adj_mat = get_adj();
		#pragma omp parallel for schedule(schedule_type,block_size) \
			firstprivate(T_vu,W_u, n_count,adj_mat) private(v,u,w,r,i,j,tr_i,e)
		for (e=0; e<E_ordered.size(); e++) {
			long long edge_id = E_ordered[e].get_id();
			long long v = e_v[edge_id], u = e_u[edge_id];
			long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];
			unsigned long long w_local_count = 0, tri_count = 0, clique4_count = 0, cycle4_count = 0;
			adj_mat[v][u] = 0;
			triangles_and_wedges_adj(v,u,T_vu,tri_count,W_u,w_local_count,adj_mat);
			tri[e] = tri_count;

			solve_graphlet_equations(thread_n_count,n_count,thread_tmp_triangles,thread_tmp_3_star,thread_tmp_4_2edge,thread_tmp_3_1edge,
					deg_v,deg_u,tri_count,w_local_count,m,n);
			cycle_adj(w_local_count, W_u, cycle4_count, v, adj_mat);
			clique_adj(tri_count, T_vu, clique4_count, v, adj_mat);
			thread_tmp_cliques[omp_get_thread_num()] += clique4_count;
			thread_tmp_cycles[omp_get_thread_num()] += cycle4_count;
			local_4_cycle[e] = cycle4_count;
			local_4_clique[e] = clique4_count;
			local_4_chordal_cycle[e] = (tri_count*(tri_count-1) / 2.0) - local_4_clique[e];

			/** @brief number of 2-star patterns centered at v */
			deg_v = deg_v-tri_count-1;
			/** @brief number of 2-star patterns centered at u */
			deg_u = deg_u-tri_count-1;

			/** @brief compute chains */
			local_4_path[e] = (deg_v * deg_u) - local_4_cycle[e];

			/** @brief stars */
			local_3_star[e] =  (deg_v * (deg_v-1) / 2.0);
			local_3_star[e] = local_3_star[e] + (deg_u * (deg_u-1) / 2.0);
			local_3_star[e] = local_3_star[e] - local_4_tailed_tris[e];
			adj_mat[v][u] = 1;
		}
	}
	else {
		if (verbose) cout << "[PGD: CSC Graph Rep.]  Using optimized sparse data structs" <<endl;
		vector<long long> ind(n, 0);

		#pragma omp parallel for schedule(schedule_type,block_size) \
			firstprivate(ind,T_vu,W_u,n_count) private(v,u,w,r,i,j,tr_i,e)
		for (e=0; e<E_ordered.size(); e++) {
			long long edge_id = E_ordered[e].get_id();
			long long v = e_v[edge_id], u = e_u[edge_id];
			long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];
			unsigned long long w_local_count = 0, tri_count = 0, star3_count = 0, disconn_nodes = 0, clique4_count = 0, cycle4_count = 0;
			unsigned long long tailed_triangles_4 = 0;
			mark_neighbors(v,u,ind);
			triangles_and_wedges(v,u,T_vu,tri_count,W_u,w_local_count,ind);
			tri[e] = tri_count;
			solve_graphlet_equations(thread_n_count,n_count,thread_tmp_triangles,thread_tmp_3_star,
					thread_tmp_4_2edge,thread_tmp_3_1edge,
					deg_v,deg_u,tri_count,w_local_count,m,n);
			cycle(w_local_count, W_u, cycle4_count, v, ind);
			clique(tri_count, T_vu, clique4_count, v, ind);
			thread_tmp_cliques[omp_get_thread_num()] += clique4_count;
			thread_tmp_cycles[omp_get_thread_num()] += cycle4_count;
			local_4_cycle[e] = cycle4_count;
			local_4_clique[e] = clique4_count;
			local_4_chordal_cycle[e] = (tri_count*(tri_count-1) / 2.0) - local_4_clique[e];
			deg_v = deg_v-tri_count-1;
			deg_u = deg_u-tri_count-1;
			local_4_path[e] = (deg_v * deg_u) - local_4_cycle[e];
			local_3_star[e] =  (deg_v * (deg_v-1) / 2.0);
			local_3_star[e] = local_3_star[e] + (deg_u * (deg_u-1) / 2.0);
			local_3_star[e] = local_3_star[e] - local_4_tailed_tris[e];
			reset_perfect_hash(v,ind);
		}
		ind.clear();
	}
	reset_graphlet_counts();

	for (int tid = 1; tid < max_num_workers; ++tid) {
		thread_tmp_triangles[0] += thread_tmp_triangles[tid];
		thread_tmp_3_star[0] += thread_tmp_3_star[tid];
		thread_tmp_3_1edge[0] += thread_tmp_3_1edge[tid];
		thread_tmp_cliques[0] += thread_tmp_cliques[tid];
		thread_tmp_cycles[0] += thread_tmp_cycles[tid];

		thread_n_count[0][1] += thread_n_count[tid][1];
		thread_n_count[0][2] += thread_n_count[tid][2];
		thread_n_count[0][3] += thread_n_count[tid][3];
		thread_n_count[0][5] += thread_n_count[tid][5];
		thread_n_count[0][4] += thread_n_count[tid][4];
		thread_n_count[0][6] += thread_n_count[tid][6];
		thread_n_count[0][7] += thread_n_count[tid][7];
		thread_n_count[0][8] += thread_n_count[tid][8];
		thread_n_count[0][9] += thread_n_count[tid][9];
		thread_n_count[0][10] += thread_n_count[tid][10];
	}
	total_2_1edge = m;
	total_2_indep = (n*(n-1) / 2.0) - m;
	total_3_tris = thread_tmp_triangles[0] / 3.0;
	total_2_star = thread_tmp_3_star[0] / 2.0;
	total_3_1edge = thread_tmp_3_1edge[0];
	total_3_indep = (n*(n-1)*(n-2)) / (3*2);
	total_3_indep = total_3_indep - (total_3_tris+total_2_star+total_3_1edge);
	total_wedges = (3*total_3_tris) + total_2_star;
	global_cc = (3*(long double)total_3_tris)/(long double)total_wedges;
	total_4_clique = thread_tmp_cliques[0] / 6.0;
	total_4_chordcycle = thread_n_count[0][3] - (6*total_4_clique);
	total_4_cycle = thread_tmp_cycles[0] / 4.0;
	total_4_path = thread_n_count[0][5] - (4*total_4_cycle);
	total_4_tailed_tris = (thread_n_count[0][4] - (4*total_4_chordcycle)) / 2.0;
	total_3_star = (thread_n_count[0][6] - total_4_tailed_tris) / 3.0;
	total_4_2edge = thread_n_count[0][7];
	total_4_2edge = total_4_2edge - (6*total_4_clique) - (4*total_4_cycle) - (4*total_4_chordcycle) - (2*total_4_path) - (2*total_4_tailed_tris);
	total_4_2edge = total_4_2edge / 2.0;
	total_4_1edge = thread_n_count[0][8] - (2*total_4_2edge);
	total_4_2star = thread_n_count[0][9] - (2*total_4_path);
	total_4_2star = total_4_2star / 2.0;
	total_4_tri = thread_n_count[0][10] - total_4_tailed_tris;
	total_4_tri = total_4_tri / 3.0;
	total_4_indep = compute_4indep(n);
	if (verbose) { test_graphlet_counts(thread_n_count);}
}

inline
void graphlet_core::get_progress(long long &e) {
	if (e % 1000 == 0) {
		#pragma omp single nowait
		{
			std::cout << (int((e/E_ordered.size())*100)/100) << "%\r";
			std::cout.flush();
		}
	}
}

/**
 * The paramters are:
 * x the current index we are on,
 * n the number of indicies to process,
 * r the number of times we want to update the display
 * (doing it every time will cause programs with large n to slow down majorly)
 * and w for the width of the bar. And here’s the code:
 */
inline void graphlet_core::load_progress_bar(long long x, long long n, long long w) {
	if ( (x != n) && (x % (n/100+1) != 0) ) return;
	float ratio  =  x/(float)n;
	int   c      =  ratio * w;
	cout << (int)(ratio*100) << "% [";
	for (int x=0; x<c; x++) cout << "=";
	for (int x=c; x<w; x++) cout << " ";
	cout << "]\r" << flush;
}

 /**
  * @brief Get the micro-level graphlet counts for each edge (or vertex) in the graph
  *
  * Example:
  * 	src,dst,triangle,2-star,3-node-1-edge,3-node-independent,4-clique,4-chordal-cycle,4-tailed-triangle, ...
  * 	1, 2, ...
  *
  * @param filename is the output filename containing the "micro-level" motif counts (per edge/vertex)
  * @param output_id
  * @param buffer_size to store before outputting to file (128000000 = 128MB).
  *
  */
 void graphlet_core::write_micro_stats(string &filename, bool output_id, unsigned int buffer_size) {
     ofstream myfile;
     char *fn = (char*)filename.c_str();
     myfile.open(fn);
     string delim = ",";
     ostringstream str_stream;
     str_stream << "% ";

     if (output_id) { str_stream << "src"<<delim<<"dst"<<delim; }
     str_stream << "triangle"<<delim<<
		     "2-star"<<delim<<
		     "4-clique"<<delim<<
		     "4-chordal-cycle"<<delim<<
		     "4-tailed-triangle"<<delim<<
		     "4-cycle"<<delim<<
		     "3-star"<<delim<<
		     "4-path"<<"\n";

     if (buffer_size==0) {
    	 myfile << str_stream.str();
    	 if (output_id) {
    		 for (long long e=0; e<E_ordered.size(); e++) {
    			 long long edge_id = E_ordered[e].get_id();
    			 long long v = e_v[edge_id], u = e_u[edge_id];
    			 long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];

    			 myfile << (v+1) << delim << (u+1) << delim << tri[e] << delim <<
    					 get_2_star(e,deg_v,deg_u)<<delim<<
						 local_4_clique[e] << delim <<
						 local_4_chordal_cycle[e] << delim <<
						 local_4_tailed_tris[e] << delim <<
						 local_4_cycle[e] << delim <<
						 local_3_star[e] << delim <<
						 local_4_path[e] << "\n";
    		 }
    	 }
    	 else {
    		 for (long long e=0; e<E_ordered.size(); e++) {
    			 long long edge_id = E_ordered[e].get_id();
    			 long long v = e_v[edge_id], u = e_u[edge_id];
    			 long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];

    			 myfile << tri[e] << delim <<
    					 get_2_star(e,deg_v,deg_u)<<delim<<
						 local_4_clique[e] << delim <<
						 local_4_chordal_cycle[e] << delim <<
						 local_4_tailed_tris[e] << delim <<
						 local_4_cycle[e] << delim <<
						 local_3_star[e] << delim <<
						 local_4_path[e] << "\n";
    		 }
    	 }
     } else {
    	 if (output_id) {
    		 for (long long e=0; e<E_ordered.size(); e++) {
    			 long long edge_id = E_ordered[e].get_id();
    			 long long v = e_v[edge_id], u = e_u[edge_id];
    			 long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];

    			 str_stream << (v+1) << delim << (u+1) << delim;
    			 str_stream << tri[e] << delim <<
    					 get_2_star(e,deg_v,deg_u)<<delim<<
						 local_4_clique[e] << delim <<
						 local_4_chordal_cycle[e] << delim <<
						 local_4_tailed_tris[e] << delim <<
						 local_4_cycle[e] << delim <<
						 local_3_star[e] << delim <<
						 local_4_path[e];
    			 str_stream << "\n";
    			 if (str_stream.str().size() > buffer_size) { write_buffer(myfile, str_stream); };
    		 }
    	 }
    	 else {
    		 for (long long e=0; e<E_ordered.size(); e++) {
    			 long long edge_id = E_ordered[e].get_id();
    			 long long v = e_v[edge_id], u = e_u[edge_id];
    			 long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];

    			 str_stream << tri[e] << delim <<
    					 get_2_star(e,deg_v,deg_u)<<delim<<
						 local_4_clique[e] << delim <<
						 local_4_chordal_cycle[e] << delim <<
						 local_4_tailed_tris[e] << delim <<
						 local_4_cycle[e] << delim <<
						 local_3_star[e] << delim <<
						 local_4_path[e];
    			 str_stream << "\n";
    			 if (str_stream.str().size() > buffer_size) { write_buffer(myfile, str_stream); };
    		 }
    	 }
    	 write_buffer(myfile, str_stream);
     }
     myfile.close();
 }

 /**
  * @brief Get string consisting of the motif counts for each edge (node) in the graph
  * Example: string micro_stat_str = get_micro_stats();
  *
  * @param filename
  * @param output_id
  * @param buffer_size
  * @return String containing the motif counts for each edge
  */
string graphlet_core::get_micro_stats(bool output_id, string delim) {
     ostringstream str_stream;
     if (output_id) { str_stream << "src"<<delim<<"dst"<<delim; }
     str_stream << "triangle"<<delim<<
		     "2-star"<<delim<<
		     "4-clique"<<delim<<
		     "4-chordal-cycle"<<delim<<
		     "4-tailed-triangle"<<delim<<
		     "4-cycle"<<delim<<
		     "3-star"<<delim<<
		     "4-path"<<"\n";
     if (output_id) {
	     for (long long e=0; e<E_ordered.size(); e++) {
		     long long edge_id = E_ordered[e].get_id();
		     long long v = e_v[edge_id], u = e_u[edge_id];
		     long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];
		     str_stream << (v+1) << delim << (u+1) << delim;
		     str_stream << tri[e] << delim <<
				     get_2_star(e,deg_v,deg_u)<<delim<<
				     local_4_clique[e] << delim <<
				     local_4_chordal_cycle[e] << delim <<
				     local_4_tailed_tris[e] << delim <<
				     local_4_cycle[e] << delim <<
				     local_3_star[e] << delim <<
				     local_4_path[e];
		     str_stream << "\n";
	     }
     }
     else {
	     for (long long e=0; e<E_ordered.size(); e++) {
		     long long edge_id = E_ordered[e].get_id();
		     long long v = e_v[edge_id], u = e_u[edge_id];
		     long long deg_v = vertices[v+1]-vertices[v], deg_u = vertices[u+1]-vertices[u];
		     str_stream << tri[e] << delim <<
				     get_2_star(e,deg_v,deg_u)<<delim<<
				     local_4_clique[e] << delim <<
				     local_4_chordal_cycle[e] << delim <<
				     local_4_tailed_tris[e] << delim <<
				     local_4_cycle[e] << delim <<
				     local_3_star[e] << delim <<
				     local_4_path[e];
		     str_stream << "\n";
	     }
     }
     return str_stream.str();
 }

/**
 * @brief Print the micro-level graphlet statistics, that is, motif counts for each edge in the graph
 *
 * @param output_id is a boolean flag indicating whether to include the edge id. Note that if true, then first two columns is the src and dst node ids
 * @param delim is the delimiter to use (default is comma)
 */
void graphlet_core::print_micro_stats(bool output_id, string delim) {
	print_line(80); cout << "\n\nMicro graphlet statistics" <<endl; print_line(80);
	cout << get_micro_stats(output_id, delim) <<endl;
}

 void graphlet_core::write_macro_stats(string filename) {
     ofstream myfile;
     char *fn = (char*)filename.c_str();
     myfile.open(fn);
     myfile << get_graphlet_counts();
     myfile.close();
 }

/**
 * @brief Sort the neighbors of each vertex (in CSC/CSR)
 *
 * Each worker sorts the edges of a vertex, and
 * afterwards, grabs the next available vertex to sort
 *
 * Note: sort neighbors by degree (largest to smallest)
 */
void graphlet_core::bucket_sort_neighbors_parallel(vector<int> &bound, bool is_small_to_large) {
	vector<int> tmp_edges(edges.size(),0);
	int v = 0;
	double sec = get_time();
	#pragma omp parallel for schedule(schedule_type,block_size) shared(tmp_edges)
	for (v = 0; v < num_vertices(); ++v) {
		int u, md=0, md_end, d, num, n = vertices[v+1] - vertices[v] + 1;
		vector<int> vert(n), pos(n), deg(n);
		for(u=1; u<n; u++) {
			deg[u] = bound[edges[vertices[v] + (u-1)]];
			if (deg[u] > md) md = deg[u];
		}
		md_end = md+1;
		vector < int > bin(md_end,0);
		for (u=1; u < n; u++)  bin[deg[u]]++;
		int start = 1;
		for (d=0; d < md_end; d++) {
			num = bin[d];
			bin[d] = start;
			start = start + num;
		}
		for (u=1; u<n; u++) {
			pos[u] = bin[deg[u]];
			vert[pos[u]] = edges[vertices[v] + (u-1)];
			bin[deg[u]]++;
		}
		if (is_small_to_large)
			for (u=n-1; u>0; --u) { tmp_edges[vertices[v] + (u-1)] = vert[u]; }
		else  for (u = 1; u < n; ++u) { tmp_edges[vertices[v] + (u-1)] = vert[n-u]; }
	}
	edges = tmp_edges;
}

/**
 * @brief Compute graph measure for ordering if it does not exist or is stale.
 *
 * Used in "order_vertex_neighbors"
 *
 * @param ordering_technique
 */
void graphlet_core::compute_graph_measure(string &o) {
	if (o == "kcore" && kcore.size()==0) { compute_cores(); }
	else if ((o == "degree" || o == "deg") && degree.size()==0) { vertex_degrees(); }
	else if (o == "degree_kcore" || o == "deg_kcore") {
		if (kcore.size()==0) 	compute_cores();
		if (degree.size()==0)	vertex_degrees();
	}
	else if ((o == "degree_vol" || o == "deg_vol") && degree.size()==0) { vertex_degrees(); }
	else if ((o == "kcore_vol" || o == "kcore_vol") && kcore.size()==0) { compute_cores(); }
	else if (o == "degree_kcore_vol" || o == "deg_kcore_vol") {
		if (kcore.size()==0) 	compute_cores();
		if (degree.size()==0)	vertex_degrees();
	}

}

/**
 * @brief Sort the neighbors of each vertex (in CSC/CSR)
 *
 * Each worker sorts the edges of a vertex, and
 * afterwards, grabs the next available vertex to sort
 *
 * Note: sort neighbors by degree (largest to smallest)
 */
void graphlet_core::order_vertex_neighbors(string &ordering_technique, bool is_small_to_large) {
	if (verbose) cout << "csc neighbor ordering strategy: " << ordering_technique <<endl;
	compute_graph_measure(ordering_technique);

	if (ordering_technique == "kcore") {
		bucket_sort_neighbors_parallel(kcore, is_small_to_large);
	}
	else if (ordering_technique == "degree" || ordering_technique == "deg") {
		bucket_sort_neighbors_parallel(degree, is_small_to_large);
	}
	else if (ordering_technique == "degree_kcore" || ordering_technique == "deg_kcore") {
		vector<int> feature(vertices.size(),0);
		#pragma omp parallel for schedule(schedule_type, block_size) shared(feature)
		for (int i=0; i<vertices.size(); i++) { feature[i] = degree[i] + kcore[i]; }
		bucket_sort_neighbors_parallel(feature, is_small_to_large);
	}
	else if (ordering_technique == "degree_vol" || ordering_technique == "deg_vol") {
		vector<int> feature(vertices.size(),0);
		#pragma omp parallel for schedule(schedule_type, block_size) shared(feature)
		for (long long v=0; v<vertices.size(); v++)
			for (long long j = vertices[v]; j < vertices[v+1]; ++j) feature[v] += degree[edges[j]];
		bucket_sort_neighbors_parallel(feature, is_small_to_large);
	}
	else if (ordering_technique == "kcore_vol" || ordering_technique == "kcore_vol") {
		vector<int> feature(vertices.size(),0);
		#pragma omp parallel for schedule(schedule_type, block_size) shared(feature)
		for (long long v=0; v<vertices.size(); v++)
			for (long long j = vertices[v]; j < vertices[v+1]; ++j) feature[v] += kcore[edges[j]];
		bucket_sort_neighbors_parallel(feature, is_small_to_large);
	}
	else if (ordering_technique == "degree_kcore_vol" || ordering_technique == "deg_kcore_vol") {
		vector<int> feature(vertices.size(),0);
		#pragma omp parallel for schedule(schedule_type, block_size) shared(feature)
		for (long long v=0; v<vertices.size(); v++)
			for (long long j = vertices[v]; j < vertices[v+1]; ++j) feature[v] += degree[edges[j]] + kcore[edges[j]];
		bucket_sort_neighbors_parallel(feature, is_small_to_large);
	}
	else if (ordering_technique == "rand") {
		vector<int> rand_perm(vertices.size(),0);
		set_custom_seed(get_time());
		for (long long i=0; i<vertices.size(); i++) { rand_perm[i] = get_cust_rand_int(); }
		bucket_sort_neighbors_parallel(rand_perm, is_small_to_large);
	}
}

string graphlet_core::compute_connected_GFD() {
	vector<unsigned long long> G4;
	G4.push_back(total_4_clique);
	G4.push_back(total_4_chordcycle);
	G4.push_back(total_4_tailed_tris);
	G4.push_back(total_4_cycle);
	G4.push_back(total_3_star);
	G4.push_back(total_4_path);

	const char* G4_name_tmp[] = {"4-clique     ", "4-chordal-cycle", "4-tailed-tri", "4-cycle       ", "3-star       ", "4-path       "};
	vector<string> G4_names(G4_name_tmp, G4_name_tmp+G4.size());
	ostringstream os;
	unsigned long long sum = 0;
	for (int i=0; i<G4.size(); i++) { sum += G4[i]; }
	vector<double> gfd_tmp(G4.size(),0);
	for (int i=0; i<G4.size(); i++) {
		if (sum>0) gfd_tmp[i] = ((double)G4[i]/(double)sum);
		else gfd_tmp[i] = 0;
		os << G4_names[i] << "\t" << gfd_tmp[i] << "\n";
	}
	os << "\n";
	connected_GFD = gfd_tmp;
	return os.str();
}

string graphlet_core::compute_disconnected_GFD() {
	vector<unsigned long long> G4;
	G4.push_back(total_4_tri);
	G4.push_back(total_4_2star);
	G4.push_back(total_4_2edge);
	G4.push_back(total_4_1edge);
	G4.push_back(total_4_indep);
	const char* G4_name_tmp[] = {"4-node-1-tri   ", "4-node-2-star", "4-node-2-edge  ", "4-node-1-edge  ", "4-node-indep  "};
	vector<string> G4_names(G4_name_tmp, G4_name_tmp+G4.size());
	ostringstream os;

	unsigned long long sum = 0;
	for (int i=0; i<G4.size(); i++) { sum += G4[i]; }

	double cum_sum=0.0;
	vector<double> gfd_tmp(G4.size(),0);
	for (int i=0; i<G4.size(); i++) {
		if (sum>0) gfd_tmp[i] = ((double)G4[i]/(double)sum);
		else gfd_tmp[i] = 0;
		cum_sum += gfd_tmp[i];
		os << G4_names[i] << "\t" << gfd_tmp[i] << "\n";
	}
	gfd_tmp[G4.size()-1] = 1.0-cum_sum;
	os << "\n";
	disconnected_GFD = gfd_tmp;
	return os.str();
}

string graphlet_core::compute_GFD() {
	// CONNECTED MOTIFS
	vector<unsigned long long> G4;
	G4.push_back(total_4_clique);
	G4.push_back(total_4_chordcycle);
	G4.push_back(total_4_tailed_tris);
	G4.push_back(total_4_cycle);
	G4.push_back(total_3_star);
	G4.push_back(total_4_path);
	// DISCONNECTED MOTIFS
	G4.push_back(total_4_tri);
	G4.push_back(total_4_2star);
	G4.push_back(total_4_2edge);
	G4.push_back(total_4_1edge);
	G4.push_back(total_4_indep);

	const char* G4_name_tmp[] = {"4-clique     ", "4-chordal-cycle", "4-tailed-tri", "4-cycle       ", "3-star       ", "4-path       ",
			"4-node-1-tri   ", "4-node-2-star  ", "4-node-2-edge", "4-node-1-edge  ", "4-node-indep  "};
	vector<string> G4_names(G4_name_tmp, G4_name_tmp+G4.size());
	ostringstream os;
	unsigned long long sum = 0;
	for (int i=0; i<G4.size(); i++) { sum += G4[i]; }
	double cum_sum = 0;
	vector<double> gfd_tmp(G4.size(),0);
	for (int i=0; i<G4.size(); i++) {
		gfd_tmp[i] = ((double)G4[i]/(double)sum);
		os << G4_names[i] << "\t" << gfd_tmp[i] << "\n";
		cum_sum += gfd_tmp[i];
	}
	gfd_tmp[G4.size()-1] = 1.0-cum_sum; // indep prob.
	os << "\n";
	GFD = gfd_tmp;
	return os.str();
}

void graphlet_core::print_connected_GFD() {
	print_line(80); cout << "Connected Motif Graphlet Frequency Distribution (GFD)\n"; print_line(80);
	cout << compute_connected_GFD() <<endl;
}

void graphlet_core::print_disconnected_GFD() {
	print_line(80); cout << "Disconnected Motif Graphlet Frequency Distribution (GFD)\n"; print_line(80);
	cout << compute_disconnected_GFD() <<endl;
}

void graphlet_core::print_GFD() {
	print_line(80); cout << "Graphlet Frequency Distribution (GFD)\n"; print_line(80);
	cout << compute_GFD() <<endl;
}

/**
 * @brief Output frequency for each graphlet motif
 */
void graphlet_core::print_graphlet_counts() {
	print_line(60,"*");
	cout << "total_2_1edge = " << total_2_1edge <<endl;
	cout << "total_2_indep = " << total_2_indep <<endl;
	print_line(40);
	cout << "total_3_tris = " << total_3_tris <<endl;
	cout << "total_2_star = " << total_2_star <<endl;
	cout << "total_3_1edge = " << total_3_1edge <<endl;
	cout << "total_3_indep = " << total_3_indep <<endl;
	print_line(40);
	cout << "total_4_clique = " << total_4_clique <<endl;
	cout << "total_4_chordcycle = " << total_4_chordcycle <<endl;
	cout << "total_4_tailed_tris = " << total_4_tailed_tris <<endl;
	cout << "total_4_cycle = " << total_4_cycle <<endl;
	cout << "total_3_star = " << total_3_star <<endl;
	cout << "total_4_path = " << total_4_path <<endl;
	print_line(40);
	cout << "total_4_1edge = " << total_4_1edge <<endl;
	cout << "total_4_2edge = " << total_4_2edge <<endl;
	cout << "total_4_2star = " << total_4_2star <<endl;
	cout << "total_4_tri = " << total_4_tri <<endl;
	cout << "total_4_indep = " << total_4_indep <<endl;
	print_line(60,"*");
}

string graphlet_core::get_graphlet_names_line(string delim) {
	ostringstream str_stream;
	str_stream << "total_2_1edge" << delim;
	str_stream << "total_2_indep" << delim;
	// connected k=3 motifs
	str_stream << "total_3_tris" << delim;
	str_stream << "total_2_star" << delim;
	// disconnected k=3 motifs
	str_stream << "total_3_1edge" << delim;
	str_stream << "total_3_indep" << delim;
	// connected k=4 motifs
	str_stream << "total_4_clique" << delim;
	str_stream << "total_4_chordcycle" << delim;
	str_stream << "total_4_tailed_tris" << delim;
	str_stream << "total_4_cycle" << delim;
	str_stream << "total_3_star" << delim;
	str_stream << "total_4_path" << delim;
	// disconnected k=4 motifs
	str_stream << "total_4_1edge" << delim;
	str_stream << "total_4_2edge" << delim;
	str_stream << "total_4_2star" << delim;
	str_stream << "total_4_tri" << delim;
	str_stream << "total_4_indep" << delim;
	return str_stream.str();
}

string graphlet_core::get_graphlet_counts_line(string delim) {
	ostringstream str_stream;
	str_stream << total_2_1edge <<delim;
	str_stream <<  total_2_indep <<delim;
	// connected k=3 motifs
	str_stream << total_3_tris <<delim;
	str_stream <<  total_2_star <<delim;
	// disconnected k=3 motifs
	str_stream << total_3_1edge <<delim;
	str_stream << total_3_indep <<delim;
	// connected k=4 motifs
	str_stream << total_4_clique <<delim;
	str_stream << total_4_chordcycle <<delim;
	str_stream << total_4_tailed_tris <<delim;
	str_stream << total_4_cycle <<delim;
	str_stream <<  total_3_star <<delim;
	str_stream << total_4_path <<delim;
	// disconnected k=4 motifs
	str_stream << total_4_1edge <<delim;
	str_stream << total_4_2edge <<delim;
	str_stream << total_4_2star <<delim;
	str_stream << total_4_tri <<delim;
	str_stream << total_4_indep <<delim;
	return str_stream.str();
}

string graphlet_core::get_graphlet_counts() {
	ostringstream str_stream;
	str_stream << "total_2_1edge = " << total_2_1edge <<"\n";
	str_stream << "total_2_indep = " << total_2_indep <<"\n";
	// connected k=3 motifs
	str_stream << "total_3_tris = " << total_3_tris <<"\n";
	str_stream << "total_2_star = " << total_2_star <<"\n";
	// disconnected k=3 motifs
	str_stream << "total_3_1edge = " << total_3_1edge <<"\n";
	str_stream << "total_3_indep = " << total_3_indep <<"\n";
	// connected k=4 motifs
	str_stream << "total_4_clique = " << total_4_clique <<"\n";
	str_stream << "total_4_chordcycle = " << total_4_chordcycle <<"\n";
	str_stream << "total_4_tailed_tris = " << total_4_tailed_tris <<"\n";
	str_stream << "total_4_cycle = " << total_4_cycle <<"\n";
	str_stream << "total_3_star = " << total_3_star <<"\n";
	str_stream << "total_4_path = " << total_4_path <<"\n";
	// disconnected k=4 motifs
	str_stream << "total_4_1edge = " << total_4_1edge <<"\n";
	str_stream << "total_4_2edge = " << total_4_2edge <<"\n";
	str_stream << "total_4_2star = " << total_4_2star <<"\n";
	str_stream << "total_4_tri = " << total_4_tri <<"\n";
	str_stream << "total_4_indep = " << total_4_indep <<"\n";
	return str_stream.str();
}

string graphlet_core::get_graphlet_frequency_distribution() {
	return "";
}
