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

#ifndef GRAPHLET_VERTEX_H_
#define GRAPHLET_VERTEX_H_

using namespace std;

namespace graphlet {
    class Vertex {
        private:
            int id, b;
        public:
            Vertex(int vertex_id, int bound): id(vertex_id), b(bound) {};

            void set_id(int vid)        { id = vid; }
            int get_id()                { return id; }

            void set_bound(int value)   { b = value; }
            int get_bound()             { return b; }
    };

    /**
     * @brief Order from largest to smallest value
     * Note that ties are broken by vertex id
     *
     * @param v is a vertex object containing a vertex id and its value
     * @param u is a vertex object containing a vertex id and its value
     * @return Returns true if the value of v is larger than the value of u
     *
     */
    static bool decr_bound(Vertex v,  Vertex u) {
        return (v.get_bound() > u.get_bound() ||
                (v.get_bound() == u.get_bound() && v.get_id() > u.get_id()));
    }
    /**
     * @brief Order from smallest to largest value
     * Note that ties are broken by vertex id
     *
     * @param v is a vertex object containing a vertex id and its value
     * @param u is a vertex object containing a vertex id and its value
     * @return Returns true if the value of v is smaller than the value of u
     */
    static bool incr_bound(Vertex v,  Vertex u) {
        return (v.get_bound() < u.get_bound() ||
                (v.get_bound() == u.get_bound() && v.get_id() > u.get_id()));
    };
};
#endif
