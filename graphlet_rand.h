#ifndef GRAPHLET_RAND_H_
#define GRAPHLET_RAND_H_

class cust_rand {
public:
	static cust_rand controller;
private:
	cust_rand() { set_cust_seed(5489UL); }
	virtual ~cust_rand() {}
public:
	void set_cust_seed(unsigned long);
	unsigned long operator()(){ return rand_int32(); }
	double next_closed() { return static_cast<double>(rand_int32()) * (1. / 4294967295.); } 		// divided by 2^32 - 1
	double next_opened() { return (static_cast<double>(rand_int32()) + .5) * (1. / 4294967296.); } 	// divided by 2^32
	double next_rand() { return (static_cast<unsigned long>(rand_int32())); } 	/// divided by 2^32
protected:
	inline unsigned long rand_int32() { /// generate 32 bit random int
		if (P == N) gen_state();
		unsigned long x = state[P++];
		x ^= (x >> 11);
		x ^= (x << 7) & 0x9D2C5680UL;
		x ^= (x << 15) & 0xEFC60000UL;
		return x ^ (x >> 18);
	}
private:
	static const int N = 624;
	static const int M = 397;
	static unsigned long state[N]; 	/// state vector array
	static int P; 					/// position in state array
	unsigned long twiddle(unsigned long u, unsigned long v) {
		return (((u & 0x80000000UL) | (v & 0x7FFFFFFFUL)) >> 1)^ ((v & 1UL) ? 0x9908B0DFUL : 0x0UL);
	}
	void gen_state();
};
inline void set_custom_seed(unsigned long s) 	{ cust_rand::controller.set_cust_seed(s); }
inline double get_cust_rand() 					{ return cust_rand::controller.next_closed(); }
inline unsigned long get_cust_rand_int() 					{ return cust_rand::controller.next_rand(); }

#endif /* GRAPHLET_RAND_H_ */
