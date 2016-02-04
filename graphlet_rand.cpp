#include "graphlet_rand.h"

cust_rand cust_rand::controller;				/** singleton */
unsigned long cust_rand::state[N] = {0x0UL}; 	/** static variables */
int cust_rand::P = 0;
void cust_rand::gen_state() { 					/** generate new state vector */
	for (int i = 0; i < (N - M); ++i) 		state[i] = state[i + M] ^ twiddle(state[i], state[i + 1]);
	for (int i = N - M; i < (N - 1); ++i) 	state[i] = state[i + M - N] ^ twiddle(state[i], state[i + 1]);
	state[N - 1] = state[M - 1] ^ twiddle(state[N - 1], state[0]);
	P = 0; // reset position
}
/** init by 32 bit seed */
void cust_rand::set_cust_seed(unsigned long s) {
	state[0] = s & 0xFFFFFFFFUL;
	for (int i = 1; i < N; ++i) {
		state[i] = 1812433253UL * (state[i - 1] ^ (state[i - 1] >> 30)) + i;
		state[i] &= 0xFFFFFFFFUL;
	}
	P = N;
}
