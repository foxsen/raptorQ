#pragma once

#include "Symbol.h"

struct TuplS {
	int d;
	int a;
	int b;
	int d1;
	int a1;
	int b1;
};

typedef struct Element {
	unsigned char val;
#ifdef SPARSE
	int rprev;
	int rnext;
	int cprev;
	int cnext;
#endif
} Elem;

struct Degree {
	int ori;
	int curr;
	int gtone;
};

class Generators {
public:
	Generators();
	virtual ~Generators();

public:
	bool gen(int _K, int _N, int _T);
	bool _0_init(int _K, int _N, int _T);
	void _1_Tuples(void);
	void _2_Matrix_GLDPC();
	void _3_Matrix_GHDPC();
	void _4_Matrix_GLT();

	bool prepare(char **source, int _N, int *esi);
	void swap_row(int i1, int i2);
	void swap_col(int j1, int j2);
	void xxor(int i1, int i2, int u);
	int gaussian_elimination(int starti, int startj);
	Symbol **generate_intermediates(void);
	Symbol **generate_repairs(int count);	
	Symbol *recover_symbol(int x);

public:
	int getL();
	int getK();

private:
	/* refer to RFC5053 for the naming */
	int K; //symbol number of a source block
	int I; // index of K1 in lookup table
	int K1; //symbol number of a extended source block
	int J_K; // systematic index of K1
	int T; //symbol size
	int X; //not use
	int S; //LDPC row number
	int H; //Half row number
	int W; //LT symbol number
	int P; //Permanent inactivated symbol number
	int P1; //smallest prime>=P
	int U; // P - H
	int B; // W - S
	int H1; //ceil(H/2)

	int L;  //K+S+H
	int N;  //received symbol number; set to K1 for encoder 
	int N1; //received symbol plus extended one
	int M;  //N+S+H
	TuplS *Tuples; //Tuples list
	int tupl_len;
	Symbol **C1; // size M: S+H zero + N symbols
	Symbol **C;  // L intermediate symbols 
	char **sources; // pointer to original source array
	Symbol **R; // repair symbols 
	Elem **A; // generator matrix
	Elem **Abak; //backup of A
#ifdef SPARSE
	int *rowh; // row list head
	int *colh; // column list head
#endif
	Degree *degree; //number of 1 in row i
	int dgh; // degree list head
	int *isi; //Encoding Symbol ID list
	int status; /* 1: para inited, 2: source filled 3: intermediate generated 4: repair generated */

public:
	const TuplS Tupl(int X);
	const unsigned int RandYim(unsigned int y, unsigned char i,unsigned int m);
	const unsigned int Deg(unsigned int v);
	Symbol* LTEnc(Symbol** C_L, TuplS tupl);
	void verify(void);
	void PrintMatrix();
	void ToString();

};

extern unsigned char octmul(unsigned char u, unsigned char v);
extern unsigned char octdiv(unsigned char u, unsigned char v);

