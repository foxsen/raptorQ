
#include "Helper.h"
#include "Symbol.h"
#include "Tables.h"
#include "stdlib.h"

using namespace std;

//#define DEBUG 1

/* initialize parameters for raptor code 
  arg:
   _K: symbol number in a source block
   _N: received symbol number; for encoder _N = _K, for decoder, let _N = _K during init, later prepare() will fixup
   _T: symbol size in bytes
  return:
    true if initialized successfully, false otherwise
 */
bool Generators::gen(int _K, int _N, int _T) {
	int i,j;

	if (!_0_init(_K, _N, _T))
		return false;

	_1_Tuples();

	_2_Matrix_GLDPC();

	_3_Matrix_GHDPC();

	_4_Matrix_GLT();

	/* save for reuse, A will be changed in calculating intermediates */
	for( i = 0; i < (M); i++)
	{
		for(j = 0; j < (L); j++)
			Abak[i][j] = A[i][j];
	}
	status = 1;

	ToString();

	//while(true);

	return true;
}

/* calculate parameters, allocate data structures */
bool Generators::_0_init(int _K, int _N, int _T) {
	int i,j;

	if (_K < 1 || _K > 56403) {
		cout << "Invalid K, should in [1,56403]\n";
		return false;
	}

	if (_N < _K) {
		cout << "N should not be smaller than K\n";
		return false;
	}

	K = _K;
	T = _T;
	N = _N;

	I=0;
	while (K > lookup_table[I][0])
		I++;

	K1 = lookup_table[I][0];
	J_K = lookup_table[I][1];
	S = lookup_table[I][2];
	H = lookup_table[I][3];
	W = lookup_table[I][4];

	N1 = K1 - K + N;

	L = K1 + S + H;
	M = N1 + S + H;

	P = L - W;
	U = P - H;
	B = W - S;
	
	H1 = (int)ceil(H/2.0);

	//P1 be the smallest prime that is greater than or equal to P
	int _P1_len = P;
	int flag = false;
	while(!flag)
	{		
		for(i = 2; i <= sqrt((double)_P1_len); i++)
		{
			if(_P1_len % i == 0)
			{
				_P1_len++;
				flag = false;
				break;
			}
			else
				flag = true;
		}
	}

	P1 = _P1_len;

	try {

		C1 = new Symbol*[M];
		for (i = 0; i < M; i++)
			C1[i] = new Symbol(T);

		C = new Symbol*[L];
		for (i = 0; i < L; i++) {
			C[i] = new Symbol(T);
			C[i]->esi = i;
		}

		R = NULL;

		A = new Elem*[M];
		for(i = 0; i < (M); i++)
			A[i] = new Elem[L];

		Abak = new Elem*[M];
		for(i = 0; i < (M); i++)
			Abak[i] = new Elem[L];


		for( i = 0; i < (M); i++)
		{
			for(j = 0; j < (L); j++)
				A[i][j].val = 0;
				//memset(&A[i][j], 0, sizeof(Elem));
		}

	

		/* isi[i]=i for encoder;
		   for decoder, isi[i] == the ESI for the ith received symbol for i in[0,N), ==i-N+K for [N,N1) */
		
		isi = new int[N1];
		for (i=0;i<N1;i++)
			isi[i] = i;

		degree = new Degree[M];
	
	} catch ( bad_alloc &e) {
		cout << "Allocation failed!" << e.what() << endl;
		abort();
	}

	return true;
}


/* calculate Tuples. 
   to avoid repeated calculation, we generate most possible ones during init
*/
#define MAX_OVERHEAD (40)
void Generators::_1_Tuples(void) {
	int i;

	try {	
		this->Tuples = new TuplS[M + MAX_OVERHEAD];
	} catch ( bad_alloc &e) {
		cout << "Allocate Tuples[] failed!" << e.what() << endl;
		abort();
	}

	for (i=0; i < M + MAX_OVERHEAD; i++) {
		Tuples[i] = Tupl(i);
	}

	tupl_len = M + MAX_OVERHEAD;

}

/* Generate the LDPC rows of matrix A */ 
void Generators::_2_Matrix_GLDPC() {
	int i;
	int a,b;

	/* G_LDPC,1 */
	for(i = 0; i < B; i++)
	{
		a = 1 + i / S;
		b = i % S;
		A[b][i].val = 1;

		b = (b + a) % S;
		A[b][i].val = 1;

    	b = (b + a) % S;
		A[b][i].val = 1;
	}

	/* identity part */
	for(i = 0; i < S; i++)
		A[i][B + i].val = 1;

	/* G_LDPC,2 */
	for (i = 0; i < S; i++) {
		a = i % P;
		b = (i + 1) % P;
		A[i][W + a].val = 1;
		A[i][W + b].val = 1;
	}
  
	return ;
}

/* Generate G_HDPC */
void Generators::_3_Matrix_GHDPC()
{
	int i,j,k;

	for (j=0; j< K1 + S - 1; j++) {
		i = RandYim(j+1, 6, H);
		A[S+i][j].val = 1;
		i = (RandYim(j+1, 6, H) + RandYim(j+1, 7, H - 1) + 1) % H;
		A[S+i][j].val = 1;
	}

	for (i=S; i < S + H; i++) {
		A[i][K1 + S - 1].val = OCT_EXP[i-S];
	}

	for (i=S; i< S + H; i++) {
		for (j=0; j < K1 + S; j++) {
			unsigned char tmp = 0;
			for (k=j; k < K1 + S; k++) 
				if (A[i][k].val) tmp ^= octmul(A[i][k].val,OCT_EXP[k - j]); 
			A[i][j].val = tmp;
		}
	}

	/* identity part */
	for(i = S; i < (S + H ) ; i++)
		A[i][K1  + i].val = 1;	
	
	return;
}

/* LT rows of maxtrix A */
void Generators::_4_Matrix_GLT()
{
	int i,j;
	int a,b,d;

	TuplS tupl;

	int flag1 = 0;	

	i = 0;
	while(i < N1)
	{ 
		if ( isi[i] < tupl_len) 
			tupl = Tuples[isi[i]];
		else
			tupl = Tupl(isi[i]);

		a = tupl.a;
		b = tupl.b;
		d = tupl.d;

		A[S + H + i][b].val = 1;

		for(j = 1; j < d; j++)
		{
			b = (b + a) % W;
			A[S + H + i][b].val = 1;
		}

		a = tupl.a1;
		b = tupl.b1;
		d = tupl.d1;

		while (b >= P)
			b = (b + a) % P1;
		A[S + H + i][W + b].val = 1;

		for(j = 1; j < d; j++)
		{
			b = (b + a) % P1;
			while(b >= P)
				b = (b + a) % P1;
			A[S + H + i][W + b].val = 1;
		}

		i++;
	}

	return;
}


/* LT encoder
arg:
   _K is K1
   C_L is the intermediate symbols
   tupl is the Tuple for giving ESI x
return:
   pointer to encoded symbol
   caller need to handle the delete of symbol space
 */
Symbol* Generators::LTEnc(Symbol** C_L ,TuplS tupl) {
	int a = tupl.a;
	int b = tupl.b;
	int d = tupl.d;
	Symbol *s;

	s = new Symbol(T);

	*s = *C_L[b];
	for (int j=1; j < d; j++) {
		b=(b + a) % W;
		s->xxor(C_L[b]);
	}

	a = tupl.a1;
	b = tupl.b1;
	d = tupl.d1;

	while (b >= P)
			b = (b + a) % P1;
	s->xxor(C_L[W + b]);

	for(int j = 1; j < d; j++)
	{
		b = (b + a) % P1;
		while(b >= P)
			b = (b + a) % P1;
		s->xxor(C_L[W + b]);
	}

	return s;
}

/* Tuple generator */
const TuplS Generators::Tupl(int X) {
	TuplS tupl;
	
	unsigned int A = (53591 + J_K*997);
	if (A % 2 == 0) A = A + 1;
	unsigned int B = 10267*(J_K + 1);
	unsigned int y = (B + X*A);

	int v = RandYim(y,0,1048576);
	tupl.d = Deg(v);
	tupl.a = 1 + RandYim(y, 1, W - 1);
	tupl.b = RandYim(y, 2, W);

	if (tupl.d < 4) tupl.d1 = 2 + RandYim(X, 3, 2); else tupl.d1 = 2; 
	tupl.a1 = 1 + RandYim(X, 4, P1 - 1);
	tupl.b1 = RandYim(X, 5, P1);

	return tupl;
}


const unsigned int Generators::RandYim(unsigned int y, unsigned char i,unsigned int m) {
	return (V0[((y & 0xff) +i) & 0xff] ^ V1[(((y>>8) & 0xff) + i) & 0xff] ^ V2[(((y>>16) & 0xff) + i) & 0xff] ^ V3[(((y>>24) & 0xff) + i) & 0xff]) % m;
}


/* find j, f[j-1] <= v <f[j]
   d = min(j, W - 2)
   v must < 2^^20
 */
const unsigned int Generators::Deg(unsigned int v) {
	int j=0;
	while( v > f[j])
		j++;
	return min(j, W - 2);
}

/* clean up leftovers of last run, fill in new source block; 
   called by the encoder/decoder, it can be called many times
        encoder->init
		while(get more source) {
		   prepare
		   decode
	    }
arg:
   source: source symbol list with Length _N
   _esi: the ESI list for the source list
return:
   success or not
 */
bool Generators::prepare(char **source, int _N, int *_esi) {

	int i,j;

	if (_N < K) {
		cout << "Invalid N in prepare!" << N << endl;
		return false;
	}

	if (status != 1) {

		/* Not the first run, recover A;
		 */
		for( i = 0; i < (M); i++)
		{
			for(j = 0; j < (L); j++)
				A[i][j] = Abak[i][j];
		}
	}

	status = 2;

	try {
		/* only decoder will provide esi(for each source block) */
		if (_esi) {
			int _N1 = _N + K1 - K;

			/* maxtrix LT parts changed, the data struct is not efficent now */
			for (i=0; i < M ;i++)
				delete []A[i];
			delete []A;

			A = new Elem*[_N1 + S + H];
			for(i = 0; i < (_N1 + S + H); i++)
				A[i] = new Elem[L];

			for(i = 0; i < (_N1 + S + H); i++)
			{
				for(j = 0; j < (L); j++)
					if (i < S + H)
						A[i][j] = Abak[i][j];
					else
						A[i][j].val = 0;
			}

			for (i=0; i < M ;i++)
				delete []Abak[i];
			delete []Abak;

			Abak = new Elem*[_N1 + S + H];
			for(i = 0; i < (_N1 + S + H); i++)
				Abak[i] = new Elem[L];

			for(i = 0; i < (_N1 + S + H); i++)
			{
				for(j = 0; j < (L); j++)
					if (i < S + H)
						Abak[i][j] = A[i][j];
					else
						Abak[i][j].val = 0;

			}


			delete[] degree;
			degree = new Degree[_N1 + S + H];
			
			for (i = 0; i < M; i++)
				delete C1[i];
			delete[] C1;

			this->C1 = new Symbol*[_N1 + S + H];
			for (i = 0; i < _N1 + S + H; i++)
				C1[i] = new Symbol(T); 

			M = _N1 + S + H;
			N = _N;
			N1 = _N1;

			delete[] isi;

			isi = new int[N1];

			for (i = 0; i < N; i++) {
				if (_esi[i] < K)
					isi[i] = _esi[i];
				else
					isi[i] = _esi[i] + K1 - K;
			}
			/* K1 - K (N1 - N) padding symbols */
			for (i=N; i < N1; i++)
				isi[i] = i - N + K;

			_4_Matrix_GLT();

			for (i = S + H; i < M; i++)
				for (j=0; j < L; j++)
					Abak[i][j] = A[i][j];

		}
	} catch ( bad_alloc &e) {
		cout << "Allocation in prepare() failed!" << e.what() << endl;
		abort();
	}


	for (i=0; i<L; i++) {
		C[i]->init(T);
		C[i]->esi = i;
	}

	for (i=S + H; i< M; i++)
		if (i < S + H + N)
			C1[i]->fillData(source[i - S - H],T);
		else //padding
			C1[i]->init(T);

	sources = source;

	return true;
}

void Generators::swap_row(int i1, int i2)
{
	if (i1 == i2) return;

	Elem *e;
	e = A[i1];
	A[i1] = A[i2];
	A[i2] = e;

	Degree d = degree[i1];
	degree[i1] = degree[i2];
	degree[i2] = d;

	
	Symbol *s;
	s = C1[i1];
	C1[i1] = C1[i2];
	C1[i2] = s;

}

void Generators::swap_col(int j1, int j2)
{
	if (j1 == j2) return;
	
	for (int i=0; i<M; i++) {
		int t = A[i][j1].val;
		A[i][j1].val = A[i][j2].val;
		A[i][j2].val = t;
	}

	Symbol *s = C[j1];
	C[j1] = C[j2];
	C[j2] = s;	
}

/* A[i1] = A[i1] ^ A[i2]; C1[i1] = C1[i1] ^ C1[i2]; */
void Generators::xxor(int i1, int i2, int U)
{
	if (i1 == i2) return;

	int d = 0;
	for (int j=i2; j<L; j++) {
		A[i1][j].val ^= A[i2][j].val;
		if (j < U && A[i1][j].val == 1) d++;
	}
	degree[i1].curr = d;
	C1[i1]->xxor(C1[i2]);
}

int Generators::gaussian_elimination(int starti, int startj)
{
	int i, k, q, jj, kk;
	int firstone;

	int* HI  = new int[L];
	int* LOW = new int[M];	
		
	for (jj=startj; jj<L; jj++)
	{
		//PrintMatrix();
		k=0;
		for (i=starti; i<=jj-1; i++)
		{
			if (A[i][jj].val)
			{
				HI[k]=i;
				k++;
			}
		}
			
		kk=0;	
		firstone = -1;
		for (i=jj; i<M; i++)
		{
			if(A[i][jj].val)
			{
				LOW[kk]=i;
				if (A[i][jj].val == 1 && firstone == -1) firstone = kk;
				kk++;
			}
		}
		
		if (kk==0){
			cout << " Encoder: due to unclear reasons the process can not continue" << endl;
			delete[] HI;
			delete[] LOW;
			return 0;
		}
		
		
		if (firstone > 0) {
			int t = LOW[0];
			LOW[0] = LOW[firstone];
			LOW[firstone] = t;
		}
		

		if (A[LOW[0]][jj].val != 1) {
			unsigned char v = A[LOW[0]][jj].val;

			C1[LOW[0]]->div(v);
			for (q=jj; q<L; q++)
				A[LOW[0]][q].val = octdiv(A[LOW[0]][q].val, v);
		}

		for (i=1; i<kk; i++)
		{
			unsigned char v = A[LOW[i]][jj].val;

			C1[LOW[i]]->muladd(C1[LOW[0]],v);
			
			for (q=jj; q<L; q++)
				A[LOW[i]][q].val ^= octmul(A[LOW[0]][q].val, v);
		}
		
		for (i=0; i<k; i++)
		{
			unsigned char v = A[HI[i]][jj].val;

			C1[HI[i]]->muladd(C1[LOW[0]], v);
			
			for (q=jj; q<L; q++ )
				A[HI[i]][q].val ^= octmul(A[LOW[0]][q].val, v);
		}
		
		if (LOW[0] != jj) {
			Elem* temp;
			Symbol *tempo;

			temp = A[jj];
			A[jj] = A[LOW[0]];
			A[LOW[0]] = temp;

			tempo = C1[jj];
			C1[jj] = C1[LOW[0]];
			C1[LOW[0]] = tempo;
		}
		//PrintMatrix();
	}

	return 1;
}

// check correctness of intermediates
void Generators::verify(void)
{
	int i,j;
	Symbol *s,*s1;
	char *p;

	s = new Symbol(T);
	s1 = new Symbol(T);

	for (i = 0; i < M; i++) {
		s->init(T);
		for (j = 0; j < L; j++) {
			if (Abak[i][j].val)
				s->muladd(C[j], Abak[i][j].val);
		}

		if (i < S + H || i >= S + H + N)
			p = (char*)s1->data;
		else
			p = (char*)sources[i - S - H];
		
		if (memcmp(s->data, p, T) != 0) {
			printf("Check fail for line %d,%x vs %x\n", i, *(int*)s->data, *(int*)p);
		}
	}
	
	delete s;
	delete s1;
}


/* core decode algorithm here. Calculate C from equation A * C = C1 
   return the intermediate symbol list C, NULL if decode fail
 */
Symbol ** Generators::generate_intermediates(void)
{
	
  /* calculate A^^-1 to get intermediate symbol list C by Gaussian elimination */

	if (status != 2) {
		cout << "Wrong call sequence! Filling the source block before generate intermediates" << endl;
		return NULL;
	}

//#define GAUSSIAN
#define RFC6330

#ifdef GAUSSIAN
	int i;
	if (!gaussian_elimination(0,0)) return NULL;
#endif

#ifdef RFC6330
	int* cols1 = new int[L];
	int* cols2 = new int[L];

	int k1,k2;

	/* init degree list */
	int i,j,d,gtone;
	for (i = 0; i < M; i++) 
	{
		d = 0;
		gtone = 0;
		for (j=0; j < L - P; j++)
			if (A[i][j].val) {
				d++;
				if (A[i][j].val > 1) gtone++;
			}
		degree[i].curr = degree[i].ori = d;
		degree[i].gtone = gtone; 
	}

	PrintMatrix();

	/* step 1 */
	int _I, _U, r;
	int gtone_start = 0;

	_I = 0;
	_U = P;

	while (_I + _U < L) {
		int index, o;

		//select minimal current degree with minimal original degree
		//gtone == 0 first
retry:
		index = M; o = L; r = L;
		for (i = _I; i < M; i ++) {
			if ((gtone_start || (gtone_start==0 && degree[i].gtone==0)) && degree[i].curr > 0 && degree[i].curr <= r) {
				index = i; 
				if (degree[i].curr < r || (degree[i].curr == r && degree[i].ori < o)) {
					o = degree[i].ori;
					r = degree[i].curr;
				}
			}
		}

		if (index == M) {
		    if (gtone_start) goto retry;
			cout << "Cannot find enough rows to decode" << endl;
			PrintMatrix();
			return NULL;
		}

		swap_row(_I, index);

#ifdef DEBUG
		cout << "swap row:" << I << "index" << index << endl;
		PrintMatrix();
#endif

		k1 = k2 = 0;
		for (j = _I; j < L - _U; j++) {
			if (j < L - _U - r + 1) {
				if (A[_I][j].val != 0 ) {
					cols1[k1++] = j;
				}
			} else {
				if (A[_I][j].val == 0)
					cols2[k2++] = j;
			}
		}

		if (k1 != k2 + 1) {
			printf("Assert fail: %d!= %d + 1, _I=%d\n", k1, k2, _I);
			return NULL;
			//PrintMatrix();
		}

		/* put one nonezero to [_I][_I], r-1 1 to [L-_U-r+1,L-_U) */
		swap_col(_I,cols1[0]);
		//cout << "swap col:" << _I << "," << cols1[0] << endl;
		//PrintMatrix();
		for (j=0; j<k2; j++) {
			swap_col(cols2[j], cols1[j+1]);
			//cout << "swap col:" << cols2[j] << "," << cols1[j+1] << endl;
			//PrintMatrix();
		}

		if (A[_I][_I].val > 1) {
			unsigned char v = A[_I][_I].val;

			C1[_I]->div(v);
			for (j=_I; j<L; j++)
				A[_I][j].val = octdiv(A[_I][j].val, v);
		}


		for (i=_I+1;i<M;i++) {
			unsigned char v = A[i][_I].val;	
			if (v) {					
				A[i][_I].val = 0; 
				degree[i].curr--; 
				if (v > 1) degree[i].gtone--;
				for (j = L - _U - (r - 1); j < L; j++) {
					int oldv = A[i][j].val;
					A[i][j].val ^= octmul(v, A[_I][j].val);
					if (j < L - _U) {
						if (A[i][j].val > 0) {
							degree[i].curr ++;
							if (A[i][j].val > 1) degree[i].gtone ++;
						}
						if (oldv > 0) { 
							degree[i].curr --;
							if (oldv > 1) degree[i].gtone --;
						}
					}
				}
				C1[i]->muladd(C1[_I], v);
			} 

			/* we calculate only the intersect part of A&V 
			   Note: lines with A[i][_I] == 0 also needs this 
			 */
			for (j=L - _U - (r - 1); j < L - _U; j++)
				if (A[i][j].val) {
					degree[i].curr--;
					if (A[i][j].val > 1) degree[i].gtone--;
				}
		}

		_I++;
		_U += r - 1;
	}
	delete[] cols1;
	delete[] cols2;

	cout << "_I=" << _I << "_U=" << _U << endl;
	//PrintMatrix();

/* step 2 */
//gaussian elimination on the (M - _I) x _U matrix
    if (!gaussian_elimination(_I, _I)) return NULL;
	//cout << "After Step 2" << endl;
	//PrintMatrix();

	/* step 3 */
	for (int jj=_I; jj<L; jj++)
		for (int i=0; i< _I; i++) {
			unsigned char v = A[i][jj].val;
			if (v) {
				A[i][jj].val = 0;
				C1[i]->muladd(C1[jj],v);
			}
		}
		
#ifdef DEBUG
	cout << "After Step 3" << endl;
	PrintMatrix();
#endif

#endif 

	/* result now in C1, copy to C; C1 is useless from now on */
	for (i=0; i<L; i++) {
		//borrow C1->esi to remember the real pos of C[i]
		C1[i]->esi = C[i]->esi;
	}

	for (i=0; i<L; i++) {
		//copy C1[i] to its real pos in C
		//to avoid data copy, just swap the pointers
		Symbol *s;
		s = C[C1[i]->esi];
		C[C1[i]->esi] = C1[i];
		C1[i] = s; 
	}

#if 0
	verify();
#endif

	status = 3;

	return C;
}

/*Generate repair symbols using LT encoder,must be called after C is calculated
  The ESI of generated symbols starts from K
arg:
  count: number of repair symbols to generate
Return:
  the list of request count of repair symbols
 */
Symbol **Generators::generate_repairs(int count)
{
	int i, isi;
	
	if (status != 3) {
		cout << "Wrong call sequence! Generate intermediates before generate repairs" << endl;
		return NULL;
	}

	/* caller needs to free R */
	R = new Symbol*[count];

	/*
	for (i=0; i < count; i++) 
	{
		R[i].init(T);
	}
	*/

	for (i = K; i < K + count; i++) 
	{
		TuplS tupl;
		isi = i + K1 - K;
		if (isi < tupl_len) 
			tupl = Tuples[isi];
		else 
			tupl = Tupl(isi);
		R[i - K] = LTEnc(C, tupl);
	}

	status = 4;

	return R;
}

Symbol *Generators::recover_symbol(int x)
{
	Symbol *s;

	if (x >= K) {
		printf("try to recover non-source symbols!\n");
		return NULL;
	}

	TuplS tupl;
	tupl = Tuples[x];

	s = LTEnc(C, tupl);

	return s;
}

int Generators::getL() {
	return this->L;
}

int Generators::getK() {
	return this->K;
}

void Generators::PrintMatrix(void)
{
#ifdef DEBUG
	int i;
	cout<<"Press a number  to print the generation matrix:" << endl;
	cin >> i;
	for (i=0; i < M ;i++) {
		for (int j=0; j < L;j++) 
			printf("%2x ", A[i][j].val);
		cout<<endl;
	}
#endif
}

void Generators::ToString() {
	cout << "K=" << K << endl; 
	cout << "T=" << T << endl; 
	cout << "H=" << H << endl; 
	cout << "S=" << S << endl; 
	cout << "L=" << L << endl; 
	cout << "N=" << N << endl; 
	cout << "M=" << M << endl; 
	cout << "K1=" << K1 << endl; 
	cout << "W=" << W << endl; 
	cout << "P=" << P << endl; 
	cout << "B=" << B << endl; 
	cout << "U=" << U << endl; 
	cout << "P1=" << P1 << endl; 


	cout <<"Tuples:" << endl;
	for (int i=0;i < L; i++) {
		cout << "Tuple " << i << " d,a,b=" << Tuples[i].d << "," << Tuples[i].a << "," << Tuples[i].b << ",";
		cout << " d1,a1,b1=" << Tuples[i].d1 << "," << Tuples[i].a1 << "," << Tuples[i].b1 << endl;

	}

	cout<<"The generation maxtrix:" << endl;
	for (int i=0; i < M ;i++) {
		for (int j=0; j < L;j++) 
			printf("%2x ", A[i][j].val);
			//cout<<A[i][j].val<<" ";
		cout<<endl;
	}
}

Generators::Generators() {
	status = 0;
}

Generators::~Generators() {
	int i;
	for (i=0; i< M; i++)
		delete C1[i];
	delete[] C1;

	for (i=0; i < M ;i++)
		delete []A[i];
	delete []A;
	for (i=0; i < M ;i++)
		delete []Abak[i];
	delete []Abak;
	delete[] Tuples;
	for (i=0; i < L; i++)
		delete C[i];
	delete[] C;
	// caller handle the delete of R
	//delete[] R;
	delete[] isi;
	delete[] degree;
}

/* section 5.7.2 */
unsigned char octmul(unsigned char u, unsigned char v)
{
	if (u == 0 || v == 0) return 0;
	if (v == 1) return u;
	if (u == 1) return v;
	return OCT_EXP[OCT_LOG[u] + OCT_LOG[v]];
}

unsigned char octdiv(unsigned char u, unsigned char v)
{
	if (u == 0) return 0;
	if (v == 1) return u;

	return OCT_EXP[OCT_LOG[u] - OCT_LOG[v] + 255];
}
