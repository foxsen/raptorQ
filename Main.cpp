#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>

#ifdef WINDOWS
#include "StdAfx.h"
#include <windows.h>
#include <WinBase.h>
#endif

#include <ctime>

#include "Helper.h"
#include "Symbol.h"

using namespace std;

#define TEST

time_t GetTickCount()
{
	struct timeval t;
	
	gettimeofday(&t, NULL);

	return t.tv_sec * 1000 + t.tv_usec/1000;
}

int main(int argc, char* argv[]) {

#ifdef TEST
	/*
	Helper helper;	

	int al = 4;
	int kmin = 128;
	int gmax = 16;
	int packagesize = 1024;
	int w = 4096000;
	int overhead = 20;
	double loss = 0.1;
	int size = 256 * 1024;

	bool ret = helper.init(al,kmin,gmax,packagesize,size,w);
	
	if (ret)
		helper.toString();
	else 
		cout<<"error"<<endl;

	PartitionS p = helper.getKtZ();

	int i;
	int t = helper.getT();
    int k = p.IS;
	int z = p.JS;
	*/

#define OVERHEAD  4 // K + OVERHEAD received symbols can have high possibility to recover

	int i,j,K,T;
	double lossrate;
	int overhead;
	
	K = 8;
	T = 4;
	lossrate = 0.1;
	overhead = (int)((K*lossrate + 10) / (1 - lossrate));

	cout << "K=" << K << " T=" << T << " Overhead=" << overhead << " lossrate=" << lossrate << endl;


    char **source=new char*[K];
	for (i=0; i<K; i++)
	{
		source[i] = new char[T];
		for (j=0; j<T; j++) {
			*(char*)(source[i] + j) = j;
		}
	}

	Encoder* encoder;

#if 0
	//test whether we can succeed for all K
	int success=0, fail=0;

	for (K=4;K<30;K++) {

		encoder = new Encoder();
		encoder->init(K,T);
		cout << "initialized with K=" << K << " and T=" << T << endl;

		if (encoder->encode(source, 2)) 
			success ++;
		else
			fail++;
			
		delete encoder;

	}
	cout << "success " << success << " fail " << fail << endl;
	while (true);
#endif

	Symbol **repairs;

	int start, end;

	/* encode */
	start = GetTickCount();
	
	encoder = new Encoder();
	encoder->init(K,T);

	repairs = encoder->encode(source, overhead);

	end = GetTickCount();

	cout << "encode bandwidth=" << K * T  /  ((end - start) * 1000.0) << "MB/s" << endl;

	/* send in packet erasure channel */
	srand((unsigned)time(NULL));

	char **received = new char*[K + overhead];
	for (i=0;i<K + overhead; i++)
		received[i] = new char[T];
	int n = 0; //received count

	int *esi = new int[K + overhead];

	int l = 0; // lost source packet count
	int *lost = new int[K];

	/* send source */
	for (i=0; i< K; i++) {
		if (rand()/(RAND_MAX + 1.0) > lossrate && i!=2) {
			memcpy(received[n], source[i], T);
			esi[n] = i;
			n++;
		} else {
			lost[l++] = i;
		}
	}

	/* send repairs */
	for (i=0; i< overhead; i++) {
		if (rand()/(RAND_MAX + 1.0) > lossrate) {
			memcpy(received[n], repairs[i]->data, T);
			esi[n] = K + i;
			n++;
		}
		delete repairs[i];
	}
	delete[] repairs;

	cout << "Received " << n << " packets out of total " << K + overhead << endl;

    /* decode */
	int no_decode = 0, success_decode = 0;
	int failed_decode1 = 0, failed_decode2 = 0, failed_decode3 = 0;
	Decoder *decoder;

	start = GetTickCount();

	if (l == 0) { 
		cout << "All source arrived!" << endl;
		no_decode ++;
	} else {
		if ( n < K ) {
			cout << "Too few packets got " << n << endl;
			failed_decode1++;
		} else {
			Symbol *s;

			decoder = new Decoder();
			decoder->init(K,T);

			if (decoder->decode(received, n, esi)) {
				for (i=0;i < l; i++) {
					s = decoder->recover(lost[i]);
#if 1
					if (memcmp(s->data, source[lost[i]], T)!=0) {
						cout << "Recoverd symbol is not correct! x=" << lost[i] << endl;
						failed_decode3++;
						//break;
					} else {
						cout << "Recoverd symbol is ok! x=" << lost[i] << endl;
					}
#endif
					delete s;
				}
				cout << "All lost symbol recovered!" << endl;
				success_decode ++;
			}else{
				cout << "Decode failed" << endl;
				failed_decode2++;
			}

			delete decoder;
		}
	}

	end = GetTickCount();

	cout << "decode bandwidth=" << K * T  /  ((end - start) * 1000.0) << "MB/s" << endl;


	for (i=0;i<K;i++)
		delete[] source[i];
	delete[] source;

	for (i=0;i<K + overhead; i++)
		delete[] received[i];
	delete[] received;
 
	/* handled in generator. the primciple is that consumer handle the free of memories */
	//delete[] esi;
	delete[] lost;	

	while(true);
	
#endif


	return 0;
}
