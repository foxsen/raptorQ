#pragma once
#include <cmath>
#include <iostream>
#include "string.h"
#include <cmath>
#include "Encoder.h"
#include "Decoder.h"
#include "Generators.h"
#include <stdio.h>

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))

struct PartitionS {
	int IL;
	int IS;
	int JL;
	int JS;
};

class Helper {
public:
	Helper();
	virtual ~Helper();

private:
	int Kmax;
	int P;
	int F;
	int W;
	int Al;
	int Kmin;
	int Gmax;
	int G;
	int T;
	int Kt;
	int Z;
	int N;
	PartitionS KtZ;
	PartitionS TAlN;

public:
	const int getKmax();
	const int getP();
	const int getF();
	const int getW();
	const int getAl();
	const int getKmin();
	const int getGmax();
	const int getG();
	const int getT();
	const int getKt();
	const int getZ();
	const int getN();
	const PartitionS getKtZ();
	const PartitionS getTAlN();

public:
	bool init(int Al, int Kmin , int Gmax, int P, int F, int W);

public:
	const PartitionS partition(int I, int J);

public:
	void toString();

};

const static long Combin(int m, int n) {
	if (n==1||n==0) 
		return m;
	if (n > m-n) {
		return Combin(m, m-n);
	}	
	return (long)((double)Combin(m-1,n-1)*m/n);
}