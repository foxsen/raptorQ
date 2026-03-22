#include "StdAfx.h"
#include <iostream>
#include <fstream>

#include "Symbol.h"
#include "Tables.h"
#include "Generators.h"

using namespace std;

Symbol::Symbol(void)
{
	nbytes = 0;
	data = NULL;
	sbn = -1;
	esi = -1;
}

Symbol::Symbol(unsigned int size)
{
	if (size % sizeof(int) != 0) {
		cout << "size must be a multiple of " << sizeof(int) << " bytes" << endl;
	}

	data = new int[size/sizeof(int)];
	memset(data, 0, size);
	nbytes = size;
}

Symbol::~Symbol(void)
{
	if (data) delete[] data;
	nbytes = 0;
}

void Symbol::init(int size)
{
	if (nbytes != size) {
		if (data) delete[] data;
		if (size % sizeof(int) != 0) {
			cout << "size must be a multiple of " << sizeof(int) << " bytes" << endl;
		}
		nbytes = size;
		data = new int[size/sizeof(int)];
	}
	memset(data, 0, size);
}

void Symbol::fillData(char *src, int size)
{
	if (nbytes != size) {
		if (data) delete[] data;
		nbytes = size;
		data = new int[size/sizeof(int)];
	}
	memcpy(data, src, size);
}

void Symbol::print(void)
{
	int i;
	
	for (i=0; i< nbytes/sizeof(int); i++)
	{
		cout << data[i] << " ";
	}
	cout << endl;
	
}

Symbol& Symbol::operator=(const Symbol &s) 
{
	if (nbytes != s.nbytes) {
	  if (data) delete[] data;
	  nbytes = s.nbytes;
	  data = new int[nbytes/sizeof(int)];
	}
	memcpy(data, s.data, s.nbytes);
	
	return *this;
}

Symbol& Symbol::operator^(const Symbol &s) 
{
	int i;

    if (nbytes != s.nbytes) 
		cout << "Error! try to xor symbols with unmatched size" << endl;

	for (i=0;i<nbytes/4; i++)
		data[i] = data[i] ^ s.data[i];

	return *this;
}

void Symbol::xxor(Symbol *s)
{
	int i;

    if (nbytes != s->nbytes) 
		cout << "Error! try to xor symbols with unmatched size" << endl;

	for (i=0;i<nbytes/4; i++)
		data[i] = data[i] ^ s->data[i];

}

void Symbol::mul(unsigned char u)
{
	int i;
	unsigned char *p = (unsigned char*)data;
	const unsigned char *mul_row = OCT_MUL[u];

	for (i=0;i<nbytes; i++)
		p[i] = mul_row[p[i]];
}

void Symbol::div(unsigned char u)
{
	int i;
	unsigned char *p = (unsigned char*)data;

	if (u == 1) return;
	/* compute multiplicative inverse: u^(-1) = exp(255 - log(u)) */
	unsigned char inv = OCT_EXP[255 - OCT_LOG[u]];
	const unsigned char *mul_row = OCT_MUL[inv];

	for (i=0;i<nbytes; i++)
		p[i] = mul_row[p[i]];
}

void Symbol::muladd(Symbol *s, unsigned char u)
{
	int i;

	if (nbytes != s->nbytes)
		cout << "Error! try to muladd symbols with unmatched size" << endl;

	if (u == 1) {
		/* XOR path: process 4 bytes at a time */
		for (i = 0; i < nbytes/4; i++)
			data[i] ^= s->data[i];
	} else {
		unsigned char *p = (unsigned char*)data;
		unsigned char *p1 = (unsigned char*)s->data;
		const unsigned char *mul_row = OCT_MUL[u];
		for (i = 0; i < nbytes; i++)
			p[i] ^= mul_row[p1[i]];
	}
}
