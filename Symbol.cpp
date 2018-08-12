#include "StdAfx.h"
#include <iostream>
#include <fstream>

#include "Symbol.h"
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
		cout << "size must be a multiple of %d bytes" << endl;
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
			cout << "size must be a multiple of %d bytes" << endl;
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
	
	for (i=0; i< nbytes; i++) 
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

	for (i=0;i<nbytes; i++)
		p[i] = octmul(p[i],u);
}

void Symbol::div(unsigned char u)
{
	int i;
	unsigned char *p = (unsigned char*)data;

	for (i=0;i<nbytes; i++)
		p[i] = octdiv(p[i],u);
}

void Symbol::muladd(Symbol *s, unsigned char u)
{
	int i;
	unsigned char *p, *p1;

	if (nbytes != s->nbytes) 
		cout << "Error! try to muladd symbols with unmatched size" << endl;

	p = (unsigned char*)data;
	p1 = (unsigned char*)s->data;

	for (i=0;i<nbytes; i++) {
		if (u==1) p[i] ^= p1[i];
		else p[i] ^= octmul(p1[i],u);
	}
}
