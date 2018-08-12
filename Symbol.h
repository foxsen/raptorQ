#pragma once 

class Symbol {

public:
	int *data;
	int nbytes;
	int sbn; /* source block number */
	int esi; /* encoding symbol id */

	Symbol(void);
	Symbol(unsigned int size); 
	~Symbol(void);

	void init(int size);
	void print(void);
	void fillData(char *src, int size);
	
	Symbol& operator=(const Symbol &s);
	Symbol& operator^(const Symbol &s);
   
	void xxor(Symbol *s);
	void mul(unsigned char u);
	void div(unsigned char u);
	void muladd(Symbol *s, unsigned char u);

	/*
	int getSBN();
	int getESI();

	void setSBN(int sbn);
	void setESI(int esi);
	*/
};
