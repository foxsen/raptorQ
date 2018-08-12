#pragma once
#include "Generators.h"

class Decoder {
public:
	Decoder();
	virtual ~Decoder();

private:
	Generators* gen;

public:
	bool init(int K, int T);
	Symbol** decode(char **source, int _N, int *esi);
	Symbol* recover(int x);
};