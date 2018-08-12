#pragma once
#include "Generators.h"
#include "Symbol.h"

class Encoder {
public:
	Encoder();
	virtual ~Encoder();

private:
	Generators* gen;

public:
	bool init(int K, int T);
	Symbol** encode(char **source, int overhead);

};