#include "Helper.h"
#include "Symbol.h"

Encoder::Encoder(){
}

Encoder::~Encoder() {
	delete gen;
}

bool Encoder::init(int K,int T) {
	gen = new Generators();
	return gen->gen(K, K, T);
}

Symbol** Encoder::encode(char **source, int overhead) {
	Symbol **s;

	gen->prepare(source,gen->getK(),NULL);
	s = gen->generate_intermediates();
	if (!s) return NULL;
	return gen->generate_repairs(overhead);
}

