#pragma once
#include "Helper.h"
#include "Symbol.h"

Decoder::Decoder() {
}

Decoder::~Decoder() {
	delete gen;
}

/* initialize general parameters as much as possible
   this is meant to be called only once.
   By default, N is set to K. It will be fixup in prepare()
 */
bool Decoder::init(int K, int T) {
	gen = new Generators();
	return gen->gen(K, K, T);
}

Symbol** Decoder::decode(char **source, int _N, int *esi)
{
	Symbol **s;

	gen->prepare(source,_N,esi);
	s = gen->generate_intermediates();
	if (!s) return NULL;
	return s;
}

Symbol* Decoder::recover(int x)
{
	return gen->recover_symbol(x);
}
