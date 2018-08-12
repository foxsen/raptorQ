#include "Helper.h"

Helper::Helper() {
}

Helper::~Helper() {
}

bool Helper::init(int Al, int Kmin, int Gmax, int P, int F, int W) {
	this->Al = Al;
	this->P = P;
	if (P % Al != 0) {
		return false;
	}
	this->Kmin = Kmin;
	this->Kmax = 8192;
	this->Gmax = Gmax;
	this->F = F;
	this->W = W;

	this->G = min((min((int)ceil((double)P * Kmin /F), (int)P/Al)), Gmax);
	this->T = (int)floor((double)P/(Al*G)) * Al;
	this->Kt = (int)ceil((double)F/T);
	this->Z = (int)ceil((double)Kt/Kmax);
	this->N = min((int)ceil(ceil((double)Kt/Z)*T/W),(int)T/Al);

	if (ceil(ceil((double)F/T)/Z) > Kmax) {
		return false;
	}

	this->KtZ = partition(Kt, Z);
	this->TAlN = partition(T/Al, N);

	return true;
};

const PartitionS Helper::partition(int I, int J) {
	PartitionS result;
	result.IL = (int)ceil((double)I/J);
	result.IS = (int)floor((double)I/J);
	result.JL = I - result.IS * J;
	result.JS = J - result.JL;

	return result;
}


const int Helper::getKmax() {
	return Kmax;
}

const int Helper::getP() {
	return P;
}

const int Helper::getF() {
	return F;
}

const int Helper::getW() {
	return W;
}

const int Helper::getAl() {
	return Al;
}

const int Helper::getKmin() {
	return Kmin;
}

const int Helper::getGmax() {
	return Gmax;
}

const int Helper::getG() {
	return G;
}

const int Helper::getT() {
	return T;
}

const int Helper::getKt() {
	return Kt;
}

const int Helper::getZ() {
	return Z;
}

const int Helper::getN() {
	return N;
}

const PartitionS Helper::getKtZ() {
	return KtZ;
}

const PartitionS Helper::getTAlN() {
	return TAlN;
}

void Helper::toString() {
	std::cout<<"Kmax="<<Kmax<<std::endl
		<<"P="<<P<<std::endl
		<<"F="<<F<<std::endl
		<<"W="<<W<<std::endl
		<<"Al="<<Al<<std::endl
		<<"Kmin="<<Kmin<<std::endl
		<<"Gmax="<<Gmax<<std::endl
		<<"G="<<G<<std::endl
		<<"T="<<T<<std::endl
		<<"Kt="<<Kt<<std::endl
		<<"Z="<<Z<<std::endl
		<<"N="<<N<<std::endl;	

	std::cout<<"KL,KS,ZL,ZS"<<std::endl<<KtZ.IL<<","<<KtZ.IS<<","<<KtZ.JL<<","<<KtZ.JS<<std::endl;
	std::cout<<"TL,TS,NL,NS"<<std::endl<<TAlN.IL<<","<<TAlN.IS<<","<<TAlN.JL<<","<<TAlN.JS<<std::endl;
}
