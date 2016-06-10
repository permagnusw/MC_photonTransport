#include "RandomGenerator.h"

using namespace photon_IS;



Generator::Generator(){
	seed = 0;
	count = 0;
	T = gsl_rng_mt19937; 
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);
}

Generator::Generator(size_t s){
	seed = s;
	count = 0;
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);
}

Generator::~Generator(){
	gsl_rng_free(r);
}

void Generator::set_seed(size_t s){
	gsl_rng_set(r, s);
	seed = s;
}

double Generator::uniform(){
	double ret = gsl_rng_uniform(r);
	count++;
	return ret;
}

double Generator::gaussian(){
	double ret = gsl_ran_ugaussian(r);
	count++;
	return ret;
}
