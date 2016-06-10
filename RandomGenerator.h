#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace photon_IS{

class Generator{
	public:
		Generator();
		Generator(size_t s);
		~Generator();
		void set_seed(size_t seed);
		size_t seed;
		size_t count;
		const gsl_rng_type * T;
		gsl_rng * r;
		double uniform();
		double gaussian();

};

}

#endif
