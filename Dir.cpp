#include "Dir.h"
#include <cmath>

using namespace photon_IS;

Dir::Dir(){
	u_x = 0.0;
	u_y = 0.0;
	u_z = 1.0;
}

Dir::Dir(double u_x0, double u_y0, double u_z0){
	u_x = u_x0;
	u_y = u_y0;
	u_z = u_z0;
}

double Dir::square(){
	return (u_x*u_x + u_y*u_y + u_z*u_z);
}

double Dir::norm(){
	return std::sqrt(this->square());
}

/*
void Dir::copy(const Dir & dir){
	u_x = dir.u_x;
	u_y = dir.u_y;
	u_z = dir.u_z;
}

*/
