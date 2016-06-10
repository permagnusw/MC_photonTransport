#include "Coord.h"
#include <cmath>
using namespace photon_IS;

Coord::Coord(){
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

Coord::Coord(double x0, double y0, double z0){
	x = x0;
	y = y0;
	z = z0;
}

double Coord::square(){
	return (x*x + y*y + z*z);
}

double Coord::norm(){
	return std::sqrt(this->square());
}

/*
Coord::Coord(const Coord &coord){
	x = coord.x;
	y = coord.y;
	z = coord.z;
}



void Coord::copy(const Coord &coord){
	x = coord.x;
	y = coord.y;
	z = coord.z;
}
*/
