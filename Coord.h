#ifndef COORD_H
#define COORD_H

namespace photon_IS{

class Coord{
	public:
		Coord();
		Coord(double x0, double y0, double z0);
		//Coord(Coord &coord);
		//void copy(Coord &coord);
		double x;
		double y;
		double z;
		double square();
		double norm();
		
	private:

};

} //End namespace photon_IS

#endif
