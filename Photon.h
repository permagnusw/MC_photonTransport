#ifndef PHOTON_H
#define PHOTON_H

#include "Coord.h"
#include "Dir.h"
#include <cstring> //for size_t

namespace photon_IS{

class Photon{
	public:
		Photon();
		Photon(const Coord &coord, const Dir &dir);

		Coord r;
		Dir u;
		double W;
		double s_tot; //Total distance travelled
		size_t N_w;   //Number of wall collisions
		size_t N_s;   //Number of interior scattering events

		bool isCollected;
		bool isAbsorbedInterior;
		bool isAbsorbedWall;
		bool isAbsorbedEntrance;
		bool isAbsorbedPort;

		double W_interior;



	private:

};

} //End namespace photon_IS

#endif
