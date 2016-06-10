#include "Photon.h"
using namespace photon_IS;

Photon::Photon(): r(), u(), W(1.0), s_tot(0.0), N_w(0), N_s(0), W_interior(0.0),
		isCollected(false), isAbsorbedInterior(false),
		isAbsorbedWall(false), isAbsorbedEntrance(false), isAbsorbedPort(false){

}


Photon::Photon(const Coord &coord, const Dir &dir): r(coord), u(dir), W(1.0), s_tot(0.0), 
		N_w(0), N_s(0), W_interior(0.0),
		isCollected(false), isAbsorbedInterior(false),
		isAbsorbedWall(false), isAbsorbedEntrance(false), isAbsorbedPort(false){

}


