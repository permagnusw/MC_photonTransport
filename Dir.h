#ifndef DIR_H
#define DIR_H

namespace photon_IS{

class Dir{
	public:
		Dir();
		Dir(double u_x0, double u_y0, double u_z0);
		double u_x;
		double u_y;
		double u_z;
		double square();
		double norm();

	private:

};


} //End photon_IS namespace

#endif
