#include "functions.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <armadillo>
#include <fstream>
#include <iomanip>

using namespace photon_IS;
using namespace std;

void diffuseLocal(Dir &u, Generator &G){
	double chi_1 = G.uniform();
	double chi_2 = G.uniform();
	double c = std::cos(2.0*M_PI*chi_1);
	double chi_2_sqrt = std::sqrt(chi_2);
	u.u_x = std::cos(2.0*M_PI*chi_1)*chi_2_sqrt;
	u.u_y = std::sin(2.0*M_PI*chi_1)*chi_2_sqrt; //make faster with sqrt and sign
	u.u_z = std::sqrt(1.0-chi_2);
}


Dir convertToGlobal(Dir &u, Coord &r, double &R){
	Dir u_global;
	double cos_theta = r.z/R;
	double sin_theta = std::sqrt(1.0-cos_theta*cos_theta);
	if(sin_theta < 0.0001){
		std::cout << "Warning, sin_theta is small" << std::endl;
	}
	double cos_phi = r.x/(R*sin_theta);  //WATCH OUT FOR LIMITS FOR THETA
	double sin_phi = r.y/(R*sin_theta);
	u_global.u_x = -u.u_x*cos_phi*cos_theta - u.u_y*sin_phi - u.u_z*cos_phi*sin_theta;
	u_global.u_y = -u.u_x*sin_phi*cos_theta + u.u_y*cos_phi - u.u_z*sin_phi*sin_theta;
	u_global.u_z =  u.u_x*sin_theta                         - u.u_z*cos_theta;

	//Force normalization? I think necessary. Small error blows up.
	double norm = std::sqrt(u_global.u_x*u_global.u_x + u_global.u_y*u_global.u_y + u_global.u_z*u_global.u_z);
	u_global.u_x = u_global.u_x/norm;
	u_global.u_y = u_global.u_y/norm;
	u_global.u_z = u_global.u_z/norm;

	return u_global;
}

double drawPathLength(double &alpha, Generator &G){
	return (-log(G.uniform())/alpha);
}

double psTheory(double &alpha, double &R){
	return ((1.0/(2.0*alpha*alpha*R*R))*(1.0 - (1.0 + 2.0*alpha*R)*exp(-2.0*alpha*R)));
}

	

void traceSingleDiffuse(Photon &p, double &R, double &rho, size_t &lim, Generator &G){
	double s0 = 0.0;
	//Move from initial position and direction
	s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
	p.s_tot += s0;
	p.N_w ++;

	while(p.N_w < lim){
		p.r.x = p.r.x + s0*p.u.u_x;
		p.r.y = p.r.y + s0*p.u.u_y;
		p.r.z = p.r.z + s0*p.u.u_z;
		if(G.uniform() < rho){
			p.N_w ++;
			diffuseLocal(p.u, G);
			p.u = convertToGlobal(p.u, p.r, R);
			s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
			p.s_tot += s0;
		}
		else{
			break;
		}
	}
}

void traceSingleDiffuseAbsorption(Photon &p, double &R, double &rho, double &alpha, size_t &lim, Generator &G){
	double s0 = 0.0;
	double path = drawPathLength(alpha, G);
	//Move from initial position and direction
	s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
	if(s0 < path){
		p.s_tot += s0;
		p.N_w ++;
		path -= s0;
	}
	else{
		p.s_tot = path;
		return;
	}

	while(p.N_w < lim){
		//Move to wall
		p.r.x = p.r.x + s0*p.u.u_x;
		p.r.y = p.r.y + s0*p.u.u_y;
		p.r.z = p.r.z + s0*p.u.u_z;
		//Reflect?
		if(G.uniform() < rho){
			//Yes
		//	p.N_w ++;
			diffuseLocal(p.u, G);
			p.u = convertToGlobal(p.u, p.r, R);
			s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
			//Can photon propagate distance s0?
			if(s0 <= path){
				p.N_w ++;
				p.s_tot += s0;
				path -= s0;
			}
			else{
				p.s_tot += path;
				break;
			}
		}
		else{
			//No
			break;
		}
	}
}

void traceSingleDiffuseInvariant(Photon &p, double &R, double &rho, size_t &lim, Generator &G){
	double s0 = 0.0;
	//Move from initial position and direction
	s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
	//diffuseLocal(p.u,G);
	//s0 = 2.0*R*p.u.u_z;
	p.s_tot += s0;
	p.N_w ++;

	while(p.N_w < lim){
		if(G.uniform() < rho){
			p.N_w ++;
			diffuseLocal(p.u, G);
			s0 = 2.0*R*p.u.u_z;
			p.s_tot += s0;
		}
		else{
			break;
		}
	}
}

void traceCollectionDiffuse(vector<Photon> &P, double &R, double &rho, size_t &lim, Generator &G){
	for(size_t i=0; i<P.size(); i++){
		traceSingleDiffuse(P[i], R, rho, lim, G);
	}
}



double meanNumberCollision(std::vector<Photon> &P){
	size_t col = 0;
	size_t len = P.size();
	for(size_t i=0; i<len; i++){
		col += P[i].N_w;
	}
	double collisions = (double) col;
	double N = (double) len;
	return (collisions/N);
}

double stddevNumberCollision(std::vector<Photon> &P, double &mean){
	double var = 0.0;
	size_t len = P.size();
	for(size_t i=0; i<len; i++){
		double col = (double) P[i].N_w;
		var += (col-mean)*(col-mean);
	}
	double N = (double) len;
	return sqrt(var/(N-1.0));
}

double meanPathLength(std::vector<Photon> &P){
	//using namespace arma; 
	double s = 0.0;
	//double s_test = 0.0;
	size_t len = P.size();
	//vec path(len);
	for(size_t i=0; i<len; i++){
		s += P[i].s_tot;
		//path(i) = P[i].s_tot;
	}
	/*
	path = sort(path);
	for(size_t i=0; i<len; i++){
		s_test += path(i);
	}
	*/
	double N = (double) len;
	return (s/N);
	//return (s_test/N);
}

double meanPathLengthCollected(std::vector<Photon> &photons){
	double s = 0.0;
	size_t N = photons.size();
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			s += photons[i].s_tot;
		}
	}
	double N_d = (double) N;
	return (s/N);
}

double stddevPathLength(std::vector<Photon> &P, double &mean){
	double var = 0.0;
	size_t len = P.size();
	for(size_t i=0; i<len; i++){
		var += (P[i].s_tot - mean)*(P[i].s_tot - mean);
	}
	double N = (double) len;
	return sqrt(var/(N-1.0));
}

double stddevPathLengthCollected(std::vector<Photon> &P, double &mean){
	double var = 0.0;
	size_t len = P.size();
	for(size_t i=0; i<len; i++){
		if(P[i].isCollected){
			var += (P[i].s_tot - mean)*(P[i].s_tot - mean);
		}
	}
	double N = (double) len;
	return sqrt(var/(N-1.0));
}

void histCollisions(vector<Photon> &photons, string fileName){

	using namespace arma;
	size_t N = photons.size();
	uvec collision_counter(N);
	for(size_t i=0; i<N; i++){
		collision_counter(i) = photons[i].N_w;
	}
	size_t min_N_w = collision_counter.min();
	size_t max_N_w = collision_counter.max();
	size_t range = max_N_w - min_N_w + 1;

	vec centers = linspace<vec>((double) min_N_w, (double) max_N_w, range);
	vec collisions = conv_to< vec >::from(collision_counter);
	uvec hist_int = hist(collisions, centers);
	vec hist_double = conv_to< vec >::from(hist_int);
	hist_double = hist_double/((double) N);

	mat collision_hist(range,2);	
	collision_hist.col(0) = centers;
	collision_hist.col(1) = hist_double;
	collision_hist.save(fileName, raw_ascii);
}

void histPathLength(vector<Photon> &photons, string fileName, double min_s, double max_s, size_t bins){
	using namespace arma;
	size_t N = photons.size();
	vec path(N);
	for(size_t i=0; i<N; i++){
		path(i) = photons[i].N_w;
	}
	double min_s0 = path.min();
	double max_s0 = path.max();

	vec centers = linspace<vec>(min_s0, max_s0, bins);
	vec paths = conv_to< vec >::from(path);
	uvec hist_int = hist(paths, centers);
	vec hist_double = conv_to< vec >::from(hist_int);
	hist_double = hist_double/((double) N);

	mat path_hist(bins,2);	
	path_hist.col(0) = centers;
	path_hist.col(1) = hist_double;
	path_hist.save(fileName, raw_ascii);
}

void savePhotons(std::vector<Photon> &photons, std::string fileName, double R, double rho, double N_photons, size_t seed){

	ofstream myFile;
	myFile.open(fileName.c_str());
	myFile << "#Radius: " << R << "\n" << "#Reflectivity: " << rho << "\n" << "#No. photons: "
	       << N_photons << "\n" << "#Seed: " << seed << "\n";
	myFile << "#N_w s_tot\n";
	for(size_t i=0; i<N_photons;i++){
		myFile << photons[i].N_w << " " << std::setprecision (17) << photons[i].s_tot << "\n";
	}
	myFile.close();
}


void initPhotons(std::vector<Photon> &photons, double R, double z_s, double cos_theta0, double a_s, double b_s, Generator &G){

	size_t N = photons.size();
	for(size_t i=0; i<N; i++){
		//random points inside unit circle
		//double chi_x = G.uniform(); //ERROR HERE, THIS IS ONLY FIRST QUADRANT
		//double chi_y = G.uniform();
		double chi_x = 2.0*G.uniform() - 1.0;
		double chi_y = 2.0*G.uniform() - 1.0;
		while((chi_x*chi_x + chi_y*chi_y) > 1.0){
			//chi_x = G.uniform();
			//chi_y = G.uniform();
			chi_x = 2.0*G.uniform() - 1.0;
			chi_y = 2.0*G.uniform() - 1.0;
		}
		//scaling
		chi_x = a_s*chi_x;
		chi_y = a_s*chi_y;
		//set position
		photons[i].r.z = z_s;
		photons[i].r.x = chi_x - b_s;
		photons[i].r.y = chi_y;
		double chi_1 = G.uniform();
		double chi_2 = G.uniform();
		double z = (1.0 - cos_theta0)*chi_2 + cos_theta0;
		//set direction
		photons[i].u.u_z = -z;
		photons[i].u.u_x = cos(2.0*M_PI*chi_1)*sqrt(1.0 - z*z);
		photons[i].u.u_y = sin(2.0*M_PI*chi_1)*sqrt(1.0 - z*z);
		//cout << "u norm: " << photons[i].u.norm() << endl;
	}
}


void tracePhotonEmptyIS(Photon &p, double &R, double &rho, double &z_s, double &cos_theta0, double &a_p, double &b_p, size_t &lim, Generator &G){
	//Move from initial position and direction
	double dot_prod = p.r.x*p.u.u_x + p.r.y*p.u.u_y + p.r.z*p.u.u_z;
//	cout << "dot_prod: " << dot_prod << endl;
	double s0 = -dot_prod + sqrt(dot_prod*dot_prod + R*R - p.r.square());
//	cout << "s0 first: " << s0 << endl;


	p.s_tot += s0;
	p.N_w ++;

	while(p.N_w < lim){
		//Move to wall
		p.r.x = p.r.x + s0*p.u.u_x;
		p.r.y = p.r.y + s0*p.u.u_y;
		p.r.z = p.r.z + s0*p.u.u_z;
		//cout << "R: " << p.r.norm() << endl;
		//Reflect?
		if(G.uniform() < rho){
			//Yes
		//	p.N_w ++;
			diffuseLocal(p.u, G);
			p.u = convertToGlobal(p.u, p.r, R);
			s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
			//Can photon propagate distance s0?
			if((p.r.z + s0*p.u.u_z) > z_s){
				double s = (z_s - p.r.z)/p.u.u_z;
				if(p.u.u_z < cos_theta0){
					p.isAbsorbedEntrance = true;
					p.s_tot += s;
					break;
				}
				else{
					double x = p.r.x + s*p.u.u_x;
					double y = p.r.y + s*p.u.u_y;
					if(((x-b_p)*(x-b_p) + y*y) < (a_p*a_p)){
						p.isCollected = true;
						//cout << "collected!" << endl;
						p.s_tot += s;
						break;
					} else{ p.isAbsorbedEntrance = true;
						p.s_tot += s;
						break;
					}
				}
			}
			else{
				p.N_w ++;
				p.s_tot += s0;
			}
		}
		else{
			//No
			p.isAbsorbedWall = true;
			break;
		}
	}
}

void getStats(std::vector<Photon> &photons, double &eps_c, double &eps_e, double &eps_w, size_t &N){
	double N_d = (double) N;
	eps_c = 0.0;
	eps_e = 0.0;
	eps_w = 0.0;
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			eps_c += 1.0;
		}else if(photons[i].isAbsorbedEntrance){
			eps_e += 1.0;
		}else{ 
			eps_w += 1.0;
		}
	}
	eps_c = eps_c/N_d;
	eps_e = eps_e/N_d;
	eps_w = eps_w/N_d;
}

void getStatsWeighted(std::vector<Photon> &photons, double &eps_c, double &eps_e, double &eps_w, double &eps_i, size_t &N){
	double N_d = (double) N;
	eps_c = 0.0;
	eps_e = 0.0;
	eps_w = 0.0;
	eps_i = 0.0;
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			eps_c += photons[i].W;
		}else if(photons[i].isAbsorbedEntrance){
			eps_e += photons[i].W;
		}else if(photons[i].isAbsorbedWall){ 
			eps_w += photons[i].W;
		}
		
		/*
		else if(photons[i].isAbsorbedInterior){
			eps_i += photons[i].W;
		}
		*/
		eps_i += photons[i].W_interior;
	}
	eps_c = eps_c/N_d;
	eps_e = eps_e/N_d;
	eps_w = eps_w/N_d;
	eps_i = eps_i/N_d;
}
		
double getMeanPathLength(std::vector<Photon> &photons, size_t &N){
	double meanPath = 0.0;
	double N_collected = 0.0;
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			meanPath += photons[i].s_tot;
			N_collected += 1.0;
		}
	}
	meanPath = (meanPath/N_collected);
	return meanPath;
}

double getMeanPathLengthWeighted(std::vector<Photon> &photons, size_t &N){
	double meanPath = 0.0;
	double collected_weight = 0.0;
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			meanPath += photons[i].W*photons[i].s_tot;
			collected_weight += photons[i].W;
		}
	}
	meanPath = (meanPath/collected_weight);
	return meanPath;
}

double getStdPathLength(std::vector<Photon> &photons, double meanPath, size_t &N){
	double stdPath = 0.0;
	double N_collected = 0.0;
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			stdPath += (photons[i].s_tot - meanPath)*(photons[i].s_tot - meanPath);
			N_collected += 1.0;
		}
	}
	if(N_collected < 1.0){
		cout << "BAD" << endl;
	}
	stdPath = std::sqrt(stdPath/(N_collected-1.0)); //WAS ERROR HERE
	return stdPath;
}

double getStdPathLengthWeighted(std::vector<Photon> &photons, double meanPath, size_t &N){
	double stdPath = 0.0;
	double N_collected = 0.0;
	double collected_weight = 0.0;
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			stdPath += photons[i].W*(photons[i].s_tot - meanPath)*(photons[i].s_tot - meanPath);
			collected_weight += photons[i].W;
			N_collected += 1.0;
		}
	}
	if(N_collected < 1.0){
		cout << "BAD" << endl;
	}
	stdPath = std::sqrt((N_collected/(collected_weight*(N_collected - 1.0)))*stdPath);
	return stdPath;
}

void savePaths(std::vector<Photon> &photons, std::string fileName, size_t &N){
	ofstream myFile;
	myFile.open(fileName.c_str());
	for(size_t i=0; i<N; i++){
		if(photons[i].isCollected){
			myFile << setprecision(17) << photons[i].s_tot << endl;
		}
	}
	myFile.close();
}

void tracePhotonAbsorbingIS(Photon &p, double &R, double &rho, double &z_s, double &cos_theta0, double &a_p, double &b_p, size_t &lim, double &alpha, Generator &G){
	double s0 = 0.0;
	double path = drawPathLength(alpha, G);
	//Move from initial position and direction
	double dot_prod = p.r.x*p.u.u_x + p.r.y*p.u.u_y + p.r.z*p.u.u_z;
//	cout << "dot_prod: " << dot_prod << endl;
	s0 = -dot_prod + sqrt(dot_prod*dot_prod + R*R - p.r.square());
//	cout << "s0 first: " << s0 << endl;

	if(s0 < path){
		p.s_tot += s0;
		p.N_w ++;
		path -= s0;
	}
	else{
		p.s_tot = path;
		p.isAbsorbedInterior = true;
		//cout << "Absorbed from launch!" << endl;
		return;
	}

	while(p.N_w < lim){
		//Move to wall
		p.r.x = p.r.x + s0*p.u.u_x;
		p.r.y = p.r.y + s0*p.u.u_y;
		p.r.z = p.r.z + s0*p.u.u_z;
		//Reflect?
		if(G.uniform() < rho){
			//Yes
		//	p.N_w ++;
			diffuseLocal(p.u, G);
			p.u = convertToGlobal(p.u, p.r, R);
			s0 = - 2.0*(p.u.u_x*p.r.x + p.u.u_y*p.r.y + p.u.u_z*p.r.z);
			//Can photon propagate distance s0?
			if(s0 < path){
				//Is photon in fiber area?
				if((p.r.z + s0*p.u.u_z) > z_s){
					double s = (z_s - p.r.z)/p.u.u_z;
					if(p.u.u_z < cos_theta0){
						p.isAbsorbedEntrance = true;
						p.s_tot += s;
						break;
					}
					else{
						double x = p.r.x + s*p.u.u_x;
						double y = p.r.y + s*p.u.u_y;
						if(((x-b_p)*(x-b_p) + y*y) < (a_p*a_p)){
							p.isCollected = true;
							//cout << "collected!" << endl;
							p.s_tot += s;
							break;
						}
						else{
							p.isAbsorbedEntrance = true;
							p.s_tot += s;
							break;
						}
					}
				}
				else{
					p.N_w ++;
					p.s_tot += s0;
					path -= s0;
				}
			}
			else{
				p.isAbsorbedInterior = true;
				p.s_tot += path;
				//cout << "Absorbed interior!" << endl;
				break;
			}
		}
		else{
			p.isAbsorbedWall = true;
			break;
		}
	}
}

void henyeyScatter(Dir &u, double &g, Generator &G){
	double u_x_temp = u.u_x;
	double u_y_temp = u.u_y;
	double u_z_temp = u.u_z;
	double chi_1 = G.uniform();
	double chi_2 = G.uniform();
	double cos_theta = (1.0/(2.0*g))*(1.0 + g*g - ((1.0-g*g)/(1.0 - g + 2.0*g*chi_1))*((1.0-g*g)/(1.0 - g + 2.0*g*chi_1)));
	double sin_theta = std::sqrt(1.0-cos_theta*cos_theta);
	double psi = 2.0*M_PI*chi_2;
	double cos_psi = std::cos(psi);
	double sin_psi = std::sin(psi);

	u.u_x = (sin_theta/(std::sqrt(1.0-u_z_temp*u_z_temp)))*(u_x_temp*u_z_temp*cos_psi - u_y_temp*sin_psi) + u_x_temp*cos_theta;
	u.u_y = (sin_theta/(std::sqrt(1.0-u_z_temp*u_z_temp)))*(u_y_temp*u_z_temp*cos_psi + u_x_temp*sin_psi) + u_y_temp*cos_theta;
	u.u_z = -sin_theta*cos_psi*std::sqrt(1.0-u_z_temp*u_z_temp) + u_z_temp*cos_theta; 

	double norm = u.norm(); 	//Force normalization
	u.u_x = u.u_x/norm;
	u.u_y = u.u_y/norm;
	u.u_z = u.u_z/norm;
}

void tracePhotonTurbidIS(Photon &p, double &R, double &rho, double &z_s, double &cos_theta0, double &a_p, double &b_p, size_t &scatter_limit, double &albedo, double &alpha_t, double &g, double &m,double &w_t,  Generator &G){

	double l = drawPathLength(alpha_t, G); 		//Choose a step size
	Coord r;

	while(p.N_s < scatter_limit){ 			//Main loop, limited by total number of scattering events
		//New position
		r.x = p.r.x + l*p.u.u_x;
		r.y = p.r.y + l*p.u.u_y;
		r.z = p.r.z + l*p.u.u_z;

		//Is position inside sphere?
		if(r.norm() < R){
			p.s_tot += l; 			//Update path length
			p.W_interior += p.W*albedo; 	//Deposit weight
			p.W -= p.W*albedo;		//Update weight	
			if(p.W < w_t){
				if(G.uniform() < (1.0/m)){
					p.W = m*p.W; 			//Survived roulette
				}
				else{
					p.isAbsorbedInterior = true;	//Did not survive roulette
					break;
				}
			}
			l = drawPathLength(alpha_t, G); 		//New step size
			henyeyScatter(p.u, g, G); 			//Get new scattering direction
			p.N_s ++; 					//Update scattering event
			p.r.x = r.x;					//Update position
			p.r.y = r.y;
			p.r.z = r.z;
		}
		else{
			//Find wall intercept
			double dot_prod = p.r.x*p.u.u_x + p.r.y*p.u.u_y + p.r.z*p.u.u_z;
			double r_square = p.r.square();
			double s0_plus = - dot_prod + std::sqrt(dot_prod*dot_prod + R*R - r_square);
			double s0_minus = - dot_prod - std::sqrt(dot_prod*dot_prod + R*R - r_square);
			double s0 = 0.0;
			if(s0_plus > 0.0){
				s0 = s0_plus;
			}else{
				s0 = s0_minus;
			}
			//Check if in port area
			double z = p.r.z + s0*p.u.u_z;
			if(z>z_s){
				//Distance to entrance port plane
				double s = (z_s - p.r.z)/p.u.u_z;
				if(p.u.u_z < cos_theta0){
					p.isAbsorbedEntrance = true; 		//Photon does not pass NA limit, is not collected
					p.s_tot += s;
					break;
				}
				else{
					double x = p.r.x + s*p.u.u_x;
					double y = p.r.y + s*p.u.u_y;
					if(((x-b_p)*(x-b_p) + y*y) < (a_p*a_p)){
						p.isCollected = true;		//Photon is collected
						p.s_tot += s;
						break;
					}
					else{
						p.isAbsorbedEntrance = true;	//Photon is lost to entrance port
						p.s_tot += s;
						break;
					}
				}
			}
			//Photon is not in entrance, move it
			p.r.x = p.r.x + s0*p.u.u_x;
			p.r.y = p.r.y + s0*p.u.u_y;
			p.r.z = p.r.z + s0*p.u.u_z;
			p.s_tot += s0; 			//Update path length
			l -= s0;
			//Photon is at wall, reflected?
			if(G.uniform() < rho){
				diffuseLocal(p.u, G);
				p.u = convertToGlobal(p.u, p.r, R); 		//Photon is reflected, calculate new direction
				p.N_w ++;
			}
			else{
				p.isAbsorbedWall = true;			//Photon is transmitted through the wall
				break;
			}
		}
	}
}



