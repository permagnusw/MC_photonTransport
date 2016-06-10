#include "Coord.h"
#include "Dir.h"
#include "Photon.h"
#include "RandomGenerator.h"
#include "functions.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>


using namespace photon_IS;
using namespace std;

int main(){

	size_t seed = 9922;			//Seed
	Generator G(seed);
	//GEOMETRY AND OPTICAL CONSTANTS
	double R = 0.5;				//Cavity radius
	double a_s = 0.105/2.0;			//Fiber geometry
	double a_p = 0.105/2.0;
	double b_s = 0.125/2.0;
	double b_p = 0.125/2.0;
	double r_e = 2.0*b_p;			//Port radius

	double alpha_s = 3.0;			//Scattering coefficient
	double m = 10.0;			//Survival rate
	double g = 0.6;				//Anisotropy factor
	double w_t = 0.00001;			//Weight threshold
	double rho = 0.95;			//Wall reflectivity
	double n = 1.0;				//Refractive index
	double NA_s = 0.22;			//NA of source fiber
	double NA = 0.22;			//NA of collection fiber
	double cos_theta0 = sqrt(1.0 - (NA*NA)/(n*n));
	double cos_theta0_s = sqrt(1.0 - (NA_s*NA_s)/(n*n));
	double z_s = sqrt(R*R - r_e*r_e);

	size_t collision_limit = 50000;		//Max number of wall collisions
	size_t scatter_limit = 1000;		//Max number of scattering evenst
	size_t N_photons = 1e6;			//Number of photons launched
	double N_d = (double) N_photons;

	
	double alpha_a0 = 0.01;			//Initial absorption coefficient
	double d_alpha_a = 0.04;		//Delta absorption coefficient
	int N_alpha_a = 25;			//Number of absorption coefficients
	int N_avg = 100;			//Number of system realizations
	double N_avg_d = (double) N_avg;

	vector<double> alpha_a_vec(N_alpha_a);	//Absorption coefficients
	vector<double> alpha_s_vec(N_alpha_a);
	vector<double> albedo_vec(N_alpha_a);
	vector<double> alpha_t_vec(N_alpha_a);

	vector<double> eps_c_mean(N_alpha_a);	//Book keeping, statistics
	vector<double> eps_e_mean(N_alpha_a);
	vector<double> eps_w_mean(N_alpha_a);
	vector<double> eps_i_mean(N_alpha_a);

	vector<double> eps_c_std(N_alpha_a);
	vector<double> eps_e_std(N_alpha_a);
	vector<double> eps_w_std(N_alpha_a);
	vector<double> eps_i_std(N_alpha_a);

	vector<double> mean_path_avg(N_alpha_a);
	vector<double> std_mean_path_avg(N_alpha_a);
	vector<double> std_path_avg(N_alpha_a);
	vector<double> std_std_path_avg(N_alpha_a);

	alpha_a_vec[0] = alpha_a0;
	alpha_t_vec[0] = alpha_a_vec[0] + alpha_s;
	albedo_vec[0] = alpha_a_vec[0]/(alpha_s + alpha_a_vec[0]);

	for(int i=1; i<N_alpha_a; i++){
		alpha_a_vec[i] = alpha_a_vec[i-1] + d_alpha_a;
		alpha_t_vec[i] = alpha_a_vec[i] + alpha_s;
		albedo_vec[i] = alpha_a_vec[i]/alpha_t_vec[i];

	}

	//REFERENCE, alpha = 0
	double eps_c_ref_mean = 0.0;
	double eps_e_ref_mean = 0.0;
	double eps_w_ref_mean = 0.0;
	double eps_i_ref_mean = 0.0;
	for(int i=0; i<N_avg; i++){
		vector<Photon> photons(N_photons);
		initPhotons(photons, R, z_s, cos_theta0_s, a_s, b_s, G);
		for(size_t j=0; j<N_photons; j++){
			tracePhotonEmptyIS(photons[j], R, rho, z_s, cos_theta0, a_p, b_p, collision_limit, G);
		}
		double eps_c_temp = 0.0;
		double eps_e_temp = 0.0;
		double eps_w_temp = 0.0;
		double eps_i_temp = 0.0;
		getStatsWeighted(photons, eps_c_temp, eps_e_temp, eps_w_temp, eps_i_temp, N_photons);
		eps_c_ref_mean += eps_c_temp;	
		eps_e_ref_mean += eps_e_temp;	
		eps_w_ref_mean += eps_w_temp;	
		eps_i_ref_mean += eps_i_temp;
	}
	eps_c_ref_mean = eps_c_ref_mean/N_avg_d;
	eps_e_ref_mean = eps_e_ref_mean/N_avg_d;
	eps_w_ref_mean = eps_w_ref_mean/N_avg_d;
	eps_i_ref_mean = eps_i_ref_mean/N_avg_d;

	//Rest of simulation
	for(int j=0; j<N_alpha_a; j++){

		vector<double> eps_c(N_avg);
		vector<double> eps_e(N_avg);
		vector<double> eps_w(N_avg);
		vector<double> eps_i(N_avg);

		vector<double> mean_path(N_avg);
		vector<double> std_path(N_avg);

		for(int k=0; k<N_avg; k++){

			vector<Photon> photons(N_photons);
			initPhotons(photons, R, z_s, cos_theta0_s, a_s, b_s, G); //CHECK NA source

			eps_c[k] = 0.0;
			eps_e[k] = 0.0;
			eps_w[k] = 0.0;
			eps_i[k] = 0.0;

			mean_path[k] = 0.0;
			std_path[k] = 0.0;

			for(size_t i=0; i<N_photons; i++){
				tracePhotonTurbidIS(photons[i], R, rho, z_s, cos_theta0, a_p, b_p, scatter_limit, albedo_vec[j], alpha_t_vec[j], g, m, w_t, G);

			}
			getStatsWeighted(photons, eps_c[k], eps_e[k], eps_w[k], eps_i[k], N_photons);
			mean_path[k] = getMeanPathLengthWeighted(photons,N_photons);
			std_path[k] = getStdPathLengthWeighted(photons, mean_path[k], N_photons);
		}

		eps_c_mean[j] = 0.0;
		eps_e_mean[j] = 0.0;
		eps_w_mean[j] = 0.0;
		eps_i_mean[j] = 0.0;

		mean_path_avg[j] = 0.0;
		std_path_avg[j] = 0.0;

		for(int i=0; i<N_avg; i++){
			eps_c_mean[j] += eps_c[i];
			eps_e_mean[j] += eps_e[i];
			eps_w_mean[j] += eps_w[i];
			eps_i_mean[j] += eps_i[i];

			mean_path_avg[j] += mean_path[i];
			std_path_avg[j] += std_path[i];
		}
		eps_c_mean[j] = eps_c_mean[j]/N_avg_d;
		eps_e_mean[j] = eps_e_mean[j]/N_avg_d;
		eps_w_mean[j] = eps_w_mean[j]/N_avg_d;
		eps_i_mean[j] = eps_i_mean[j]/N_avg_d;

		mean_path_avg[j] = mean_path_avg[j]/N_avg_d;
		std_path_avg[j] = std_path_avg[j]/N_avg_d;

		eps_c_std[j] = 0.0;
		eps_e_std[j] = 0.0;
		eps_w_std[j] = 0.0;
		eps_i_std[j] = 0.0;

		std_mean_path_avg[j] = 0.0;
		std_std_path_avg[j] = 0.0;

		for(int i=0; i<N_avg; i++){
			eps_c_std[j] += (eps_c[i] - eps_c_mean[j])*(eps_c[i] - eps_c_mean[j]);
			eps_e_std[j] += (eps_e[i] - eps_e_mean[j])*(eps_e[i] - eps_e_mean[j]);
			eps_w_std[j] += (eps_w[i] - eps_w_mean[j])*(eps_w[i] - eps_w_mean[j]);
			eps_i_std[j] += (eps_i[i] - eps_i_mean[j])*(eps_i[i] - eps_i_mean[j]);


			std_mean_path_avg[j] += (mean_path[i] - mean_path_avg[j])*(mean_path[i] - mean_path_avg[j]);
			std_std_path_avg[j] += (std_path[i] - std_path_avg[j])*(std_path[i] - std_path_avg[j]);

		}

		eps_c_std[j] = sqrt(eps_c_std[j]/(N_avg_d - 1.0));
		eps_e_std[j] = sqrt(eps_e_std[j]/(N_avg_d - 1.0));
		eps_w_std[j] = sqrt(eps_w_std[j]/(N_avg_d - 1.0));
		eps_i_std[j] = sqrt(eps_i_std[j]/(N_avg_d - 1.0));

		std_mean_path_avg[j] = sqrt(std_mean_path_avg[j]/(N_avg_d - 1.0));
		std_std_path_avg[j] = sqrt(std_std_path_avg[j]/(N_avg_d - 1.0));

	}

	//Save all data
	ofstream myFile;
	myFile.open("data/turbidIS_final/full_rho95_g0.6_alphaS_3_25.txt");
	myFile << "#Empty integrating sphere. Lambertian reflectance. No ports for fluid." << endl
		<< "#Length units in [mm]" << endl
		<< "#Reflectivity of interior wall: rho = " << rho << endl
		<< "#Scattering coefficient: alpha_s = " << alpha_s << endl
		<< "#Anisotropy: g " << g << endl
		<< "#Refractive index interior medium: n = " << n << endl << "#" << endl
		<< "#Source fiber:" << endl
		<< "# - core diameter: 2*a_s = " << 2.0*a_s << endl
		<< "# - clad diameter: 2*b_s = " << 2.0*b_s << endl
		<< "# - numerical aperture: " << "NA_s = " << NA_s << endl << "#" << endl
		<< "#Pick-up fiber:" << endl
		<< "# - core diameter: 2*a_p = " << 2.0*a_p << endl
		<< "# - clad diameter: 2*b_p = " << 2.0*b_p << endl
		<< "# - numerical aperture: " << "NA_p = " << NA << endl << "#" << endl
		<< "#Entrance radius: r_e = 2*b_p = " << r_e << endl
		<< "#Number of photons: N = " << N_photons << endl
		<< "#Number of averages: N_avg = " << N_avg << endl
		<< "#Number of alpha_s-values: N_alpha_s = " << N_alpha_a << endl
		<< "#Delta alpha_a: d_alpha_a = " << d_alpha_a << endl
		<< "#Seed: " << seed << endl << "#" << endl
		<< "#alpha_a <eps_c> <eps_e> <eps_w> <eps_i> std_eps_c std_eps_e std_eps_w std_eps_i <<path>> std<path> <std_path> <<std_path>>" << endl;
	for(int i=0; i<N_alpha_a; i++){
		myFile << setprecision(5) << fixed << alpha_a_vec[i] << " " << setprecision(10) << fixed << eps_c_mean[i] << " " << eps_e_mean[i] << " " << eps_w_mean[i]
		       	<< " " << eps_i_mean[i] << " " << eps_c_std[i] << " " << eps_e_std[i] << " " << eps_w_std[i] << " " << eps_i_std[i] << " " << mean_path_avg[i] 
			<< " " << std_mean_path_avg[i] << " " << std_path_avg[i] << " " << std_std_path_avg[i] << endl;
	}
	myFile.close();

	//Compute absorbance and transmission
	vector<double> A(N_alpha_a);
	vector<double> T(N_alpha_a);
	vector<double> L_eff(N_alpha_a);
	for(size_t i=0; i<N_alpha_a; i++){
		T[i] = eps_c_mean[i]/eps_c_ref_mean;
		A[i] = std::log(1/T[i]);
		L_eff[i] = A[i]/alpha_a_vec[i];
	}

	//Save absorbance and transmission data
	ofstream myFile2;
	myFile2.open("data/turbidIS_final/A_rho95_g0.6_alphaS_3_25.txt");
	myFile2 << "#Empty integrating sphere. Lambertian reflectance. No ports for fluid." << endl
		<< "#Length units in [mm]" << endl
		<< "#Reflectivity of interior wall: rho = " << rho << endl
		<< "#Scattering coefficient: alpha_s = " << alpha_s << endl
		<< "#Anisotropy: g = " << g << endl
		<< "#Radius of sphere: R = " << R << endl
		<< "#Refractive index interior medium: n = " << n << endl << "#" << endl
		<< "#Source fiber:" << endl
		<< "# - core diameter: 2*a_s = " << 2.0*a_s << endl
		<< "# - clad diameter: 2*b_s = " << 2.0*b_s << endl
		<< "# - numerical aperture: " << "NA_s = " << NA_s << endl << "#" << endl
		<< "#Pick-up fiber:" << endl
		<< "# - core diameter: 2*a_p = " << 2.0*a_p << endl
		<< "# - clad diameter: 2*b_p = " << 2.0*b_p << endl
		<< "# - numerical aperture: " << "NA_p = " << NA << endl << "#" << endl
		<< "#Entrance radius: r_e = 2*b_p = " << r_e << endl
		<< "#Number of photons: N = " << N_photons << endl
		<< "#Number of averages: N_avg = " << N_avg << endl
		<< "#Number of alpha-values: N_alpha_s = " << N_alpha_a << endl
		<< "#Delta alpha_a: d_alpha_a = " << d_alpha_a << endl
		<< "#Seed: " << seed << endl << "#" << endl
		<< "#alpha T A L_eff" << endl;
	for(int i=0; i<N_alpha_a; i++){
		myFile2 << setprecision(5) << fixed << alpha_a_vec[i] << " " << setprecision(10) << fixed << T[i] << " " << A[i] << " " << L_eff[i] << endl; 
	}
	myFile2.close();

	return 0;
}
	
