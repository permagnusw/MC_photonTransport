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

	size_t seed = 5001;
	Generator G(seed);
	//GEOMETRY AND OPTICAL CONSTANTS
	double R = 0.5;
	double a_s = 0.105/2.0;
	double a_p = 0.105/2.0;
	double b_s = 0.125/2.0;
	double b_p = 0.125/2.0;
	double r_e = 2.0*b_p;
	double r_port = 0.0;
	double alpha = 1.2786;
	double rho = 0.99;
	double n = 1.0;
	double NA_s = 0.22;
	double NA = 0.22;
	//double cos_theta0 = sqrt(1.0 - (NA*NA)/(n*n));
	double cos_theta0_s = sqrt(1.0 - (NA_s*NA_s)/(n*n));
	//double z_s = sqrt(R*R - r_e*r_e);

	size_t collision_limit = 50000;
	size_t N_photons = 1e6;
	double N_d = (double) N_photons;


	
	double NA_0 = 0.05;
	double d_NA = 0.025;
	int N_R = 25;
	int N_avg = 100;
	double N_avg_d = (double) N_avg;
	vector<double> NA_vec(N_R);

	vector<double> eps_c_mean(N_R);
	vector<double> eps_e_mean(N_R);
	vector<double> eps_w_mean(N_R);

	vector<double> eps_c_std(N_R);
	vector<double> eps_e_std(N_R);
	vector<double> eps_w_std(N_R);

	vector<double> mean_path_avg(N_R);
	vector<double> std_mean_path_avg(N_R);
	vector<double> std_path_avg(N_R);
	vector<double> std_std_path_avg(N_R);

	NA_vec[0] = NA_0;
	for(int i=1; i<N_R; i++){
		NA_vec[i] = NA_vec[i-1] + d_NA;
	}

	for(int j=0; j<N_R; j++){
		double z_s = sqrt(R*R - r_e*r_e);
		double cos_theta0 = sqrt(1.0 - (NA_vec[j]*NA_vec[j])/(n*n));
		cout << "cos theta0 " << cos_theta0 << endl;


		vector<double> eps_c(N_avg);
		vector<double> eps_e(N_avg);
		vector<double> eps_w(N_avg);

		vector<double> mean_path(N_avg);
		vector<double> std_path(N_avg);

		for(int k=0; k<N_avg; k++){

			vector<Photon> photons(N_photons);
			initPhotons(photons, R, z_s, cos_theta0_s, a_s, b_s, G);

			eps_c[k] = 0.0;
			eps_e[k] = 0.0;
			eps_w[k] = 0.0;

			mean_path[k] = 0.0;
			std_path[k] = 0.0;

			for(size_t i=0; i<N_photons; i++){
				tracePhotonEmptyIS(photons[i], R, rho, z_s, cos_theta0, a_p, b_p, collision_limit, G);
			}
			getStats(photons, eps_c[k], eps_e[k], eps_w[k], N_photons);
			//savePaths(photons, "data/emptyIS-pathLength/paths.txt",N_photons);
			mean_path[k] = getMeanPathLength(photons,N_photons);
			std_path[k] = getStdPathLength(photons, mean_path[k], N_photons);
			//cout << "mean path: " << meanPath << " std path: " << stdPath << endl;
			//cout << "stats: " << eps_c << " " << eps_e << " " << eps_w << endl;
		}

		eps_c_mean[j] = 0.0;
		eps_e_mean[j] = 0.0;
		eps_w_mean[j] = 0.0;

		mean_path_avg[j] = 0.0;
		std_path_avg[j] = 0.0;

		for(int i=0; i<N_avg; i++){
			eps_c_mean[j] += eps_c[i];
			eps_e_mean[j] += eps_e[i];
			eps_w_mean[j] += eps_w[i];

			mean_path_avg[j] += mean_path[i];
			std_path_avg[j] += std_path[i];
		}
		eps_c_mean[j] = eps_c_mean[j]/N_avg_d;
		eps_e_mean[j] = eps_e_mean[j]/N_avg_d;
		eps_w_mean[j] = eps_w_mean[j]/N_avg_d;

		mean_path_avg[j] = mean_path_avg[j]/N_avg_d;
		std_path_avg[j] = std_path_avg[j]/N_avg_d;



		eps_c_std[j] = 0.0;
		eps_e_std[j] = 0.0;
		eps_w_std[j] = 0.0;

		std_mean_path_avg[j] = 0.0;
		std_std_path_avg[j] = 0.0;

		for(int i=0; i<N_avg; i++){
			eps_c_std[j] += (eps_c[i] - eps_c_mean[j])*(eps_c[i] - eps_c_mean[j]);
			eps_e_std[j] += (eps_e[i] - eps_e_mean[j])*(eps_e[i] - eps_e_mean[j]);
			eps_w_std[j] += (eps_w[i] - eps_w_mean[j])*(eps_w[i] - eps_w_mean[j]);

			std_mean_path_avg[j] += (mean_path[i] - mean_path_avg[j])*(mean_path[i] - mean_path_avg[j]);
			std_std_path_avg[j] += (std_path[i] - std_path_avg[j])*(std_path[i] - std_path_avg[j]);

		}

		eps_c_std[j] = sqrt(eps_c_std[j]/(N_avg_d - 1.0));
		eps_e_std[j] = sqrt(eps_e_std[j]/(N_avg_d - 1.0));
		eps_w_std[j] = sqrt(eps_w_std[j]/(N_avg_d - 1.0));

		std_mean_path_avg[j] = sqrt(std_mean_path_avg[j]/(N_avg_d - 1.0));
		std_std_path_avg[j] = sqrt(std_std_path_avg[j]/(N_avg_d - 1.0));


		cout << "Finished number: " << j << endl;
	}

	/*
	for(int i=0; i<N_R; i++){
		cout << "R: " << R_vec[i] << endl;
		cout << "eps_c: " << eps_c_mean[i] << "  eps_e: " << eps_e_mean[i] << "  eps_w: " << eps_w_mean[i] << endl;
		cout << "eps_c_std: " << eps_c_std[i] << "  eps_e_std: " << eps_e_std[i] << "  eps_w_std: " << eps_w_std[i] << endl << endl;

	}
	*/
	cout << "numbers generated: " << G.count << endl;

	ofstream myFile;
	myFile.open("data/transparentIS-NA_final/R0.5rho99_t.txt");
	myFile << "#Empty integrating sphere. Lambertian reflectance. No ports for fluid." << endl
		<< "#Length units in [mm]" << endl
		<< "#Reflectivity of interior wall: rho = " << rho << endl
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
		<< "#Sphere radius: R = " << R << endl
		<< "#Number of photons: N = " << N_photons << endl
		<< "#Number of averages: N_avg = " << N_avg << endl
		<< "#Number of NA-values: N_R = " << N_R << endl
		<< "#Delta NA: d_NA = " << d_NA << endl
		<< "#Seed: " << seed << endl << "#" << endl
		<< "#NA <eps_c> <eps_e> <eps_w> std_eps_c std_eps_e std_eps_w <<path>> std<path> <std_path> <<std_path>>" << endl;
	for(int i=0; i<N_R; i++){
		myFile << setprecision(5) << fixed << NA_vec[i] << " " << setprecision(10) << fixed << eps_c_mean[i] << " " << eps_e_mean[i] << " " << eps_w_mean[i]
		       	<< " " << eps_c_std[i] << " " << eps_e_std[i] << " " << eps_w_std[i] << " " << mean_path_avg[i] << " " << std_mean_path_avg[i]
		        << " " << std_path_avg[i] << " " << std_std_path_avg[i] << endl;
	}
	myFile.close();

	/*
	cout << "entrance fraction: " << (r_e*r_e/(4.0*R*R)) << endl;
	cout << "A_core/A_entrance: " << (a_p*a_p/(r_e*r_e)) << endl;
	cout << "Rho->inf: " << (a_p*a_p/(r_e*r_e))*(1.0-cos_theta0) << endl;
	cout << "Collected: " << collected/N_d << endl;
	cout << "entrance: " << absorbedEntrance/N_d << endl;
	cout << "wall: " << absorbedWall/N_d << endl;
	cout << "sum: " << collected/N_d + absorbedEntrance/N_d + absorbedWall/N_d << endl;
	

	//double Ps = psTheory(alpha,R);
	//double meanCollision = meanNumberCollision(photons);
	//cout << "mean collision: " << meanCollision << endl;
	//cout << "mean collision theoretical: " << (Ps/(1.0-Ps*rho)) << endl;
	//cout << "stddev collision: " << stddevNumberCollision(photons, meanCollision) << endl; 
	double meanPath = meanPathLengthCollected(photons);
	cout << "mean path: " << meanPath << endl;
	//cout << "mean path theoretical: " << ((1.0/alpha)*((1.0-Ps)/(1.0-rho*Ps))) << endl;
	cout << "stddev path: " << stddevPathLengthCollected(photons, meanPath) << endl;
	*/

	return 0;
}
	
