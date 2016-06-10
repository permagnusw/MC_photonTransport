#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "Coord.h"
#include "Dir.h"
#include "Photon.h"
#include "RandomGenerator.h"
#include <vector>
#include <string>

using namespace photon_IS;


void diffuseLocal(Dir &u, Generator &G);

Dir convertToGlobal(Dir &u, Coord &r, double &R);

double drawPathLength(double &alpha, Generator &G);

double psTheory(double &alpha, double &R);

void traceSingleDiffuse(Photon &p, double &R, double &rho, size_t &lim, Generator &G);

void traceSingleDiffuseAbsorption(Photon &p, double &R, double &rho, double &alpha, size_t &lim, Generator &G);

void traceCollectionDiffuse(std::vector<Photon> &P, double &R, double &rho, size_t &lim, Generator &G);

double meanNumberCollision(std::vector<Photon> &P);

double stddevNumberCollision(std::vector<Photon> &P, double &mean);

double meanPathLength(std::vector<Photon> &P);

double meanPathLengthCollected(std::vector<Photon> &photons);

double stddevPathLength(std::vector<Photon> &P, double &mean);

double stddevPathLengthCollected(std::vector<Photon> &P, double &mean);

void histCollisions(std::vector<Photon> &photons, std::string fileName);

void histPathLength(std::vector<Photon> &photons, std::string fileName, double min_s, double max_s, size_t bins);

void savePhotons(std::vector<Photon> &photons, std::string fileName, double R, double rho, double N_photons, size_t seed);

void traceSingleDiffuseInvariant(Photon &p, double &R, double &rho, size_t &lim, Generator &G);

void initPhotons(std::vector<Photon> &photons, double R, double z_s, double cos_theta0, double a_s, double b_s, Generator &G);

void tracePhotonEmptyIS(Photon &p, double &R, double &rho, double &z_s, double &cos_theta0, double &a_p, double &b_p, size_t &lim, Generator &G);

void tracePhotonAbsorbingIS(Photon &p, double &R, double &rho, double &z_s, double &cos_theta0, double &a_p, double &b_p, size_t &lim, double &alpha, Generator &G);

void getStats(std::vector<Photon> &photons, double &eps_c, double &eps_e, double &eps_w, size_t &N);

void getStatsWeighted(std::vector<Photon> &photons, double &eps_c, double &eps_e, double &eps_w, double &eps_i, size_t &N);


double getMeanPathLength(std::vector<Photon> &photons, size_t &N);

double getMeanPathLengthWeighted(std::vector<Photon> &photons, size_t &N);

double getStdPathLength(std::vector<Photon> &photons, double meanPath, size_t &N);

double getStdPathLengthWeighted(std::vector<Photon> &photons, double meanPath, size_t &N);

void savePaths(std::vector<Photon> &photons, std::string fileName, size_t &N);

void henveyScatter(Dir &u, double &g, Generator &G);

void tracePhotonTurbidIS(Photon &p, double &R, double &rho, double &z_s, double &cos_theta0, double &a_p, double &b_p, size_t &scatter_limit, double &albedo, double &alpha_t, double &g, double &m,double &w_t,  Generator &G);


#endif

