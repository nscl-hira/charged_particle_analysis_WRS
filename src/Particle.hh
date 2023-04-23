#ifndef Particle_hh
#define Particle_hh

#include "TMath.h"
#include "Physics.hh"

struct Particle
{
    int N, Z;
    double px, py, pz;
    double pmag_trans, pmag_lab, kinergy_lab, theta_lab, phi, rapidity_lab, rapidity_lab_normed;
    void initialize(const double &betacms, const double &beam_rapidity, const double &m = 0.);
};

#endif