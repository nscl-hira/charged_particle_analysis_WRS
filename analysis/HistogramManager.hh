#ifndef HisotogramManager_hh
#define HisotogramManager_hh

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TH2D.h"

#include "AME.hh"
#include "Particle.hh"

class HistogramManager
{
public:
    HistogramManager(const std::string &str);
    ~HistogramManager();
    void _Initialize();
    void SetAcceptedParticles(const std::vector<std::string> &particles);
    void Fill(const Particle &particle, const double &weight = 1.);
    void Normalize(const double &norm);
    void Write();

private:
    std::string mStr;
    std::vector<std::string> mParticles;
    std::map<std::string, TH2D *> h2_pta_rapidity_lab;
    std::map<std::string, TH2D *> h2_theta_kinery_lab;

protected:
    std::vector<std::string> ACCEPTED_PARTICLES = {
        "p",
        "d",
        "t",
        "3He",
        "4He",
    };
};

#endif