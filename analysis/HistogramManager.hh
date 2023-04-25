#ifndef HisotogramManager_hh
#define HisotogramManager_hh

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TH2D.h"

#include "Particle.hh"

class HistogramManager
{
public:
    HistogramManager(const std::string &str);
    ~HistogramManager();
    void _Initialize();
    void Fill(const Particle &particle, const double &weight = 1.);
    void Normalize(const double &norm);
    void Write();

private:
    std::string mStr;
    std::vector<std::string> mParticlesNames;
    std::map<std::string, TH2D *> h2_pta_rapidity_lab;

protected:
    std::vector<std::string> PARTICLES = {
        "p",
        "d",
        "t",
        "3He",
        "4He",
    };
};

#endif