#ifndef ReactionLost_hh
#define ReactionLost_hh

#include <map>
#include <array>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "TF1.h"
#include "TMath.h"

#include "AME.hh"

class ReactionLost
{
public:
    ReactionLost(const std::string &pth, const std::string &fcn_name = "f1_ReactionLost_CorEff");
    ~ReactionLost() { ; }
    double Get_ReactionLost_CorEff(const int &z, const int &a, const double &Ekin_Lab);

    void _inititalize(const std::string &pth);
    bool _accepted(const std::string &name);
    void _define_acceptable_particles();

    void set_heavy_particle_substitute(const std::string &name);
    void set_efficiency_function(const std::string &fcn) { this->efficiency_function = fcn; }
    void set_efficiency_range(const double &min, const double &max) { this->efficiency_range = {min, max}; }

private:
    std::string function_name;
    std::string efficiency_function;
    std::array<double, 2> efficiency_range;

    std::map<std::string, TF1 *> f1_ReactionLost_CorEff;
    std::vector<std::string> acceptable_particles;
    std::string heavy_particle_substitute;

protected:
    const double unphysical_effeciency = 1.;
    const std::array<double, 2> default_efficiency_range = {0, 800};
    const std::string default_efficiency_function = "TMath::Exp([0] * x * x + [1] * x)";
    const std::vector<std::string> default_acceptable_particles = {"p", "d", "t", "3He", "4He"};
    // currently, Sean only check the reaction lost correction up to 4He, the heavier particle will be taken as 4He.
    const std::string default_heavy_particle_substitute = "4He";
};

#endif
