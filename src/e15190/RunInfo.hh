#ifndef RunInfo_hh
#define RunInfo_hh

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <map>
#include <sstream>

#include "TMath.h"
#include "TString.h"

class RunInfo
{
public:
    RunInfo(const std::string &path);
    ~RunInfo();

    void AddRunInfo(const std::string &system, const unsigned int &runID_first, const unsigned int &runID_last, const std::string &badmap_ver, const unsigned int &shadowbar, const std::string &trigger);

    std::map<std::string, unsigned int> RunNumbers;
    std::map<std::string, std::vector<unsigned int>> RunID;
    std::map<std::string, std::vector<std::string>> BadMapVersion;
    std::map<std::string, std::vector<std::string>> TriggerCondition;
    std::map<std::string, std::vector<unsigned int>> ShadowBars;

    std::map<std::string, double> BetaCMS;
    std::map<std::string, double> BeamRapidity;

    /**
     * @brief Get total number of runs (files) for a given reaction system
     *
     * @param system name of the reaction system
     * @return int
     */
    int GetRunNumber(const std::string &system) { return this->RunNumbers[system]; }
    unsigned int GetRunIndex(const std::string &system, const unsigned int &iRun);
    std::string GetBadMapVersion(const std::string &system, const unsigned int &iRun);
    std::string GetTriggerCondition(const std::string &system, const unsigned int &iRun);
    unsigned int GetShadowBar(const std::string &system, const unsigned int &iRun);

    double GetBetaCMS(const std::string &system) { return this->BetaCMS[system]; }
    double GetBeamRapidity(const std::string &system) { return this->BeamRapidity[system]; }

protected:
    // unit = MeV/c^2
    double mass_1u = 931.49410242;
    std::vector<std::string> beam_name = {"Ca40", "Ca48"};
    std::vector<std::string> target_name = {"Ni58", "Ni64", "Sn112", "Sn124"};

    // unit = mass_1u
    std::vector<double> beam_mass = {39.962590866, 47.95252290};
    std::vector<double> target_mass = {57.935343, 63.9279660, 111.904818, 123.9052739};

    std::vector<int> beam_energy = {56, 140};
    std::vector<int> beam_A = {40, 48};
};
#endif