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
};
#endif