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
    RunInfo(const std::string &path)
    {
        auto remove_space = [](std::string &line) -> void
        {
            while (line[0] == ' ')
            {
                line = line.substr(1, line.length() - 1);
            }
        };

        std::ifstream file(path.c_str());
        for (std::string line; std::getline(file, line);)
        {
            remove_space(line);
            if (line[0] == '#' || line.length() == 0)
            {
                continue;
            }
            std::istringstream ss(line);
            std::string system, badmap_ver, trigger;
            unsigned int runID_first, runID_last, shadowbar;
            ss >> system >> runID_first >> runID_last >> badmap_ver >> shadowbar >> trigger;
            this->AddRunInfo(system, runID_first, runID_last, badmap_ver, shadowbar, trigger);
        }

        for (auto &[system, indices] : this->RunID)
        {
            this->RunNumbers[system] = indices.size();
        }
    }
    ~RunInfo(){};

    void AddRunInfo(const std::string &system, const unsigned int &runID_first, const unsigned int &runID_last, const std::string &badmap_ver, const unsigned int &shadowbar, const std::string &trigger)
    {
        for (unsigned int id = runID_first; id <= runID_last; id++)
        {
            this->RunID[system].push_back(id);
            this->BadMapVersion[system].push_back(badmap_ver);
            this->TriggerCondition[system].push_back(trigger);
            this->ShadowBars[system].push_back(shadowbar);
        }
        return;
    }

    std::map<std::string, unsigned int> RunNumbers;
    std::map<std::string, std::vector<unsigned int>> RunID;
    std::map<std::string, std::vector<std::string>> BadMapVersion;
    std::map<std::string, std::vector<std::string>> TriggerCondition;
    std::map<std::string, std::vector<unsigned int>> ShadowBars;
    std::map<std::string, double> BetaCMS;
    std::map<std::string, double> BeamRapidity;

    int GetRunNumber(const std::string &system) { return this->RunNumbers[system]; }

    unsigned int GetRunIndex(const std::string &system, const unsigned int &iRun) { return (iRun < this->RunNumbers[system]) ? this->RunID[system][iRun] : 9999; }

    std::string GetBadMapVersion(const std::string &system, const unsigned int &iRun) { return (iRun < this->RunNumbers[system]) ? this->BadMapVersion[system][iRun] : ""; }

    std::string GetTriggerCondition(const std::string &system, const unsigned int &iRun) { return (iRun < this->RunNumbers[system]) ? this->TriggerCondition[system][iRun] : ""; }

    unsigned int GetShadowBar(const std::string &system, const unsigned int &iRun) { return (iRun < this->RunNumbers[system]) ? this->ShadowBars[system][iRun] : 9999; }
};
#endif