#include "RunInfo.hh"

RunInfo::RunInfo(const std::string &path)
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

RunInfo::~RunInfo(){};

/**
 * @brief add run information for each batch
 *
 * @param system name of reaction system
 * @param runID_first first ID of the batch, e.g. CalibratedData_${4023}.root
 * @param runID_last last ID of the batch
 * @param badmap_ver badmap version of the batch
 * @param shadowbar
 * @param trigger
 */
void RunInfo::AddRunInfo(const std::string &system, const unsigned int &runID_first, const unsigned int &runID_last, const std::string &badmap_ver, const unsigned int &shadowbar, const std::string &trigger)
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

/**
 * @brief Get the ID of the i-th file. e.g. CalibratedData_${ID}.root
 *
 * @param system name of the reaction system
 * @param iRun i-th file
 * @return unsigned int
 */
unsigned int RunInfo::GetRunIndex(const std::string &system, const unsigned int &iRun)
{
    return (iRun < this->RunNumbers[system]) ? this->RunID[system][iRun] : 9999;
}

/**
 * @brief Get the bad map version of the i-th file for a given system
 *
 * @param system
 * @param iRun
 * @return std::string
 */
std::string RunInfo::GetBadMapVersion(const std::string &system, const unsigned int &iRun)
{
    return (iRun < this->RunNumbers[system]) ? this->BadMapVersion[system][iRun] : "";
}

/**
 * @brief Get the trigger condition of the i-th file for a given system
 *
 * @param system
 * @param iRun
 * @return std::string
 */
std::string RunInfo::GetTriggerCondition(const std::string &system, const unsigned int &iRun)
{
    return (iRun < this->RunNumbers[system]) ? this->TriggerCondition[system][iRun] : "";
}

unsigned int RunInfo::GetShadowBar(const std::string &system, const unsigned int &iRun)
{
    return (iRun < this->RunNumbers[system]) ? this->ShadowBars[system][iRun] : 9999;
}