
#include "runinfo.hh"

runinfo::runinfo(const std::string &path)
{
    for (unsigned int i = 0; i < this->beam_mass.size(); i++)
    {
        this->beam_mass[i] = this->beam_mass[i] * this->mass_1u;
    }
    for (unsigned int i = 0; i < this->target_mass.size(); i++)
    {
        this->target_mass[i] = this->target_mass[i] * this->mass_1u;
    }

    for (unsigned int beam_id = 0; beam_id < this->beam_name.size(); beam_id++)
    {
        for (unsigned int target_id = 0; target_id < this->target_name.size(); target_id++)
        {
            for (unsigned int en_id = 0; en_id < this->beam_energy.size(); en_id++)
            {

                double beam_ke = this->beam_energy[en_id] * this->beam_A[beam_id];
                double beam_energy_tot = beam_ke + this->beam_mass[beam_id];
                double mom_beam = TMath::Sqrt(pow(beam_ke, 2.) + 2. * beam_ke * this->beam_mass[beam_id]);
                double gamma = beam_energy_tot / this->beam_mass[beam_id];

                std::string system = Form("%s%sE%i", this->beam_name[beam_id].c_str(), this->target_name[target_id].c_str(), this->beam_energy[en_id]);

                //                 std::cout << system << "\t" << this->beam_mass[beam_id] << "\t" << this->target_mass[target_id] << "\t" << this->beam_energy[en_id] << std::endl;
                this->BetaCMS[system] = mom_beam / (gamma * this->beam_mass[beam_id] + this->target_mass[target_id]);
                this->BeamRapidity[system] = 0.5 * TMath::Log((beam_energy_tot + mom_beam) / (beam_energy_tot - mom_beam));
            }
        }
    }

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
        // std::cout << system << "\t " << runID_first << "\t " << runID_last << "\t " << badmap_ver << "\t " << shadowbar << "\t " << trigger << std::endl;
        this->AddRunInfo(system, runID_first, runID_last, badmap_ver, shadowbar, trigger);
    }

    for (auto &[system, indices] : this->RunID)
    {
        this->RunNumbers[system] = indices.size();
    }

    //     for (auto& [system, number] : this->RunNumbers) {
    //         std::cout << system << "\t" << number << std::endl;
    //     }

    //     std::cout << "Initializing RunInfo..." << std::endl;
}

runinfo::~runinfo(){};

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
void runinfo::AddRunInfo(const std::string &system, const unsigned int &runID_first, const unsigned int &runID_last, const std::string &badmap_ver, const unsigned int &shadowbar, const std::string &trigger)
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
unsigned int runinfo::GetRunIndex(const std::string &system, const unsigned int &iRun)
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
std::string runinfo::GetBadMapVersion(const std::string &system, const unsigned int &iRun)
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
std::string runinfo::GetTriggerCondition(const std::string &system, const unsigned int &iRun)
{
    return (iRun < this->RunNumbers[system]) ? this->TriggerCondition[system][iRun] : "";
}

unsigned int runinfo::GetShadowBar(const std::string &system, const unsigned int &iRun)
{
    return (iRun < this->RunNumbers[system]) ? this->ShadowBars[system][iRun] : 9999;
}