#include "TFile.h"
#include "TChain.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <regex>
#include <filesystem>
namespace fs = std::filesystem;

fs::path PROJECT_DIR = std::getenv("PROJECT_DIR");
constexpr int UBALL_MULTI_MIN = 5;
constexpr double TDC_TRIGGER_THRESHOLD = -9990;

fs::path DIR_DATA = "/data/HiRA_Cali";
fs::path DIR_GEOEFF = PROJECT_DIR / "database/e15190/hira/GeoEff";
fs::path PATH_PIXEL_ANGLES = PROJECT_DIR / "database/e15190/hira/Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat";
fs::path PATH_RUNINFO = PROJECT_DIR / "database/e15190/RunInfo.dat";
fs::path PATH_BMULTI_MAP = PROJECT_DIR / "database/e15190/microball/bimp_mapping";


void read_reaction(const std::string &reaction, int &beamA, int &beamZ, int &targetA, int &targetZ, int &beam_energy)
{
    {
        std::regex pattern("[A-Z][a-z]");
        std::vector<std::string> tokens;
        std::sregex_iterator iter(reaction.begin(), reaction.end(), pattern);
        std::sregex_iterator end;

        while (iter != end)
        {
            std::smatch match = *iter;
            tokens.push_back(match.str());
            ++iter;
        }

        std::string beam = tokens[0];
        std::string target = tokens[1];

        auto GetZ = [](const std::string &nuclei) -> double
        {
            if (nuclei == "Ca")
                return 20;
            else if (nuclei == "Ni")
                return 28;
            else if (nuclei == "Sn")
                return 50;
            else
                return 0;
        };
        beamZ = GetZ(tokens[0]);
        targetZ = GetZ(tokens[1]);
    }

    {
        std::regex pattern("[0-9]+");
        std::vector<int> tokens;

        std::sregex_iterator iter(reaction.begin(), reaction.end(), pattern);
        std::sregex_iterator end;

        while (iter != end)
        {
            std::smatch match = *iter;
            tokens.push_back(std::stoi(match.str()));
            ++iter;
        }
        beamA = tokens[0];
        targetA = tokens[1];
        beam_energy = tokens[2];
    }
}

std::vector<std::string> particles = {
    "p",
    "d",
    "t",
    "3He",
    "4He",
};

std::map<std::array<int,2>, std::string> nz2particle = {
    {{0,1}, "p"},
    {{1,1}, "d"},
    {{2,1}, "t"},
    {{1,2}, "3He"},
    {{2,2}, "4He"},
};

std::vector<std::string> modes = {"raw", "geo", "georeaction"};
std::map<std::string, std::map<std::string, TH2D *>> h2_kinergy_theta_lab;
std::map<std::string, std::map<std::string, TH2D *>> h2_pt_rapidity_lab;

void define_histograms()
{
    for (auto &mode : modes)
    {
        for (auto &particle : particles)
        {
            h2_kinergy_theta_lab[mode][particle] = new TH2D(
                ("h2_kinergy_theta_lab_" + mode + "_" + particle).c_str(),
                "",
                250, 0, 250,
                600, 20, 80);

            h2_pt_rapidity_lab[mode][particle] = new TH2D(
                ("h2_pt_rapidity_lab_" + mode + "_" + particle).c_str(),
                "",
                150, 0, 1.5,
                800, 0, 800);

            h2_kinergy_theta_lab[mode][particle]->Sumw2();
            h2_pt_rapidity_lab[mode][particle]->Sumw2();
            
        }
    }
    return;
}