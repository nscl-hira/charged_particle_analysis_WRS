

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <map>
#include <filesystem>
namespace fs = std::filesystem;

#include "../sources/root_io.cpp"
#include "../sources/hira_particle.cpp"
#include "../sources/runinfo.cpp"
#include "../sources/hira.cpp"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

fs::path PROJECT_DIR = "/home/kin/Desktop/charged_particle_analysis_WRS";

struct event
{
    int multi_hira, multi_uball;
    double bimp;
    double tdc_uball_ds;
    std::vector<hira_particle> particles;
};

struct eventcut
{
    std::array<int, 2> multicut_hira;
    std::array<int, 2> multicut_uball;
    std::array<double, 2> impact_parameter_cut;

    bool pass(const event &event);
    bool pass_normalization(const event &event);
};

bool eventcut::pass_normalization(const event &event)
{
    int multi_uball = event.multi_uball;
    double bimp = event.bimp;
    double tdc_uball_ds = event.tdc_uball_ds;
    return (
        multi_uball >= this->multicut_uball[0] &&
        multi_uball <= this->multicut_uball[1] &&
        bimp >= this->impact_parameter_cut[0] &&
        bimp < this->impact_parameter_cut[1] &&
        tdc_uball_ds > -9990);
}

bool eventcut::pass(const event &event)
{
    int multi_hira = event.multi_hira;
    int multi_uball = event.multi_uball;
    double bimp = event.bimp;
    return (multi_hira >= this->multicut_hira[0] && multi_hira <= this->multicut_hira[1] && multi_uball >= this->multicut_uball[0] && multi_uball <= this->multicut_uball[1] && bimp >= this->impact_parameter_cut[0] && bimp < this->impact_parameter_cut[1]);
}

struct histograms
{
    std::string mode;
    double betacms, beam_rapidity;
    TH1D *h1_invariant_mass_8Be;

    void init();
    void fill(const event &event);
    void normalize(const double &norm);
    void write();
};

void histograms::init()
{
    std::string hname = "h1_invariant_mass_" + this->mode + "_8Be";
    this->h1_invariant_mass_8Be = new TH1D(hname.c_str(), "", 600, -100, 500);
    this->h1_invariant_mass_8Be->Sumw2();
    this->h1_invariant_mass_8Be->SetDirectory(0);
}

void histograms::fill(const event &event)
{
    std::vector<hira_particle> par_collection;
    for (auto &par : event.particles)
    {
        if (par.zid == 2 && par.aid == 4)
        {
            par_collection.push_back(par);
        }
    }

    if (par_collection.size() > 1)
    {
        for (unsigned int i = 0; i < par_collection.size(); i++)
        {
            for (unsigned int j = i + 1; j < par_collection.size(); j++)
            {

                // pair cut
                int fnumtel1 = par_collection[i].fnumtel;
                int fnumtel2 = par_collection[j].fnumtel;

                int fnumcsi1 = par_collection[i].fnumcsi;
                int fnumcsi2 = par_collection[j].fnumcsi;

                int fnumstripf1 = par_collection[i].fnumstripf;
                int fnumstripf2 = par_collection[j].fnumstripf;

                int fnumstripb1 = par_collection[i].fnumstripb;
                int fnumstripb2 = par_collection[j].fnumstripb;

                if ((fnumtel1 == fnumtel2 && fnumcsi1 == fnumcsi2) ||
                    (fnumtel1 == fnumtel2 && fnumstripf1 == fnumstripf2) ||
                    (fnumtel1 == fnumtel2 && fnumstripb1 == fnumstripb2))
                {
                    continue;
                }

                double px1 = par_collection[i].px;
                double px2 = par_collection[j].px;
                double py1 = par_collection[i].py;
                double py2 = par_collection[j].py;
                double pz1 = par_collection[i].pz;
                double pz2 = par_collection[j].pz;

                double m1 = (pow(par_collection[i].plab, 2.) - pow(par_collection[i].ekinlab, 2.)) / 2. / par_collection[i].ekinlab;

                double m2 = (pow(par_collection[j].plab, 2.) - pow(par_collection[j].ekinlab, 2.)) / 2. / par_collection[j].ekinlab;

                double e1 = par_collection[i].ekinlab + m1;
                double e2 = par_collection[j].ekinlab + m2;
                double eff1 = par_collection[i].geoeff + par_collection[i].reaction_eff;
                double eff2 = par_collection[j].geoeff + par_collection[j].reaction_eff;

                double p1p2 = px1 * px2 + py1 * py2 + pz1 * pz2;
                double M2 = pow(m1, 2.) + pow(m2, 2.) + 2 * (e1 * e2 - p1p2);
                // std::cout << m1 << "\t" << m2 << "\t" << TMath::Sqrt(M2) - 7452.9 << std::endl;

                if (eff1 > 0 && eff2 > 0)
                {
                    double weight = 1. / eff1 / eff2;
                    this->h1_invariant_mass_8Be->Fill(TMath::Sqrt(M2) - 7452.9, weight);
                }
            }
        }
    }
}

void histograms::normalize(const double &norm)
{
    this->h1_invariant_mass_8Be->Scale(norm);
}

void histograms::write()
{
    h1_invariant_mass_8Be->Write();
}
