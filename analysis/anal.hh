

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

    std::vector<std::string> particlenames = {"p", "d", "t", "3He", "4He"};
    std::vector<std::string> atom_name = {"H", "He"};

    std::map<std::string, TH2D *> h2_pta_rapidity_lab;
    std::map<std::string, TH1D *> h1_multi;

    void init();
    void fill(const hira_particle &particle, const double &weight = 1.);
    void normalize(const double &norm);
    void write();
};

void histograms::init()
{

    for (auto &name : this->particlenames)
    {
        this->h2_pta_rapidity_lab[name] = new TH2D(("h2_pta_rapidity_lab_" + this->mode + "_" + name).c_str(), "", 100, 0, 1, 800, 0, 800);
        this->h2_pta_rapidity_lab[name]->Sumw2();
        this->h2_pta_rapidity_lab[name]->SetDirectory(0);
    }
    for (auto &name : this->atom_name)
    {
        this->h1_multi[name] = new TH1D(("h1_multi_" + this->mode + "_" + name).c_str(), "", 10, 0, 10);
        this->h1_multi[name]->Sumw2();
        this->h1_multi[name]->SetDirectory(0);
    }
}

void histograms::fill(const hira_particle &particle, const double &weight)
{

    if (std::find(this->particlenames.begin(), this->particlenames.end(), particle.name) - this->particlenames.begin() == this->particlenames.size())
    {
        return;
    }

    if (particle.zid == 1)
    {
        this->h1_multi["H"]->Fill(particle.aid, weight);
    }
    if (particle.zid == 2)
    {
        this->h1_multi["He"]->Fill(particle.aid, weight);
    }

    this->h2_pta_rapidity_lab[particle.name]->Fill(particle.rapidity_lab / this->beam_rapidity, particle.pt / particle.aid, weight);
}

void histograms::normalize(const double &norm)
{
    for (auto &[name, h2] : this->h2_pta_rapidity_lab)
    {
        h2->Scale(norm);
    }
    for (auto &[name, h1] : this->h1_multi)
    {
        h1->Scale(norm);
    }
}

void histograms::write()
{
    for (auto &[name, h2] : this->h2_pta_rapidity_lab)
    {
        h2->Write();
    }
    for (auto &[name, h1] : this->h1_multi)
    {
        h1->Write();
    }
}
