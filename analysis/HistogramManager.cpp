#include "HistogramManager.hh"

HistogramManager::HistogramManager(const std::string &str)
{
    mStr = str;
    mParticlesNames = this->PARTICLES;
    this->_Initialize();
}

HistogramManager::~HistogramManager()
{
    for (auto &[str, hist] : h2_pta_rapidity_lab)
    {
        if (hist)
        {
            delete hist;
        }
    }
}

void HistogramManager::Fill(const Particle &particle, const double &weight)
{
    std::string name = "";
    if (particle.Z == 1 && particle.N == 0)
    {
        name = "p";
    }
    else if (particle.Z == 1 && particle.N == 1)
    {
        name = "d";
    }
    else if (particle.Z == 1 && particle.N == 2)
    {
        name = "t";
    }
    else if (particle.Z == 2 && particle.N == 1)
    {
        name = "3He";
    }
    else if (particle.Z == 2 && particle.N == 2)
    {
        name = "4He";
    }
    else
    {
        return;
    }
    if (this->h2_pta_rapidity_lab.count(name) == 0)
    {
        return;
    }
    this->h2_pta_rapidity_lab[name]->Fill(particle.rapidity_lab_normed, particle.pmag_trans / (particle.N + particle.Z), weight);
    return;
}

void HistogramManager::Write()
{
    for (auto &p : this->mParticlesNames)
    {
        this->h2_pta_rapidity_lab[p]->Write();
    }
}

void HistogramManager::Normalize(const double &norm)
{
    for (auto &p : this->mParticlesNames)
    {
        this->h2_pta_rapidity_lab[p]->Scale(1. / norm);
    }
}

void HistogramManager::_Initialize()
{
    for (auto &p : this->mParticlesNames)
    {
        this->h2_pta_rapidity_lab[p] = new TH2D(("h2_pta_rapidity_lab_" + mStr + "_" + p).c_str(), "", 100, 0, 1, 800, 0, 800);
        this->h2_pta_rapidity_lab[p]->Sumw2();
        this->h2_pta_rapidity_lab[p]->SetDirectory(0);
    }
}
