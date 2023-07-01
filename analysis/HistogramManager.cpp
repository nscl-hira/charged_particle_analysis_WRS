#include "HistogramManager.hh"

HistogramManager::HistogramManager(const std::string &str)
{
    mStr = str;
    mParticles = this->ACCEPTED_PARTICLES;
    this->_Initialize();
}

HistogramManager::~HistogramManager()
{
    for (auto &[key, hist] : this->h2_pta_rapidity_lab)
    {
        delete hist;
    }

    for (auto &[key, hist] : this->h2_theta_kinery_lab)
    {
        delete hist;
    }
}

void HistogramManager::Fill(const Particle &particle, const double &weight)
{
    int A = particle.N + particle.Z;
    std::string ame_name = AME::get_instance()->GetSymbol(particle.Z, A).value_or("");
    if (this->h2_pta_rapidity_lab.count(ame_name) == 0)
    {
        return;
    }
    this->h2_pta_rapidity_lab[ame_name]->Fill(particle.rapidity_lab_normed, particle.pmag_trans / A, weight);
    this->h2_theta_kinery_lab[ame_name]->Fill(particle.kinergy_lab / A, particle.theta_lab * TMath::RadToDeg(), weight);
    return;
}

void HistogramManager::Write()
{
    for (auto &[name, hist] : this->h2_pta_rapidity_lab)
    {
        hist->Write();
    }

    for (auto &[name, hist] : this->h2_theta_kinery_lab)
    {
        hist->Write();
    }
}

void HistogramManager::Normalize(const double &norm)
{
    for (auto &[name, hist] : this->h2_pta_rapidity_lab)
    {
        hist->Scale(1. / norm);
    }

    for (auto &[name, hist] : this->h2_theta_kinery_lab)
    {
        hist->Scale(1. / norm);
    }
}

void HistogramManager::_Initialize()
{
    auto ame = AME::get_instance();
    for (auto &p : this->mParticles)
    {
        if (!ame->IsPhysical(p))
        {
            continue;
        }
        std::string ame_name = AME::get_instance()->GetSymbol(p).value_or("unknown");
        this->h2_pta_rapidity_lab[ame_name] = new TH2D(("h2_pta_rapidity_lab_" + mStr + "_" + p).c_str(), "", 150, 0, 1.5, 800, 0, 800);

        this->h2_theta_kinery_lab[ame_name] = new TH2D(("h2_theta_kinery_lab_" + mStr + "_" + p).c_str(), "", 250, 0, 250, 600, 20, 80);

        this->h2_pta_rapidity_lab[ame_name]->Sumw2();
        this->h2_theta_kinery_lab[ame_name]->Sumw2();
        this->h2_pta_rapidity_lab[ame_name]->SetDirectory(0);
        this->h2_theta_kinery_lab[ame_name]->SetDirectory(0);
    }
}

void HistogramManager::SetAcceptedParticles(const std::vector<std::string> &particles)
{
    for (auto &p : particles)
    {
        if (std::find(this->ACCEPTED_PARTICLES.begin(), this->ACCEPTED_PARTICLES.end(), p) != this->ACCEPTED_PARTICLES.end())
        {
            this->mParticles.push_back(p);
        }
    }
}
