#include "Particle.hh"

void Particle::initialize(const double &betacms, const double &beam_rapidity, const double &m)
{
    double A = this->N + this->Z;
    double mass = m;
    if (mass == 0.)
    {
        mass = 938.272 * A;
    }

    this->pmag_trans = TMath::Sqrt(pow(this->px, 2.) + pow(this->py, 2.));
    this->pmag_lab = TMath::Sqrt(pow(pmag_trans, 2.) + pow(this->pz, 2.));
    this->kinergy_lab = TMath::Sqrt(pow(pmag_lab, 2.) + pow(mass, 2.)) - mass;

    this->theta_lab = TMath::ATan2(this->pmag_trans, this->pz) * TMath::RadToDeg();
    this->phi = TMath::ATan2(this->py, this->px) * TMath::RadToDeg();
    this->rapidity_lab = Physics::GetRapidity(this->kinergy_lab, this->pz, mass);
    this->rapidity_lab_normed = this->rapidity_lab / beam_rapidity;
}