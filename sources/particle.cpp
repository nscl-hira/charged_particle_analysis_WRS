#include "particle.hh"

void particle::init(const double &betacms)
{
    std::map<std::pair<int, int>, double> mass_table;
    std::map<std::pair<int, int>, std::string> name_table;

    name_table[{1, 0}] = "n";
    name_table[{1, 1}] = "p";
    name_table[{2, 1}] = "d";
    name_table[{3, 1}] = "t";
    name_table[{3, 2}] = "3He";
    name_table[{4, 2}] = "4He";

    mass_table[{1, 0}] = 939.2;     // n
    mass_table[{1, 1}] = 938.272;   // p
    mass_table[{2, 1}] = 1875.61;   // d
    mass_table[{3, 1}] = 2808.9182; // t
    mass_table[{3, 2}] = 2808.387;  // 3He
    mass_table[{4, 2}] = 3727.374;  // 4He

    this->aid = this->nid + this->zid;
    double mass = mass_table[std::make_pair(this->aid, this->zid)];
    this->name = name_table[std::make_pair(this->aid, this->zid)];
    double gamma = 1. / TMath::Sqrt(1 - pow(betacms, 2.));

    this->pt = TMath::Sqrt(pow(this->px, 2.) + pow(this->py, 2.));
    this->plab = TMath::Sqrt(pow(this->pt, 2.) + pow(this->pz, 2.));
    this->ekinlab = TMath::Sqrt(pow(this->plab, 2.) + pow(mass, 2.)) - mass;

    this->pzcms = gamma * (this->pz - betacms * (this->ekinlab + mass));
    this->pcms = TMath::Sqrt(pow(this->pt, 2.) + pow(this->pzcms, 2.));
    this->ekincms = TMath::Sqrt(pow(this->pcms, 2.) + pow(mass, 2.)) - mass;

    this->thetacms = TMath::ATan2(this->pt, this->pzcms) * TMath::RadToDeg();
    this->thetalab = TMath::ATan2(this->pt, this->pz) * TMath::RadToDeg();
    this->phi = TMath::ATan2(this->py, this->px) * TMath::RadToDeg();

    if (this->coordinate == "hira")
    {
        if (this->phi < 0)
        {
            this->phi += 360.;
        }
        if (this->phi > 360)
        {
            this->phi -= 360.;
        }
    }
    else if (this->coordinate == "uball")
    {
        if (this->phi < -18.)
        {
            this->phi += 360.;
        }
        if (this->phi > 342.)
        {
            this->phi -= 360.;
        }
    }
    this->etrans = TMath::Sqrt(pow(this->pt, 2.) + pow(mass, 2.)) - mass;

    this->rapidity_lab = 0.5 * TMath::Log((this->ekinlab + this->pz + mass) / (this->ekinlab - this->pz + mass));

    this->rapidity_cms = 0.5 * TMath::Log((this->ekincms + this->pzcms + mass) / (this->ekincms - this->pzcms + mass));

    this->pseudo_rapidity_lab = 0.5 * TMath::Log((this->plab + this->pz) / (this->plab - this->pz));

    this->pseudo_rapidity_cms = 0.5 * TMath::Log((this->pcms + this->pzcms) / (this->pcms - this->pzcms));
}