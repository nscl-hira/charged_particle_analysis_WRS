#include "hira.hh"

void hira::init()
{

    this->m_reactionlost = new reactionlost();
    this->m_angles = new angles(this->path_angles);

    this->ekinlabcut["p"] = {20.0, 198.0};
    this->ekinlabcut["d"] = {15.0, 263. / 2};
    this->ekinlabcut["t"] = {12.0, 312. / 3};
    this->ekinlabcut["3He"] = {20.0, 200.0};
    this->ekinlabcut["4He"] = {18.0, 200.0};
    for (auto &[name, cut] : this->ekinlabcut)
    {
        this->thetalabcut[name] = {30., 75.};
    }
}

bool hira::pass_badmap(const int &tel, const int &ef, const int &eb)
{
    if (!this->m_badmap)
    {
        throw std::invalid_argument("badmap is not initiated.");
    }
    bool IsBad_Hira = this->m_badmap->IsBad_Hira(tel);
    bool IsBad_CsI = this->m_badmap->IsBad_CsI(tel, ef, eb);
    bool IsBad_StripX = this->m_badmap->IsBad_StripX(tel, ef);
    bool IsBad_StripY = this->m_badmap->IsBad_StripY(tel, eb);
    return (IsBad_Hira == 0 && IsBad_CsI == 0 && IsBad_StripX == 0 && IsBad_StripY == 0);
}

bool hira::pass(const particle &particle)
{
    std::string name = particle.name;
    if (this->ekinlabcut.count(name) == 0 || this->thetalabcut.count(name) == 0)
    {
        return false;
    }
    double ekinlab = particle.ekinlab / particle.aid;
    double thetalab = particle.thetalab;

    return (
        ekinlab >= this->ekinlabcut[name][0] && ekinlab <= this->ekinlabcut[name][1] &&
        thetalab >= this->thetalabcut[name][0] && thetalab <= this->thetalabcut[name][1]);
}
