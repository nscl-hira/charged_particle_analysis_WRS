#include "hira.hh"

Hira::Hira()
{
    m_badmap = 0;
    m_geoeff = 0;
    m_angles = 0;
    m_reactionlost = new reactionlost();

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
Hira::~Hira()
{
    delete m_badmap;
    delete m_reactionlost;
    delete m_geoeff;
    delete m_angles;
}

void Hira::Init_Angle(const std::string &path)
{
    if (m_angles)
    {
        delete m_angles;
    }

    if (!fs::exists(path))
    {
        std::string msg = path + " dost not exist.";
        throw std::invalid_argument(msg.c_str());
    }
    path_angles = path;
    m_angles = new angles(path);
}
void Hira::Init_Badmap(const std::string &path, const std::string &version)
{
    if (m_badmap)
    {
        delete m_badmap;
    }
    dir_badmap = path;
    if (!fs::is_directory(dir_badmap))
    {
        std::string msg = path + " is not a directory.";
        throw std::invalid_argument(msg.c_str());
    }
    m_badmap = new badmap(path, version);
}

void Hira::Init_GeoEff(const std::string &path)
{
    if (m_geoeff)
    {
        delete m_geoeff;
    }
    m_geoeff = new geoeff();
    m_geoeff->ReadGeoEffHistogram(path.c_str());
}

double Hira::GetTheta(const int &tel, const int &ef, const int &eb)
{
    return m_angles->GetTheta(tel, ef, eb);
}

double Hira::GetPhi(const int &tel, const int &ef, const int &eb)
{
    return m_angles->GetPhi(tel, ef, eb);
}

double Hira::Get_GeoEff(const double &theta_deg)
{
    return m_geoeff->Get_GeoEff(theta_deg);
}

double Hira::Get_ReactionLost_CorEff(const int &z, const int &a, const double &ekinlab)
{
    return this->m_reactionlost->Get_ReactionLost_CorEff(z, a, ekinlab);
}

bool Hira::pass_badmap(const int &tel, const int &ef, const int &eb)
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

// bool Hira::pass(const particle &particle)
bool Hira::pass(const particle &particle)
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

// void Hira::reset()
// {
//     delete m_geoeff;
//     delete m_badmap;
// }