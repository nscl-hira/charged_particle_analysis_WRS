#ifndef hira_h
#define hira_h

#include "hira/badmap.cpp"
#include "hira/geoeff.cpp"
#include "hira/reactionlost.cpp"
#include "hira/angles.cpp"

#include "hira_particle.hh"
/**
 * @brief A wrapping class for handling all HiRA-related data.
 *
 */
class Hira
{
public:
    Hira();
    ~Hira();

    void Init_Angle(const std::string &path);
    void Init_Badmap(const std::string &path, const std::string &version);
    void Init_GeoEff(const std::string &path);

    double GetTheta(const int &tel, const int &ef, const int &eb);
    double GetPhi(const int &tel, const int &ef, const int &eb);
    double Get_GeoEff(const double &theta);
    double Get_ReactionLost_CorEff(const int &z, const int &a, const double &ekinlab);

    bool pass(const hira_particle &particle);
    bool pass_badmap(const int &tel, const int &ef, const int &eb);

private:
    fs::path dir_badmap, path_angles;
    badmap *m_badmap;
    geoeff *m_geoeff;
    reactionlost *m_reactionlost;
    angles *m_angles;

    std::map<std::string, std::array<double, 2>> ekinlabcut;
    std::map<std::string, std::array<double, 2>> thetalabcut;
};

#endif