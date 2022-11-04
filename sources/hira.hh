#ifndef hira_h
#define hira_h

#include "hira/badmap.cpp"
#include "hira/geoeff.cpp"
#include "hira/reactionlost.cpp"
#include "hira/angles.cpp"

struct hira
{
    fs::path dir_badmap, path_angles;

    badmap *m_badmap;
    geoeff *m_geoeff;
    reactionlost *m_reactionlost;
    angles *m_angles;

    std::map<std::string, std::array<double, 2>> ekinlabcut;
    std::map<std::string, std::array<double, 2>> thetalabcut;

    void init();
    bool pass(const particle &particle);
    bool pass_badmap(const int &tel, const int &ef, const int &eb);
};
#endif