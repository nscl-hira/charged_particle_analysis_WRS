#ifndef HiRA_hh
#define HiRA_hh

#include "Particle.hh"
#include "HiRA/HiRAStripMap.hh"
#include "HiRA/GeometricEfficiency.hh"
#include "HiRA/ReactionLost.hh"
#include "HiRA/LaserAngles.hh"

class HiRA
{
public:
    HiRA();
    ~HiRA();

    void Initialize_LaserAngles(const std::string &path);
    void Initialize_StripMap(const std::string &path, const std::string &version);
    void Initialize_GeometricEfficiency(const std::string &path);

    double GetTheta(const int &tel, const int &ef, const int &eb);
    double GetPhi(const int &tel, const int &ef, const int &eb);
    double Get_GeometricEfficiency(const double &theta);
    double Get_ReactionLost_CorEff(const int &z, const int &a, const double &ekinlab);

    bool Pass(const Particle &particle);
    bool IsGoodStrip(const int &tel, const int &ef, const int &eb);

private:
    fs::path mDir_HiRAStripMap;
    HiRAStripMap *mHiRAStripMapReader;
    GeometricEfficiency *mGeoEffReader;
    ReactionLost *mReactionLostReader;
    LaserAngles *mLaserAnglesReader;

    std::map<std::string, std::array<double, 2>> mKinergyLabCut;
    std::map<std::string, std::array<double, 2>> mThetaLabCut;

protected:
    std::map<std::string, std::array<double, 2>> mDefaultKinergyLabCut = {
        {"p", {20.0, 198.0}},
        {"d", {15.0, 263. / 2}},
        {"t", {12.0, 312. / 3}},
        {"3He", {20.0, 200.0}},
        {"4He", {18.0, 200.0}},
    };

    std::map<std::string, std::array<double, 2>> mDefaultThetaLabCut = {
        {"p", {30., 75.}},
        {"d", {30., 75.}},
        {"t", {30., 75.}},
        {"3He", {30., 75.}},
        {"4He", {30., 75.}},
    };
};

#endif