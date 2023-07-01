#ifndef HiRA_hh
#define HiRA_hh

#include "AME.hh"
#include "Particle.hh"
#include "HiRA/HiRAStripMap.hh"
#include "HiRA/GeometricEfficiency.hh"
#include "HiRA/ReactionLost.hh"
#include "HiRA/LaserAngles.hh"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

class HiRA
{
public:
    HiRA();
    ~HiRA();

    void Initalize_Coverage(const std::string &path);
    void Initialize_LaserAngles(const std::string &path);
    void Initialize_StripMap(const std::string &path, const std::string &version);
    void Initialize_GeometricEfficiency(const std::string &path);
    void Initialize_ReactionLost(const std::string &path);

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
    std::vector<std::string> ACCEPTED_PARTICLES = {
        "p",
        "d",
        "t",
        "3He",
        "4He",
    };
};

#endif