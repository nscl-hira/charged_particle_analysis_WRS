#ifndef HiRA_hh
#define HiRA_hh

#include <hira10/geometric_efficiency.hpp>
#include <hira10/sillicon_strips.hpp>
#include <hira10/pixel_angles.hpp>
#include <hira10/reactionlost.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

class HiRA
{
public:
    HiRA();
    ~HiRA();

    void Initialize_ReactionLost();
    void Initalize_Coverage(const std::string &path);
    void Initialize_LaserAngles(const std::string &path);
    void Initialize_StripMap(const std::string &path, const std::string &version);
    void Initialize_GeometricEfficiency(const std::string &path);

    double GetTheta(const int &tel, const int &ef, const int &eb);
    double GetPhi(const int &tel, const int &ef, const int &eb);
    double Get_GeometricEfficiency(const double &theta);
    double Get_ReactionLost_CorEff(const int &z, const int &a, const double &Ekin_Lab);

    bool IsGoodStrip(const int &tel, const int &ef, const int &eb);
    bool AcceptTheta(const std::pair<int,int> & za, const double& theta_lab_deg) {return theta_lab_deg >= this->mThetaLabCut[za][0] && theta_lab_deg <= this->mThetaLabCut[za][1];}
    bool AcceptKinergy(const std::pair<int,int> & za, const double& elab_per_nuc) {return elab_per_nuc >= this->mKinergyLabCut[za][0] && elab_per_nuc <= this->mKinergyLabCut[za][1];}

private:
    fs::path mDir_HiRAStripMap;
    HiRAStripMap *mHiRAStripMapReader;
    GeometricEfficiency *mGeoEffReader;
    LaserAngles *mLaserAnglesReader;
    ReactionLost *mReactionLostReader;

    std::map<std::pair<int,int>, std::array<double, 2>> mKinergyLabCut;
    std::map<std::pair<int,int>, std::array<double, 2>> mThetaLabCut;

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