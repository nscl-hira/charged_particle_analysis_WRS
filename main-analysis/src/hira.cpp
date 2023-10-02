#include "hira.hpp"

HiRA::HiRA()
{
    mHiRAStripMapReader = 0;
    mGeoEffReader = 0;
    mLaserAnglesReader = 0;
    mReactionLostReader = 0;

    mKinergyLabCut = {
        {{1,1},{20,198.0}},
        {{1,2},{15,131.5}},
        {{1,3},{12,104.0}},
        {{2,3},{20,200.0}},
        {{2,4},{18,200.0}},
        {{2,6},{13,200.0}},
    };

    mThetaLabCut = {
        {{1,1},{30,75.0}},
        {{1,2},{30,75.0}},
        {{1,3},{30,75.0}},
        {{2,3},{30,75.0}},
        {{2,4},{30,75.0}},
        {{2,6},{30,75.0}},
    };
}

HiRA::~HiRA()
{
    delete mHiRAStripMapReader;
    delete mGeoEffReader;
    delete mLaserAnglesReader;
    delete mReactionLostReader;
}

void HiRA::Initialize_ReactionLost() {
    if (mReactionLostReader)
    {
        delete mReactionLostReader;
    }
    mReactionLostReader = new ReactionLost();
}

void HiRA::Initialize_LaserAngles(const std::string &path)
{
    if (mLaserAnglesReader)
    {
        delete mLaserAnglesReader;
    }

    if (!fs::exists(path))
    {
        std::cerr << "HiRA::Initialize_LaserAngles: file not found: " << path << std::endl;
        std::exit(1);
    }
    mLaserAnglesReader = new LaserAngles(path);
}

void HiRA::Initialize_StripMap(const std::string &path, const std::string &version)
{
    if (mHiRAStripMapReader)
    {
        delete mHiRAStripMapReader;
    }
    if (!fs::is_directory(path))
    {
        std::cerr << "HiRA::Initialize_StripMap: directory not found: " << path << std::endl;
        std::exit(1);
    }
    mHiRAStripMapReader = new HiRAStripMap(path, version);
}


void HiRA::Initialize_GeometricEfficiency(const std::string &path)
{
    if (mGeoEffReader)
    {
        delete mGeoEffReader;
    }
    mGeoEffReader = new GeometricEfficiency();
    mGeoEffReader->ReadGeometricEfficiencyHistogram(path.c_str());
}

double HiRA::GetTheta(const int &tel, const int &ef, const int &eb)
{
    return mLaserAnglesReader->GetTheta(tel, ef, eb);
}

double HiRA::GetPhi(const int &tel, const int &ef, const int &eb)
{
    return mLaserAnglesReader->GetPhi(tel, ef, eb);
}

double HiRA::Get_GeometricEfficiency(const double &theta_deg)
{
    return mGeoEffReader->Get_GeometricEfficiency(theta_deg);
}

bool HiRA::IsGoodStrip(const int &tel, const int &ef, const int &eb)
{
    if (!mHiRAStripMapReader)
    {
        throw std::invalid_argument("badmap is not initiated.");
    }
    bool IsBad_HiRA = mHiRAStripMapReader->IsBad_Hira(tel);
    bool IsBad_CsI = mHiRAStripMapReader->IsBad_CsI(tel, ef, eb);
    bool IsBad_StripX = mHiRAStripMapReader->IsBad_StripX(tel, ef);
    bool IsBad_StripY = mHiRAStripMapReader->IsBad_StripY(tel, eb);
    return (IsBad_HiRA == 0 && IsBad_CsI == 0 && IsBad_StripX == 0 && IsBad_StripY == 0);
}

double HiRA::Get_ReactionLost_CorEff(const int &z, const int &a, const double &Ekin_Lab)
{
    if (!mReactionLostReader)
    {
        throw std::invalid_argument("reaction lost correction is not initiated.");
    }
    return mReactionLostReader->Get_ReactionLost_CorEff(z, a, Ekin_Lab);
}
