#include "HiRA.hh"

HiRA::HiRA()
{
    mHiRAStripMapReader = 0;
    mGeoEffReader = 0;
    mLaserAnglesReader = 0;
    mReactionLostReader = new ReactionLost();
    mKinergyLabCut = mDefaultKinergyLabCut;
    mThetaLabCut = mDefaultThetaLabCut;
}
HiRA::~HiRA()
{
    delete mHiRAStripMapReader;
    delete mReactionLostReader;
    delete mGeoEffReader;
    delete mLaserAnglesReader;
}

void HiRA::Initialize_LaserAngles(const std::string &path)
{
    if (mLaserAnglesReader)
    {
        delete mLaserAnglesReader;
    }

    if (!fs::exists(path))
    {
        std::string msg = path + " does not exist.";
        throw std::invalid_argument(msg.c_str());
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
        std::string msg = path + " is not a directory.";
        throw std::invalid_argument(msg.c_str());
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

double HiRA::Get_ReactionLost_CorEff(const int &z, const int &a, const double &ekinlab)
{
    return this->mReactionLostReader->Get_ReactionLost_CorEff(z, a, ekinlab);
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

bool HiRA::Pass(const Particle &particle)
{
    std::string name = "";
    if (particle.Z == 1 && particle.N == 0)
    {
        name = "p";
    }
    else if (particle.Z == 1 && particle.N == 1)
    {
        name = "d";
    }
    else if (particle.Z == 1 && particle.N == 2)
    {
        name = "t";
    }
    else if (particle.Z == 2 && particle.N == 1)
    {
        name = "3He";
    }
    else if (particle.Z == 2 && particle.N == 2)
    {
        name = "4He";
    }
    else
    {
        return false;
    }

    if (this->mKinergyLabCut.count(name) == 0 || this->mThetaLabCut.count(name) == 0)
    {
        return false;
    }
    double ekinlab = particle.kinergy_lab / (particle.N + particle.Z);
    double thetalab = particle.theta_lab * TMath::RadToDeg();

    return (
        ekinlab >= this->mKinergyLabCut[name][0] && ekinlab <= this->mKinergyLabCut[name][1] &&
        thetalab >= this->mThetaLabCut[name][0] && thetalab <= this->mThetaLabCut[name][1]);
}
