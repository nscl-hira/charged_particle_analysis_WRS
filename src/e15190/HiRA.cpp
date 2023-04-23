#include "HiRA.hh"

HiRA::HiRA()
{
    mReader = 0;
    mGeoEffReader = 0;
    mLaserAnglesReader = 0;
    mReactionLostReader = new ReactionLost();
    mKinergyLabCut = mDefaultKinergyLabCut;
    mThetaLabCut = mDefaultThetaLabCut;
}
HiRA::~HiRA()
{
    delete mReader;
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
void HiRA::Initialize_(const std::string &path, const std::string &version)
{
    if (mReader)
    {
        delete mReader;
    }
    mDir_ = path;
    if (!fs::is_directory(mDir_))
    {
        std::string msg = path + " is not a directory.";
        throw std::invalid_argument(msg.c_str());
    }
    mReader = new (path, version);
}

void HiRA::Initialize_GeometricEfficiency(const std::string &path)
{
    if (mGeoEffReader)
    {
        delete mGeoEffReader;
    }
    mGeoEffReader = new GeometricEfficiency();
    mGeoEffReader->ReadGeoEffHistogram(path.c_str());
}

double HiRA::GetTheta(const int &tel, const int &ef, const int &eb)
{
    return mLaserAnglesReader->GetTheta(tel, ef, eb);
}

double HiRA::GetPhi(const int &tel, const int &ef, const int &eb)
{
    return mLaserAnglesReader->GetPhi(tel, ef, eb);
}

double HiRA::Get_GeoEff(const double &theta_deg)
{
    return mGeoEffReader->Get_GeoEff(theta_deg);
}

double HiRA::Get_ReactionLost_CorEff(const int &z, const int &a, const double &ekinlab)
{
    return this->mReactionLostReader->Get_ReactionLost_CorEff(z, a, ekinlab);
}

bool HiRA::IsGoodStrip(const int &tel, const int &ef, const int &eb)
{
    if (!mReader)
    {
        throw std::invalid_argument("badmap is not initiated.");
    }
    bool IsBad_HiRA = mReader->IsBad_HiRA(tel);
    bool IsBad_CsI = mReader->IsBad_CsI(tel, ef, eb);
    bool IsBad_StripX = mReader->IsBad_StripX(tel, ef);
    bool IsBad_StripY = mReader->IsBad_StripY(tel, eb);
    return (IsBad_HiRA == 0 && IsBad_CsI == 0 && IsBad_StripX == 0 && IsBad_StripY == 0);
}

bool HiRA::Pass(const Particle &particle)
{

    if (this->mKinergyLabCut.count(name) == 0 || this->mThetaLabCut.count(name) == 0)
    {
        return false;
    }
    double ekinlab = particle.ekinlab / particle.aid;
    double thetalab = particle.thetalab;

    return (
        ekinlab >= this->mKinergyLabCut[name][0] && ekinlab <= this->mKinergyLabCut[name][1] &&
        thetalab >= this->mThetaLabCut[name][0] && thetalab <= this->mThetaLabCut[name][1] &&
        this->pass_badmap(particle.fnumtel, particle.fnumstripf, particle.fnumstripb));
}
