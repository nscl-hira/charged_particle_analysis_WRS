#include "HiRA.hh"

HiRA::HiRA()
{
    mHiRAStripMapReader = 0;
    mGeoEffReader = 0;
    mLaserAnglesReader = 0;
    mReactionLostReader = 0;
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

void HiRA::Initalize_Coverage(const std::string &path)
{
    json hira_kinematic_cut;
    std::ifstream json_file(path.c_str());
    json_file >> hira_kinematic_cut;

    auto is_inside = [this](const std::string &key) -> bool
    {
        return std::find(this->ACCEPTED_PARTICLES.begin(), this->ACCEPTED_PARTICLES.end(), key) != this->ACCEPTED_PARTICLES.end();
    };

    for (auto &[name, array] : hira_kinematic_cut["thetalab"].items())
    {
        if (!is_inside(name))
            continue;
        std::string ame_name = AME::get_instance()->GetSymbol(name).value();
        mThetaLabCut[ame_name] = array;
    }

    for (auto &[name, array] : hira_kinematic_cut["kinergy_per_nucleon"].items())
    {
        if (!is_inside(name))
            continue;
        std::string ame_name = AME::get_instance()->GetSymbol(name).value();
        mThetaLabCut[ame_name] = array;
    }

    json_file.close();
    return;
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

void HiRA::Initialize_ReactionLost(const std::string &path)
{
    if (mReactionLostReader)
    {
        delete mReactionLostReader;
    }
    mReactionLostReader = new ReactionLost(path);
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
    auto ame = AME::get_instance();
    int A = particle.Z + particle.N;
    if (!ame->IsPhysical(particle.Z, particle.A))
    {
        return false;
    }
    std::string ame_name = AME::get_instance()->GetSymbol(particle.Z, particle.A).value();

    if (this->mKinergyLabCut.count(ame_name) == 0 || this->mThetaLabCut.count(ame_name) == 0)
    {
        return false;
    }
    double ekinlab = particle.kinergy_lab / particle.A;
    double thetalab = particle.theta_lab * TMath::RadToDeg();

    auto accepted = [](const double &v, const std::array<double, 2> &arr) -> bool
    {
        return v >= arr[0] && v <= arr[1];
    };

    return (
        accepted(ekinlab, this->mKinergyLabCut[ame_name]) &&
        accepted(thetalab, this->mThetaLabCut[ame_name])
        //
    );
}
