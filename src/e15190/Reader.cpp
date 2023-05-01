#include "Reader.hh"
Reader::Reader(const std::string &reaction, const std::string &dir_data, const std::string &path_runinfo)
{
    this->_Initialize_ReactionSystem(reaction);

    mDirData = dir_data;
    mCurrentRun = 0;

    mRunInfo = new RunInfo(path_runinfo.c_str());

    int NumberOfRun = mRunInfo->GetRunNumber(mReaction);
    for (int iRun = 0; iRun < NumberOfRun; iRun++)
    {
        int RunIndex = mRunInfo->GetRunIndex(mReaction, iRun);

        std::string path_file = Form(
            "%s/%i%s%i%s_%.0fMeVu/CalibratedData_%d.root",
            dir_data.c_str(),
            mBeamA,
            mBeamNuclei.c_str(),
            mTargetA,
            mTargetNuclei.c_str(),
            mBeamEnergy,
            RunIndex);

        mChains[iRun] = new TChain("E15190");
        mChains[iRun]->AddFile(path_file.c_str());
        mRunIndex[iRun] = RunIndex;
        mBadMapVersion[iRun] = mRunInfo->GetBadMapVersion(mReaction, iRun);
        mTriggerCondition[iRun] = mRunInfo->GetTriggerCondition(mReaction, iRun);
    }
}

Reader::~Reader()
{
    delete mRunInfo;
    for (auto &[iRun, chain] : mChains)
    {
        if (chain)
        {
            delete chain;
        }
    }
}

void Reader::_Initialize_ReactionSystem(const std::string &reaction)
{
    mReaction = reaction;
    {
        std::regex pattern("[A-Z][a-z]");
        std::vector<std::string> tokens;
        std::sregex_iterator iter(reaction.begin(), reaction.end(), pattern);
        std::sregex_iterator end;

        while (iter != end)
        {
            std::smatch match = *iter;
            tokens.push_back(match.str());
            ++iter;
        }

        mBeamNuclei = tokens[0];
        mTargetNuclei = tokens[1];

        auto GetZ = [](const std::string &nuclei) -> double
        {
            if (nuclei == "Ca")
                return 20;
            else if (nuclei == "Ni")
                return 28;
            else if (nuclei == "Sn")
                return 50;
            else
                return 0;
        };
        mBeamZ = GetZ(tokens[0]);
        mTargetZ = GetZ(tokens[1]);
    }

    {
        std::regex pattern("[0-9]+");
        std::vector<int> tokens;

        std::sregex_iterator iter(reaction.begin(), reaction.end(), pattern);
        std::sregex_iterator end;

        while (iter != end)
        {
            std::smatch match = *iter;
            tokens.push_back(std::stoi(match.str()));
            ++iter;
        }
        mBeamA = tokens[0];
        mTargetA = tokens[1];
        mBeamEnergy = tokens[2];
    }
}

void Reader::Initialize_Chain(const int &iRun)
{
    TChain *chain = mChains[iRun];
    mCurrentRun = iRun;
    chain->SetBranchAddress("TDCTriggers.uBall_DS_TRG", &this->tdc_trigger_uball_ds);
    chain->SetBranchAddress("TDCTriggers.uBallHiRA_TRG", &this->tdc_trigger_uball_hira);
    chain->SetBranchAddress("TDCTriggers.uBallNW_TRG", &this->tdc_trigger_uball_nw);

    chain->SetBranchAddress("uBall.fmulti", &this->uball_multi);
    chain->SetBranchAddress("HiRA.fmulti", &this->hira_multi);
    chain->SetBranchAddress("HiRA.fAId", &this->hira_A[0]);
    chain->SetBranchAddress("HiRA.fZId", &this->hira_Z[0]);
    chain->SetBranchAddress("HiRA.fMomentum", &this->hira_pmag[0]);

    chain->SetBranchAddress("HiRA.fnumtel", &this->hira_numtel[0]);
    chain->SetBranchAddress("HiRA.fnumstripf", &this->hira_numstripf[0]);
    chain->SetBranchAddress("HiRA.fnumstripb", &this->hira_numstripb[0]);
    chain->SetBranchAddress("HiRA.fnumcsi", &this->hira_numcsi[0]);
    return;
}

void Reader::Clear_Chain(const int &iRun)
{
    delete mChains[iRun];
    mChains[iRun] = 0;
}

int Reader::GetEntry(const int &ievt)
{
    return mChains[mCurrentRun]->GetEntry(ievt);
}
