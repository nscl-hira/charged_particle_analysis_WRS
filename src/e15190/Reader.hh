#ifndef Reader_hh
#define Reader_hh

#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <any>
#include <map>
#include <regex>
#include <sstream>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <filesystem>
namespace fs = std::filesystem;

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"

#include "e15190/RunInfo.hh"

class Reader
{
public:
    Reader(const std::string &reaction,
           const std::string &dir_data,
           const std::string &path_runinfo);
    ~Reader();

    void _Initialize_ReactionSystem(const std::string &reaction);
    void Initialize_Chain(const int &iRun);
    void Clear_Chain(const int &iRun);

    std::string GetBeamNuclei() { return mBeamNuclei; }
    std::string GetTargetNuclei() { return mTargetNuclei; }
    int GetBeamZ() { return mBeamZ; }
    int GetBeamA() { return mBeamA; }
    int GetTargetZ() { return mTargetZ; }
    int GetTargetA() { return mTargetA; }
    double GetBeamEnergy() { return mBeamEnergy; }

    long GetEntries(const int &iRun) { return mChains[iRun]->GetEntries(); }
    int GetRuns() { return mChains.size(); }

    int GetEntry(const int &ievt);

    std::map<int, TChain *> mChains;
    std::map<int, int> mRunIndex;
    std::map<int, std::string> mBadMapVersion;
    std::map<int, std::string> mTriggerCondition;

    RunInfo *mRunInfo;

    // Getters
    double GetTDC_Trig_Uball_DS() { return tdc_trigger_uball_ds; }
    double GetTDC_Trig_Uball_HiRA() { return tdc_trigger_uball_hira; }
    double GetTDC_Trig_Uball_NW() { return tdc_trigger_uball_nw; }

    int GetUball_Multi() { return uball_multi; }
    int GetHiRA_Multi() { return hira_multi; }

    int *GetHiRA_A() { return &hira_A[0]; }
    int *GetHiRA_Z() { return &hira_Z[0]; }
    double *GetHiRA_Kinergy() { return &hira_kinergy[0]; }

    int *GetHiRA_NumTel() { return &hira_numtel[0]; }
    int *GetHiRA_NumStripF() { return &hira_numstripf[0]; }
    int *GetHiRA_NumStripB() { return &hira_numstripb[0]; }
    int *GetHiRA_NumCSI() { return &hira_numcsi[0]; }

private:
    fs::path mDirData;

    std::string mReaction;
    std::string mBeamNuclei;
    std::string mTargetNuclei;
    int mBeamZ;
    int mBeamA;
    int mTargetZ;
    int mTargetA;
    double mBeamEnergy;

    int mCurrentRun;

protected:
    const static int MAX_MULTI = 128;

    // tdc trigger
    double tdc_trigger_uball_ds;
    double tdc_trigger_uball_hira;
    double tdc_trigger_uball_nw;

    // microball
    int uball_multi;

    // hira
    int hira_multi;
    std::array<int, MAX_MULTI> hira_A;
    std::array<int, MAX_MULTI> hira_Z;
    std::array<double, MAX_MULTI> hira_kinergy;

    std::array<int, MAX_MULTI> hira_numtel;
    std::array<int, MAX_MULTI> hira_numstripf;
    std::array<int, MAX_MULTI> hira_numstripb;
    std::array<int, MAX_MULTI> hira_numcsi;

    // veto wall
    // neutron wall
};
#endif