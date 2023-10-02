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

#include "runinfo.hpp"

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

    std::string GetBeamNuclei() const { return mBeamNuclei; }
    std::string GetTargetNuclei() const { return mTargetNuclei; }
    int GetBeamZ() const { return mBeamZ; }
    int GetBeamA() const { return mBeamA; }
    int GetTargetZ() const { return mTargetZ; }
    int GetTargetA() const { return mTargetA; }
    double GetBeamEnergy() const { return mBeamEnergy; }

    long GetEntries();
    long GetEntries(const int &iRun) { return mChains[iRun]->GetEntries(); }
    int GetRuns() const { return mChains.size(); }

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

    double *GetHiRA_Kinergy_SiF_Cal() { return &hira_kinergy_sif_cal[0]; }
    double *GetHiRA_Kinergy_SiB_Cal() { return &hira_kinergy_sib_cal[0]; }
    double *GetHiRA_Kinergy_CsI_Cal() { return &hira_kinergy_csi_cal[0]; }
    double* GetHiRA_Kinergy_SiF_Matched() { return &hira_kinergy_sif_matched[0]; }
    double* GetHiRA_Kinergy_SiB_Matched() { return &hira_kinergy_sib_matched[0]; }

    unsigned short* GetHiRA_Kinergy_SiF_High() { return &hira_kinergy_sif_high[0]; }
    unsigned short* GetHiRA_Kinergy_SiF_Low() { return &hira_kinergy_sif_low[0]; }
    unsigned short* GetHiRA_Kinergy_SiB_High() { return &hira_kinergy_sib_high[0]; }
    unsigned short* GetHiRA_Kinergy_SiB_Low() { return &hira_kinergy_sib_low[0]; }
    unsigned short* GetHiRA_Kinergy_CsI() { return &hira_kinergy_csi[0]; }

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

    std::array<unsigned short, MAX_MULTI> hira_kinergy_sif_high;
    std::array<unsigned short, MAX_MULTI> hira_kinergy_sif_low;
    std::array<unsigned short, MAX_MULTI> hira_kinergy_sib_high;
    std::array<unsigned short, MAX_MULTI> hira_kinergy_sib_low;
    std::array<double, MAX_MULTI> hira_kinergy_sif_matched;
    std::array<double, MAX_MULTI> hira_kinergy_sib_matched;
    std::array<double, MAX_MULTI> hira_kinergy_sif_cal;
    std::array<double, MAX_MULTI> hira_kinergy_sib_cal;
    std::array<unsigned short, MAX_MULTI> hira_kinergy_csi;
    std::array<double, MAX_MULTI> hira_kinergy_csi_cal;

};
#endif