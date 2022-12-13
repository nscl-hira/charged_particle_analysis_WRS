#include <filesystem>
namespace fs = std::filesystem;
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <array>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"

#include "../sources/runinfo.hh"
#include "../sources/hira.hh"
#include "../sources/root_io.hh"
#include "../sources/hira_particle.hh"

struct DataStructure
{
    static constexpr int max_multi = 128;
    double TDC_hira_ds_nwtdc;
    double TDC_hira_live;
    double TDC_master;
    double TDC_master_nw;
    double TDC_master_vw;
    double TDC_nw_ds;
    double TDC_nw_ds_nwtdc;
    double TDC_rf_nwtdc;
    double TDC_mb_hira;
    double TDC_mb_hira_nwtdc;
    double TDC_mb_nw;
    double TDC_mb_nw_nwtdc;
    double TDC_mb_ds;

    // Microball
    int MB_multi;
    std::array<int, max_multi> MB_ring; // ring number
    std::array<int, max_multi> MB_det;  // detector (crystal) number
    std::array<short, max_multi> MB_tail;
    std::array<short, max_multi> MB_fast;
    std::array<short, max_multi> MB_time;

    // Hira
    int Hira_multi;
    std::array<int, max_multi> Hira_Z;
    std::array<int, max_multi> Hira_A;
    std::array<short, max_multi> Hira_telescope;
    std::array<short, max_multi> Hira_csi;
    std::array<short, max_multi> Hira_stripf;
    std::array<short, max_multi> Hira_stripb;
    std::array<double, max_multi> Hira_momentum;
    std::array<double, max_multi> Hira_theta;
    std::array<double, max_multi> Hira_phi;
};

DataStructure structure;

TChain *input_chain(
    const std::string &path, const std::string tr_name)
{
    TChain *chain = new TChain(tr_name.c_str());
    chain->Add(path.c_str());

    chain->SetBranchAddress("TDCTriggers.HiRA_DS_TRG_NWTDC", &structure.TDC_hira_ds_nwtdc);
    chain->SetBranchAddress("TDCTriggers.HiRA_LIVE", &structure.TDC_hira_live);
    chain->SetBranchAddress("TDCTriggers.MASTER_TRG", &structure.TDC_master);
    chain->SetBranchAddress("TDCTriggers.MASTER_TRG_NWTDC", &structure.TDC_master_nw);
    chain->SetBranchAddress("TDCTriggers.MASTER_TRG_VWTDC", &structure.TDC_master_vw);
    chain->SetBranchAddress("TDCTriggers.NW_DS_TRG", &structure.TDC_nw_ds);
    chain->SetBranchAddress("TDCTriggers.NW_DS_TRG_NWTDC", &structure.TDC_nw_ds_nwtdc);
    chain->SetBranchAddress("TDCTriggers.RF_TRG_NWTDC", &structure.TDC_rf_nwtdc);
    chain->SetBranchAddress("TDCTriggers.uBallHiRA_TRG", &structure.TDC_mb_hira);
    chain->SetBranchAddress("TDCTriggers.uBallHiRA_TRG_NWTDC", &structure.TDC_mb_hira_nwtdc);
    chain->SetBranchAddress("TDCTriggers.uBallNW_TRG", &structure.TDC_mb_nw);
    chain->SetBranchAddress("TDCTriggers.uBallNW_TRG_NWTDC", &structure.TDC_mb_nw_nwtdc);
    chain->SetBranchAddress("TDCTriggers.uBall_DS_TRG", &structure.TDC_mb_ds);
    // Microball
    chain->SetBranchAddress("uBall.fmulti", &structure.MB_multi);
    chain->SetBranchAddress("uBall.fnumring", &structure.MB_ring[0]);
    chain->SetBranchAddress("uBall.fnumdet", &structure.MB_det[0]);
    chain->SetBranchAddress("uBall.fTail", &structure.MB_tail[0]);
    chain->SetBranchAddress("uBall.fFast", &structure.MB_fast[0]);
    chain->SetBranchAddress("uBall.fTime", &structure.MB_time[0]);

    // Hira
    chain->SetBranchAddress("HiRA.fmulti", &structure.Hira_multi);
    chain->SetBranchAddress("HiRA.fZId", &structure.Hira_Z[0]);
    chain->SetBranchAddress("HiRA.fAId", &structure.Hira_A[0]);
    chain->SetBranchAddress("HiRA.fnumtel", &structure.Hira_telescope[0]);
    chain->SetBranchAddress("HiRA.fnumcsi", &structure.Hira_csi[0]);
    chain->SetBranchAddress("HiRA.fnumstripf", &structure.Hira_stripf[0]);
    chain->SetBranchAddress("HiRA.fnumstripb", &structure.Hira_stripb[0]);
    chain->SetBranchAddress("HiRA.fMomentum", &structure.Hira_momentum[0]);

    // enable class objects
    chain->SetMakeClass(1);

    // set branch status
    chain->SetBranchStatus("*", false);
    // TDC triggers
    chain->SetBranchStatus("TDCTriggers.HiRA_DS_TRG_NWTDC", true);
    chain->SetBranchStatus("TDCTriggers.HiRA_LIVE", true);
    chain->SetBranchStatus("TDCTriggers.MASTER_TRG", true);
    chain->SetBranchStatus("TDCTriggers.MASTER_TRG_NWTDC", true);
    chain->SetBranchStatus("TDCTriggers.MASTER_TRG_VWTDC", true);
    chain->SetBranchStatus("TDCTriggers.NW_DS_TRG", true);
    chain->SetBranchStatus("TDCTriggers.NW_DS_TRG_NWTDC", true);
    chain->SetBranchStatus("TDCTriggers.RF_TRG_NWTDC", true);
    chain->SetBranchStatus("TDCTriggers.uBallHiRA_TRG", true);
    chain->SetBranchStatus("TDCTriggers.uBallHiRA_TRG_NWTDC", true);
    chain->SetBranchStatus("TDCTriggers.uBallNW_TRG", true);
    chain->SetBranchStatus("TDCTriggers.uBallNW_TRG_NWTDC", true);
    chain->SetBranchStatus("TDCTriggers.uBall_DS_TRG", true);
    // Microball
    chain->SetBranchStatus("uBall.fmulti", true);
    chain->SetBranchStatus("uBall.fnumring", true);
    chain->SetBranchStatus("uBall.fnumdet", true);
    chain->SetBranchStatus("uBall.fTail", true);
    chain->SetBranchStatus("uBall.fFast", true);
    chain->SetBranchStatus("uBall.fTime", true);
    // Hira
    chain->SetBranchStatus("HiRA.fmulti", true);
    chain->SetBranchStatus("HiRA.fZId", true);
    chain->SetBranchStatus("HiRA.fAId", true);
    chain->SetBranchStatus("HiRA.fnumtel", true);
    chain->SetBranchStatus("HiRA.fnumcsi", true);
    chain->SetBranchStatus("HiRA.fnumstripf", true);
    chain->SetBranchStatus("HiRA.fnumstripb", true);
    chain->SetBranchStatus("HiRA.fMomentum", true);

    return chain;
}

TTree *get_tree(const std::string tr_name = "E15190")
{
    TTree *tr = new TTree(tr_name.c_str(), tr_name.c_str());
    //   TDC triggers
    tr->Branch("TDC_hira_ds_nwtdc", &structure.TDC_hira_ds_nwtdc, "TDC_hira_ds_nwtdc/D");
    tr->Branch("TDC_hira_live", &structure.TDC_hira_live, "TDC_hira_live/D");
    tr->Branch("TDC_master", &structure.TDC_master, "TDC_master/D");
    tr->Branch("TDC_master_nw", &structure.TDC_master_nw, "TDC_master_nw/D");
    tr->Branch("TDC_master_vw", &structure.TDC_master_vw, "TDC_master_vw/D");
    tr->Branch("TDC_nw_ds", &structure.TDC_nw_ds, "TDC_nw_ds/D");
    tr->Branch("TDC_nw_ds_nwtdc", &structure.TDC_nw_ds_nwtdc, "TDC_nw_ds_nwtdc/D");
    tr->Branch("TDC_rf_nwtdc", &structure.TDC_rf_nwtdc, "TDC_rf_nwtdc/D");
    tr->Branch("TDC_mb_hira", &structure.TDC_mb_hira, "TDC_mb_hira/D");
    tr->Branch("TDC_mb_hira_nwtdc", &structure.TDC_mb_hira_nwtdc, "TDC_mb_hira_nwtdc/D");
    tr->Branch("TDC_mb_nw", &structure.TDC_mb_nw, "TDC_mb_nw/D");
    tr->Branch("TDC_mb_nw_nwtdc", &structure.TDC_mb_nw_nwtdc, "TDC_mb_nw_nwtdc/D");
    tr->Branch("TDC_mb_ds", &structure.TDC_mb_ds, "TDC_mb_ds/D");

    // Microball
    tr->Branch("MB_multi", &structure.MB_multi, "MB_multi/I");
    tr->Branch("MB_ring", &structure.MB_ring[0], "MB_ring[MB_multi]/I");
    tr->Branch("MB_det", &structure.MB_det[0], "MB_det[MB_multi]/I");
    tr->Branch("MB_tail", &structure.MB_tail[0], "MB_tail[MB_multi]/S");
    tr->Branch("MB_fast", &structure.MB_fast[0], "MB_fast[MB_multi]/S");
    tr->Branch("MB_time", &structure.MB_time[0], "MB_time[MB_multi]/S");

    // Hira
    tr->Branch("Hira_multi", &structure.Hira_multi, "Hira_multi/I");
    tr->Branch("Hira_Z", &structure.Hira_Z[0], "Hira_Z[Hira_multi]/I");
    tr->Branch("Hira_A", &structure.Hira_A[0], "Hira_A[Hira_multi]/I");
    tr->Branch("Hira_telescope", &structure.Hira_telescope[0], "Hira_telescope[Hira_multi]/S");
    tr->Branch("Hira_csi", &structure.Hira_csi[0], "Hira_csi[Hira_multi]/S");
    tr->Branch("Hira_stripb", &structure.Hira_stripb[0], "Hira_stripb[Hira_multi]/S");
    tr->Branch("Hira_stripf", &structure.Hira_stripf[0], "Hira_stripf[Hira_multi]/S");
    tr->Branch("Hira_momentum", &structure.Hira_momentum[0], "Hira_momentum[Hira_multi]/D");
    tr->Branch("Hira_theta", &structure.Hira_theta[0], "Hira_theta[Hira_multi]/D");
    tr->Branch("Hira_phi", &structure.Hira_phi[0], "Hira_phi[Hira_multi]/D");

    return tr;
}