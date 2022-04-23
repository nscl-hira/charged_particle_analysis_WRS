#ifndef HiRA_GeoEff_h
#define HiRA_GeoEff_h 1

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "iostream"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "fstream"
#include "TCutG.h"
#include "HiRA_BadMap.hh"
#include "HiRA_PosCali.hh"

//#define NPTOOL 1

#ifdef NPTOOL
#include "THiraData.h"
#include "TInitialConditions.h"
#endif

using namespace std;
using namespace TMath;

class HiRA_GeoEff : public TObject
{
public:
  HiRA_GeoEff();
  ~HiRA_GeoEff();
  void printProgress (double percentage);
  
  double ThetaRange[2];
  double ThetaBinSize;
  void SetTheta(double* Range,double BinSize) { ThetaRange[0] = Range[0]; ThetaRange[1] = Range[1]; ThetaBinSize = BinSize; }
  
  HiRA_BadMap* Hira_BadMapper;
  void Initial_Hira_BadMapper(HiRA_BadMap* tem) { Hira_BadMapper = tem; }
  void Cal_GeoEff(TH1D* h1_Count,TH1D* h1_GeoEff,int TotalSimEvtNum);
  
  HiRA_PosCali* HiraPos;
  void Initial_Hira_PosCali(HiRA_PosCali* tem) { HiraPos = tem; }
  bool IsApplyPixelAngle;
  void Set_IsApplyPixelAngle(bool tem) { IsApplyPixelAngle = tem; }
  
  TFile* f1_Results;
  TH2D* h2_WholeHira_Theta_Phi_Lab;
  TH2D* h2_noBadMap_Theta_Phi_Lab;
  TH2D* h2_BadMap_Theta_Phi_Lab;
  TH1D* h1_BadMap_Theta_Lab_HitCount;
  TH1D* h1_BadMap_Theta_Lab_Eff;
  TH1D* h1_noBadMap_Theta_Lab_HitCount;
  TH1D* h1_noBadMap_Theta_Lab_Eff;
  TH1D* h1_HiraHit_Multi;

  //the below is for producing the GeoEff from the NPTool simulation.
#ifdef NPTOOL
  void ReadSimData(string SimFileName,string RootFileNameForStore);
#endif
  void Draw_Info();
  //the below is for get the GeoEff from the Histogram.
  void ReadGeoEffHistogram(string RootFileName);
  double Get_GeoEff(double Theta_Lab); //unit: Deg.
  TH1D* Get_ThetaEff_Histogram() { return h1_BadMap_Theta_Lab_Eff; }
  
  ClassDef(HiRA_GeoEff,1)
};

#endif
