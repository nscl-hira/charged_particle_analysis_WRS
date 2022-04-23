// The correction function : TF1* f1_ReactionLost_CorEff comes from Sean.

#ifndef HiRA_ReactionLost_h
#define HiRA_ReactionLost_h 1

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
#include "TLegend.h"

using namespace std;
using namespace TMath;

class HiRA_ReactionLost : public TObject
{
public:
  HiRA_ReactionLost();
  ~HiRA_ReactionLost();
  
  TF1* f1_ReactionLost_CorEff[20];
  int ReactionLost_Cor_ParticleNum;
  string ReactionLost_Cor_ParticleName[20];
  int ReactionLost_Cor_Z[20]; int ReactionLost_Cor_A[20];
  void Initial_ReactionLost_CorEff();
  double Get_ReactionLost_CorEff(int Z, int A, double Ekin_Lab);
  TH1D* CorrESpec_perA(string HistoName, TH1D* h1_ESpec_perA_Origin, int Z, int A);
  
  //Draw the reaction lost correction lost correction function in 1 canvas.
  void Draw_ReactionEff();
  
  ClassDef(HiRA_ReactionLost,1)
};

#endif
