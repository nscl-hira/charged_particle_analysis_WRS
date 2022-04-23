#ifndef HiRA_Acceptance_h
#define HiRA_Acceptance_h 1

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

class HiRA_Acceptance : public TObject
{
public:
    HiRA_Acceptance();
    ~HiRA_Acceptance();
    
    double GetThetaMin(int pid) {return mThetaRange[pid][0];}   
    double GetThetaMax(int pid) {return mThetaRange[pid][1];}   
    double GetPunchThroughLow(int pid) {return mPunchthrough[pid][0];}
    double GetPunchThroughHigh(int pid) {return mPunchthrough[pid][1];}
    int GetPid(int zid, int nid);
    string GetParticleName(int pid) {return mParticleName[pid];}
    
private:    
        
    string mParticleName[12] = {"p","d","t","3He","4He","6He","6Li","7Li","8Li","7Be","8Be","10Be"};
    int mParticleZ[12] = {1,1,1,2,2,2,3,3,3,4,4,4};
    int mParticleN[12] = {0,1,2,1,2,4,3,4,5,3,5,6};
    
    double mPunchthrough[12][2] = {{20.,198.},{15.,263./2},{12.,312./3},{20.,200.},{18,200},{13,200},{22,200},{22,200},{22,200},{22,200},{22,200},{22,200}};
    
    double mThetaRange[12][2] = {{30.,75.},{30.,75.},{30.,75.},{30.,75.},{30.,75.},{30.,75.},{30.,75.},{30.,75.},{30.,75.},{30.,75.},{30.,75.},{30.,75.}};
  
ClassDef(HiRA_Acceptance,1)
};

#endif
