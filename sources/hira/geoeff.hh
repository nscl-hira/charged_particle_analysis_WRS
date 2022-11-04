#ifndef geoeff_h
#define geoeff_h

#include <iostream>
#include <fstream>

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"

#include "badmap.hh"
#include "angles.hh"

class geoeff
{
public:
    geoeff();
    ~geoeff();

    double ThetaRange[2];
    double ThetaBinSize;
    void SetTheta(double *Range, double BinSize)
    {
        ThetaRange[0] = Range[0];
        ThetaRange[1] = Range[1];
        ThetaBinSize = BinSize;
    }

    badmap *Hira_BadMapper;
    void Initial_Hira_BadMapper(badmap *tem) { Hira_BadMapper = tem; }
    void Cal_GeoEff(TH1D *h1_Count, TH1D *h1_GeoEff, int TotalSimEvtNum);

    angles *HiraPos;
    void Initial_Hira_PosCali(angles *tem) { HiraPos = tem; }
    bool IsApplyPixelAngle;
    void Set_IsApplyPixelAngle(bool tem) { IsApplyPixelAngle = tem; }

    TFile *f1_Results;
    TH2D *h2_WholeHira_Theta_Phi_Lab;
    TH2D *h2_noBadMap_Theta_Phi_Lab;
    TH2D *h2_BadMap_Theta_Phi_Lab;
    TH1D *h1_BadMap_Theta_Lab_HitCount;
    TH1D *h1_BadMap_Theta_Lab_Eff;
    TH1D *h1_noBadMap_Theta_Lab_HitCount;
    TH1D *h1_noBadMap_Theta_Lab_Eff;
    TH1D *h1_HiraHit_Multi;

    // the below is for get the GeoEff from the Histogram.
    void ReadGeoEffHistogram(const std::string &RootFileName);
    double Get_GeoEff(const double &Theta_Lab); // unit: Deg.
    TH1D *Get_ThetaEff_Histogram() { return h1_BadMap_Theta_Lab_Eff; }
};

#endif
